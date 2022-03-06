
clear 
close all
clc;
digits(50)

% call parameters
params=PARAMS();
v2struct(params.envr);
v2struct(params.modl);
v2struct(params.prod);

RUN_NLPF_HAT_SS = 1; 
RUN_NLPF_HAT    = 1; 
RUN_DGP         = 1; 
RUN_RECUR       = 1;

%% Obtain non-linear level outcome in HAT (Obtain initial Steady State)
if RUN_NLPF_HAT_SS ==1
disp('#################')
disp('Running NLPF_HAT_SS')
    load('DATA/BASE_FOURSECTOR.mat','mu0','L0')    

    %Change to Biannual basis
    [VV,D] = eig(mu0(:,:));
    mu0 = real(VV * (D)^8 * inv(VV));
    L00 =L0(:);
    % start from steady state
    for i=1:500
        L00 =  mu0'*L00;
    end
    L0=reshape(L00,J,R);

    v_td=ones(R*(J),TIME_SS); %Initial guess for the Ys (exp(Vt+1-V)^1/NU)
    %load('DATA/NLPF_HAT_SS.mat', 'eqm_nlpf_HAT_SS'); %one-shot convergence
    %v_td(:,1:length(eqm_nlpf_HAT_SS.v_td)) = eqm_nlpf_HAT_SS.v_td(:,1:length(eqm_nlpf_HAT_SS.v_td));
    
    initial_SS.L0 = L0;
    initial_SS.mu0 = mu0;
    initial_SS.v_td = v_td;
    initial_SS.T_HAT = T_HAT_SS;
    initial_SS.TIME = TIME_SS;
    SS =1;
    [eqm_nlpf_HAT_SS, approx_nlpf_HAT_SS] = NLPF_HAT(params, initial_SS, SS);
else
    %loading the equilibrium paths in the baseline economy
    load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS'); 
end

%% Obtain non-linear level outcome in HAT
if RUN_NLPF_HAT ==1
    disp('#################')
    disp('Running NLPF_HAT')
    L0 = eqm_nlpf_HAT_SS.Ldyn(:,:,TIME_SS);
    mu0 = approx_nlpf_HAT_SS.mu(:,:,TIME_SS);
    %clear eqm_nlpf_HAT_SS
    %clear approx_nlpf_HAT_SS

    v_td=ones(R*(J),TIME); %Initial guess for the Ys (exp(Vt+1-V)^1/NU)
    load('DATA/NLPF_HAT.mat', 'eqm_nlpf_HAT'); %one-shot convergence
    v_td(:,1:TIME) = eqm_nlpf_HAT.v_td(:,1:TIME);
    
    initial.L0 = L0;
    initial.mu0 = mu0;
    initial.v_td = v_td;
    initial.T_HAT = T_HAT;
    SS =0;
    [eqm_nlpf_HAT, approx_nlpf_HAT] = NLPF_HAT(params, initial, SS);
else
    %loading the equilibrium paths in the baseline economy
    load('DATA/NLPF_HAT.mat','eqm_nlpf_HAT','approx_nlpf_HAT'); %loading the equilibrium values in the baseline economy
end



%% Obtain DGP path
if RUN_DGP ==1
disp('#################')
disp('Running DGP')
    [eqm_dgp, approx_dgp] = DGP(params, W_TRUE, eqm_nlpf_HAT, approx_nlpf_HAT);
else
    load('DATA/DGP.mat', 'eqm_dgp','approx_dgp'); 
end

%{
%% Test (nonlinear solution to the belief productivity in the first period)
T_HAT_belief = ones(J,N,TIME);
for t=1:TIME-1
    T_HAT_belief(:,:,t+1)=T_belief(:,:,t+1,1)./T_belief(:,:,t,1); %time difference in belief
end
load('DATA/NLPF_HAT_BELIEF.mat','eqm_nlpf_HAT_belief');
initial_belief.L0 = eqm_nlpf_HAT_SS.Ldyn(:,:,TIME_SS);
initial_belief.mu0 = approx_nlpf_HAT_SS.mu(:,:,TIME_SS);
initial_belief.T_HAT = T_HAT_belief;
initial_belief.v_td = eqm_nlpf_HAT_belief.v_td;
SS=2;
[eqm_nlpf_HAT_belief, approx_nlpf_HAT_belief] = NLPF_HAT(params, initial_belief, SS);

Ldynamic = permute(sum(eqm_nlpf_HAT.Ldyn,1),[2,3,1]);
Ldynamic_belief = permute(sum(eqm_nlpf_HAT_belief.Ldyn,1),[2,3,1]);
L_belief_dgp = eqm_dgp.L_belief_dgp;
L_belief_agg_dgp = permute(sum(L_belief_dgp,1),[2,3,4,1]);

figure
hold on
title('Labor California (level)')
plot(1:TIME-1,Ldynamic(5,1:TIME-1,1))
plot(1:TIME-1,Ldynamic_belief(5,1:TIME-1,1),'--')
plot(1:TIME-1,L_belief_agg_dgp(5,1:TIME-1,1,1),':')
legend('Nonlinear PF','Nonlinear Belief','Linear Belief','location','best')
saveas(gcf,'figures/NLPF_TEST.png')
%}


%% Obtain Period by period DGP & PF deviation
if RUN_RECUR ==1
disp('#################')
disp('Running RECURSIVE')
    initial_recur.v_hat = eqm_dgp.v_hat;
    initial_recur.w_hat = eqm_dgp.w_hat;
    [eqm_recur] = RECURSIVE(params, W_TRUE, initial_recur, eqm_dgp, approx_dgp);
else
    load('DATA/RECURSIVE.mat', 'eqm_recur'); 
end

%% figures
FIGURES(params, eqm_nlpf_HAT_SS, eqm_nlpf_HAT, eqm_dgp, eqm_recur)

