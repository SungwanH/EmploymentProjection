clear 
close all
clc;
digits(50)

% Call parameters
params=PARAMS();
v2struct(params.envr);
%v2struct(params.modl);
%v2struct(params.prod);

RUN_NLPF_HAT_SS = 1; 
RUN_NLPF_HAT    = 1; 
RUN_DGP         = 1; 
RUN_RECUR       = 1;

% Load initial data
data_4_sector=load('DATA/BASE_FOURSECTOR.mat', 'VALjn00', 'Din00','mu0','L0');


%%Using fake data
 data_4_sector = FAKE_DATA(params);

%% Obtain non-linear level outcome in HAT (Obtain initial Steady State)
params_ss=rmfield(params,'prod');
hat_fundamental_ss.T_HAT=params.prod.T_HAT_SS;
hat_fundamental_ss.TIME=TIME_SS;

if RUN_NLPF_HAT_SS ==1
    disp('#################')
    disp('Running NLPF_HAT_SS')

    % Generate the initial value of the model: run the temporary
    % equilibrium once for consistency; everything about it will be stored
    % in the structure named 'temporary_struct'
    temporary_struct. w_guess   = ones(J,N); %initial guess for wage
    temporary_struct.p_guess   = ones(J,N); %initial guess for good prices
    temporary_struct.kappa_hat = ones(J*N,N); % relative change in trade cost
    
    temporary_struct.VALjn00   = data_4_sector.VALjn00;  
    temporary_struct.Ljn_hat00 = ones(J,N);
    temporary_struct.T_hat00   = ones(J,N);
    temporary_struct.Din00=data_4_sector.Din00;
    temporary_struct.Din0=temporary_struct.Din00./sum(temporary_struct.Din00,2);
    [~, ~, ~, temporary_struct.Din00_matched, temporary_struct.X_matched, temporary_struct.VALjn00_matched] =...
        NLPF_TEMP_HAT(params, temporary_struct.VALjn00, temporary_struct.Din00, temporary_struct.kappa_hat, temporary_struct.T_hat00, temporary_struct.Ljn_hat00, ...
        temporary_struct. w_guess, temporary_struct.p_guess);
     
    % Scale initial labor allocation it to biannually
    % for quicker convergence
    %Change to Annual basis
    %[temporary_struct.VV,temporary_struct.D] = eig(data_4_sector.mu0(:,:));   
    %temporary_struct.mu0 = real(temporary_struct.VV * (temporary_struct.D)^8 * inv(temporary_struct.VV));
    temporary_struct.mu0=data_4_sector.mu0(:,:);
    temporary_struct.L00 =data_4_sector.L0(:);

    % start from steady state
    for i=1:500
        temporary_struct.L00 =  temporary_struct.mu0'*temporary_struct.L00;
    end
    temporary_struct.L0=reshape(temporary_struct.L00,J,R);  

    starting_point_ss.VALjn0    = temporary_struct.VALjn00_matched; %Labor compensation (w*L)
    starting_point_ss.Din0      = temporary_struct.Din00_matched; %bilateral trade shares
    starting_point_ss.X0        = temporary_struct.X_matched; %Total expenditure 
    starting_point_ss.L0=temporary_struct.L0;
    starting_point_ss.mu0=temporary_struct.mu0;
    
    initial_guess_ss.v_td=ones(R*(J),TIME_SS);
%    load('DATA/NLPF_HAT_SS.mat', 'eqm_nlpf_HAT_SS'); %one-shot convergence
 %   initial_guess_ss.v_td= eqm_nlpf_HAT_SS.v_td(:,1:length(eqm_nlpf_HAT_SS.v_td));
    
         
    [eqm_nlpf_HAT_SS, approx_nlpf_HAT_SS] = NLPF_HAT(params, starting_point_ss,hat_fundamental_ss,initial_guess_ss);
     save('DATA/NLPF_HAT_SS.mat', 'eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS');     
     clear temporary_struct eqm_nlpf_HAT_SS approx_nlpf_HAT_SS
else
    %loading the equilibrium paths in the baseline economy
    load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS'); 
end





%% Obtain non-linear level outcome in HAT
params_NLPF=rmfield(params,'prod');
hat_fundamental_NLPF.T_HAT=params.prod.T_HAT;
hat_fundamental_NLPF.TIME=TIME;

if RUN_NLPF_HAT ==1
    disp('#################')
    disp('Running NLPF_HAT')
    
     v_td=ones(R*(J),TIME); %Initial guess for the Ys (exp(Vt+1-V)^1/NU)
    load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS')

    starting_point_nlpf.VALjn0    = eqm_nlpf_HAT_SS.VALjn00(:,:,TIME_SS); %Labor compensation (w*L)
    starting_point_nlpf.X0        = eqm_nlpf_HAT_SS.X(:,:,TIME_SS); %Total expenditure     
    starting_point_nlpf.L0=eqm_nlpf_HAT_SS.Ldyn(:,:,TIME_SS);
    starting_point_nlpf.mu0=approx_nlpf_HAT_SS.mu(:,:,TIME_SS);
    starting_point_nlpf.Din0      = approx_nlpf_HAT_SS.pi(:,:,TIME_SS); %bilateral trade shares

    initial_guess_nlpf.v_td=ones(R*(J),TIME); %Initial guess for the Ys (exp(Vt+1-V)^1/NU)
    
%    load('DATA/NLPF_HAT.mat', 'eqm_nlpf_HAT'); %one-shot convergence
%    initial_guess_nlpf.v_td= eqm_nlpf_HAT.v_td(:,1:TIME);
    
    [eqm_nlpf_HAT, approx_nlpf_HAT] = NLPF_HAT(params_NLPF, starting_point_nlpf,hat_fundamental_NLPF,initial_guess_nlpf);

    save('DATA/NLPF_HAT.mat', 'eqm_nlpf_HAT','approx_nlpf_HAT');     
    clear  eqm_nlpf_HAT_SS approx_nlpf_HAT_SS   
else
    %loading the equilibrium paths in the baseline economy
    load('DATA/NLPF_HAT.mat','eqm_nlpf_HAT','approx_nlpf_HAT'); %loading the equilibrium values in the baseline economy
end



%% Obtain DGP path
if RUN_DGP ==1
disp('#################')
disp('Running DGP')
    [eqm_dgp, approx_dgp] = DGP(params, eqm_nlpf_HAT, approx_nlpf_HAT);
else
    load('DATA/DGP.mat', 'eqm_dgp','approx_dgp'); 
end



%% DGP linearize around SS
load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS'); 

params_lin_around_ss=params;
params_lin_around_ss.prod.E_T_hat=params.prod.E_T_hat;
eqm_nlpf_HAT_SS_TIME=eqm_nlpf_HAT_SS;
approx_nlpf_HAT_SS_TIME=approx_nlpf_HAT_SS;

eqm_nlpf_HAT_SS_TIME.v_td=eqm_nlpf_HAT_SS_TIME.v_td(:,end-TIME+1:end);
eqm_nlpf_HAT_SS_TIME.Ldyn=eqm_nlpf_HAT_SS_TIME.Ldyn(:,:,end-TIME+1:end);
eqm_nlpf_HAT_SS_TIME.realwages=eqm_nlpf_HAT_SS_TIME.realwages(:,:,end-TIME+1:end);
eqm_nlpf_HAT_SS_TIME.wf00=eqm_nlpf_HAT_SS_TIME.wf00(:,:,end-TIME+1:end);
eqm_nlpf_HAT_SS_TIME.pf00=eqm_nlpf_HAT_SS_TIME.pf00(:,:,end-TIME+1:end);
eqm_nlpf_HAT_SS_TIME.VALjn00=eqm_nlpf_HAT_SS_TIME.VALjn00(:,:,end-TIME+1:end);
eqm_nlpf_HAT_SS_TIME.X=eqm_nlpf_HAT_SS_TIME.X(:,:,end-TIME+1:end);

approx_nlpf_HAT_SS_TIME.mu=approx_nlpf_HAT_SS_TIME.mu(:,:,end-TIME+1:end);
approx_nlpf_HAT_SS_TIME.pi=approx_nlpf_HAT_SS_TIME.pi(:,:,end-TIME+1:end);
approx_nlpf_HAT_SS_TIME.varrho=approx_nlpf_HAT_SS_TIME.varrho(:,:,end-TIME+1:end);
approx_nlpf_HAT_SS_TIME.chi=approx_nlpf_HAT_SS_TIME.chi(:,:,end-TIME+1:end);
approx_nlpf_HAT_SS_TIME.zeta=approx_nlpf_HAT_SS_TIME.zeta(:,:,end-TIME+1:end);
approx_nlpf_HAT_SS_TIME.lambda=approx_nlpf_HAT_SS_TIME.lambda(:,:,end-TIME+1:end);




[eqm_lin_around_ss, approx_lin_around_ss] = DGP(params_lin_around_ss, eqm_nlpf_HAT_SS_TIME, approx_nlpf_HAT_SS_TIME);



%% Test (nonlinear solution to the belief productivity in the first period)'
params_NLPF_belief=rmfield(params,'prod');
hat_fundamental_NLPF_belief.TIME=TIME;
hat_fundamental_NLPF_belief.T_HAT=ones(J,N,TIME);
for t=1:TIME-1
    hat_fundamental_NLPF_belief.T_HAT(:,:,t+1)=params.prod.T_belief(:,:,t+1,1)./params.prod.T_belief(:,:,t,1); 
end

disp('#################')
disp('Running NLPF_HAT for belief in the first period')
load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS')
starting_point_nlpf_belief.VALjn0    = eqm_nlpf_HAT_SS.VALjn00(:,:,TIME_SS); %Labor compensation (w*L)
starting_point_nlpf_belief.X0        = eqm_nlpf_HAT_SS.X(:,:,TIME_SS); %Total expenditure     
starting_point_nlpf_belief.L0=eqm_nlpf_HAT_SS.Ldyn(:,:,TIME_SS);
starting_point_nlpf_belief.mu0=approx_nlpf_HAT_SS.mu(:,:,TIME_SS);
starting_point_nlpf_belief.Din0      = approx_nlpf_HAT_SS.pi(:,:,TIME_SS); %bilateral trade shares

initial_guess_nlpf_belief.v_td=ones(R*(J),TIME); %Initial guess for the Ys (exp(Vt+1-V)^1/NU)

[eqm_nlpf_HAT_belief, approx_nlpf_HAT_belief] = NLPF_HAT(params_NLPF, starting_point_nlpf_belief,hat_fundamental_NLPF_belief,initial_guess_nlpf_belief);

Ldynamic = permute(sum(eqm_nlpf_HAT.Ldyn,1),[2,3,1]);
Ldynamic_belief = permute(sum(eqm_nlpf_HAT_belief.Ldyn,1),[2,3,1]);
L_belief_dgp = eqm_dgp.L_belief_dgp;
L_belief_agg_dgp = permute(sum(L_belief_dgp,1),[2,3,4,1]);

figure
hold on
title('Labor California (level)')
plot(1:TIME-1,Ldynamic(5,1:TIME-1,1))
plot(1:TIME-1,Ldynamic_belief(5,1:TIME-1,1),'--')
plot(1:TIME-1,L_belief_agg_dgp(5,1:TIME-1,2),':')
legend('Nonlinear PF','Nonlinear Belief','Linear Belief','location','best')
saveas(gcf,'figures/NLPF_TEST.png')


% Temporary: manufacturing employment share
temp1=reshape(sum(eqm_nlpf_HAT.Ldyn(1,CHINA,:),1),1,100);
temp3=reshape(sum(eqm_lin_around_ss.L_belief_dgp(1,CHINA,:,1),1),1,100);

plot(1:100,log(temp1)-log(0.1))
hold on
plot(1:100,reshape(eqm_lin_around_ss.L_hat(1,1,:,1),1,100))

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

