%% Test code for the simulation

clear all
close all
clc;
digits(50)

% Size Test
NN=61;
JJ=3;
TT=50;
%IF T=100, N>51 doesn't work
%IF T=50, N>62 doesn't work
F = zeros(NN,JJ,NN,NN,NN,TT);
%params=PARAMS_SIM(1);
DRAW_SD = 0.01;
NUM_SIM = 50;
params=PARAMS_MPS(1,DRAW_SD,NUM_SIM);
v2struct(params.envr);


%% Switchers for this program; loading data
% If the parameter values are changed, check which one to run!
RUN_NLPF_HAT_SS = 1; 
RUN_NLPF_HAT    = 1; 
RUN_NLPF_DD     = 1; 
%RUN_NLPF_DD_BELIEF = 0; 


% Load initial data
%load('Base_year_three_sectors_US_only.mat', 'VA','Din00','mu0','L0','rw0')
load('Base_year_three_sectors.mat', 'VA','Din00','mu0','L0','rw0')
VALjn00=VA;
rw0 = rw0*100;
data_3_sector=v2struct(VALjn00,Din00,mu0,L0,rw0);
%% Obtain non-linear level outcome in HAT (obtain the path representating initial steady state)
params_ss=rmfield(params,'prod');
hat_fundamental_ss.T_HAT=params.prod.T_HAT_SS;
hat_fundamental_ss.TIME=TIME_SS;

if RUN_NLPF_HAT_SS ==1
    disp('#################')
    disp('Running NLPF_HAT_SS')

    % Generate the initial value of the model: run the temporary
    % equilibrium once for consistency; everything about it will be stored
    % in the structure named 'temporary_struct'
    temporary_struct.w_guess   = ones(J,N); %initial guess for wage
    temporary_struct.p_guess   = ones(J,N); %initial guess for good prices
    temporary_struct.kappa_hat = ones(J*N,N); % relative change in trade cost
    
    temporary_struct.VALjn00   = data_3_sector.VALjn00;  
    temporary_struct.rw0       = data_3_sector.rw0;  
    temporary_struct.Ljn_hat00 = ones(J,N);
    temporary_struct.T_hat00   = ones(J,N);
    temporary_struct.Din00=data_3_sector.Din00;
    temporary_struct.Din0=temporary_struct.Din00./sum(temporary_struct.Din00,2);

    [temporary_struct.wf0_matched, temporary_struct.pf0_matched, temporary_struct.Pf0_matched, temporary_struct.Din00_matched, temporary_struct.X_matched, temporary_struct.VALjn00_matched] =...
        NLPF_TEMP_HAT(params, temporary_struct.VALjn00, temporary_struct.Din00, temporary_struct.kappa_hat, temporary_struct.T_hat00, temporary_struct.Ljn_hat00, ...
        temporary_struct.w_guess, temporary_struct.p_guess);
     
    temporary_struct.mu0=data_3_sector.mu0(:,:);
    temporary_struct.L00 =data_3_sector.L0(:);

    % start from steady state
    for i=1:500
        temporary_struct.L00 =  temporary_struct.mu0'*temporary_struct.L00;
    end
    
    %for t=1:TIME-1
    %    starting_point_ss.rw_lev_base(:,:,t+1) = starting_point_ss.rw_lev_base(:,:,t) .* temporary_struct.wf0_matched(:,:,t+1) ./ temporary_struct.Pf0_matched(:,:,t+1);
    %end

    temporary_struct.L0           =  reshape(temporary_struct.L00,J,R);  
    starting_point_ss.rw0         = data_3_sector.rw0 .* temporary_struct.wf0_matched(:,1:R,1) ./ repmat(temporary_struct.Pf0_matched(:,1:R,1),[J,1,1]);
    starting_point_ss.VALjn0      = temporary_struct.VALjn00_matched; %Labor compensation (w*L)
    starting_point_ss.Din0        = temporary_struct.Din00_matched; %bilateral trade shares
    starting_point_ss.X0          = temporary_struct.X_matched; %Total expenditure 
    starting_point_ss.L0          = temporary_struct.L0;
    starting_point_ss.mu0         = temporary_struct.mu0;
    
    initial_guess_ss.v_td=ones(R*(J),TIME_SS);
    
    [eqm_nlpf_HAT_SS, approx_nlpf_HAT_SS] = NLPF_HAT(params, starting_point_ss,hat_fundamental_ss,initial_guess_ss);
     save('DATA/NLPF_HAT_SS.mat', 'eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS');     
     clear temporary_struct % eqm_nlpf_HAT_SS approx_nlpf_HAT_SS

else

    load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS'); 
end

%% Obtain non-linear level outcome in HAT
params_NLPF=rmfield(params,'prod');
hat_fundamental_NLPF.T_HAT=params.prod.T_HAT_SS; 
hat_fundamental_NLPF.TIME=TIME;
if RUN_NLPF_HAT ==1
    disp('#################')
    disp('Running NLPF_HAT')
    
     v_td=ones(R*(J),TIME);  
    load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS')

    starting_point_nlpf.VALjn0    = eqm_nlpf_HAT_SS.VALjn00(:,:,TIME_SS); %Labor compensation (w*L)
    starting_point_nlpf.rw0       = eqm_nlpf_HAT_SS.rw_lev_base(:,:,TIME_SS); %Labor compensation (w*L)
    starting_point_nlpf.X0        = eqm_nlpf_HAT_SS.X(:,:,TIME_SS); %Total expenditure     
    starting_point_nlpf.L0        = eqm_nlpf_HAT_SS.Ldyn(:,:,TIME_SS);
    starting_point_nlpf.mu0       = approx_nlpf_HAT_SS.mu(:,:,TIME_SS);
    starting_point_nlpf.Din0      = approx_nlpf_HAT_SS.pi(:,:,TIME_SS); %bilateral trade shares

    initial_guess_nlpf.v_td=ones(R*(J),TIME); %Initial guess for the Ys (exp(Vt+1-V)^1/NU)
    
%    load('DATA/NLPF_HAT.mat', 'eqm_nlpf_HAT'); %one-shot convergence
%    initial_guess_nlpf.v_td= eqm_nlpf_HAT.v_td(:,1:TIME);
    
    [eqm_nlpf_HAT, approx_nlpf_HAT] = NLPF_HAT(params_NLPF, starting_point_nlpf,hat_fundamental_NLPF,initial_guess_nlpf);

    save('DATA/NLPF_HAT.mat', 'eqm_nlpf_HAT','approx_nlpf_HAT');     
%    clear  eqm_nlpf_HAT_SS approx_nlpf_HAT_SS   
else
    load('DATA/NLPF_HAT.mat','eqm_nlpf_HAT','approx_nlpf_HAT'); %loading the equilibrium values in the baseline economy
end
%% Nonlinear eqm
params_NLPF=rmfield(params,'prod');
hat_fundamental_cross.T_HAT = ones(J,N,TIME);
for t=1:TIME-1
    hat_fundamental_cross.T_HAT(:,:,t+1)= (params.prod.T(:,:,t+1)./params.prod.T(:,:,t))./ones(J,N,1); 
end
%    hat_fundamental_cross.T_HAT(:,:,2) = params.belief.E_T_hat_test(:,:,2,1)+1;
hat_fundamental_cross.TIME=TIME;
if RUN_NLPF_DD ==1
   disp('#################')
    disp('Running NLPF_DD (With Objective Productivity)')
    
    v_cd=ones(R*(J),TIME); %Initial guess for the Ys (exp(Vt+1-V)^1/NU)
    load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS')
    
    % starting points for counterfactual path
    starting_point_dd.VALjn0    = eqm_nlpf_HAT_SS.VALjn00(:,:,TIME_SS); %Labor compensation (w*L)
    starting_point_dd.X0        = eqm_nlpf_HAT_SS.X(:,:,TIME_SS); %Total expenditure     
    starting_point_dd.L0        = eqm_nlpf_HAT_SS.Ldyn(:,:,TIME_SS);
    starting_point_dd.mu0       = approx_nlpf_HAT_SS.mu(:,:,TIME_SS);
    starting_point_dd.Din0      = approx_nlpf_HAT_SS.pi(:,:,TIME_SS); %bilateral trade shares
    
    %compute dot values from the baseline path level variables
    base_point_dd.mu_base_dot(:,:,1)    = approx_nlpf_HAT.mu(:,:,1)./approx_nlpf_HAT_SS.mu(:,:,TIME_SS);
    for t=1:TIME-1
        base_point_dd.VALjn_base_dot(:,:,t+1) = eqm_nlpf_HAT.VALjn00(:,:,t+1)./eqm_nlpf_HAT.VALjn00(:,:,t);
        base_point_dd.Ldyn_base_dot(:,:,t+1)  = eqm_nlpf_HAT.Ldyn(:,:,t+1)./eqm_nlpf_HAT.Ldyn(:,:,t);
        base_point_dd.Din_base_dot(:,:,t+1)   = approx_nlpf_HAT.pi(:,:,t+1)./approx_nlpf_HAT.pi(:,:,t);
        base_point_dd.mu_base_dot(:,:,t+1)    = approx_nlpf_HAT.mu(:,:,t+1)./approx_nlpf_HAT.mu(:,:,t);
    end
    base_point_dd.mu_base_dot(isnan(base_point_dd.mu_base_dot))   = 0;
    base_point_dd.Din_base_dot(isnan(base_point_dd.Din_base_dot)) = 0;
    base_point_dd.VALjn_base    = eqm_nlpf_HAT.VALjn00;
    base_point_dd.Ldyn_base     = eqm_nlpf_HAT.Ldyn;
    base_point_dd.Din_base      = approx_nlpf_HAT.pi;
    base_point_dd.mu_base       = approx_nlpf_HAT.mu;
    base_point_dd.wf00_base     = eqm_nlpf_HAT.wf00;
    base_point_dd.pf00_base     = eqm_nlpf_HAT.pf00;
    base_point_dd.rw_lev_base   = eqm_nlpf_HAT.rw_lev_base;

    initial_guess_cross.v_cd       = ones(R*(J),TIME); %Initial guess for the Ys exp((V1'-V0')-(V1-V0))^1/NU 

    [eqm_nlpf_dd, approx_nlpf_dd] = NLPF_DD_NEW(params_NLPF, starting_point_dd, base_point_dd, hat_fundamental_cross, initial_guess_cross);

    save('DATA/NLPF_DD.mat', 'eqm_nlpf_dd','approx_nlpf_dd');     
%    clear  eqm_nlpf_HAT_SS approx_nlpf_HAT_SS   
else
    load('DATA/NLPF_DD.mat','eqm_nlpf_dd','approx_nlpf_dd'); %loading the equilibrium values in the counterfactual economy
end


%% Compute the impact of uncertainty using second-order approximation
tic
mat_pbp = MAT_CMEX(params, approx_nlpf_dd);
toc
%mat_pbp = MAT_CMEX(params, approx_nlpf_HAT_SS);
E_T_hat_sim = params.belief.E_T_hat_sim;
L_hat_FO=zeros(J,N,TIME,NUM_SIM);
Ldyn_FO=zeros(J,N,TIME,NUM_SIM);
%{
kappa_hat = zeros(N*J,N,TIME);
L = zeros(R*J,1);
V = zeros(R*J,TIME);
[eqm_temp] = PBP_DYN(params, 1, 1, E_T_hat_sim(:,:,:,2,1,1), kappa_hat, L, V, eqm_nlpf_dd, approx_nlpf_dd, mat_pbp);
figure
hold on
plot(1:TIME,(eqm_temp.L(5,1:TIME)))
[eqm_temp] = PBP_DYN(params, 1, 1, E_T_hat_sim(:,:,:,2,2,1), kappa_hat, L, V, eqm_nlpf_dd, approx_nlpf_dd, mat_pbp);
plot(1:TIME,(eqm_temp.L(5,1:TIME)))
%}
for t1=1:ENDT+1
    if t1==1
        L = zeros(R*J,1); %initial deviation from the path around which we are linearizing
        V = zeros(R*J,TIME);
    else
        L = eqm_SO.L(:,t1);
        V = eqm_SO.v;
    end    
    [eqm_sim, approx_sim] = DGP_SIM(params, t1, L, eqm_nlpf_dd, approx_nlpf_dd, mat_pbp);
    v_check(:,:,:,t1) = eqm_sim.v_sim;
    rw_check(:,:,:,:,t1) = eqm_sim.w_sim - repmat(eqm_sim.P_sim,[J,1,1,1]);
    [mom_dyn, mom_temp] = SO_MOMENTS(params, t1, eqm_sim);
    TT(:,:,:,:,t1) = mom_temp.TT;
    vv(:,:,:,t1) = mom_dyn.vv;
    [SO_temp, SO_dyn] = SO_TERMS(params, t1, eqm_sim, eqm_nlpf_dd, approx_nlpf_dd);
    kappa_hat = zeros(N*J,N,TIME);
    %V = zeros(R*J,TIME);
    [eqm_SO] = PBP_DYN_SO_SIM(params, t1, t1, zeros(J,N,TIME,ENDT+1), kappa_hat, L, V, eqm_nlpf_dd, approx_nlpf_dd, SO_temp, SO_dyn);
    v_FO(:,:,t1) = eqm_SO.FO;
    v_SO(:,:,t1) = eqm_SO.SO;
    w_hat_SO = eqm_SO.w;
    p_hat_SO = eqm_SO.p;
    P_hat_SO = eqm_SO.P;
    v_hat_SO = eqm_SO.v;
    if t1<ENDT+1
        L_hat_FO(:,:,t1+1,:) = (eqm_sim.L_sim(:,:,t1+1,:)); 
        L_hat_SO(:,:,t1+1) = reshape(eqm_SO.L(:,t1+1),J,R,1);
        Ldyn_SO(:,:,t1+1) = exp(L_hat_SO(:,:,t1+1))  .* eqm_nlpf_dd.Ldyn(:,:,t1+1);
        Ldyn_FO(:,:,t1+1,1:NUM_SIM) = exp(eqm_sim.L_sim(:,:,t1+1,1:NUM_SIM))  .* eqm_nlpf_dd.Ldyn(:,:,t1+1);
    else
        L_hat_FO(:,:,t1+1:TIME,:) = (eqm_sim.L_sim(:,:,t1+1:TIME,:)); 
        L_hat_SO(:,:,t1+1:TIME) = reshape(eqm_SO.L(:,t1+1:TIME),J,R,TIME-t1);
        Ldyn_SO(:,:,t1+1:TIME) = exp(L_hat_SO(:,:,t1+1:TIME))  .* eqm_nlpf_dd.Ldyn(:,:,t1+1:TIME);
        Ldyn_FO(:,:,t1+1:TIME,1:NUM_SIM) = exp(eqm_sim.L_sim(:,:,t1+1:TIME,1:NUM_SIM))  .* eqm_nlpf_dd.Ldyn(:,:,t1+1:TIME);
    end
end
mean(v_check(5,2,:,1))
L_hat_SO(:,:,2)
Ldyn_SO(:,:,1) = eqm_nlpf_dd.Ldyn(:,:,1);
for s=1:NUM_SIM
    Ldyn_FO(:,:,1,s) = eqm_nlpf_dd.Ldyn(:,:,1);
end

figure
hold on
title('Second order Labor (SD=0.5) (#sim=50) CRRA:5/ Nonlinear Labor (one sector)')
%plot(1:TIME,permute(Ldyn_FO(1,1,1:TIME,3),[4,3,1,2]),':')
%plot(1:TIME,permute(Ldyn_FO(1,1,1:TIME,2),[4,3,1,2]),':')
plot(1:TIME,permute(Ldyn_FO_mean(1,1,1:TIME),[3,1,2]),':')
plot(1:TIME,permute(Ldyn_SO(1,1,1:TIME),[3,1,2]),'--')
plot(1:TIME,permute(eqm_nlpf_dd.Ldyn(1,1,1:TIME),[3,1,2]))
%plot(1:TIME,permute(Ldyn_FO(1,5,1:TIME,:),[4,3,1,2]),':')
%plot(1:TIME,permute(Ldyn_FO(1,5,1:TIME,2),[4,3,1,2]),':')
plot(1:TIME,permute(Ldyn_FO_mean(1,5,1:TIME),[3,1,2]),':')
plot(1:TIME,permute(Ldyn_SO(1,5,1:TIME),[3,1,2]),'--')
plot(1:TIME,permute(eqm_nlpf_dd.Ldyn(1,5,1:TIME),[3,1,2]))
legend('Region 1 First-order','Region 1 Second-order','Region 1 Nonlinear','Region 5 First-order','Region 5 Second-order','Region 5 Nonlinear','location','best')
%legend('Region 1 First-order (SD=0.5)','Region 1 Second-order (SD=0.5)','Region 1 Nonlinear','Region 5 First-order (SD=0.5)','Region 5 Second-order (SD=0.5)','Region 5 Nonlinear','location','best')
saveas(gcf,'FIGURES/Labor_05_onesect_CRRA_PAR_5.png')
saveas(gcf,'FIGURES/Labor_05_onesect_CRRA_PAR_5.fig')
%{
figure
hold on
title('Second order Labor (SD=0.5) (#sim=50) CRRA/ Nonlinear Labor (one sector)')
%plot(1:TIME, permute(eqm_sim.Ldyn(1,5,1:TIME,1),[3,4,2,1]),'--')
%plot(2:TIME,permute(Ldyn_SO_01(1,5,2:TIME),[3,1,2]),'--')
%plot(1:TIME,permute(Ldyn_FO(1,1,1:TIME,:),[3,4,1,2]),':')
%plot(1:TIME,permute(Ldyn_SO(1,1,1:TIME),[3,1,2]),'--')
%plot(1:TIME,permute(eqm_nlpf_dd.Ldyn(1,1,1:TIME),[3,1,2]))
plot(1:TIME,permute(Ldyn_FO(1,5,1:TIME,1:NUM_SIM),[3,4,1,2]),':')
plot(1:TIME,permute(Ldyn_SO(1,5,1:TIME),[3,1,2]),'--')
plot(1:TIME,permute(eqm_nlpf_dd.Ldyn(1,5,1:TIME),[3,1,2]))
%saveas(gcf,'FIGURES/Labor_06_onesect_CRRA_PAR_5_FO.png')
%saveas(gcf,'FIGURES/Labor_06_onesect_CRRA_PAR_5_FO.fig')
%}
stop

figure
hold on
plot(1:TIME,permute(Ldyn_FO(1,1,1:TIME,:),[3,4,1,2]),':')
figure
hold on
title('Second order Labor (SD=0.5) / Nonlinear Labor (one sector)')
%plot(1:TIME, permute(eqm_sim.Ldyn(1,5,1:TIME,1),[3,4,2,1]),'--')
%plot(2:TIME,permute(Ldyn_SO_01(1,5,2:TIME),[3,1,2]),'--')
plot(1:TIME,permute(Ldyn_SO(1,1,1:TIME),[3,1,2]),'--')
plot(1:TIME,permute(eqm_nlpf_dd.Ldyn(1,1,1:TIME),[3,1,2]))
plot(1:TIME,permute(Ldyn_SO(1,5,1:TIME),[3,1,2]),'--')
plot(1:TIME,permute(eqm_nlpf_dd.Ldyn(1,5,1:TIME),[3,1,2]))
legend('Region 1 Second-order (SD=0.5)','Region 1 Nonlinear','Region 5 Second-order (SD=0.5)','Region 5 Nonlinear','location','best')
saveas(gcf,'FIGURES/Labor_sd05_onesect.png')
saveas(gcf,'FIGURES/Labor_sd05_onesect.fig')

save('DATA/SO_Dev_from_05_onesect.mat', 'Ldyn_SO'); 
