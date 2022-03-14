clear 
close all
clc;
digits(50)

% Call parameters
params=PARAMS();
v2struct(params.envr);
%v2struct(params.modl);
v2struct(params.prod);

RUN_NLPF_HAT_SS     = 1; 
RUN_NLPF_HAT        = 1; 
RUN_NLPF_CROSS_TIME = 1; 
RUN_DGP             = 1; 
RUN_RECUR           = 1;


figure
hold on
title('Objective Productivity and Belief (sector 1)')
plot(1:TIME, permute(T(1,CHINA,1:TIME),[2,3,1]))
plot(1:TIME, permute(T_belief(1,CHINA,1:TIME,1),[2,3,4,1]),'--')
plot(2:TIME, permute(T_belief(1,CHINA,2:TIME,2),[2,3,4,1]),'--')
plot(5:TIME, permute(T_belief(1,CHINA,5:TIME,5),[2,3,4,1]),'--')
%plot(7:TIME, permute(T_belief(1,CHINA,7:TIME,7),[2,3,4,1]),'--')
plot(10:TIME, permute(T_belief(1,CHINA,10:TIME,10),[2,3,4,1]),'--')
legend('Productivity','Belief at t=1','Belief at t=2','location','best')

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
    temporary_struct.L0      =  reshape(temporary_struct.L00,J,R);  

    starting_point_ss.VALjn0    = temporary_struct.VALjn00_matched; %Labor compensation (w*L)
    starting_point_ss.Din0      = temporary_struct.Din00_matched; %bilateral trade shares
    starting_point_ss.X0        = temporary_struct.X_matched; %Total expenditure 
    starting_point_ss.L0        = temporary_struct.L0;
    starting_point_ss.mu0       = temporary_struct.mu0;
    
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
%hat_fundamental_NLPF.T_HAT=params.prod.T_HAT;
hat_fundamental_NLPF.T_HAT=params.prod.T_HAT_SS;
hat_fundamental_NLPF.TIME=TIME;

if RUN_NLPF_HAT ==1
    disp('#################')
    disp('Running NLPF_HAT')
    
     v_td=ones(R*(J),TIME); %Initial guess for the Ys (exp(Vt+1-V)^1/NU)
    load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS')

    starting_point_nlpf.VALjn0    = eqm_nlpf_HAT_SS.VALjn00(:,:,TIME_SS); %Labor compensation (w*L)
    starting_point_nlpf.X0        = eqm_nlpf_HAT_SS.X(:,:,TIME_SS); %Total expenditure     
    starting_point_nlpf.L0        = eqm_nlpf_HAT_SS.Ldyn(:,:,TIME_SS);
    starting_point_nlpf.mu0       = approx_nlpf_HAT_SS.mu(:,:,TIME_SS);
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

%% Obtain non-linear Cross-Time difference (Counterfactual Objective productivity)
% This part derives counterfactual impact of productivity shock
params_NLPF=rmfield(params,'prod');
hat_fundamental_cross.T_HAT = ones(J,N,TIME);
for t=1:TIME-1
    hat_fundamental_cross.T_HAT(:,:,t+1)=(params.prod.T(:,:,t+1)./params.prod.T(:,:,t))./ones(J,N,1); 
end
hat_fundamental_cross.TIME=TIME;

if RUN_NLPF_CROSS_TIME ==1
    disp('#################')
    disp('Running NLPF_CROSS (With Objective Productivity)')
    
    v_cd=ones(R*(J),TIME); %Initial guess for the Ys (exp(Vt+1-V)^1/NU)
    load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS')

    starting_point_cross.VALjn0    = eqm_nlpf_HAT_SS.VALjn00(:,:,TIME_SS); %Labor compensation (w*L)
    starting_point_cross.X0        = eqm_nlpf_HAT_SS.X(:,:,TIME_SS); %Total expenditure     
    starting_point_cross.L0        = eqm_nlpf_HAT_SS.Ldyn(:,:,TIME_SS);
    starting_point_cross.mu0       = approx_nlpf_HAT_SS.mu(:,:,TIME_SS);
    starting_point_cross.Din0      = approx_nlpf_HAT_SS.pi(:,:,TIME_SS); %bilateral trade shares
    
    base_point_cross.VALjn_base    = eqm_nlpf_HAT.VALjn00;
    base_point_cross.Ldyn_base     = eqm_nlpf_HAT.Ldyn;
    base_point_cross.Din_base      = approx_nlpf_HAT.pi;
    base_point_cross.mu_base       = approx_nlpf_HAT.mu;
    base_point_cross.wf00_base     = eqm_nlpf_HAT.wf00;
    base_point_cross.pf00_base     = eqm_nlpf_HAT.pf00;

    initial_guess_cross.v_cd       = ones(R*(J),TIME); %Initial guess for the Ys exp((V1'-V0')-(V1-V0))^1/NU 
    
 %   load('DATA/NLPF_HAT.mat', 'eqm_nlpf_HAT'); %one-shot convergence
 %   initial_guess_nlpf.v_td= eqm_nlpf_HAT.v_td(:,1:TIME);
    
    [eqm_nlpf_cross_time, approx_nlpf_cross_time] = NLPF_CROSS_TIME(params_NLPF, starting_point_cross, base_point_cross, hat_fundamental_cross, initial_guess_cross);

    save('DATA/NLPF_CROSS.mat', 'eqm_nlpf_cross_time','approx_nlpf_cross_time');     
    clear  eqm_nlpf_HAT_SS approx_nlpf_HAT_SS   
else
    %loading the equilibrium paths in the baseline economy
    load('DATA/NLPF_CROSS.mat','eqm_nlpf_cross_time','approx_nlpf_cross_time'); %loading the equilibrium values in the baseline economy
end

%% Obtain non-linear Cross-Time difference (Counterfactual 'Belief')
% This part derives counterfactual impact of productivity shock
params_NLPF=rmfield(params,'prod');
hat_fundamental_cross_belief.T_HAT = ones(J,N,TIME);
for t=1:TIME-1
%    hat_fundamental_cross.T_HAT(:,:,t+1)=(params.prod.T(:,:,t+1)./params.prod.T(:,:,t))./ones(J,N,1); 
%    hat_fundamental_cross_belief.T_HAT(:,:,t+1)=(params.prod.T_belief(:,:,t+1,1)./params.prod.T_belief(:,:,t,1))./(params.prod.T(:,:,t+1)./params.prod.T(:,:,t)); 
    hat_fundamental_cross_belief.T_HAT(:,:,t+1)=(params.prod.T_belief(:,:,t+1,1)./params.prod.T_belief(:,:,t,1))./ones(J,N,1); 
end
hat_fundamental_cross_belief.TIME=TIME;

if RUN_NLPF_CROSS_TIME ==1
    disp('#################')
    disp('Running NLPF_CROSS (With Belief)')
    
    load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS')

    starting_point_cross_belief.VALjn0    = eqm_nlpf_HAT_SS.VALjn00(:,:,TIME_SS); %Labor compensation (w*L)
    starting_point_cross_belief.X0        = eqm_nlpf_HAT_SS.X(:,:,TIME_SS); %Total expenditure     
    starting_point_cross_belief.L0        = eqm_nlpf_HAT_SS.Ldyn(:,:,TIME_SS);
    starting_point_cross_belief.mu0       = approx_nlpf_HAT_SS.mu(:,:,TIME_SS);
    starting_point_cross_belief.Din0      = approx_nlpf_HAT_SS.pi(:,:,TIME_SS); %bilateral trade shares
    
    base_point_cross_belief.VALjn_base    = eqm_nlpf_HAT.VALjn00;
    base_point_cross_belief.Ldyn_base     = eqm_nlpf_HAT.Ldyn;
    base_point_cross_belief.Din_base      = approx_nlpf_HAT.pi;
    base_point_cross_belief.mu_base       = approx_nlpf_HAT.mu;
    base_point_cross_belief.wf00_base     = eqm_nlpf_HAT.wf00;
    base_point_cross_belief.pf00_base     = eqm_nlpf_HAT.pf00;

    initial_guess_cross_belief.v_cd       = ones(R*(J),TIME); %Initial guess for the Ys exp((V1'-V0')-(V1-V0))^1/NU 
    
 %   load('DATA/NLPF_HAT.mat', 'eqm_nlpf_HAT'); %one-shot convergence
 %   initial_guess_nlpf.v_td= eqm_nlpf_HAT.v_td(:,1:TIME);
    
    [eqm_nlpf_cross_time_belief, approx_nlpf_cross_time_belief] = NLPF_CROSS_TIME(params_NLPF, starting_point_cross_belief, base_point_cross_belief, hat_fundamental_cross_belief, initial_guess_cross_belief);

    save('DATA/NLPF_CROSS_BELIEF.mat', 'eqm_nlpf_cross_time_belief','approx_nlpf_cross_time_belief');     
    clear  eqm_nlpf_HAT_SS approx_nlpf_HAT_SS   
else
    %loading the equilibrium paths in the baseline economy
    load('DATA/NLPF_CROSS_BELIEF.mat','eqm_nlpf_cross_time_belief','approx_nlpf_cross_time_belief'); %loading the equilibrium values in the baseline economy
end


%% Obtain DGP path
if RUN_DGP ==1
disp('#################')
disp('Running DGP')
%    [eqm_dgp, approx_dgp] = DGP(params, eqm_nlpf_HAT, approx_nlpf_HAT);
    [eqm_dgp, approx_dgp] = DGP(params, eqm_nlpf_cross_time, approx_nlpf_cross_time);
else
    load('DATA/DGP.mat', 'eqm_dgp','approx_dgp'); 
end


%{
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

%}

Ldynamic = permute(sum(eqm_nlpf_HAT.Ldyn,1),[2,3,1]);
Ldynamic_sec = permute(eqm_nlpf_HAT.Ldyn(:,CHINA,:),[1,3,2]);
Ldynamic_sec_cross = permute(eqm_nlpf_cross_time.Ldyn(:,CHINA,:),[1,3,2]);
Ldynamic_sec_cross_belief = permute(eqm_nlpf_cross_time_belief.Ldyn(:,CHINA,:),[1,3,2]);
Ldynamic_dgp = permute(eqm_dgp.Ldyn,[2,3,1]);
%Ldynamic_sec_belief = permute(eqm_nlpf_HAT_belief.Ldyn(:,CHINA,:),[1,3,2]);
LdynamicManu= reshape(sum(eqm_nlpf_HAT.Ldyn(1,:,:),2),TIME,1);
LdynamicManu_cross= reshape(sum(eqm_nlpf_cross_time.Ldyn(1,:,:),2),TIME,1);
%LdynamicManu_belief= reshape(sum(eqm_nlpf_HAT_belief.Ldyn(1,:,:),2),TIME,1);
L_belief_dgp = eqm_dgp.L_belief_dgp;
L_belief_agg_dgp = permute(sum(L_belief_dgp,1),[2,3,4,1]);
LdynamicManu_dgp= reshape(sum(L_belief_dgp(1,:,:,1),2),TIME,1);

figure
hold on
title('DATA: Labor in region 1')
plot(1:TIME,Ldynamic_sec_cross(1,1:TIME))
%plot(1:TIME,Ldynamic_sec(1,1:TIME))
plot(1:TIME,Ldynamic_dgp(1,1:TIME),'--')
plot(2:TIME,L_belief_agg_dgp(1,2:TIME,2),':')
plot(3:TIME,L_belief_agg_dgp(1,3:TIME,3),':')
plot(4:TIME,L_belief_agg_dgp(1,4:TIME,4),':')
%plot(10:TIME,L_belief_agg_dgp(1,10:TIME,10),':')
legend('Nonlin PF','DATA','Belief','location','best')
saveas(gcf,'figures/dgp_labor_region1.png')

figure
hold on
title('CHINA Employment')
plot(1:TIME,Ldynamic_sec(1,1:TIME))
plot(1:TIME,Ldynamic_sec_cross(1,1:TIME),'--')
plot(1:TIME,Ldynamic_sec_cross_belief(1,1:TIME),'--')
%plot(1:TIME,Ldynamic_sec_belief(1,1:TIME),'--')
plot(1:TIME,L_belief_agg_dgp(1,1:TIME,2),':')
%plot(1:TIME,L_belief_agg_dgp(5,1:TIME,1),':')
legend('Baseline(Nonlinear TimeDiff)','Nonlinear Cross-Time(Objective T)','Nonlinear Cross-Time(Belief)','Linear Belief','location','best')
saveas(gcf,'figures/NLPF_All_Region1.png')

%{
logdiff = -(log(Ldynamic_sec(1,:)) - log(Ldynamic_sec_cross)); 
logdiff_belief = -(log(Ldynamic_sec(1,:)) - log(L_belief_agg_dgp(:,:,2))); 

figure
hold on
title('change in log percentage CHINA Employment - NLPF and Linear Approx')
%plot(1:TIME,Ldynamic_sec(1,1:TIME))
plot(1:TIME,logdiff(1,1:TIME))
plot(1:TIME,logdiff_belief(1,1:TIME),'--')
%plot(1:TIME,L_belief_agg_dgp(5,1:TIME,1),':')
legend('Nonlinear Time-Cross(Objective T)','Linear Approximation','location','best')
saveas(gcf,'figures/Comparison_Approx_logdiff.png')
saveas(gcf,'figures/Comparison_Approx_logdiff.fig')

figure
hold on
title('CHINA Employment - NLPF and Linear Approx')
%plot(1:TIME,Ldynamic_sec(1,1:TIME))
plot(1:TIME,Ldynamic_sec_cross(1,1:TIME))
plot(1:TIME,L_belief_agg_dgp(1,1:TIME,2),'--')
%plot(1:TIME,L_belief_agg_dgp(5,1:TIME,1),':')
legend('Nonlinear Time-Cross(Objective T)','Linear Approximation','location','best')
saveas(gcf,'figures/Comparison_Approx_Sym.png')
saveas(gcf,'figures/Comparison_Approx_Sym.fig')


% Temporary: manufacturing employment share
temp1=reshape(sum(eqm_nlpf_HAT.Ldyn(1,CHINA,:),1),1,100);
temp3=reshape(sum(eqm_lin_around_ss.L_belief_dgp(1,CHINA,:,1),1),1,100);

figure
plot(1:100,log(temp1)-log(0.1))
hold on
plot(1:100,reshape(eqm_lin_around_ss.L_hat(1,1,:,1),1,100),'--')
plot(1:TIME,eqm_dgp.L_hat(1,1:TIME),'--')
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
Ldynamic_sec_cross = permute(eqm_nlpf_cross_time.Ldyn(1,:,:),[1,3,2]);
Ldyn_dgp = eqm_dgp.Ldyn;
Lhat_dgp = eqm_dgp.L_hat;
Ldynamic_dgp = permute(sum(Ldyn_dgp,1),[2,3,1]);
L_belief_dgp = eqm_dgp.L_belief_dgp;
L_belief_agg_dgp = permute(sum(L_belief_dgp,1),[2,3,4,1]);

L_pf_recur = eqm_recur.L_pf_lev;
L_belief_recur = eqm_recur.L_belief_lev;
L_pf_agg_recur = permute(sum(L_pf_recur,1),[2,3,4,1]);
L_belief_agg_recur = permute(sum(L_belief_recur,1),[2,3,4,1]);

figure
hold on
title('RECURSIVE: Labor in region 1')
plot(1:TIME,Ldynamic_sec_cross(1,1:TIME,1))
plot(1:TIME,Ldynamic_dgp(1,1:TIME))
plot(2:TIME,L_belief_agg_recur(1,2:TIME,2),':')
plot(2:TIME,L_belief_agg_dgp(1,2:TIME,2))
plot(2:TIME,L_pf_agg_recur(1,2:TIME,2),'--')
%plot(10:TIME,L_belief_agg_recur(1,10:TIME,10),':')
%plot(10:TIME,L_belief_agg_dgp(1,10:TIME,10))
%plot(10:TIME,L_pf_agg_recur(1,10:TIME,10),'--')
legend('Nonlin PF','DATA','Belief (RECOVERED)','Belief (DGP)','PF (RECOVERED)','location','best')
figure
hold on
title('RECURSIVE: Labor in region 5')
plot(1:TIME,Ldynamic_sec_cross(1,1:TIME,5))
plot(1:TIME,Ldynamic_dgp(5,1:TIME))
plot(1:TIME,L_belief_agg_recur(5,1:TIME,2),':')
plot(1:TIME,L_belief_agg_dgp(5,1:TIME,2))
plot(1:TIME,L_pf_agg_recur(5,1:TIME,2),'--')
legend('Nonlin PF','DATA','Belief (RECOVERED)','Belief (DGP)','PF (RECOVERED)','location','best')
%% figures
FIGURES(params, eqm_nlpf_HAT_SS, eqm_nlpf_HAT, eqm_dgp, eqm_recur)

