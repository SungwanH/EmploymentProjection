clear 
close all
clc;
digits(50)

% Prepare params under different assumption on the true belief updating
% rule
params=PARAMS(0.5);
params_w_01=PARAMS(0.1);
params_w_09=PARAMS(0.9);
params_w_re=PARAMS(1); % rational expectation
v2struct(params.envr);

%% Inspecting productivity and belief.
Figure1=figure;
hold on
title('Actual and belief productivity in region 1')
h(1)=plot(1:TIME, permute(params.prod.T(1,CHINA,1:TIME),[2,3,1]),'LineWidth',3);
h(2)=plot(1:TIME, permute(params.belief.T_belief(1,CHINA,1:TIME,1),[2,3,4,1]),'--','LineWidth',2);
h(2)=plot(15:TIME, permute(params.belief.T_belief(1,CHINA,15:TIME,15),[2,3,4,1]),'--','LineWidth',2);
h(2)=plot(30:TIME, permute(params.belief.T_belief(1,CHINA,30:TIME,30),[2,3,4,1]),'--','LineWidth',2);
h(3)=plot(1:TIME, permute(params_w_re.belief.T_belief(1,CHINA,1:TIME,1),[2,3,4,1]),'.-','LineWidth',1.2);
h(3)=plot(15:TIME, permute(params_w_re.belief.T_belief(1,CHINA,15:TIME,15),[2,3,4,1]),'.-','LineWidth',1.2);
h(3)=plot(30:TIME, permute(params_w_re.belief.T_belief(1,CHINA,30:TIME,30),[2,3,4,1]),'.-','LineWidth',1.2);
h(4)=plot(1:TIME, permute(params_w_01.belief.T_belief(1,CHINA,1:TIME,1),[2,3,4,1]),':','LineWidth',2);
h(4)=plot(15:TIME, permute(params_w_01.belief.T_belief(1,CHINA,15:TIME,15),[2,3,4,1]),':','LineWidth',2);
h(4)=plot(30:TIME, permute(params_w_01.belief.T_belief(1,CHINA,30:TIME,30),[2,3,4,1]),':','LineWidth',2);
h(5)=xline(30);
legend(h([1,2 ,3,4,5]),{'Actual','Belief (W=0.5)', 'Belief (W=1)','Belief(W=0.1)','Shocks stop'},'location','southeast');
print(Figure1,'figures/productivity_and_belief.png','-dpng','-r600');



%% Switchers for this program; loading data
RUN_NLPF_HAT_SS = 1; 
RUN_NLPF_HAT    = 1; 
RUN_NLPF_DD     = 1; 
RUN_DGP         = 1; 
RUN_RECUR       = 1;

% Load initial data
data_4_sector=load('DATA/BASE_FOURSECTOR.mat', 'VALjn00', 'Din00','mu0','L0');

%%Using fake data
data_4_sector = FAKE_DATA(params);

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
        temporary_struct.w_guess, temporary_struct.p_guess);
     
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



%% Obtain non-linear Double (Cross & Time) difference
% This part derives counterfactual impact of productivity shock
params_NLPF=rmfield(params,'prod');
hat_fundamental_cross.T_HAT = ones(J,N,TIME);
for t=1:TIME-1
    hat_fundamental_cross.T_HAT(:,:,t+1)= (params.prod.T(:,:,t+1)./params.prod.T(:,:,t))./ones(J,N,1); 
end
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
    base_point_dd.mu_base_dot(isnan(base_point_dd.mu_base_dot)) = 0;
    base_point_dd.Din_base_dot(isnan(base_point_dd.Din_base_dot)) = 0;
    base_point_dd.VALjn_base    = eqm_nlpf_HAT.VALjn00;
    base_point_dd.Ldyn_base     = eqm_nlpf_HAT.Ldyn;
    base_point_dd.Din_base      = approx_nlpf_HAT.pi;
    base_point_dd.mu_base       = approx_nlpf_HAT.mu;
    base_point_dd.wf00_base     = eqm_nlpf_HAT.wf00;
    base_point_dd.pf00_base     = eqm_nlpf_HAT.pf00;

    initial_guess_cross.v_cd       = ones(R*(J),TIME); %Initial guess for the Ys exp((V1'-V0')-(V1-V0))^1/NU 

    [eqm_nlpf_dd, approx_nlpf_dd] = NLPF_DD_NEW(params_NLPF, starting_point_dd, base_point_dd, hat_fundamental_cross, initial_guess_cross);

    save('DATA/NLPF_DD.mat', 'eqm_nlpf_dd','approx_nlpf_dd');     
    clear  eqm_nlpf_HAT_SS approx_nlpf_HAT_SS   
else
    load('DATA/NLPF_DD.mat','eqm_nlpf_dd','approx_nlpf_dd'); %loading the equilibrium values in the counterfactual economy
end

%% Obtain non-linear Double (Cross & Time) difference (Counterfactual 'Belief')
% This part derives counterfactual impact of productivity shock according
% to the belief
params_NLPF=rmfield(params,'prod');
hat_fundamental_cross_belief.T_HAT = ones(J,N,TIME);
for t=1:TIME-1
    hat_fundamental_cross_belief.T_HAT(:,:,t+1)=(params.belief.T_belief(:,:,t+1,1)./params.belief.T_belief(:,:,t,1))./ones(J,N,1); 
end
hat_fundamental_cross_belief.TIME=TIME;

if RUN_NLPF_DD ==1
    disp('#################')
    disp('Running NLPF_CROSS (With Belief)')
    
    load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS')

    % starting points for counterfactual path
    starting_point_dd_belief.VALjn0    = eqm_nlpf_HAT_SS.VALjn00(:,:,TIME_SS); %Labor compensation (w*L)
    starting_point_dd_belief.X0        = eqm_nlpf_HAT_SS.X(:,:,TIME_SS); %Total expenditure     
    starting_point_dd_belief.L0        = eqm_nlpf_HAT_SS.Ldyn(:,:,TIME_SS);
    starting_point_dd_belief.mu0       = approx_nlpf_HAT_SS.mu(:,:,TIME_SS);
    starting_point_dd_belief.Din0      = approx_nlpf_HAT_SS.pi(:,:,TIME_SS); %bilateral trade shares
    
    %compute dot values from the baseline path level variables
    base_point_dd_belief.mu_base_dot(:,:,1)    = approx_nlpf_HAT.mu(:,:,1)./approx_nlpf_HAT_SS.mu(:,:,TIME_SS);
    for t=1:TIME-1
        base_point_dd_belief.VALjn_base_dot(:,:,t+1) = eqm_nlpf_HAT.VALjn00(:,:,t+1)./eqm_nlpf_HAT.VALjn00(:,:,t);
        base_point_dd_belief.Ldyn_base_dot(:,:,t+1)  = eqm_nlpf_HAT.Ldyn(:,:,t+1)./eqm_nlpf_HAT.Ldyn(:,:,t);
        base_point_dd_belief.Din_base_dot(:,:,t+1)   = approx_nlpf_HAT.pi(:,:,t+1)./approx_nlpf_HAT.pi(:,:,t);
        base_point_dd_belief.mu_base_dot(:,:,t+1)    = approx_nlpf_HAT.mu(:,:,t+1)./approx_nlpf_HAT.mu(:,:,t);
    end
    base_point_dd_belief.mu_base_dot(isnan(base_point_dd_belief.mu_base_dot)) = 0;
    base_point_dd_belief.Din_base_dot(isnan(base_point_dd_belief.Din_base_dot)) = 0;
    base_point_dd_belief.VALjn_base    = eqm_nlpf_HAT.VALjn00;
    base_point_dd_belief.Ldyn_base     = eqm_nlpf_HAT.Ldyn;
    base_point_dd_belief.Din_base      = approx_nlpf_HAT.pi;
    base_point_dd_belief.mu_base       = approx_nlpf_HAT.mu;
    base_point_dd_belief.wf00_base     = eqm_nlpf_HAT.wf00;
    base_point_dd_belief.pf00_base     = eqm_nlpf_HAT.pf00;

    initial_guess_dd_belief.v_cd       = ones(R*(J),TIME); %Initial guess for the Ys exp((V1'-V0')-(V1-V0))^1/NU 
    
    [eqm_nlpf_dd_belief, approx_nlpf_dd_belief] = NLPF_DD_NEW(params_NLPF, starting_point_dd_belief, base_point_dd_belief, hat_fundamental_cross_belief, initial_guess_dd_belief);

    save('DATA/NLPF_DD_BELIEF.mat', 'eqm_nlpf_dd_belief','approx_nlpf_dd_belief');     
%    clear  eqm_nlpf_HAT_SS approx_nlpf_HAT_SS   
else
    load('DATA/NLPF_DD_BELIEF.mat','eqm_nlpf_dd_belief','approx_nlpf_dd_belief'); %loading the equilibrium values in the counterfactual economy
end



%% Obtain DGP path
if RUN_DGP ==1
disp('#################')
disp('Running DGP')
    [eqm_dgp, approx_dgp] = DGP(params, eqm_nlpf_dd, approx_nlpf_dd);
    [eqm_dgp_w_01, approx_dgp_w_01] = DGP(params_w_01, eqm_nlpf_dd, approx_nlpf_dd);
    [eqm_dgp_w_09, approx_dgp_w_09] = DGP(params_w_09, eqm_nlpf_dd, approx_nlpf_dd);        
    [eqm_dgp_w_re, approx_dgp_w_re] = DGP(params_w_re, eqm_nlpf_dd, approx_nlpf_dd);            
else
    load('DATA/DGP.mat', 'eqm_dgp','approx_dgp'); 
end




%% Obtain Period by period DGP & PF deviation
if RUN_RECUR ==1
disp('#################')
disp('Running RECURSIVE')
    initial_recur.v_hat = eqm_dgp.v_hat;
    initial_recur.w_hat = eqm_dgp.w_hat;
    [eqm_recur_w_01] = RECURSIVE(params, 0.1, initial_recur, eqm_dgp, approx_dgp);
    [eqm_recur] = RECURSIVE(params, 0.5, initial_recur, eqm_dgp, approx_dgp);
    [eqm_recur_w_09] = RECURSIVE(params, 0.9, initial_recur, eqm_dgp, approx_dgp);    
    [eqm_recur_w_re] = RECURSIVE(params, 1, initial_recur, eqm_dgp, approx_dgp);        
else
    load('DATA/RECURSIVE.mat', 'eqm_recur');     
end


%% Make figures
FIGURES_SCRIPT


%% To do the saving function is not done yet.
