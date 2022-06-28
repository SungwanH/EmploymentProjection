clear all
close all
clc;
digits(50)


params=PARAMS_TEST(0.5);
%params=PARAMS(0.5);
v2struct(params.envr);

%% Switchers for this program; loading data
RUN_NLPF_HAT_SS   = 1; 
RUN_NLPF_HAT      = 1; 
RUN_NLPF_DD       = 0; 
%RUN_NLPF_DD_BELIEF= 0; 
RUN_DGP_SO        = 0; 
RUN_DGP           = 0; 
RUN_DGP_HETERO    = 0; 
RUN_RECUR         = 0;



% Load initial data
%data_4_sector=load('DATA/BASE_FOURSECTOR.mat', 'VALjn00', 'Din00','mu0','L0');

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
    temporary_struct.w_guess   = ones(J,N); %initial guess for wage
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
%    clear  eqm_nlpf_HAT_SS approx_nlpf_HAT_SS   
else
    load('DATA/NLPF_DD.mat','eqm_nlpf_dd','approx_nlpf_dd'); %loading the equilibrium values in the counterfactual economy
end
%{
%% Obtain non-linear Double (Cross & Time) difference (Counterfactual 'Belief')
% This part derives counterfactual impact of productivity shock according
% to the belief
params_NLPF=rmfield(params,'prod');
hat_fundamental_cross_belief.T_HAT = ones(J,N,TIME);
for t=1:TIME-1
%    hat_fundamental_cross_belief.T_HAT(:,:,t+1)=(params.belief.T_belief(:,:,t+1,1)./params.belief.T_belief(:,:,t,1))./ones(J,N,1); 
    hat_fundamental_cross_belief.T_HAT(:,:,t+1)=(exp(params.belief.E_T_hat(:,:,t+1,1))./exp(params.belief.E_T_hat(:,:,t,1)))./ones(J,N,1); 
end
hat_fundamental_cross_belief.TIME=TIME;

if RUN_NLPF_DD_BELIEF ==1
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
%}
%}
%{
    load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS')
w_guess   = ones(J,N); %initial guess for wage
p_guess   = ones(J,N); %initial guess for good prices
kappa_hat = ones(J*N,N); % relative change in trade cost
Ljn_hat          = ones(J,N);
Ljn_hat(1,1)          = 1.05;
Ljn_hat(J,N)          = 0.95;
T_hat = hat_fundamental_cross.T_HAT(:,:,2);
VALjn0 = eqm_nlpf_HAT_SS.VALjn00(:,:,TIME_SS); 
Din = approx_nlpf_HAT_SS.pi(:,:,TIME_SS); 
[wf0, pf0,~,~,~,~] = NLPF_TEMP_HAT(params, VALjn0, Din, kappa_hat, T_hat, Ljn_hat, w_guess, p_guess);

L_hat = zeros(R*J,TIME);
L_hat(1,2) =log(1.05)-log(1);
L_hat(R*J,2) =log(0.95)-log(1);
T_hat = params.belief.E_T_hat_test;
kappa_hat=zeros(J*N,N,TIME);
W = zeros(J,N,TIME);
approx = approx_nlpf_HAT_SS;
[w_hat, p_hat, P_hat, pi_hat, X_hat] = PBP_TEMP(params, 2, T_hat, kappa_hat, W, L_hat, approx);

[w_hat_SO, p_hat_SO, P_hat_SO, pi_hat_SO, X_hat_SO] = PBP_TEMP_SO(params, 2, T_hat, kappa_hat, W, L_hat, approx);
pf0_FO = exp(p_hat(:,:,2));
pf0_SO = exp(p_hat_SO(:,:,2));
wf0_FO = exp(w_hat(:,:,2));
wf0_SO = exp(w_hat_SO(:,:,2));

p = [reshape(pf0,N*J,1) reshape(pf0_FO,N*J,1) reshape(pf0_SO,N*J,1)];
w = [reshape(wf0,N*J,1) reshape(wf0_FO,N*J,1) reshape(wf0_SO,N*J,1)];
rw = w./p;
%[reshape(pf0-wf0,N*J,1) reshape(pf0_FO(:,:,2)-wf0_FO(:,:,2),N*J,1) reshape(pf0_SO(:,:,2)-wf0_SO(:,:,2),N*J,1)];

sum(abs(p(:,1)-p(:,2))) % error of the first order approx.
sum(abs(p(:,1)-p(:,3))) % error of the second order approx.
sum(abs(w(:,1)-w(:,2))) % error of the first order approx.
sum(abs(w(:,1)-w(:,3))) % error of the second order approx.
sum(abs(rw(:,1)-rw(:,2))) % error of the first order approx.
sum(abs(rw(:,1)-rw(:,3))) % error of the second order approx.

x = linspace(min(min(rw)),max(max(rw)));
y = x;
figure
hold on
plot(x,y,'black')
scatter(rw(:,1),rw(:,2),'red')
scatter(rw(:,1),rw(:,3),'blue','d')
xlabel('Nonlinear hat (realwage)') 
ylabel('Approximation (realwage)') 
legend('45 degree line', 'First Order', 'Second Order', 'location', 'best')
title('Level / FO / SO Comparison')

%}

tic
%    mat_pbp = MAT_CMEX(params, approx_nlpf_dd);
    mat_pbp = MAT_CMEX(params, approx_nlpf_HAT_SS);
toc
[eqm_dgp_so] = SO_OUT(params, eqm_nlpf_HAT_SS, approx_nlpf_HAT_SS);
[eqm_dgp] = FO_OUT(params, eqm_nlpf_HAT_SS, approx_nlpf_HAT_SS, mat_pbp);




%% Obtain DGP path using second order
if RUN_DGP_SO ==1
disp('#################')
disp('Running DGP_Second Order')
%tic
    [eqm_dgp_so, approx_dgp_so] = DGP_SO(params, eqm_nlpf_HAT_SS, approx_nlpf_HAT_SS);
%    [eqm_dgp_so, approx_dgp_so] = DGP_SO(params, eqm_nlpf_dd, approx_nlpf_dd);
%toc
else
    load('DATA/DGP_SO.mat', 'eqm_dgp_so','approx_dgp_so'); 
end

%% Generate matrices for temporary equilibrium
tic
    mat_pbp = MAT(params, approx_nlpf_dd);
toc


tic
%    mat_pbp = MAT_CMEX(params, approx_nlpf_dd);
    mat_pbp = MAT_CMEX(params, approx_nlpf_HAT_SS);
toc
%{
tic
    mat_pbp2 = MAT_CMEXB(params, approx_nlpf_dd);
toc
%}
%% Obtain DGP path
if RUN_DGP ==1
disp('#################')
disp('Running DGP')
tic
%     [eqm_dgp, approx_dgp] = DGP(params, eqm_nlpf_dd, approx_nlpf_dd,mat_pbp);
    [eqm_dgp, approx_dgp] = DGP(params, eqm_nlpf_HAT_SS, approx_nlpf_HAT_SS, mat_pbp);
toc
else
    load('DATA/DGP.mat', 'eqm_dgp','approx_dgp'); 
end


rw_dd = eqm_nlpf_dd.w_lev ./ eqm_nlpf_dd.p_lev;
rw_dgp = eqm_dgp.wf00 ./ eqm_dgp.pf00;
rw_dgp_so = eqm_dgp_so.wf00 ./ eqm_dgp_so.pf00;
p_dd = eqm_nlpf_dd.p_lev;
p_dgp = eqm_dgp.pf00;
p_dgp_so = eqm_dgp_so.pf00;
for t=1:TIME
    Labor(:,:,t) =[reshape(eqm_nlpf_dd.Ldyn(:,:,t),N*J,1) reshape(eqm_dgp.Ldyn(:,:,t),N*J,1) reshape(eqm_dgp_so.Ldyn(:,:,t),N*J,1)];
end
disp('########################################')
disp('Approximation error for Labor')
sum(abs(Labor(:,1,2)-Labor(:,2,2))) % error of the first order on Labor allocation
sum(abs(Labor(:,1,2)-Labor(:,3,2))) % error of the second order on Labor allocation 
for t=1:TIME
    Labor_err(t) = sum(abs(Labor(:,1,t)-Labor(:,2,t)))/sum(abs(Labor(:,1,t)-Labor(:,3,t))); % difference in error in Labor allocation
    rw(:,:,t) = [reshape(rw_dd(:,:,t),N*J,1) reshape(rw_dgp(:,:,t),N*J,1) reshape(rw_dgp_so(:,:,t),N*J,1)];
    p(:,:,t) = [reshape(p_dd(:,:,t),N*J,1) reshape(p_dgp(:,:,t),N*J,1) reshape(p_dgp_so(:,:,t),N*J,1)];
    rw_err(t) = sum(abs(rw(:,1,t)-rw(:,2,t)))/sum(abs(rw(:,1,t)-rw(:,3,t))); % difference in error in realwgae
    p_err(t) = sum(abs(p(:,1,t)-p(:,2,t)))/sum(abs(p(:,1,t)-p(:,3,t))); % difference in error in price
end
disp('########################################')
disp('Approximation error for price')
sum(abs(p(:,1,2)-p(:,2,2))) % error of the first order
sum(abs(p(:,1,2)-p(:,3,2))) % error of the second order

(p(:,1,2)-p(:,2,2))' % error of the first order
(p(:,1,2)-p(:,3,2))' % error of the second order
disp('Approximation error for realwage')
sum(abs(rw(:,1,2)-rw(:,2,2))) % error of the first order
sum(abs(rw(:,1,2)-rw(:,3,2))) % error of the second order

(abs(rw(:,1,2)-rw(:,2,2)))' % error of the first order
(abs(rw(:,1,2)-rw(:,3,2)))' % error of the second order

sum(sum(abs(Labor(:,1,:)-Labor(:,2,:)))) % error of the second order
sum(sum(abs(Labor(:,1,:)-Labor(:,3,:)))) % error of the second order
x = linspace(min(min(rw(:,:,2))),max(max(rw(:,:,2))));
y = x;
figure
hold on
plot(x,y,'black')
scatter(rw(:,1,2),rw(:,2,2),'red')
scatter(rw(:,1,2),rw(:,3,2),'blue','d')
xlabel('Nonlinear hat (realwage)') 
ylabel('Approximation (realwage)') 
legend('45 degree line', 'First Order', 'Second Order', 'location', 'best')
title('Level / FO / SO Comparison(Temp only)')
saveas(gcf,'FIGURES/Temp_Dyn_rwcomparison.png')

figure
hold on
title('Labor in region 5 sector 1 ')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp_so.Ldyn(1,5,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_Labor_region5.png')
saveas(gcf,'FIGURES/SO_Labor_region5.fig')

figure
hold on
title('Labor in region 3 sector 1 ')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(1,3,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,3,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp_so.Ldyn(1,3,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_Labor_region3.png')
saveas(gcf,'FIGURES/SO_Labor_region3.fig')
figure
hold on
title('Labor in region 1 sector 1 ')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(1,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp_so.Ldyn(1,1,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_Labor_region1.png')
saveas(gcf,'FIGURES/SO_Labor_region1.fig')


figure
hold on
title('Real wage in region 5 sector 1 ')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dd(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dgp(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dgp_so(1,5,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_realwage_region5.png')
saveas(gcf,'FIGURES/SO_realwage_region5.fig')

figure
hold on
title('Real wage in region 3 sector 1 ')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dd(1,3,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dgp(1,3,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dgp_so(1,3,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_realwage_region3.png')
saveas(gcf,'FIGURES/SO_realwage_region3.fig')

figure
hold on
title('Real wage in region 1 sector 1 ')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dd(1,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dgp(1,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dgp_so(1,1,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_realwage_region1.png')
saveas(gcf,'FIGURES/SO_realwage_region1.fig')

%{
%% Heterogeneous Belief
if RUN_DGP_HETERO == 1
disp('#################')
disp('Running DGP: HETEROGENEOUS AGENT')
tic
    [eqm_dgp_hetero, approx_dgp_hetero] = DGP_HETERO(params, eqm_nlpf_dd, approx_nlpf_dd, mat_pbp);
toc
else
    load('DATA/DGP_HETERO.mat', 'eqm_dgp_hetero','approx_dgp_hetero'); 
end


figure
hold on
title('Realized labor in sector 2 region 5')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp.L_dgp(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp_hetero.L_dgp(1,5,1:TIME-1),[1,3,2]),'.')
plot(1:TIME-1,permute(eqm_dgp_hetero.L_A_dgp(1,5,1:TIME-1,1),[2,3,4,1]),'--')
plot(1:TIME-1,permute(eqm_dgp_hetero.L_B_dgp(1,5,1:TIME-1,1),[2,3,4,1]),'--')
legend('Homogeneous Agent (A only)','HETERO:A+B','HETERO:A','HETERO:B','location','best')
saveas(gcf,'FIGURES/hetero_L_realized.png')

figure
hold on
title('Belief in labor path for sector 1 region 5 at t=1')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp.L_hat(1,5,1:TIME-1,1),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp_hetero.L_A_hat(1,5,1:TIME-1,1),[2,3,4,1]),'--')
plot(1:TIME-1,permute(eqm_dgp_hetero.L_B_hat(1,5,1:TIME-1,1),[2,3,4,1]),':')
legend('Homogeneous Agent (A only)','HETERO:A','HETERO:B','location','best')
saveas(gcf,'FIGURES/hetero_L_belief_1.png')
figure
hold on
title('Belief in labor path for sector 1 region 5')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp.L_hat(1,5,1:TIME-1,1),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp_hetero.L_A_hat(1,5,1:TIME-1,1),[2,3,4,1]),'--')
plot(1:TIME-1,permute(eqm_dgp_hetero.L_B_hat(1,5,1:TIME-1,1),[2,3,4,1]),':')
plot(1:TIME-1,permute(eqm_dgp.L_hat(1,5,1:TIME-1,5),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp_hetero.L_A_hat(1,5,1:TIME-1,5),[2,3,4,1]),'--')
plot(1:TIME-1,permute(eqm_dgp_hetero.L_B_hat(1,5,1:TIME-1,5),[2,3,4,1]),':')
plot(1:TIME-1,permute(eqm_dgp.L_hat(1,5,1:TIME-1,10),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp_hetero.L_A_hat(1,5,1:TIME-1,10),[2,3,4,1]),'--')
plot(1:TIME-1,permute(eqm_dgp_hetero.L_B_hat(1,5,1:TIME-1,10),[2,3,4,1]),':')
plot(1:TIME-1,permute(eqm_dgp.L_hat(1,5,1:TIME-1,20),[1,3,2]))
plot(1:TIME-1,permute(eqm_dgp_hetero.L_A_hat(1,5,1:TIME-1,20),[2,3,4,1]),'--')
plot(1:TIME-1,permute(eqm_dgp_hetero.L_B_hat(1,5,1:TIME-1,20),[2,3,4,1]),':')
legend('Homogeneous Agent (A only)','HETERO:A','HETERO:B','location','best')
saveas(gcf,'FIGURES/hetero_L_belief_2.png')
%}
%{
E_A_T_hat = params.belief.E_A_T_hat;
E_B_T_hat = params.belief.E_B_T_hat;
kappa_A_hat = zeros(N*J,N,TIME);
kappa_B_hat = zeros(N*J,N,TIME);
V = zeros(R*J,TIME);

L_A = zeros(N*J,1);
L_B = zeros(N*J,1);
L_A_init = zeros(N*J,1);
L_B_init = zeros(N*J,1);
L_A_2 = zeros(N*J,1);
L_B_2 = zeros(N*J,1);   
MAXIT       = 1E+8;
TOLFP       = 1E-7;
for t1=1:ENDT+1
%L_B_2 = zeros(N*J,1);
ITER_FP = 0;
maxerror=1;
    while (ITER_FP <= MAXIT) && (maxerror > TOLFP)
    %for t1=1:ENDT+1
    %L_B_2 = L_B(:,t1+1);
        [eqm_A] = PBP_DYN_HETERO(params, t1, t1, E_A_T_hat, kappa_A_hat, L_A_init, L_B_2, approx_nlpf_dd, mat_pbp, 1);
        L_A_new = eqm_A.L_own;
        L_A_2_new = eqm_A.L_own(:,t1+1);
        [eqm_B] = PBP_DYN_HETERO(params, t1, t1, E_B_T_hat, kappa_B_hat, L_B_init, L_A_2_new, approx_nlpf_dd, mat_pbp, 2);
        L_B_new = eqm_B.L_own;
        L_B_2_new = eqm_B.L_own(:,t1+1);
        
%        for t=t1:TIME
        maxerror=max(abs(L_B_2_new-L_B_2))
%        end
        if ITER_FP >3000
            maxerror(t1:TIME)
            t1
            disp('Fixed Point loop err')
            ITER_FP
            stop
        end
        
%        L_A_2 = L_A_2_new;
        L_B_2 = L_B_2_new;

        ITER_FP=ITER_FP+1
    end
    L_A_init = L_A_2_new;
    L_B_init = L_B_2_new;
    p_A_hat(:,:,t1) = eqm_A.p(:,:,t1);
    p_B_hat(:,:,t1) = eqm_B.p(:,:,t1);
    L_A_hat(:,t1+1)=L_A_2_new;
    L_B_hat(:,t1+1)=L_B_2_new;
t1
end
L_A_hat
L_B_hat
p_A_dgp
p_B_dgp
%}