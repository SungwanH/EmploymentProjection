%% Test code for the  second-order approximation (Dynamic part) by comparing with nonlinear and first-order approximation

clear all
close all
clc;
digits(50)

params=PARAMS_TEST(0.5);
v2struct(params.envr);

%% Switchers for this program; loading data
% If the parameter values are changed, check which one to run!
RUN_NLPF_HAT_SS = 1; 
RUN_NLPF_HAT    = 1; 
RUN_NLPF_DD     = 1; 
%RUN_NLPF_DD_BELIEF = 0; 


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
    
    [eqm_nlpf_HAT_SS, approx_nlpf_HAT_SS] = NLPF_HAT(params, starting_point_ss,hat_fundamental_ss,initial_guess_ss);
     save('DATA/NLPF_HAT_SS.mat', 'eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS');     
     clear temporary_struct % eqm_nlpf_HAT_SS approx_nlpf_HAT_SS


params_NLPF=rmfield(params,'prod');
hat_fundamental_cross.T_HAT = ones(J,N,TIME);
for t=1:TIME-1
    hat_fundamental_cross.T_HAT(:,:,t+1)= (params.prod.T(:,:,t+1)./params.prod.T(:,:,t))./ones(J,N,1); 
end
%    hat_fundamental_cross.T_HAT(:,:,2) = params.belief.E_T_hat_test(:,:,2,1)+1;
hat_fundamental_cross.TIME=TIME;
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
%% Nonlinear eqm

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


%% Compute first-order and second-order approx. to the temporary equilibrium

t1=1;
mat_pbp = MAT_CMEX(params, approx_nlpf_HAT_SS);
kappa_hat = zeros(N*J,N,TIME);
V = zeros(R*J,TIME);
L = zeros(R*J,1); %initial deviation from the path around which we are linearizing
E_T_hat_test = params.belief.E_T_hat_test;
[eqm_FO] = PBP_DYN(params, t1, t1, E_T_hat_test, kappa_hat, L, V, approx_nlpf_HAT_SS,mat_pbp);
w_hat_FO = eqm_FO.w;
p_hat_FO = eqm_FO.p;
P_hat_FO = eqm_FO.P;
v_hat_FO = eqm_FO.v;
L_hat_FO = reshape(eqm_FO.L,J,R,TIME);
% generate level values
Ldyn_FO = exp(L_hat_FO)  .* eqm_nlpf_HAT_SS.Ldyn;
wf00_FO = exp(w_hat_FO) .* eqm_nlpf_HAT_SS.wf00;
pf00_FO = exp(p_hat_FO) .* eqm_nlpf_HAT_SS.pf00;

% Second-order approximation
%{
[eqm_SO] = PBP_DYN_SO(params, t1, t1, E_T_hat_test, kappa_hat, L, V, approx_nlpf_HAT_SS);
w_hat_SO = eqm_SO.w;
p_hat_SO = eqm_SO.p;
P_hat_SO = eqm_SO.P;
v_hat_SO = eqm_SO.v;
L_hat_SO = reshape(eqm_SO.L,J,R,TIME);
Ldyn_SO = exp(L_hat_SO)  .* eqm_nlpf_HAT_SS.Ldyn;
wf00_SO = exp(w_hat_SO) .* eqm_nlpf_HAT_SS.wf00;
pf00_SO = exp(p_hat_SO) .* eqm_nlpf_HAT_SS.pf00;

rw_dd = eqm_nlpf_dd.w_lev ./ eqm_nlpf_dd.p_lev;
rw_FO = wf00_FO ./ pf00_FO;
rw_SO = wf00_SO ./ pf00_SO;
p_dd = eqm_nlpf_dd.p_lev;
p_FO = pf00_FO;
p_SO = pf00_SO;
L_hat_NL = zeros(J,N,TIME);
%}


% using constants for second order terms
[eqm_SO_CONS] = PBP_DYN_SO_CONS(params, t1, t1, E_T_hat_test, kappa_hat, L, V, approx_nlpf_HAT_SS, eqm_FO);
w_hat_SO_CONS = eqm_SO_CONS.w;
p_hat_SO_CONS = eqm_SO_CONS.p;
P_hat_SO_CONS = eqm_SO_CONS.P;
v_hat_SO_CONS = eqm_SO_CONS.v;
L_hat_SO_CONS = reshape(eqm_SO_CONS.L,J,R,TIME);
Ldyn_SO_CONS = exp(L_hat_SO_CONS)  .* eqm_nlpf_HAT_SS.Ldyn;
wf00_SO_CONS = exp(w_hat_SO_CONS) .* eqm_nlpf_HAT_SS.wf00;
pf00_SO_CONS = exp(p_hat_SO_CONS) .* eqm_nlpf_HAT_SS.pf00;

for t=1:TIME-1
    L_hat_NL(:,:,t) = log(eqm_nlpf_dd.Ldyn(:,:,t)./eqm_nlpf_dd.Ldyn(:,:,1));
end
for t=1:TIME
%    Labor(:,:,t) =[reshape(eqm_nlpf_dd.Ldyn(:,:,t),N*J,1) reshape(eqm_dgp.Ldyn(:,:,t),N*J,1) reshape(eqm_SO.Ldyn(:,:,t),N*J,1)];
%    Labor(:,:,t) =[reshape(L_hat_NL(:,:,t),N*J,1) reshape((L_hat_FO(:,:,t,1)),N*J,1) reshape((L_hat_SO(:,:,t,1)),N*J,1) reshape((L_hat_SO_CONS(:,:,t,1)),N*J,1)];
    Labor(:,:,t) =[reshape(L_hat_NL(:,:,t),N*J,1) reshape((L_hat_FO(:,:,t,1)),N*J,1) reshape((L_hat_SO(:,:,t,1)),N*J,1)];
end
disp('########################################')

figure
hold on
title('Labor in region 5 sector 1(one-sector)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO_CONS(1,5,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_Labor_region5_onesect.png')
saveas(gcf,'FIGURES/SO_Labor_region5_onesect.fig')

for t=1:TIME
    Labor_err(t) = sum(abs(Labor(:,1,t)-Labor(:,2,t)))/sum(abs(Labor(:,1,t)-Labor(:,3,t))); % difference in error in Labor allocation
    rw(:,:,t) = [reshape(rw_dd(:,:,t),N*J,1) reshape(rw_FO(:,:,t),N*J,1) reshape(rw_SO(:,:,t),N*J,1)];
    p(:,:,t) = [reshape(p_dd(:,:,t),N*J,1) reshape(p_FO(:,:,t),N*J,1) reshape(p_SO(:,:,t),N*J,1)];
    rw_err(t) = sum(abs(rw(:,1,t)-rw(:,2,t)))/sum(abs(rw(:,1,t)-rw(:,3,t))); % difference in error in realwgae
    p_err(t) = sum(abs(p(:,1,t)-p(:,2,t)))/sum(abs(p(:,1,t)-p(:,3,t))); % difference in error in price
end
%price = [p_prime_lev(:) p_FO(:)];
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

disp('Approximation error for Labor')
sum(sum(abs(Labor(:,1,:)-Labor(:,2,:)))) % error of the first order
sum(sum(abs(Labor(:,1,:)-Labor(:,3,:)))) % error of the second order

error_FO=(((Labor(:,1,2)-Labor(:,2,2)))./abs(Labor(:,1,2)))' % pct. diff. error of the first order
error_SO=(((Labor(:,1,2)-Labor(:,3,2)))./abs(Labor(:,1,2)))' % pct. diff. error of the second order
(((Labor(:,1,49)-Labor(:,2,49)))./(Labor(:,1,49)))' % pct. diff. error of the first order
(((Labor(:,1,49)-Labor(:,3,49)))./(Labor(:,1,49)))' % pct. diff. error of the second order
abs(error_SO)./abs(error_FO)
mean(abs(error_SO)./abs(error_FO))
error_FO_T=(((Labor(:,1,TIME-1)-Labor(:,2,TIME-1)))./abs(Labor(:,1,TIME-1)))' % pct. diff. error of the first order
error_SO_T=(((Labor(:,1,TIME-1)-Labor(:,3,TIME-1)))./abs(Labor(:,1,TIME-1)))' % pct. diff. error of the second order
abs(error_SO_T)./abs(error_FO_T)


x = linspace(min(min(rw(:,:,TIME-1))),max(max(rw(:,:,TIME-1))));
y = x;
figure
hold on
plot(x,y,'black')
scatter(rw(:,1,TIME-1),rw(:,2,TIME-1),'red')
scatter(rw(:,1,TIME-1),rw(:,3,TIME-1),'blue','d')
xlabel('Nonlinear hat (realwage)') 
ylabel('Approximation (realwage)') 
legend('45 degree line', 'First Order', 'Second Order', 'location', 'best')
title('Level / FO / SO Comparison (Real Wage, New SS)')
saveas(gcf,'FIGURES/rwcomparison_ext.png')

x = linspace(min(min(Labor(:,:,TIME-1))),max(max(Labor(:,:,TIME-1))));
y = x;
figure
hold on
plot(x,y,'black')
scatter(Labor(:,1,TIME-1),Labor(:,2,TIME-1),'red')
scatter(Labor(:,1,TIME-1),Labor(:,3,TIME-1),'blue','d')
xlabel('Nonlinear hat (Labor)') 
ylabel('Approximation (Labor)') 
legend('45 degree line', 'First Order', 'Second Order', 'location', 'best')
title('Level / FO / SO Comparison (Labor, New SS)')
saveas(gcf,'FIGURES/laborcomparison_ext.png')


figure
hold on
title('Productivity in level (extended)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(params.prod.T(1,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(params.prod.T(1,6,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(params.prod.T(1,10,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(params.prod.T(1,15,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(params.prod.T(1,20,1:TIME-1),[1,3,2]))
legend('Region 1 Sector 1','Region 6 Sector 1','Region 10 Sector 1','Region 15 Sector 1','Region 20 Sector 1','location','best')
saveas(gcf,'FIGURES/Productivity_ext.png')
saveas(gcf,'FIGURES/Productivity_ext.fig')

figure
hold on
title('Labor in region 1 sector 1 (extended)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(1,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(1,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO(1,1,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_Labor_region1_ext.png')
saveas(gcf,'FIGURES/SO_Labor_region1_ext.fig')

figure
hold on
title('Labor in region 1 sector 2 (extended)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(2,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(2,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO(2,1,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_Labor_region1_sec2_ext.png')
saveas(gcf,'FIGURES/SO_Labor_region1_sec2_ext.fig')

figure
hold on
title('Labor in region 10 sector 1 (extended)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(1,10,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(1,10,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO(1,10,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_Labor_region10_ext.png')
saveas(gcf,'FIGURES/SO_Labor_region10_ext.fig')

figure
hold on
title('Labor in region 18 sector 2 (extended)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(2,18,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(2,18,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO(2,18,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_Labor_region18_ext.png')
saveas(gcf,'FIGURES/SO_Labor_region18_ext.fig')

figure
hold on
title('Labor in region 10 sector 2 (extended)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(2,10,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(2,10,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO(2,10,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_Labor_region10_sec2_ext.png')
saveas(gcf,'FIGURES/SO_Labor_region10_sec2_ext.fig')
figure
hold on
title('Labor in region 20 sector 1 (extended)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(1,20,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(1,20,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO(1,20,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_Labor_region20_ext.png')
saveas(gcf,'FIGURES/SO_Labor_region20_ext.fig')

figure
hold on
title('Labor in region 20 sector 2 (extended)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(2,20,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(2,20,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO(2,20,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_Labor_region20_sec2_ext.png')
saveas(gcf,'FIGURES/SO_Labor_region20_sec2_ext.fig')

%{
%% Figures for SO_CONS
x = linspace(min(min(Labor(:,:,2))),max(max(Labor(:,:,2))));
y = x;
figure
hold on
plot(x,y,'black')
scatter(Labor(:,1,TIME-1),Labor(:,2,TIME-1),'red')
scatter(Labor(:,1,TIME-1),Labor(:,3,TIME-1),'blue','d')
scatter(Labor(:,1,TIME-1),Labor(:,4,TIME-1),'green','d')
xlabel('Nonlinear hat (Labor)') 
ylabel('Approximation (Labor)') 
legend('45 degree line', 'First Order', 'Second Order', 'Second Order:Constant', 'location', 'best')
title('Level / FO / SO Comparison (Labor)')
saveas(gcf,'FIGURES/laborcomparison_SO_CONS.png')


figure
hold on
title('Labor in region 1 sector 1 (with SO CONS)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(1,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(1,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO(1,1,1:TIME-1),[1,3,2]),'--')
plot(1:TIME-1,permute(Ldyn_SO_CONS(1,1,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','Second-order:Cons','location','best')
saveas(gcf,'FIGURES/SO_Labor_region1_Cons.png')
saveas(gcf,'FIGURES/SO_Labor_region1_Cons.fig')

figure
hold on
title('Labor in region 1 sector 2 (with SO CONS)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(2,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(2,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO(2,1,1:TIME-1),[1,3,2]),'--')
plot(1:TIME-1,permute(Ldyn_SO_CONS(2,1,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','Second-order:Cons','location','best')
saveas(gcf,'FIGURES/SO_Labor_region1_sec2_Cons.png')
saveas(gcf,'FIGURES/SO_Labor_region1_sec2_Cons.fig')

figure
hold on
title('Labor in region 3 sector 1 (with SO CONS)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(1,3,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(1,3,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO(1,3,1:TIME-1),[1,3,2]),'--')
plot(1:TIME-1,permute(Ldyn_SO_CONS(1,3,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','Second-order:Cons','location','best')
saveas(gcf,'FIGURES/SO_Labor_region3_Cons.png')
saveas(gcf,'FIGURES/SO_Labor_region3_Cons.fig')

figure
hold on
title('Labor in region 3 sector 2 (with SO CONS)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(2,3,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(2,3,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO(2,3,1:TIME-1),[1,3,2]),'--')
plot(1:TIME-1,permute(Ldyn_SO_CONS(2,3,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','Second-order:Cons','location','best')
saveas(gcf,'FIGURES/SO_Labor_region3_Sec2_Cons.png')
saveas(gcf,'FIGURES/SO_Labor_region3_Sec2_Cons.fig')

figure
hold on
title('Labor in region 5 sector 1 (with SO CONS)')
%plot(1:TIME-1,permute(eqm_dgp.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_FO(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_SO(1,5,1:TIME-1),[1,3,2]),'--')
plot(1:TIME-1,permute(Ldyn_SO_CONS(1,5,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','Second-order:Cons','location','best')
saveas(gcf,'FIGURES/SO_Labor_region5_Cons.png')
saveas(gcf,'FIGURES/SO_Labor_region5_Cons.fig')

%}

figure
hold on
title('Real wage in region 1 sector 1   ')
%plot(1:TIME-1,permute(eqm_FO.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dd(1,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_FO(1,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_SO(1,1,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_realwage_region1.png')
saveas(gcf,'FIGURES/SO_realwage_region1.fig')


figure
hold on
title('Real wage in region 1 sector 2   ')
%plot(1:TIME-1,permute(eqm_FO.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dd(2,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_FO(2,1,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_SO(2,1,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_realwage_region1_sec2.png')
saveas(gcf,'FIGURES/SO_realwage_region1_sec2.fig')

figure
hold on
title('Real wage in region 3 sector 1   ')
%plot(1:TIME-1,permute(eqm_FO.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dd(1,3,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_FO(1,3,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_SO(1,3,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_realwage_region3.png')
saveas(gcf,'FIGURES/SO_realwage_region3.fig')

figure
hold on
title('Real wage in region 5 sector 1   ')
%plot(1:TIME-1,permute(eqm_FO.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dd(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_FO(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_SO(1,5,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_realwage_region5.png')
saveas(gcf,'FIGURES/SO_realwage_region5.fig')


figure
hold on
title('Real wage in region 10 sector 1   ')
%plot(1:TIME-1,permute(eqm_FO.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dd(1,10,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_FO(1,10,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_SO(1,10,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_realwage_region10.png')
saveas(gcf,'FIGURES/SO_realwage_region10.fig')

figure
hold on
title('Real wage in region 10 sector 2   ')
%plot(1:TIME-1,permute(eqm_FO.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dd(2,10,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_FO(2,10,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_SO(2,10,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_realwage_region10_sec2.png')
saveas(gcf,'FIGURES/SO_realwage_region10_sec2.fig')

figure
hold on
title('Real wage in region 20 sector 1   ')
%plot(1:TIME-1,permute(eqm_FO.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dd(1,20,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_FO(1,20,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_SO(1,20,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_realwage_region20.png')
saveas(gcf,'FIGURES/SO_realwage_region20.fig')


figure
hold on
title('Real wage in region 20 sector 2   ')
%plot(1:TIME-1,permute(eqm_FO.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_dd(2,20,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_FO(2,20,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(rw_SO(2,20,1:TIME-1),[1,3,2]),'--')
legend('Nonlinear','First-order','Second-order','location','best')
saveas(gcf,'FIGURES/SO_realwage_region20_sec2.png')
saveas(gcf,'FIGURES/SO_realwage_region20_sec2.fig')

