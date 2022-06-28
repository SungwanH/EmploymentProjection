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
RUN_NLPF_DD     = 0; 
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

mat_pbp = MAT_CMEX(params, approx_nlpf_HAT_SS);
[eqm_dgp] = FO_OUT(params, eqm_nlpf_HAT_SS, approx_nlpf_HAT_SS, mat_pbp);  % First-order approx.
[eqm_dgp_so] = SO_OUT(params, eqm_nlpf_HAT_SS, approx_nlpf_HAT_SS);        % Second-order approx.


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

for t=1:TIME
    Labor_err(t) = sum(abs(Labor(:,1,t)-Labor(:,2,t)))/sum(abs(Labor(:,1,t)-Labor(:,3,t))); % difference in error in Labor allocation
    rw(:,:,t) = [reshape(rw_dd(:,:,t),N*J,1) reshape(rw_dgp(:,:,t),N*J,1) reshape(rw_dgp_so(:,:,t),N*J,1)];
    p(:,:,t) = [reshape(p_dd(:,:,t),N*J,1) reshape(p_dgp(:,:,t),N*J,1) reshape(p_dgp_so(:,:,t),N*J,1)];
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
title('Level / FO / SO Comparison')
saveas(gcf,'FIGURES/rwcomparison.png')

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