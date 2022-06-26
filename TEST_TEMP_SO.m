%% Test code for the  second-order approximation (temporary part) by comparing with nonlinear and first-order approximation

clear all
close all
clc;
digits(50)


params=PARAMS_TEST(0.5);
v2struct(params.envr);

%% Switchers for this program; loading data
RUN_NLPF_HAT_SS = 1; 


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
else
    load('DATA/NLPF_HAT_SS.mat','eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS'); 
end

%% Comptute non-linear temporary outcomes
w_guess   = ones(J,N); %initial guess for wage
p_guess   = ones(J,N); %initial guess for good prices
kappa_hat = ones(J*N,N); % relative change in trade cost
Ljn_hat          = ones(J,N);
%Ljn_hat(1,1)          = 1.05;
%Ljn_hat(J,N)          = 0.95;
t1 = 2; % The initial period of a productivity shock
T_hat = params.prod.T_HAT(:,:,t1); %This is double diff T_HAT
VALjn0 = eqm_nlpf_HAT_SS.VALjn00(:,:,TIME_SS); %value added
Din = approx_nlpf_HAT_SS.pi(:,:,TIME_SS);  % initial trade share
[wf0, pf0,~,~,~,~] = NLPF_TEMP_HAT(params, VALjn0, Din, kappa_hat, T_hat, Ljn_hat, w_guess, p_guess);

%% Compute first-order and second-order approx. to the temporary equilibrium
L_hat = zeros(R*J,TIME);
%L_hat(1,2) =log(1.05)-log(1);
%L_hat(R*J,2) =log(0.95)-log(1);
T_hat = params.belief.E_T_hat_test; % log difference btw initial steady state
kappa_hat=zeros(J*N,N,TIME);
W = zeros(J,N,TIME); %initial guess
approx = approx_nlpf_HAT_SS;
[w_hat, p_hat, P_hat, pi_hat, X_hat] = PBP_TEMP(params, t1, T_hat, kappa_hat, W, L_hat, approx); %log-linear

[w_hat_SO, p_hat_SO, P_hat_SO, pi_hat_SO, X_hat_SO] = PBP_TEMP_SO(params, t1, T_hat, kappa_hat, W, L_hat, approx); %second-order approx

%% Comparing the approximation errors
pf0_FO = exp(p_hat(:,:,t1)); % price from first-order approx.
pf0_SO = exp(p_hat_SO(:,:,t1));% price from second-order approx.
wf0_FO = exp(w_hat(:,:,t1)); % wage from first-order approx.
wf0_SO = exp(w_hat_SO(:,:,t1)); % wage from second-order approx.

p = [reshape(pf0,N*J,1) reshape(pf0_FO,N*J,1) reshape(pf0_SO,N*J,1)];
w = [reshape(wf0,N*J,1) reshape(wf0_FO,N*J,1) reshape(wf0_SO,N*J,1)];
rw = w/p;

disp('Approximation Errors')
disp('approx. error on price')
disp(sum(abs(p(:,1)-p(:,2)))) % approx. error on price (first order approx.)
disp(sum(abs(p(:,1)-p(:,3)))) % approx. error on price (second order approx.)
disp('approx. error on wage')
disp(sum(abs(w(:,1)-w(:,2)))) % approx. error on wage (first order approx.)
disp(sum(abs(w(:,1)-w(:,3)))) % approx. error on wage (second order approx.)
disp('approx. error on realwage')
disp(sum(abs(rw(:,1)-rw(:,2)))) % approx. error on realwage (first order approx.)
disp(sum(abs(rw(:,1)-rw(:,3)))) % approx. error on realwage (second order approx.)

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
