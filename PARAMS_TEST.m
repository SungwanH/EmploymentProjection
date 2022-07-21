function [params] = PARAMS_TEST(W_TRUE)
TIME    =50; 
TIME_SS = 50;
N       = 87; 
R       = 50; %number of US regions
C       = 37; %number of non-US countries (To be used later)
J       = 4; %number of sectors
US      = 1;
CHINA   = 5; %region number for now assume the shock is to California

% Elasticities
load('DATA/BASE_FOURSECTOR.mat', 'alphas','T','gamma')
ALPHAS = alphas; %consumption share
THETA  = 1./T; 
GAMMA      = gamma; % share of labor in production
clear alphas
clear T

NU      = 5.3436; %dispersion of taste shocks (In CDP: 5.3436)
BETA    = 0.9227; %discount rate (0.99^8: transformed from quarterly to bi-annual)


% Baseline:  
SIGMA       = 0.0; % s.d. of catch-up
ENDT        = 15; %Last period of learning (PF starts from here)
ENDT_SAMPLE = 1;
ENDT_DGP    = 1;
EPST        = 15; %Length of epsilon draw


%% FAKE DATA
%N       = 20; 
%R       = 20; %number of US regions
N       = 5; 
R       = 5; %number of US regions
C       = 0; %number of non-US countries (To be used later)
J       = 1; %number of sectors
US      = 2;
CHINA   = 5; %region number for now assume the shock is to California

rng(20220626)
ALPHAS=1/J*ones(J,N);
%ALPHAS(1,:)=0.6;
%ALPHAS(2,:)=0.4;
%ALPHAS(1,1:3)=0.4;
%ALPHAS(2,1:3)=0.6;
%ALPHAS(1,4:5)=0.7;
%ALPHAS(2,4:5)=0.3;
THETA=7*ones(J,1);
GAMMA=0.5*ones(J,N);
%GAMMA=max(rand(J,N),0.3);%0.5*ones(J,N);
%GAMMA(1,1:3)=0.7;
%GAMMA(2,1:3)=0.4;
%GAMMA(2,4:N)=0.6;

BETA = 0.8;
NU = 1.3;
%BETA =0.85;



%% Technical parameters
ESTM_BOTH   = 0; % if ESTM_BOTH =0, estimate MU only. if =1, estimate both MU and RHO
UPDT_V      = 0.4; %update speed for value loop (lower value->conservative)
UPDT_W      = 0.2; %update speed for wage loop (lower value->conservative)
UPDT_V_NL   = 0.2; %update speed for nonlinear value loop (lower value->conservative)
UPDT_W_NL   = 0.1; %update speed for nonlinear wage loop (lower value->conservative)
TOL_NL      = 1E-8;  %tolerance rate for nonlinear dynamic equilibrium (outerloop)
TOL_NL_TEMP = 1E-8;  %tolerance rate for nonlinear temporary equilibrium (inner loop)
TOLDYN      = 1E-8;  %tolerance rate for linear dynamic equilibrium
TOLTEMP     = 1E-8;  % tolerance rate for linear temporary equilibrium
TOLFP       = 1E-7;  % tolerance rate for fixed point iteration
MAXIT       = 1E+8; %maximum number of iterations


NUM_SIM =5; % number of simulations
DRAW_SD = 0.3; % standard error of draws

params.envr = v2struct(TIME, TIME_SS, N, R, C, J, US, CHINA, ENDT, ENDT_SAMPLE, ENDT_DGP, EPST, NUM_SIM, DRAW_SD); %remoevd TAU from this structure
params.modl = v2struct(THETA, NU, BETA, ALPHAS, THETA, GAMMA);
params.tech = v2struct(ESTM_BOTH, UPDT_V, UPDT_W, UPDT_V_NL, UPDT_W_NL, TOL_NL, TOL_NL_TEMP, TOLDYN, TOLTEMP, TOLFP, MAXIT);


%% Productivity
% Baseline productivity (Used in deriving initial steady state)
T_BASE = ones(J,N,TIME*3)*2; % US (& other countries except China) productivity is constant for all period
T_BASE(:,CHINA,1:TIME*3) = 1.3;

% Learning
RHO     = 0.7; % catch-up speed (Higher RHO -> slower catch-up)
MU      = 0.0;

%Previous productivity of China (used in estimating RHO and MU)
t_bef=10; % # of observations from the previous period
T_PREV=NaN(t_bef,1);
T_PREV(1)=0.6;
MU_PREV=0.45;

for tt=2:t_bef
    T_PREV(tt)=exp(log(T_BASE(1,US,tt))-(1-RHO)*MU_PREV-RHO*(log(T_BASE(1,US,tt))-log(T_PREV(tt-1))));
end

params.prod = v2struct(T_BASE, T_PREV, MU, SIGMA, RHO);

T = PRODUCTIVITY_DGP(params); %objective productivity from the first period and on (not used in testing second order)
%T(2,5,1:TIME) = T_BASE(2,5,1);
params.prod = v2struct(T_BASE, T_PREV, MU, SIGMA, RHO, T);
%The following productivity path is used in the testcode for second-order:
%{
T(1:1,1,1:20) = repmat(reshape(linspace(T_BASE(1,1,1),T_BASE(1,1,20).*1.3,20),1,1,20),J,1,1);
T(1:J,1,21:TIME) = T_BASE(1,1,20)*1.3;
T(1:J,2,1:20) = repmat(reshape(linspace(T_BASE(1,2,1),T_BASE(1,2,20).*1,20),1,1,20),J,1,1);
T(1:J,2,21:TIME) = T_BASE(1,2,20)*1;
T(1:J,3,1:20) = repmat(reshape(linspace(T_BASE(1,3,1),T_BASE(1,3,20).*1,20),1,1,20),J,1,1);
T(1:J,3,21:TIME) = T_BASE(1,3,20)*1;
T(1:J,4,1:20) = repmat(reshape(linspace(T_BASE(1,4,1),T_BASE(1,4,20).*1,20),1,1,20),J,1,1);
T(1:J,4,21:TIME) = T_BASE(1,4,20)*1;
T(1:J,5,1:20) = repmat(reshape(linspace(T_BASE(1,5,1),T_BASE(1,5,20).*0.7,20),1,1,20),J,1,1);
T(1:J,5,21:TIME) = T_BASE(1,5,20)*0.7;
%}
%{
prod_draw(1,:) = linspace(0.5,1.5,N)';% exp(rand(J,N))-0.7;
prod_draw(2:J,:) = ones(1,N);% exp(rand(J,N))-0.7;
for n=1:N
    for j=1:J
    T(j,n,1:ENDT) = reshape(linspace(T_BASE(j,n,1),T_BASE(j,n,ENDT).*prod_draw(j,n),ENDT),1,1,ENDT);
    T(j,n,ENDT:TIME) = T_BASE(1,n,ENDT)*prod_draw(j,n);
    end
end
%}
% one time shock
%{
T = T_BASE(:,:,1:TIME);
T(1:J,1,2) = T_BASE(1,1,1).*1.3;
T(1:J,2,2) = T_BASE(1,2,1).*1.0;
T(1:J,3,2) = T_BASE(1,3,1).*1.0;
T(1:J,4,2) = T_BASE(1,4,1).*1.0;
T(1:J,5,2) = T_BASE(1,5,1).*1.0;
%}
T_SIM = zeros(J,N,TIME,NUM_SIM);

for s=1:NUM_SIM
    T_SIM(:,:,:,s) = T;
end
for s=1:NUM_SIM
%    T_SIM(1,5,2:TIME,:,s) = DRAW_SD.*randn(1,1,TIME-1,ENDT) + T(1,5,2:TIME);
    T_SIM(1,5,2:ENDT,s) = DRAW_SD.*randn(1,1,ENDT-1) + T(1,5,2:ENDT);
end
for s=1:NUM_SIM
    params_sim = params;
    params_sim.prod.T = T_SIM(:,:,:,s);
    T_belief(:,:,:,:,s) = BELIEF(params_sim, 1);
end

%T_belief = BELIEF(params, W_TRUE);
E_T_hat_sim = zeros(J,N,TIME,ENDT+1,NUM_SIM); %Deviation from the pefect foresight
for s=1:NUM_SIM
    for t1=1:TIME
        for t2=1:ENDT+1
%            E_T_hat_sim(:,:,t1,t2,s) = log(T(:,:,t1))-log(T_SIM(:,:,t1,t2,s));
            E_T_hat_sim(:,:,t1,t2,s) = log(T_SIM(:,:,t1,s))-log(T_belief(:,:,t1,t2,s));
        end
    end
end

%for s=1:NUM_SIM
%    for t=2:ENDT
%        T_SIM(1,1,t,t,s) = T_SIM(1,1,t,t-1,s);
%        T_SIM(1,5,t,s) = T(1,5,t);
%    end
%end
figure
hold on
plot(1:TIME,permute(T_SIM(1,5,1:TIME,:),[4,3,2,1]),'--')
%plot(1:TIME,permute(T_SIM(1,5,1:TIME,NUM_SIM),[4,3,2,1]),'--')
plot(1:TIME,permute(T(1,5,1:TIME),[1,3,2]))

figure
hold on
%plot(1:TIME,permute(E_T_hat_sim(1,5,1:TIME,1,1),[4,3,2,1]),'--')
%plot(1:TIME,permute(E_T_hat_sim(1,5,1:TIME,2,1),[4,3,2,1]),'--')
%plot(1:TIME,permute(E_T_hat_sim(1,5,1:TIME,3,1),[4,3,2,1]),'--')
plot(1:TIME,permute(T_belief(1,5,1:TIME,1,1),[4,3,2,1]),'--')
plot(1:TIME,permute(T_belief(1,5,1:TIME,2,1),[4,3,2,1]),'--')
plot(1:TIME,permute(T_belief(1,5,1:TIME,3,1),[4,3,2,1]),'--')
plot(1:TIME,permute(T_SIM(1,5,1:TIME,1),[1,3,2]))

% Derive Time difference productivity
T_HAT_SS = ones(J,N,TIME_SS);
T_HAT    = rand(J,N,TIME);
for t=1:TIME_SS-1
    T_HAT_SS(:,:,t+1) = T_BASE(:,:,t+1)./T_BASE(:,:,t); %relative change in productivity (US: 2 for all period, CHINA: 1 for all period except the first period)
end
for t=1:TIME-1
    T_HAT(:,:,t+1) = T(:,:,t+1)./T(:,:,t); %relative change in productivity (CHINA is catching up here)
end
params.prod = v2struct(MU, SIGMA, RHO, T_BASE, T_PREV, T, T_HAT, T_SIM, T_HAT_SS,t_bef);

figure
hold on
plot(1:TIME,permute(E_T_hat_sim(1,5,1:TIME,1,1),[1,3,2,4,5]))
plot(1:TIME,permute(E_T_hat_sim(1,5,1:TIME,1,2),[1,3,2,4,5]))
plot(1:TIME,permute(E_T_hat_sim(1,5,1:TIME,1,3),[1,3,2,4,5]))

figure
hold on
plot(1:TIME,permute(E_T_hat_sim(1,5,1:TIME,1,1),[1,3,2,4,5]))
plot(2:TIME,permute(E_T_hat_sim(1,5,2:TIME,2,1),[1,3,2,4,5]))
plot(3:TIME,permute(E_T_hat_sim(1,5,3:TIME,3,1),[1,3,2,4,5]))


% This is used in testing second order approximation error
E_T_hat_test = zeros(J,N,TIME,ENDT+1); %Deviation from the baseline level
for t1=1:TIME
    for t2=1:ENDT+1
        E_T_hat_test(:,:,t1,t2) = log(T(:,:,t1))-log(T_BASE(:,:,t1));
    end
end


%% belief for productivity given the weight
%T_belief = BELIEF(params, W_TRUE);
E_T_hat  = zeros(J,N,TIME,ENDT+1); 
E_T_hat_pf = zeros(J,N,TIME,ENDT+1);
for t1=2:TIME
    for t2=1:5
        %E_T_hat(:,:,t1,t2) = rand(J,N);
        E_T_hat(:,1:2,t1,t2) = 0.3;
        E_T_hat(:,3:N,t1,t2) = -0.3;
    end
%    for t2=6:20
%        E_T_hat(:,1:2,t1,t2) = 0.1;
%        E_T_hat(:,3:N,t1,t2) = -0.1;
%    end
end

% For heterogeneous agent case:
E_A_T_hat  = E_T_hat;
E_B_T_hat  = E_T_hat;
for t1=2:TIME
    for t2=1:5
        E_B_T_hat(:,1:2,t1,t2) = -0.3;
        E_B_T_hat(:,3:N,t1,t2) = 0.3;
    end
    for t2=6:20
        E_B_T_hat(:,1:2,t1,t2) = -0.1;
        E_B_T_hat(:,3:N,t1,t2) = 0.1;
    end
end

share_A = 0.5*ones(R*J,1); % share of type A
params.hetero=v2struct(E_A_T_hat,E_B_T_hat,share_A);

%or tt=1:ENDT+1
%       E_T_hat(:,CHINA,:,tt) = log(T_belief(:,CHINA,:,tt)) - log(T(:,CHINA,:));
%       E_T_hat_pf(:,CHINA,:,tt) = -log(T_belief(:,CHINA,:,tt)) + log(T(:,CHINA,:));        
%end
params.belief=v2struct(W_TRUE,E_T_hat,E_T_hat_sim,E_T_hat_test,E_A_T_hat,E_B_T_hat,E_T_hat_pf);

end
