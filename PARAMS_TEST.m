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
SIGMA       = 0.01; % s.d. of catch-up
ENDT        = 1; %Last period of learning (PF starts from here)
ENDT_SAMPLE = 1;
ENDT_DGP    = 1;
EPST        = 30; %Length of epsilon draw


%% FAKE DATA
N       = 5; 
R       = 5; %number of US regions
C       = 0; %number of non-US countries (To be used later)
J       = 2; %number of sectors
US      = 2;
CHINA   = 5; %region number for now assume the shock is to California

rng(20220625)
ALPHAS=1/J*ones(J,N);
%ALPHAS(1,:)=0.6;
%ALPHAS(2,:)=0.4;
%ALPHAS(1,1:3)=0.4;
%ALPHAS(2,1:3)=0.6;
%ALPHAS(1,4:5)=0.7;
%ALPHAS(2,4:5)=0.3;
THETA=11*ones(J,1);
%GAMMA=0.5*ones(J,N);
GAMMA=max(rand(J,N),0.3);%0.5*ones(J,N);
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


params.envr = v2struct(TIME, TIME_SS, N, R, C, J, US, CHINA, ENDT, ENDT_SAMPLE, ENDT_DGP, EPST); %remoevd TAU from this structure
params.modl = v2struct(THETA, NU, BETA, ALPHAS, THETA, GAMMA);
params.tech = v2struct(ESTM_BOTH, UPDT_V, UPDT_W, UPDT_V_NL, UPDT_W_NL, TOL_NL, TOL_NL_TEMP, TOLDYN, TOLTEMP, TOLFP, MAXIT);


%% Productivity
% Baseline productivity (Used in deriving initial steady state)
T_BASE = ones(J,N,TIME*3)*2; % US (& other countries except China) productivity is constant for all period
T_BASE(:,CHINA,1:TIME*3) = 1.3;

% Learning
RHO     = 0.9; % catch-up speed (Higher RHO -> slower catch-up)
MU      = 0.3;

%Previous productivity of China (used in estimating RHO and MU)
t_bef=10; % # of observations from the previous period
T_PREV=NaN(t_bef,1);
T_PREV(1)=0.6;
MU_PREV=0.45;

for tt=2:t_bef
    T_PREV(tt)=exp(log(T_BASE(1,US,tt))-(1-RHO)*MU_PREV-RHO*(log(T_BASE(1,US,tt))-log(T_PREV(tt-1))));
end

params.prod = v2struct(T_BASE, T_PREV, MU, SIGMA, RHO);

T = PRODUCTIVITY_DGP(params); %objective productivity from the first period and on
T(1:J,1,1:30) = repmat(reshape(linspace(T_BASE(1,1,1),T_BASE(1,1,30).*2,30),1,1,30),J,1,1);
T(1:J,1,31:TIME) = T_BASE(1,1,30)*2;
T(1:J,2,1:30) = repmat(reshape(linspace(T_BASE(1,2,1),T_BASE(1,2,30).*1.5,30),1,1,30),J,1,1);
T(1:J,2,31:TIME) = T_BASE(1,2,30)*1.5;
T(1:J,3,1:30) = repmat(reshape(linspace(T_BASE(1,3,1),T_BASE(1,3,30).*1.3,30),1,1,30),J,1,1);
T(1:J,3,31:TIME) = T_BASE(1,3,30)*1.3;
T(1:J,4,1:30) = repmat(reshape(linspace(T_BASE(1,4,1),T_BASE(1,4,30).*1,30),1,1,30),J,1,1);
T(1:J,4,31:TIME) = T_BASE(1,4,30)*1;
T(1:J,5,1:30) = repmat(reshape(linspace(T_BASE(1,5,1),T_BASE(1,5,30).*0.8,30),1,1,30),J,1,1);
T(1:J,5,31:TIME) = T_BASE(1,5,30)*0.8;
% Derive Time difference productivity
T_HAT_SS = ones(J,N,TIME_SS);
T_HAT    = rand(J,N,TIME);
for t=1:TIME_SS-1
    T_HAT_SS(:,:,t+1) = T_BASE(:,:,t+1)./T_BASE(:,:,t); %relative change in productivity (US: 2 for all period, CHINA: 1 for all period except the first period)
end
for t=1:TIME-1
    T_HAT(:,:,t+1) = T(:,:,t+1)./T(:,:,t); %relative change in productivity (CHINA is catching up here)
end
params.prod = v2struct(MU, SIGMA, RHO, T_BASE, T_PREV, T, T_HAT, T_HAT_SS,t_bef);

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
    for t2=6:30
        E_T_hat(:,1:2,t1,t2) = 0.1;
        E_T_hat(:,3:N,t1,t2) = -0.1;
    end
end

% For heterogeneous agent case:
E_A_T_hat  = E_T_hat;
E_B_T_hat  = E_T_hat;
for t1=2:TIME
    for t2=1:5
        E_B_T_hat(:,1:2,t1,t2) = -0.3;
        E_B_T_hat(:,3:N,t1,t2) = 0.3;
    end
    for t2=6:30
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
params.belief=v2struct(W_TRUE,E_T_hat,E_T_hat_test,E_A_T_hat,E_B_T_hat,E_T_hat_pf);

end
