function [params] = PARAMS_MPS(W_TRUE, DRAW_SD, NUM_SIM)
TIME    =50; 
TIME_SS = 50;
C=28;% # non-US of countries
R=50;% # of US states
N=R+C; % # of US states & other countries
J=3; % # of sectors
US      = 1;
CHINA   = 59; %region number for now assume the shock is to California

% Elasticities
load('Base_year_three_sectors.mat', 'X','B')
ALPHAS = X; % Final consumption share
T = 0.22*ones(J,1);
THETA  = 1./T; 
GAMMA      = B; % share of labor in production
clear X
clear B

NU      = 5.3436; %dispersion of taste shocks (In CDP: 5.3436)
BETA    = 0.9227; %discount rate (0.99^8: transformed from quarterly to bi-annual)


% Baseline:  
SIGMA       = 0.0; % s.d. of catch-up
ENDT        = 10; %Last period of learning (PF starts from here)
ENDT_SAMPLE = 1;
ENDT_DGP    = 1;
EPST        = 10; %Length of epsilon draw

CRRA_PAR = 5; %CRRA risk aversion parameter
%BETA = 0.8;
%NU = 1.3;



%% Technical parameters
ESTM_BOTH   = 0; % if ESTM_BOTH =0, estimate MU only. if =1, estimate both MU and RHO
UPDT_V      = 0.4; %update speed for value loop (lower value->conservative)
UPDT_W      = 0.2; %update speed for wage loop (lower value->conservative)
UPDT_V_NL   = 0.1; %update speed for nonlinear value loop (lower value->conservative)
UPDT_W_NL   = 0.1; %update speed for nonlinear wage loop (lower value->conservative)
TOL_NL      = 1E-6;  %tolerance rate for nonlinear dynamic equilibrium (outerloop)
TOL_NL_TEMP = 1E-6;  %tolerance rate for nonlinear temporary equilibrium (inner loop)
TOLDYN      = 1E-7;  %tolerance rate for linear dynamic equilibrium
TOLTEMP     = 1E-7;  % tolerance rate for linear temporary equilibrium
TOLFP       = 1E-7;  % tolerance rate for fixed point iteration
MAXIT       = 1E+8; %maximum number of iterations


%NUM_SIM = 50; % number of simulations
%DRAW_SD = 0.1; % standard error of draws
AR_GR   = 0;  % Whether growth rate follows AR1 or not
params.envr = v2struct(TIME, TIME_SS, N, R, C, J, US, CHINA, ENDT, ENDT_SAMPLE, ENDT_DGP, EPST, NUM_SIM, DRAW_SD); %remoevd TAU from this structure
params.modl = v2struct(THETA, NU, BETA, ALPHAS, THETA, GAMMA, CRRA_PAR);
params.tech = v2struct(ESTM_BOTH, UPDT_V, UPDT_W, UPDT_V_NL, UPDT_W_NL, TOL_NL, TOL_NL_TEMP, TOLDYN, TOLTEMP, TOLFP, MAXIT);


%% Productivity
% Baseline productivity (Used in deriving initial steady state)
T_BASE = ones(J,N,TIME*3)*2; % US (& other countries except China) productivity is constant for all period
T_BASE(:,CHINA,1:TIME*3) = 2; % We start from same level

% Learning
RHO     = 0.85; % catch-up speed (Higher RHO -> slower catch-up)
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

%% Baseline

T = PRODUCTIVITY_DGP(params); %objective productivity from the first period and on (not used in testing second order)

params.prod = v2struct(T_BASE, T_PREV, MU, SIGMA, RHO, T);

T_SIM = zeros(J,N,TIME,NUM_SIM,ENDT+1);

for s=1:NUM_SIM
    for t=1:ENDT+1
        T_SIM(:,:,:,s,t) = T;
    end
end
for s=1:NUM_SIM
    for t1=1:ENDT+1
        T_SIM(1,5,t1+1:ENDT,s,t1) = DRAW_SD.*randn(1,1,ENDT-t1) + T(1,5,t1+1:ENDT);
%        T_SIM(1,5,t1+1:ENDT,1:NUM_SIM,t1) = DRAW_SD.*randn(1,1,ENDT-t1,NUM_SIM) + T(1,5,t1+1:ENDT);
    end
end

for s=1:NUM_SIM
    for t1=1:ENDT+1
        params_sim = params;
        params_sim.prod.T = T_SIM(:,:,:,s,t1);
        T_belief(:,:,:,:,s,t1) = BELIEF(params_sim, 1);
    end
end

%T_belief = BELIEF(params, W_TRUE);
E_T_hat_sim = zeros(J,N,TIME,ENDT+1,NUM_SIM,ENDT+1); % given the economy is shifted into the t3, belief deviation from simulated productivity at period t1 when the belief is made at t2 
for t3=1:ENDT+1
    for s=1:NUM_SIM
        for t2=1:ENDT+1
            for t1=1:TIME
%                E_T_hat_sim(:,:,t1,t2,s,t3) = log(T_belief(:,:,t1,t2,s,t3))- log(T(:,:,t1,s,t3));
                E_T_hat_sim(:,:,t1,t2,s,t3) = (log(T_belief(:,:,t1,t2,s,t3))-log(T(:,:,t1)));
            end
        end
    end
end
E_T_hat_sim(:,:,ENDT+1:TIME,:,:,:) = zeros(J,N,TIME-ENDT,ENDT+1,NUM_SIM,ENDT+1);

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

params.belief=v2struct(W_TRUE,E_T_hat_sim);

end
