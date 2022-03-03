function [params] = PARAMS()

TIME    = 100; 
TIME_SS = 200;
N       = 87; 
R       = 50; %number of US regions
C       = 37; %number of non-US countries (To be used later)
J       = 4; %number of sectors
US      = 1;
%CHINA   = 57; %region number
CHINA   = 57; %region number

% Elasticities
load('DATA/BASE_FOURSECTOR.mat', 'alphas','T','gamma')
ALPHAS = alphas; %consumption share
THETA  = 1./T; 
GAMMA      = gamma; % share of labor in production
clear alphas
clear T
%THETA   = 5; %trade elasticity
NU      = 5.3436; %dispersion of taste shocks (In CDP: 5.3436)
%NU      = 3; %dispersion of taste shocks (In CDP: 5.3436)
%BETA    = 0.9605; %discount rate (0.99^4: transformed from quarterly to yearly)
BETA    = 0.9227; %discount rate (0.99^8: transformed from quarterly to bi-annual)

TAU     = ones(N*J,N); %iceberg trade cost
%TAU     = (ones(N,N)-eye(N)).*(1+rand(N))+eye(N);


% Baseline: RHO 0.85 ENDT 20
SIGMA   = 0.01; % s.d. of catch-up
W_TRUE = 0.5; % True weight on RHO when deriving RHO_HAT
ENDT    = 10; %Last period of learning (PF starts from here)
ENDT_SAMPLE = 10;
ENDT_DGP    = 10;
EPST    = ENDT; %Length of epsilon draw

% Learning
RHO     = 0.80; % catch-up speed (Higher RHO -> slower catch-up)
MU      = 0.15;

% Technical parameters
ESTM_BOTH   = 1; % if ESTM_BOTH =0, estimate MU only. if =1, estimate both MU and RHO
UPDT_V      = 0.3; %update speed for value loop (lower value->conservative)
UPDT_W      = 0.3; %update speed for wage loop (lower value->conservative)
UPDT_V_NL   = 0.5; %update speed for nonlinear value loop (lower value->conservative)
UPDT_W_NL   = 0.3; %update speed for nonlinear wage loop (lower value->conservative)
TOL_NL      = 1E-7;  %tolerance rate for nonlinear dynamic equilibrium (outerloop)
TOL_NL_TEMP = 1E-7;  %tolerance rate for nonlinear temporary equilibrium (inner loop)
TOLDYN      = 1E-7;  %tolerance rate for linear dynamic equilibrium
TOLTEMP     = 1E-7;  % tolerance rate for linear temporary equilibrium
MAXIT       = 1E+8; %maximum number of iterations


params.envr   = v2struct(TIME, TIME_SS, N, R, C, J, US, CHINA, TAU, ENDT, ENDT_SAMPLE, ENDT_DGP, EPST);
params.modl   = v2struct(THETA, NU, BETA, ALPHAS, THETA, GAMMA);
params.tech   = v2struct(ESTM_BOTH, UPDT_V, UPDT_W, UPDT_V_NL, UPDT_W_NL, TOL_NL, TOL_NL_TEMP, TOLDYN, TOLTEMP, MAXIT);


%% Productivity
% Baseline productivity (Used in deriving initial steady state)
T_BASE = ones(J,N,TIME*3)*2; % US (& other countries except China) productivity is constant for all period
T_BASE(1:J,CHINA,2:TIME*3) = 1;
%Previous productivity of China (used in estimating RHO and MU)
T_PREV = repmat(linspace(0.97,0.99,20),4,1); %Assume previous 20 periods' CHINA productivity was from 0.97 to 0.99
params.prod   = v2struct(T_BASE, T_PREV, MU, SIGMA, RHO, W_TRUE);

T = PRODUCTIVITY_DGP(params); %objective productivity (in level)


T_HAT_SS = ones(J,N,TIME_SS);
T_HAT = ones(J,N,TIME);
for t=1:TIME_SS-1
    T_HAT_SS(:,:,t+1)=T_BASE(:,:,t+1)./T_BASE(:,:,t); %relative change in technology (US: 2 for all period, CHINA: 1 for all period except the first period)
end
for t=1:TIME-1
    T_HAT(:,:,t+1)=T(:,:,t+1)./T(:,:,t); %relative change in technology (CHINA is catching up here)
end


params.prod   = v2struct(T_BASE, T_PREV, T, T_HAT, T_HAT_SS, MU, SIGMA, RHO, W_TRUE);

end
