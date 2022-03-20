function [params] = PARAMS()

TIME    = 100; 
TIME_SS = 200;
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


% Baseline: RHO 0.85 ENDT 20
SIGMA       = 0.01; % s.d. of catch-up
W_TRUE      = 0.1; % True weight on RHO when deriving RHO_HAT
ENDT        = 10; %Last period of learning (PF starts from here)
ENDT_SAMPLE = 1;
ENDT_DGP    = 1;
EPST        = ENDT; %Length of epsilon draw

% Learning
RHO     = 0.80; % catch-up speed (Higher RHO -> slower catch-up)
MU      = 0.15;



%% FAKE DATA
N       = 10; 
R       = 10; %number of US regions
C       = 0; %number of non-US countries (To be used later)
J       = 1; %number of sectors
US      = 2;
CHINA   = 1; %region number for now assume the shock is to California

ALPHAS=1/J*ones(J,N);
THETA=4*ones(J,1);
GAMMA=0.5*ones(J,N);

BETA = 0.5;
NU = 3*BETA;
%BETA =0.85;



%% FAKE DATA ENDS HERE


% Technical parameters
ESTM_BOTH   = 1; % if ESTM_BOTH =0, estimate MU only. if =1, estimate both MU and RHO
UPDT_V      = 0.5; %update speed for value loop (lower value->conservative)
UPDT_W      = 0.3; %update speed for wage loop (lower value->conservative)
UPDT_V_NL   = 0.5; %update speed for nonlinear value loop (lower value->conservative)
UPDT_W_NL   = 0.3; %update speed for nonlinear wage loop (lower value->conservative)
TOL_NL      = 1E-7;  %tolerance rate for nonlinear dynamic equilibrium (outerloop)
TOL_NL_TEMP = 1E-7;  %tolerance rate for nonlinear temporary equilibrium (inner loop)
TOLDYN      = 1E-7;  %tolerance rate for linear dynamic equilibrium
TOLTEMP     = 1E-7;  % tolerance rate for linear temporary equilibrium
MAXIT       = 1E+8; %maximum number of iterations


params.envr = v2struct(TIME, TIME_SS, N, R, C, J, US, CHINA, ENDT, ENDT_SAMPLE, ENDT_DGP, EPST); %remoevd TAU from this structure
params.modl = v2struct(THETA, NU, BETA, ALPHAS, THETA, GAMMA);
params.tech = v2struct(ESTM_BOTH, UPDT_V, UPDT_W, UPDT_V_NL, UPDT_W_NL, TOL_NL, TOL_NL_TEMP, TOLDYN, TOLTEMP, MAXIT);


%% Productivity
% Baseline productivity (Used in deriving initial steady state)
T_BASE = ones(J,N,TIME*3)*2; % US (& other countries except China) productivity is constant for all period
T_BASE(:,CHINA,1:TIME*3) = 1;
%Previous productivity of China (used in estimating RHO and MU)
T_PREV = repmat(linspace(0.9,1,20),J,1); %Assume previous 20 periods' CHINA productivity was from 0.0to 0.1
params.prod = v2struct(T_BASE, T_PREV, MU, SIGMA, RHO, W_TRUE);

T = PRODUCTIVITY_DGP(params); %objective productivity (in level)
% Derive Time difference productivity
T_HAT_SS = ones(J,N,TIME_SS);
T_HAT    = ones(J,N,TIME);
for t=1:TIME_SS-1
    T_HAT_SS(:,:,t+1) = T_BASE(:,:,t+1)./T_BASE(:,:,t); %relative change in productivity (US: 2 for all period, CHINA: 1 for all period except the first period)
end
for t=1:TIME-1
    T_HAT(:,:,t+1) = T(:,:,t+1)./T(:,:,t); %relative change in productivity (CHINA is catching up here)
end

params.prod = v2struct(MU, SIGMA, RHO, W_TRUE, T_BASE, T_PREV, T, T_HAT, T_HAT_SS);

% Derive productivity belief log deviation from objective productivity
T_belief = BELIEF(params, W_TRUE);
E_T_hat  = zeros(J,N,TIME,ENDT+1); % Except CHINA, productivity is constant
for tt=1:ENDT+1
    for j=1:J
        E_T_hat(j,CHINA,:,tt) = log(T_belief(j,CHINA,:,tt)) - log(T(j,CHINA,:));
    end
end
% perfect foresight deviation from belief
E_T_hat_pf = zeros(J,N,TIME,ENDT+1);
for tt=1:ENDT+1
    for j=1:J
       E_T_hat_pf(j,CHINA,:,tt) = -log(T_belief(j,CHINA,:,tt)) + log(T(j,CHINA,:));
    end
end
% Percentage deviation
%E_T_hat = zeros(J,N,TIME,ENDT+1); % Except CHINA, productivity is constant
%for tt=1:ENDT+1
%    for j=1:J
%         E_T_hat(j,CHINA,:,tt) = (T_belief(j,CHINA,:,tt)-T(j,CHINA,:))./(T(j,CHINA,:));
%    end
%end
%E_T_hat_pf = zeros(J,N,TIME,ENDT+1);
%for tt=1:ENDT+1
%    for j=1:J
%         E_T_hat_pf(j,CHINA,:,tt) = (T(j,CHINA,:)-T_belief(j,CHINA,:,tt))./(T_belief(j,CHINA,:,tt));
%    end
%end
% add deviation variables 

%{
%% Productivity: one time productivity shock in China (manufacturing only) 
T_BASE = ones(J,N,TIME*3)*2; % 
T=T_BASE(:,:,1:TIME);
T(1,CHINA,2:end)=0.5;

params.prod = v2struct(MU, SIGMA, RHO, W_TRUE, T_BASE, T_PREV, T, T_HAT, T_HAT_SS);

% Derive productivity belief log deviation from objective productivity
T_belief = BELIEF(params, W_TRUE);
for t=1:TIME_SS-1
    T_HAT_SS(:,:,t+1) = T_BASE(:,:,t+1)./T_BASE(:,:,t); %relative change in productivity (US: 2 for all period, CHINA: 1 for all period except the first period)
end
for t=1:TIME-1
    T_HAT(:,:,t+1) = T(:,:,t+1)./T(:,:,t); %relative change in productivity (CHINA is catching up here)
end


%E_T_hat  = zeros(J,N,TIME,ENDT+1); % Except CHINA, productivity is constant
for tt=1:ENDT+1
    E_T_hat(:,CHINA,:,tt) = -log(T_BASE(:,CHINA,1:TIME)) + log(T(:,CHINA,:));
end
%}

%{
for tt=1:ENDT+1
    for j=1:J
        E_T_hat(j,CHINA,:,tt) = log(T_belief(j,CHINA,:,tt)) - log(T(j,CHINA,:));
    end
end
% perfect foresight deviation from belief
E_T_hat_pf = zeros(J,N,TIME,ENDT+1);
for tt=1:ENDT+1
    for j=1:J
       E_T_hat_pf(j,CHINA,:,tt) = -log(T_belief(j,CHINA,:,tt)) + log(T(j,CHINA,:));
    end
end
%for tt=1:ENDT+1
%    for j=1:J
%        E_T_hat(j,CHINA,:,tt) = log(T_belief(j,CHINA,:,tt)) - log(T(j,CHINA,:));
%    end
%end
%}
params.prod = v2struct(MU, SIGMA, RHO, W_TRUE, T_BASE, T_PREV, T, T_HAT, T_HAT_SS, T_belief, E_T_hat, E_T_hat_pf);

end
