function [eqm_dgp_so] = FO_OUT(params, eqm_nlpf, approx_nlpf,mat_pbp)
% Outer loop for First-order approximation (Used in test code)
% Same with DGP.m
%% Roll down parameters
v2struct(params.envr)
v2struct(params.belief)


%% Roll down Baseline outcomes
v2struct(eqm_nlpf)
v2struct(approx_nlpf)

%mat_pbp =ones(N,J,TIME);
% period by period equilibrium
% Initialize
v_hat     = zeros(R*J,TIME); 
w_hat     = NaN(J,N,TIME,ENDT+1); 
p_hat     = NaN(J,N,TIME,ENDT+1); 
X_hat     = NaN(J,N,TIME,ENDT+1); 
pi_hat    = NaN(N*J,N,TIME,ENDT+1); 
mu_hat    = NaN(R*J,R*J,TIME,ENDT+1); 
L_hat     = NaN(J,R,TIME,ENDT+1); 
L_dgp     = zeros(J,R,TIME);
kappa_hat = zeros(N*J,N,TIME);
L = zeros(R*J,1); %initial deviation from the path around which we are linearizing
for t1=1:ENDT+1
    tic
    if t1 ==1 
         V = zeros(R*J,TIME); %initial value for v_hat
         W = zeros(J,N,TIME);
    end
    [eqm_temp] = PBP_DYN(params, t1, t1, E_T_hat_test, kappa_hat, L, V, W, approx_nlpf,mat_pbp);
%    eqm_temp.L(:,t1+1)
%    [eqm_temp_iter] = PBP_DYN(params, t1, t1, E_T_hat_test, kappa_hat, L, V, W, approx_nlpf);
%    eqm_temp_iter.L(:,t1+1)
                      
    L = eqm_temp.L(:,t1+1); % update the initial value of labor 
    V = eqm_temp.v; % update the initial value of V
    W = eqm_temp.w; % update the initial value of wage
    %save the belief path
    w_hat(:,:,:,t1)  = eqm_temp.w;
    p_hat(:,:,:,t1)  = eqm_temp.p;
    P_hat(:,:,:,t1)  = eqm_temp.P;
    v_hat(:,:,t1)    = eqm_temp.v;
    L_hat(:,:,:,t1)  = reshape(eqm_temp.L,J,R,TIME);
    pi_hat(:,:,:,t1) = eqm_temp.pi;
    mu_hat(:,:,:,t1) = eqm_temp.mu;
    X_hat(:,:,:,t1)  = eqm_temp.X;
    %save the realized path
    if t1<ENDT+1
        L_dgp(:,:,t1+1) = reshape(eqm_temp.L(:,t1+1),J,R);
        w_dgp(:,:,t1)   = eqm_temp.w(:,:,t1);
        p_dgp(:,:,t1)   = eqm_temp.p(:,:,t1);
        P_dgp(:,:,t1)   = eqm_temp.P(:,:,t1);
        pi_dgp(:,:,t1)  = eqm_temp.pi(:,:,t1);
        mu_dgp(:,:,t1)  = eqm_temp.mu(:,:,t1);
        X_dgp(:,:,t1)   = eqm_temp.X(:,:,t1);
    else
        L_dgp(:,:,t1+1:TIME) = reshape(eqm_temp.L(:,t1+1:TIME),J,R,TIME-t1);
        w_dgp(:,:,t1:TIME)   = eqm_temp.w(:,:,t1:TIME);
        p_dgp(:,:,t1:TIME)   = eqm_temp.p(:,:,t1:TIME);
        P_dgp(:,:,t1:TIME)   = eqm_temp.P(:,:,t1:TIME);
        pi_dgp(:,:,t1:TIME)  = eqm_temp.pi(:,:,t1:TIME);
        mu_dgp(:,:,t1:TIME)  = eqm_temp.mu(:,:,t1:TIME);
        X_dgp(:,:,t1:TIME)   = eqm_temp.X(:,:,t1:TIME);
    end
    t1 %display t1
    toc
end

% generate level values
Ldyn = exp(L_dgp)  .* eqm_nlpf.Ldyn;
pi   = exp(pi_dgp) .* approx_nlpf.pi;
mu   = exp(mu_dgp) .* approx_nlpf.mu;

wf00 = exp(w_dgp) .* eqm_nlpf.wf00;
pf00 = exp(p_dgp) .* eqm_nlpf.pf00;
eqm_dgp_so   = v2struct(Ldyn, wf00, pf00);
%approx_dgp_so = v2struct(mu, pi, varrho, chi, zeta, lambda);