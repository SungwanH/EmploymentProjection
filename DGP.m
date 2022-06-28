function [eqm_dgp, approx_dgp] = DGP(params, eqm_nlpf, approx_nlpf,mat_pbp)
% This code gives DGP+PF path
% From t=ENDT, people gets perfect foresight
% Until ENDT, people imperfectly guess the future path of productivity
% input:
% params: parameters (structure)
% W: Weight in learning parameter (RHO_HAT); the current set of code does
% not use, in fact, W is used below as wages...so I commented out this
% input
% eqm_base: equilibrium outcomes obtained from BASE_LEVEL (perefect foresight level outcome)
% approx_base: approximation points obtained from BASE_LEVEL (perefect foresight level outcome)

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
for t1=1:1%ENDT+1
    tic
    if t1 ==1 
         V = zeros(R*J,TIME); %initial value for v_hat
         W = zeros(J,N,TIME);
    end
    [eqm_temp] = PBP_DYN(params, t1, t1, E_T_hat_test, kappa_hat, L, V, W, approx_nlpf,mat_pbp);
%    [eqm_temp] = PBP_DYN_SO(params, t1, t1, E_T_hat_test, kappa_hat, L, V, W, eqm_nlpf, approx_nlpf);
                      
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

%normalize
%{
for t=1:TIME
%    Ldyn(:,:,t)       =    Ldyn(:,:,t)./sum(sum(Ldyn(:,:,t)));
    for i=1:N*J
        pi(i,:,t)     =    pi(i,:,t)./sum(pi(i,:,t));
    end
    for i=1:R*J
        mu(i,:,t)     =    mu(i,:,t)./sum(mu(i,:,t));
    end
end
%}
% recover belief path in level
L_belief_dgp = zeros(J,R,TIME,ENDT+1);
for t1=1:ENDT+1
L_belief_dgp(:,:,:,t1) = exp(L_hat(:,:,:,t1)) .* eqm_nlpf.Ldyn;
end
X    = exp(X_dgp) .* eqm_nlpf.X;
wf00 = exp(w_dgp) .* eqm_nlpf.wf00;
pf00 = exp(p_dgp) .* eqm_nlpf.pf00;
VALjn00(1:J,1:R,1:TIME)   = exp(w_dgp(1:J,1:R,1:TIME))  .* exp(L_dgp(1:J,1:R,1:TIME)) .* eqm_nlpf.VALjn00(1:J,1:R,1:TIME);
VALjn00(1:J,R+1:N,1:TIME) = exp(w_dgp(1:J,R+1:N,1:TIME)).* eqm_nlpf.VALjn00(1:J,R+1:N,1:TIME); % L_hat for non-US is zero

varrho = zeros(N*J,N,TIME);
zeta   = zeros(N*J,J,TIME);
chi    = zeros(N*J,N,TIME);
lambda = zeros(R*J,R*J, TIME);
for t=1:TIME
    for n=1:N
        for j=1:J
            for ii=1:N
                varrho(n+(j-1)*N,ii,t) = pi(n+(j-1)*N,ii,t) * X(j,n,t) / (X(j,ii,t));
                chi(n+(j-1)*N,ii,t)    = pi(n+(j-1)*N,ii,t) * X(j,n,t) / VALjn00(j,ii,t);
            end
        end
    end
end
for t=1:TIME
    for k=1:J
        for ii=1:N
            for j=1:J
                zeta(ii+(k-1)*N,j,t) = VALjn00(k,ii,t)./ X(j,ii,t);
            end
        end
    end
end
for t=1:TIME
    for k=1:J
        for ii=1:R
            for j=1:J
                for n=1:R
                    if t==1
                        lambda(k+(ii-1)*J,j+(n-1)*J,1) = mu(k+(ii-1)*J,j+(n-1)*J,1)*Ldyn(k,ii,1)/Ldyn(j,n,1);
                    else
                        lambda(k+(ii-1)*J,j+(n-1)*J,t) = mu(k+(ii-1)*J,j+(n-1)*J,t-1)*Ldyn(k,ii,t-1)/Ldyn(j,n,t);
                    end
                end
            end
        end
    end
end

%save the outputs
eqm_dgp   = v2struct(Ldyn, L_belief_dgp, L_hat, v_hat, w_hat, p_hat, P_hat, pi_hat, mu_hat, X_hat, E_T_hat, ...
    L_dgp, w_dgp, p_dgp, P_dgp, pi_dgp, mu_dgp, X_dgp, X, VALjn00, wf00, pf00);
approx_dgp = v2struct(mu, pi, varrho, chi, zeta, lambda);

save('DATA/DGP.mat', 'eqm_dgp', 'approx_dgp'); 

toc

