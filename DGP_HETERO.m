function [eqm_dgp_hetero, approx_dgp_hetero] = DGP_HETERO(params, eqm_nlpf, approx_nlpf, mat_pbp)
% This code gives DGP+PF path
% From t=ENDT+1, people gets perfect foresight
% Until ENDT, people imperfectly guess the future path of productivity
% input:
% params: parameters (structure)
% W: Weight in learning parameter (RHO_HAT); the current set of code does
% not use, in fact, W is used below as wages...so I commented out this
% input
% eqm_base: equilibrium outcomes obtained from BASE_LEVEL (perefect foresight level outcome)
% approx_base: approximation points obtained from BASE_LEVEL (perefect foresight level outcome)
% mat_pbp: matrices used in computing temporary equilibrium (inversion)

%% Roll down parameters
v2struct(params.envr)
v2struct(params.belief)
v2struct(params.modl)
v2struct(params.hetero)
v2struct(params.tech);


%% Roll down Baseline outcomes
v2struct(eqm_nlpf)
v2struct(approx_nlpf)

% period by period equilibrium
% Initialize
E_A_T_hat   = params.belief.E_A_T_hat;
E_B_T_hat   = params.belief.E_B_T_hat;
kappa_A_hat = zeros(N*J,N,TIME);
kappa_B_hat = zeros(N*J,N,TIME);


L_A_hat    = NaN(J,R,TIME,ENDT+1); 
w_A_hat    = NaN(J,N,TIME,ENDT+1); 
p_A_hat    = NaN(J,N,TIME,ENDT+1); 
P_A_hat    = NaN(1,N,TIME,ENDT+1); 
X_A_hat    = NaN(J,N,TIME,ENDT+1); 
pi_A_hat   = NaN(N*J,N,TIME,ENDT+1); 
mu_A_hat   = NaN(R*J,R*J,TIME,ENDT+1); 
L_B_hat    = NaN(J,R,TIME,ENDT+1); 
w_B_hat    = NaN(J,N,TIME,ENDT+1); 
p_B_hat    = NaN(J,N,TIME,ENDT+1); 
P_B_hat    = NaN(1,N,TIME,ENDT+1); 
X_B_hat    = NaN(J,N,TIME,ENDT+1); 
pi_B_hat   = NaN(N*J,N,TIME,ENDT+1); 
mu_B_hat   = NaN(R*J,R*J,TIME,ENDT+1); 
L_dgp      = zeros(J,R,TIME);
L_A_dgp    = zeros(J,R,TIME);
L_B_dgp    = zeros(J,R,TIME);
w_dgp      = NaN(J,R,TIME);
p_dgp      = NaN(J,R,TIME);
P_dgp      = NaN(1,R,TIME);
X_dgp      = NaN(J,R,TIME);
pi_dgp     = NaN(N*J,N,TIME);
mu_dgp     = NaN(R*J,R*J,TIME);

L_A_init = zeros(N*J,1); %initial deviation from the path around which we are linearizing
L_B_init = zeros(N*J,1); %initial deviation from the path around which we are linearizing
L_A_2 = zeros(N*J,1);
L_B_2 = zeros(N*J,1);   

for t1=1:ENDT+1
L_B_2 = zeros(N*J,1);
ITER_FP = 0;
maxerror=1;
    while (ITER_FP <= MAXIT) && (maxerror > TOLFP)
        % solve for type A's problem given type B's labor at t1+1
        [eqm_A] = PBP_DYN_HETERO(params, t1, t1, E_A_T_hat, kappa_A_hat, L_A_init, L_B_2, approx_nlpf, mat_pbp, 1);
        L_A_2_new = eqm_A.L_own(:,t1+1);
        % solve for type B's problem given type A's labor (L_A_2_new) at t1+1
        [eqm_B] = PBP_DYN_HETERO(params, t1, t1, E_B_T_hat, kappa_B_hat, L_B_init, L_A_2_new, approx_nlpf, mat_pbp, 2);
        L_B_2_new = eqm_B.L_own(:,t1+1);

        maxerror=max(abs(L_B_2_new-L_B_2))
        if ITER_FP >5000
            maxerror
            t1
            disp('Fixed Point loop err')
            ITER_FP
            stop
        end

        L_B_2 = L_B_2_new; %update
        ITER_FP = ITER_FP+1
    end
    L_A_init = L_A_2_new; % update the initial value of labor of type A
    L_B_init = L_B_2_new; % update the initial value of labor of type B
    L_A_hat(:,:,:,t1)  = reshape(eqm_A.L_own,J,R,TIME);
    w_A_hat(:,:,:,t1)  = eqm_A.w;
    p_A_hat(:,:,:,t1)  = eqm_A.p;
    P_A_hat(:,:,:,t1)  = eqm_A.P;
    pi_A_hat(:,:,:,t1) = eqm_A.pi;
    mu_A_hat(:,:,:,t1) = eqm_A.mu;
    X_A_hat(:,:,:,t1)  = eqm_A.X;
    L_B_hat(:,:,:,t1)  = reshape(eqm_B.L_own,J,R,TIME);
    w_B_hat(:,:,:,t1)  = eqm_B.w;
    p_B_hat(:,:,:,t1)  = eqm_B.p;
    P_B_hat(:,:,:,t1)  = eqm_B.P;
    pi_B_hat(:,:,:,t1) = eqm_A.pi;
    mu_B_hat(:,:,:,t1) = eqm_A.mu;
    X_B_hat(:,:,:,t1)  = eqm_A.X;
    if t1<ENDT+1
%        L_A_dgp(:,:,t1+1)=reshape(L_A_2_new,J,R);
%        L_B_dgp(:,:,t1+1)=reshape(L_B_2_new,J,R);
        L_dgp(:,:,t1+1)   = reshape((share_A.*eqm_A.L_own(:,t1+1) + (1-share_A).*eqm_B.L_own(:,t1+1)),J,R);
        L_A_dgp(:,:,t1+1) = reshape(eqm_A.L_own(:,t1+1),J,R);
        L_B_dgp(:,:,t1+1) = reshape(eqm_B.L_own(:,t1+1),J,R);
        w_dgp(:,:,t1)     = eqm_A.w(:,:,t1);
        p_dgp(:,:,t1)     = eqm_A.p(:,:,t1);
        P_dgp(:,:,t1)     = eqm_A.P(:,:,t1);
        pi_dgp(:,:,t1)    = eqm_A.pi(:,:,t1);
        mu_dgp(:,:,t1)    = eqm_A.mu(:,:,t1);
        X_dgp(:,:,t1)     = eqm_A.X(:,:,t1);
    else
        L_A_dgp(:,:,t1+1:TIME) = reshape(eqm_A.L_own(:,t1+1:TIME),J,R,TIME-t1);
        L_B_dgp(:,:,t1+1:TIME) = reshape(eqm_B.L_own(:,t1+1:TIME),J,R,TIME-t1);
        L_dgp(:,:,t1+1:TIME)   = reshape((share_A.*eqm_A.L_own(:,t1+1:TIME) + (1-share_A).*eqm_B.L_own(:,t1+1:TIME)),J,R,TIME-t1);
        w_dgp(:,:,t1:TIME)     = eqm_A.w(:,:,t1:TIME);
        p_dgp(:,:,t1:TIME)     = eqm_A.p(:,:,t1:TIME);
        P_dgp(:,:,t1:TIME)     = eqm_A.P(:,:,t1:TIME);
        pi_dgp(:,:,t1:TIME)    = eqm_A.pi(:,:,t1:TIME);
        mu_dgp(:,:,t1:TIME)    = eqm_A.mu(:,:,t1:TIME);
        X_dgp(:,:,t1:TIME)     = eqm_A.X(:,:,t1:TIME);
    end
t1
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
L_A_belief_dgp = zeros(J,R,TIME,ENDT+1);
L_B_belief_dgp = zeros(J,R,TIME,ENDT+1);
for t1=1:ENDT+1
L_A_belief_dgp(:,:,:,t1) = exp(L_A_hat(:,:,:,t1)) .* eqm_nlpf.Ldyn;
L_B_belief_dgp(:,:,:,t1) = exp(L_B_hat(:,:,:,t1)) .* eqm_nlpf.Ldyn;
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
eqm_dgp_hetero   = v2struct(Ldyn, L_A_belief_dgp, L_B_belief_dgp, L_A_hat, L_B_hat, w_A_hat, w_B_hat, ...
    p_A_hat, p_B_hat, P_A_hat, P_B_hat, pi_A_hat, pi_B_hat, mu_A_hat, mu_B_hat, ...
    X_A_hat, X_B_hat, E_A_T_hat, E_B_T_hat, L_dgp, L_A_dgp, L_B_dgp, w_dgp, p_dgp, P_dgp, pi_dgp, mu_dgp, X_dgp, X, VALjn00, wf00, pf00);
approx_dgp_hetero = v2struct(mu, pi, varrho, chi, zeta, lambda);

save('DATA/DGP_HETERO.mat', 'eqm_dgp_hetero', 'approx_dgp_hetero'); 

toc

