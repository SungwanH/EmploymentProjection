function [eqm_sim, approx_sim] = DGP_SIM(params, t, L_R, eqm_nlpf, approx_nlpf, mat_pbp)
% This code solves the first-order approximation of the equilibrium given
% simulated productivity path and correspopnding belief from period 1

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
w_hat     = NaN(J,N,TIME,ENDT+1,NUM_SIM); 
p_hat     = NaN(J,N,TIME,ENDT+1,NUM_SIM); 
X_hat     = NaN(J,N,TIME,ENDT+1,NUM_SIM); 
pi_hat    = NaN(N*J,N,TIME,ENDT+1,NUM_SIM); 
mu_hat    = NaN(R*J,R*J,TIME,ENDT+1,NUM_SIM); 
L_hat     = NaN(J,R,TIME,ENDT+1,NUM_SIM); 
L_sim     = zeros(J,R,TIME,NUM_SIM);
kappa_hat = zeros(N*J,N,TIME);
L_belief_sim = zeros(J,R,TIME,ENDT+1,NUM_SIM);
%L = zeros(R*J,1); %initial deviation from the path around which we are linearizing
L = L_R; %initial deviation from the path around which we are linearizing

for s=1:NUM_SIM
    L = L_R;
    V = zeros(R*J,TIME);
    for t1=t:ENDT+1
        %tic
        if t1 ==t
             V = zeros(R*J,TIME); %initial value for v_hat
%             W = zeros(J,N,TIME);
        end
        [eqm_temp] = PBP_DYN(params, t1, t1, E_T_hat_sim(:,:,:,:,s,t), kappa_hat, L, V, eqm_nlpf, approx_nlpf, mat_pbp);
    %    [eqm_temp] = PBP_DYN_SO(params, t1, t1, E_T_hat_test, kappa_hat, L, V, W, eqm_nlpf, approx_nlpf);
                          
        L = eqm_temp.L(:,t1+1); % update the initial value of labor 
        V = eqm_temp.v; % update the initial value of V
        W = eqm_temp.w; % update the initial value of wage
        %save the belief path
        w_hat(:,:,:,t1,s)  = eqm_temp.w;
        p_hat(:,:,:,t1,s)  = eqm_temp.p;
        P_hat(:,:,:,t1,s)  = eqm_temp.P;
        v_hat(:,:,t1,s)    = eqm_temp.v;
        L_hat(:,:,:,t1,s)  = reshape(eqm_temp.L,J,R,TIME);
        pi_hat(:,:,:,t1,s) = eqm_temp.pi;
        mu_hat(:,:,:,t1,s) = eqm_temp.mu;
        X_hat(:,:,:,t1,s)  = eqm_temp.X;
        %save the realized path
        if t1<ENDT+1
            L_sim(:,:,t1+1,s) = reshape(eqm_temp.L(:,t1+1),J,R);
            v_sim(:,t1+1,s)   = eqm_temp.v(:,t1+1); %check whether this should be t1 or t1+1 
            w_sim(:,:,t1,s)   = eqm_temp.w(:,:,t1);
            p_sim(:,:,t1,s)   = eqm_temp.p(:,:,t1);
            P_sim(:,:,t1,s)   = eqm_temp.P(:,:,t1);
            pi_sim(:,:,t1,s)  = eqm_temp.pi(:,:,t1);
            mu_sim(:,:,t1,s)  = eqm_temp.mu(:,:,t1);
            X_sim(:,:,t1,s)   = eqm_temp.X(:,:,t1);
            T_sim(:,:,t1,s)   = E_T_hat_sim(:,:,t1,t1,s);
        else
            L_sim(:,:,t1+1:TIME,s) = reshape(eqm_temp.L(:,t1+1:TIME),J,R,TIME-t1);
            v_sim(:,t1+1:TIME,s)   = eqm_temp.v(:,t1+1:TIME); %check whether this should be t1 or t1+1
            w_sim(:,:,t1:TIME,s)   = eqm_temp.w(:,:,t1:TIME);
            p_sim(:,:,t1:TIME,s)   = eqm_temp.p(:,:,t1:TIME);
            P_sim(:,:,t1:TIME,s)   = eqm_temp.P(:,:,t1:TIME);
            pi_sim(:,:,t1:TIME,s)  = eqm_temp.pi(:,:,t1:TIME);
            mu_sim(:,:,t1:TIME,s)  = eqm_temp.mu(:,:,t1:TIME);
            X_sim(:,:,t1:TIME,s)   = eqm_temp.X(:,:,t1:TIME);
            T_sim(:,:,t1:TIME,s)   = E_T_hat_sim(:,:,t1:TIME,t1,s);
        end
        t1 %display t1
     %   toc
    end
    % generate level values
    Ldyn(:,:,:,s) = exp(L_sim(:,:,:,s))  .* eqm_nlpf.Ldyn;
    pi(:,:,:,s)   = exp(pi_sim(:,:,:,s)) .* approx_nlpf.pi;
    mu(:,:,:,s)   = exp(mu_sim(:,:,:,s)) .* approx_nlpf.mu;
    for t1=1:ENDT+1
        L_belief_sim(:,:,:,t1,s) = exp(L_hat(:,:,:,t1,s)) .* eqm_nlpf.Ldyn;
    end
end

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
%{
% recover belief path in level
for s=1:NUM_SIM
end
X    = exp(X_sim) .* eqm_nlpf.X;
wf00 = exp(w_sim) .* eqm_nlpf.wf00;
pf00 = exp(p_sim) .* eqm_nlpf.pf00;
VALjn00(1:J,1:R,1:TIME)   = exp(w_sim(1:J,1:R,1:TIME))  .* exp(L_sim(1:J,1:R,1:TIME)) .* eqm_nlpf.VALjn00(1:J,1:R,1:TIME);
VALjn00(1:J,R+1:N,1:TIME) = exp(w_sim(1:J,R+1:N,1:TIME)).* eqm_nlpf.VALjn00(1:J,R+1:N,1:TIME); % L_hat for non-US is zero

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
%}
%save the outputs
%eqm_sim   = v2struct(Ldyn, L_belief_sim, L_hat, v_hat, w_hat, p_hat, P_hat, pi_hat, mu_hat, X_hat, ...
%    L_sim, w_sim, p_sim, P_sim, pi_sim, mu_sim, X_sim, X, VALjn00, wf00, pf00);
%approx_sim = v2struct(mu, pi, varrho, chi, zeta, lambda);

eqm_sim    = v2struct(Ldyn, L_belief_sim, L_hat, v_hat, w_hat, p_hat, P_hat, pi_hat, mu_hat, X_hat, ...
                      L_sim, v_sim, w_sim, p_sim, P_sim, pi_sim, mu_sim, X_sim, T_sim);
approx_sim = v2struct(mu, pi);
save('DATA/SIM.mat', 'eqm_sim', 'approx_sim'); 

%toc
