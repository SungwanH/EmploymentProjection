function eqm_recur = RECURSIVE(params, W, initial, eqm, approx)
% This code gives (linearized) period by period deviation from DGP+PF path
% input:
% params: parameters
% W: weight
% initial: initial v_hat, w_hat
% eqm: equilibrium outcomes obtained from linearized Data path (DGP)
% approx: approximation points obtained from linearized Data path (DGP)


%% Roll down parameters

v2struct(params.envr);
v2struct(params.tech);
E_T_hat = params.prod.E_T_hat;      % Belief on productivity (deviation from T)
E_T_hat_pf = params.prod.E_T_hat_pf;% perfect foresight belief (deviation from T_belief)

%Initial approximation point
Ldyn   = eqm.Ldyn;
wf00   = eqm.wf00;
pf00   = eqm.pf00;
X      = eqm.X;
VALjn00= eqm.VALjn00;
mu     = approx.mu;
pi     = approx.pi;
varrho = approx.varrho;
chi    = approx.chi;
zeta   = approx.zeta;
lambda = approx.lambda;

kappa_hat = zeros(N*J,N,TIME);
%% Roll up approximation points
eqm_belief    = v2struct(Ldyn, wf00, pf00, X, VALjn00);
approx_belief = v2struct(mu, pi, varrho, chi, zeta, lambda);
eqm_pf        = v2struct(Ldyn, wf00, pf00, X, VALjn00);
approx_pf     = v2struct(mu, pi, varrho, chi, zeta, lambda);

L_belief  = zeros(J,R,TIME);
w_belief  = zeros(J,N,TIME);
X_belief  = zeros(J,N,TIME);
p_belief  = zeros(J,N,TIME);
%v_belief  = zeros(R*J,TIME);
pi_belief = zeros(N*J,N,TIME);
mu_belief = zeros(R*J,R*J,TIME);

L_pf = zeros(J,R,TIME);
w_pf = zeros(J,N,TIME);
X_pf = zeros(J,N,TIME);
p_pf = zeros(J,N,TIME);
%v_pf = zeros(R*J,TIME);
pi_pf = zeros(N*J,N,TIME);
mu_pf = zeros(R*J,R*J,TIME);

for t2 = ENDT:-1:1 %Perfect Foresight begins from ENDT+1
    %% Recover BELIEF at period ENDT (recursive)
    disp('Current T:')
    disp(t2)
    disp('########################################')
    disp('Recover BELIEF')
    disp('########################################')

    L = zeros(R*J,1); %initial deviation is zero
    
    t1 = t2+1;
    V = initial.v_hat(:,:,t2); %initial value for v_hat
    W = initial.w_hat(:,:,:,t2); %initial value for w_hat

%    if t2==ENDT
%        V = zeros(R*J,TIME); %initial value for v_hat
%        V = initial.v_hat(:,:,ENDT);
%        W = initial.w_hat(:,:,:,ENDT);
%    else 
%        V = zeros(R*J,TIME);
%        W = zeros(J,N,TIME);
%    end

    % approximation points are coming from PF
    [eqm_temp_belief] = PBP_DYN(params, t1, t2, E_T_hat, kappa_hat, L, V, W, approx_pf);
    
    L_belief(:,:,:,t1) = reshape(eqm_temp_belief.L,J,R,TIME);
    w_belief(:,:,:,t1) = eqm_temp_belief.w;
    p_belief(:,:,:,t1) = eqm_temp_belief.p;
    X_belief(:,:,:,t1) = eqm_temp_belief.X;
    pi_belief(:,:,:,t1) = eqm_temp_belief.pi;
    mu_belief(:,:,:,t1) = eqm_temp_belief.mu;

    % Append New approximation points from DGP
    % update
    for t=t1:TIME
        pi(:,:,t)      =  approx_pf.pi(:,:,t)   .* exp(pi_belief(:,:,t,t1)); 
        mu(:,:,t)      =  approx_pf.mu(:,:,t)   .* exp(mu_belief(:,:,t,t1));
        X(:,:,t)       =  eqm_pf.X(:,:,t)       .* exp(X_belief(:,:,t,t1));
        Ldyn(:,:,t)    =  eqm_pf.Ldyn(:,:,t)    .* exp(L_belief(:,:,t,t1));
        wf00(:,:,t)    =  eqm_pf.wf00(:,:,t)    .* exp(w_belief(:,:,t,t1));
        pf00(:,:,t)    =  eqm_pf.pf00(:,:,t)    .* exp(p_belief(:,:,t,t1));
        VALjn00(:,1:R,t) =  eqm_pf.VALjn00(:,1:R,t) .* exp(w_belief(:,1:R,t,t1) + L_belief(:,1:R,t,t1));
        VALjn00(:,R+1:N,t) =  eqm_pf.VALjn00(:,R+1:N,t) .* exp(w_belief(:,R+1:N,t,t1));
        
        %normalize
        %{
        Ldyn(:,:,t)     =    Ldyn(:,:,t)./ sum(sum(Ldyn(:,:,t)));
        for i=1:N*J
            pi(i,:,t)     =    pi(i,:,t)./sum(pi(i,:,t));
        end
        for i=1:R*J
            mu(i,:,t)     =    mu(i,:,t)./sum(mu(i,:,t));
        end
        %}
        % compute approximation points (in Fernando's note)
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
        for n=1:N
             for j=1:J
                 for ii=1:N
                     varrho(n+(j-1)*N,ii,t) = pi(n+(j-1)*N,ii,t) * X(j,n,t) / X(j,ii,t);
                     chi(n+(j-1)*N,ii,t) = pi(n+(j-1)*N,ii,t) * X(j,n,t) / VALjn00(j,ii,t);
                 end
             end
         end
         for k=1:J
             for ii=1:N
                 for j=1:J
                     zeta(ii+(k-1)*N,j,t) = VALjn00(k,ii,t)./ X(j,ii,t);
                 end
             end
         end
     end
    %% save the updated values
    approx_belief.mu(:,:,t1:TIME)     = mu(:,:,t1:TIME);
    approx_belief.lambda(:,:,t1:TIME) = lambda(:,:,t1:TIME);
    approx_belief.pi(:,:,t1:TIME)     = pi(:,:,t1:TIME);
    approx_belief.varrho(:,:,t1:TIME) = varrho(:,:,t1:TIME);
    approx_belief.chi(:,:,t1:TIME)    = chi(:,:,t1:TIME);
    approx_belief.zeta(:,:,t1:TIME)   = zeta(:,:,t1:TIME);
    eqm_belief.wf00(:,:,t1:TIME)      = wf00(:,:,t1:TIME);
    eqm_belief.pf00(:,:,t1:TIME)      = pf00(:,:,t1:TIME);
    eqm_belief.Ldyn(:,:,t1:TIME)      = Ldyn(:,:,t1:TIME);
    eqm_belief.VALjn00(:,:,t1:TIME)   = VALjn00(:,:,t1:TIME);
    
    L_belief_lev(:,:,:,t2) = eqm_belief.Ldyn; 
    wf_belief(:,:,:,t2) = eqm_belief.wf00; 
    pf_belief(:,:,:,t2) = eqm_belief.pf00; 

    %% Recover PERFECT FORESIGHT
    disp('########################################')
    disp('Recover PF')
    disp('########################################')
    L = zeros(R*J,1);
    if t2==ENDT
        %V = zeros(R*J,TIME); %initial value for v_hat
        V = initial.v_hat(:,:,ENDT+1);
        W = initial.w_hat(:,:,:,ENDT+1);
    else
        V = zeros(R*J,TIME); %initial value for v_tildehat
        W = zeros(J,N,TIME);
    end
    % input for belief here is pefect foresight (E_T_hat_pf)
    % approximation points come from the belief path generated right before
    [eqm_temp_pf] = PBP_DYN(params, t2, t2, E_T_hat_pf, kappa_hat, L, V, W, approx_belief);
    
    % save deviation path for each period(t2)
    L_pf(:,:,:,t2) = reshape(eqm_temp_pf.L,J,R,TIME);
    w_pf(:,:,:,t2) = eqm_temp_pf.w;
    p_pf(:,:,:,t2) = eqm_temp_pf.p;
    X_pf(:,:,:,t2) = eqm_temp_pf.X;
    pi_pf(:,:,:,t2) = eqm_temp_pf.pi;
    mu_pf(:,:,:,t2) = eqm_temp_pf.mu;

    % Append New approximation points from DGP (in level)
    for t=t2:TIME
        pi(:,:,t)      =  approx_belief.pi(:,:,t)  .* exp(pi_pf(:,:,t,t2)); 
        mu(:,:,t)      =  approx_belief.mu(:,:,t)  .* exp(mu_pf(:,:,t,t2));
        X(:,:,t)       =  eqm_belief.X(:,:,t)      .* exp(X_pf(:,:,t,t2));
        Ldyn(:,:,t)    =  eqm_belief.Ldyn(:,:,t)   .* exp(L_pf(:,:,t,t2));
        wf00(:,:,t)    =  eqm_belief.wf00(:,:,t)   .* exp(w_pf(:,:,t,t2));
        pf00(:,:,t)    =  eqm_belief.pf00(:,:,t)   .* exp(p_pf(:,:,t,t2));
        VALjn00(:,1:R,t) =  eqm_belief.VALjn00(:,1:R,t)     .* exp(w_pf(:,1:R,t,t2) + L_pf(:,1:R,t,t2));
        VALjn00(:,R+1:N,t) =  eqm_belief.VALjn00(:,R+1:N,t) .* exp(w_pf(:,R+1:N,t,t2));

        %normalize
        %{
        Ldyn(:,:,t)     =    Ldyn(:,:,t)./ sum(sum(Ldyn(:,:,t)));
        for i=1:N*J
            pi(i,:,t)     =    pi(i,:,t)./sum(pi(i,:,t));
        end
        for i=1:R*J
            mu(i,:,t)     =    mu(i,:,t)./sum(mu(i,:,t));
        end
        %}
        % compute approximation points (in Fernando's note)
        for k=1:J
            for ii=1:R
                for j=1:J
                    for n=1:R
%                        lambda(k+(ii-1)*J,j+(n-1)*J,t) = mu(k+(ii-1)*J,j+(n-1)*J,t)*Ldyn(k,ii,t)/Ldyn(j,n,t-1);
                         if t==1
                            lambda(k+(ii-1)*J,j+(n-1)*J,1) = mu(k+(ii-1)*J,j+(n-1)*J,1)*Ldyn(k,ii,1)/Ldyn(j,n,1);
                         else
                            lambda(k+(ii-1)*J,j+(n-1)*J,t) = mu(k+(ii-1)*J,j+(n-1)*J,t-1)*Ldyn(k,ii,t-1)/Ldyn(j,n,t);
                         end
                    end
                end
            end
        end
        for n=1:N
            for j=1:J
                for ii=1:N
                    varrho(n+(j-1)*N,ii,t) = pi(n+(j-1)*N,ii,t) * X(j,n,t) / X(j,ii,t);
                    chi(n+(j-1)*N,ii,t) = pi(n+(j-1)*N,ii,t) * X(j,n,t) / VALjn00(j,ii,t);
                end
            end
        end
        for k=1:J
            for ii=1:N
                for j=1:J
                    zeta(ii+(k-1)*N,j,t) = VALjn00(k,ii,t)./ X(j,ii,t);
                end
            end
        end
    end
  
    %% save the updated value
    approx_pf.mu(:,:,t2:TIME)     = mu(:,:,t2:TIME);
    approx_pf.lambda(:,:,t2:TIME) = lambda(:,:,t2:TIME);
    approx_pf.pi(:,:,t2:TIME)     = pi(:,:,t2:TIME);
    approx_pf.varrho(:,:,t2:TIME) = varrho(:,:,t2:TIME);
    approx_pf.chi(:,:,t2:TIME)    = chi(:,:,t2:TIME);
    approx_pf.zeta(:,:,t2:TIME)   = zeta(:,:,t2:TIME);
    eqm_pf.wf00(:,:,t2:TIME)      = wf00(:,:,t2:TIME);
    eqm_pf.pf00(:,:,t2:TIME)      = pf00(:,:,t2:TIME);
    eqm_pf.Ldyn(:,:,t2:TIME)      = Ldyn(:,:,t2:TIME);
    eqm_pf.VALjn00(:,:,t2:TIME)   = VALjn00(:,:,t2:TIME);
    
    
    L_pf_lev(:,:,:,t2) = eqm_pf.Ldyn; 
    wf_pf(:,:,:,t2) = eqm_pf.wf00; 
    pf_pf(:,:,:,t2) = eqm_pf.pf00; 
end

eqm_recur = v2struct(L_pf_lev, L_belief_lev, wf_pf, pf_pf, L_pf, w_pf, p_pf, pi_pf, mu_pf,...
     wf_belief, pf_belief, L_belief, w_belief, p_belief, pi_belief, mu_belief);

save('DATA/RECURSIVE.mat','eqm_recur')
end
