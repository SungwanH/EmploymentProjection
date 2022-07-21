function [SO_temp, SO_dyn] = SO_TERMS(params, t1, eqm_sim, eqm, approx)
%% Roll down parameters
v2struct(params.envr);
v2struct(params.tech);
v2struct(params.modl);

%% Roll down approximation points
w_bar   =   eqm.w_lev;
P_bar   =   eqm.P_lev;
rw_bar  =   reshape(w_bar(:,1:R,:)./P_bar(:,1:R,:),[R*J,TIME]); 
mu      =   approx.mu;
lambda  =   approx.lambda;
varrho  =   approx.varrho;
zeta    =   approx.zeta;
chi     =   approx.chi;
pi      =   approx.pi;

[mom_dyn, mom_temp] = SO_MOMENTS(params, t1, eqm_sim);

v2struct(mom_dyn);
v2struct(mom_temp);
%% SO terms for Dynamic eqm 
SO_L = zeros(R*J,TIME); SO_v = zeros(R*J,TIME);
for t = t1:TIME-1
    for k=1:R*J
        for g=1:R*J
            for m=1:R*J
                L_LL(g,k,m,t)  = ((k==m)*lambda(k,g,t+1) - ( lambda(k,g,t+1)*lambda(m,g,t+1) )) * ll(k,m,t);
                L_vL(g,k,m,t)  = ((BETA/NU) * (lambda(m,g,t+1)*((k==g)-mu(m,k,t))) - (BETA/NU) * ( lambda(m,g,t+1)*( lambda(:,g,t+1)'*((k==g)-mu(:,k,t)) ) )) * vl(k,m,t);
                L_vv1(g,k,m,t) = (BETA/NU)^2 * (-(mu(:,k,t) .* lambda(:,g,t+1))' *((m==g) + (m==k) - 2*mu(:,m,t)));
                L_vv2(g,k,m,t) = (BETA/NU)^2 * ((k==g) * lambda(:,g,t+1)'*((m==g)-mu(:,m,t)));
                L_vv3(g,k,m,t) = -(BETA/NU)^2 * (lambda(:,g,t+1)'*((k==g)-mu(:,k,t)))*(lambda(:,g,t+1)'*((m==g)-mu(:,m,t)));
                L_vv(g,k,m,t)  = (L_vv1(g,k,m,t) + L_vv2(g,k,m,t) + L_vv3(g,k,m,t)) * vv(k,m,t);
            end
        end
    end
    SO_L(:,t+1) = (1/2) * sum(sum(L_LL(:,:,:,t),3),2) + (1/2) * sum(sum(L_vv(:,:,:,t),3),2) + sum(sum(L_vL(:,:,:,t),3),2);
end

for t=TIME-1:-1:t1 
    mu_aux = mu(:,:,t);
    for m=1:R*J
        for k=1:R*J
            for n=1:R*J
                v_vv(n,k,m) = (BETA^2 / NU) * mu_aux(n,m) * ((k==m) - mu_aux(n,k)) * vv(m,k,t);
            end
        end
    end
    SO_v(:,t) = (1/2) * sum(sum(v_vv,3),2); % This is R*J by 1 
end

v_ww(:,t1:TIME)  = (1-CRRA_PAR) .* rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* ww_dyn(:,t1:TIME);
v_PP(:,t1:TIME)  = (1-CRRA_PAR) .* rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* PP_dyn(:,t1:TIME);
v_wP(:,t1:TIME)  = -(1-CRRA_PAR) .* rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* wP_dyn(:,t1:TIME);
SO_rw = (1/2) * v_ww + (1/2) * v_PP + v_wP; % This is J by R 

%for t=t1:TIME
 %   for n=1:R
 %       for k=1:J
%            v_ww(k,n) = (w_bar(k,n,t)^(-CRRA_PAR-1)/P_bar(1,n,t)^(1-CRRA_PAR)) * ww_dyn(k,n,t);
%            v_PP(k,n) = (w_bar(k,n,t)^(-CRRA_PAR-1)/P_bar(1,n,t)^(3-CRRA_PAR)) * PP_dyn(k,n,t);
%            v_wP(k,n) = -(1-CRRA_PAR) * (w_bar(k,n,t)^(-CRRA_PAR)/P_bar(1,n,t)^(2-CRRA_PAR)) * wP_dyn(k,n,t);
%        end
%    end
%    SO_rw(:,t) = reshape((1/2) * v_ww + (1/2) * v_PP + v_wP, [R*J,1]); % This is J by R 
%end

%% SO terms for Temporary eqm 
SO_p = zeros(J,N,TIME); SO_X = zeros(J,N,TIME); SO_w = zeros(J,N,TIME);
for t = t1:TIME
    for j=1:J
        for n=1:N
            for k=1:N        
                for m=1:N
                    %diff w.r.t w and w
                    P_ww(j,n,k,m) = ((m==k) *  GAMMA(j,k) * GAMMA(j,k) *(-THETA(j)) * pi(n+(j-1)*N,k,t) -  GAMMA(j,k) * GAMMA(j,m) * (-THETA(j)) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t)) * ww_price(j,k,m,t);
                   %diff w.r.t p and p
                    P_pp(j,n,k,m) = ((m==k) * (-THETA(j)) * (1-GAMMA(j,k)) * (1-GAMMA(j,k)) * pi(n+(j-1)*N,k,t) - (-THETA(j)) * (1-GAMMA(j,k)) * (1-GAMMA(j,m)) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t)) * pp(j,k,m,t);
                    %diff w.r.t T and T
                    P_TT(j,n,k,m) = ((m==k) *(-1/THETA(j)) * pi(n+(j-1)*N,k,t) - (-1/THETA(j)) * pi(n+(j-1)*N,k,t)* pi(n+(j-1)*N,m,t)) * TT(j,k,m,t);
                    %diff w.r.t w and T
                    P_wT(j,n,k,m) = ((m==k) * GAMMA(j,k) * pi(n+(j-1)*N,k,t) - GAMMA(j,k) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t)) * wT(j,k,m,t);
                    %diff w.r.t p and T
                    P_pT(j,n,k,m) = ((m==k) * (1-GAMMA(j,k)) * pi(n+(j-1)*N,k,t) - (1-GAMMA(j,k)) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t) ) *pT(j,k,m,t);
                    %diff w.r.t w and P
                    P_wp(j,n,k,m) = ((m==k) * GAMMA(j,k) * (1-GAMMA(j,m))* (-THETA(j)) * pi(n+(j-1)*N,m,t) - (GAMMA(j,k)) * (1-GAMMA(j,m)) * (-THETA(j)) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t)) * wp(j,k,m,t);
                end
            end
        end
    end
SO_p(:,:,t) = (1/2) * sum(sum(P_ww+P_pp+P_TT,4),3) + sum(sum(P_wT+P_pT+P_wp,4),3);
%        p_hat(:,:,t) = sum(P_w+P_p+P_T,3); % first-order approximation
end  

% Step 4d. Solve for total expenditure
for t = t1:TIME
    for j=1:J
        for i=1:N
            for n=1:N
                for m=1:N
                    %diff w.r.t. X and X
                    X_XX(j,i,n,m)   = ((n==m) * (1-GAMMA(j,i)) * varrho(m+(j-1)*N,i,t) - ((1-GAMMA(j,i))^2) * varrho(n+(j-1)*N,i,t)* varrho(m+(j-1)*N,i,t)) * XX(j,n,m,t);
                    %diff w.r.t. pi and pi
                    X_pipi(j,i,n,m) = ((n==m) * (1-GAMMA(j,i)) * varrho(m+(j-1)*N,i,t) - ((1-GAMMA(j,i))^2) * varrho(n+(j-1)*N,i,t)* varrho(m+(j-1)*N,i,t)) * pipi(n+(j-1)*N,m+(i-1)*N,t);
                    %diff w.r.t. pi and X
                    X_piX(j,i,n,m)  = ((n==m) * (1-GAMMA(j,i)) * varrho(m+(j-1)*N,i,t) - ((1-GAMMA(j,i))^2) * varrho(n+(j-1)*N,i,t)* varrho(m+(j-1)*N,i,t)) * piX(n+(j-1)*N,m+(i-1)*N,t);
                end
                for s=1:J
                    X_piw(j,i,n,s) = -(1-GAMMA(j,i)) * ALPHAS(j,i) * varrho(n+(j-1)*N,i,t) * zeta(i+(s-1)*N,j,t) * piw(n+(j-1)*N,i+(s-1)*N,t);
                    X_piL(j,i,n,s) = -(1-GAMMA(j,i)) * ALPHAS(j,i) * varrho(n+(j-1)*N,i,t) * zeta(i+(s-1)*N,j,t) * pil_ca(n+(j-1)*N,i+(s-1)*N,t);
                    X_Xw(j,i,n,s)  = -(1-GAMMA(j,i)) * ALPHAS(j,i) * varrho(n+(j-1)*N,i,t) * zeta(i+(s-1)*N,j,t) * Xw(j,n,s,i,t);
                    X_XL(j,i,n,s)  = -(1-GAMMA(j,i)) * ALPHAS(j,i) * varrho(n+(j-1)*N,i,t) * zeta(i+(s-1)*N,j,t) * Xl(j,n,s,i,t);
                end
            end
            for k=1:J
                for s=1:J
                    %diff w.r.t. w and w
                    X_ww(j,i,k,s) = ((k==s) * ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) - (ALPHAS(j,i)^2) * zeta(i+(k-1)*N,j,t) * zeta(i+(s-1)*N,j,t)) * ww_ca(k,i,s,t);
                    %diff w.r.t. L and L
                    X_LL(j,i,k,s) = ((k==s) * ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) - (ALPHAS(j,i)^2) * zeta(i+(k-1)*N,j,t) * zeta(i+(s-1)*N,j,t)) * ll_ca(k,i,s,t);
                    %diff w.r.t. w and L
                    X_wL(j,i,k,s) = ((k==s) * ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) - (ALPHAS(j,i)^2) * zeta(i+(k-1)*N,j,t) * zeta(i+(s-1)*N,j,t)) * wl_ca(k,i,s,t);
                end
            end
        end
    end
    SO_X(:,:,t) = (1/2) * (sum(sum(X_XX+X_pipi,4),3) + sum(sum(X_ww+X_LL,4),3)) + sum(sum(X_piX,4),3) + sum(sum(X_piw+X_piL+X_Xw+X_XL,4),3) +  sum(sum(X_wL,4),3);
%         X_hat(:,:,t) = sum(X_pi+X_X,3) + sum(X_w+X_L,3); %first-order approx
end  
% Step 4e. Solve for an updated w_hat using the labor market clearing

for t = t1:TIME
    for n=1:N
        for i=1:N
            for j=1:J
                for m=1:N
                    w_pipi(j,i,n,m) =  ((m==n) * GAMMA(j,i) * chi(n+(j-1)*N,i,t) - (GAMMA(j,i)^2) * chi(n+(j-1)*N,i,t) * chi(m+(j-1)*N,i,t)) * pipi(n+(j-1)*N,m+(i-1)*N,t);
                    w_XX(j,i,n,m)   =  ((m==n) * GAMMA(j,i) * chi(n+(j-1)*N,i,t) - (GAMMA(j,i)^2) * chi(n+(j-1)*N,i,t) * chi(m+(j-1)*N,i,t)) * XX(j,n,m,t);
                    w_piX(j,i,n,m)  =  ((m==n) * GAMMA(j,i) * chi(n+(j-1)*N,i,t) - (GAMMA(j,i)^2) * chi(n+(j-1)*N,i,t) * chi(m+(j-1)*N,i,t)) * piX(n+(j-1)*N,m+(i-1)*N,t);
                end
            end
        end
    end
    SO_w(:,:,t) = (1/2) * sum(sum(w_pipi+w_XX,4),3) + +sum(sum(w_piX,4),3);
end

SO_dyn  = v2struct(SO_L, SO_v, SO_rw);
SO_temp = v2struct(SO_p, SO_X, SO_w);
end