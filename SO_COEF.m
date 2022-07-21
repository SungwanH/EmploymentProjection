function [SO_temp, SO_dyn] = SO_COEF(params, eqm, approx)
% This function computes coefficients of second order terms 
% using the approximation points obtained from ex-ante eqm outcomes

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

%% SO terms for Dynamic eqm 
SO_L = zeros(R*J,TIME); SO_v = zeros(R*J,TIME);
for t = 1:TIME-1
    for k=1:R*J
        for g=1:R*J
            for m=1:R*J
                L_LL(g,k,m,t)  = ((k==m)*lambda(k,g,t+1) - ( lambda(k,g,t+1)*lambda(m,g,t+1) ));
                L_vL(g,k,m,t)  = ((BETA/NU) * (lambda(m,g,t+1)*((k==g)-mu(m,k,t))) - (BETA/NU) * ( lambda(m,g,t+1)*( lambda(:,g,t+1)'*((k==g)-mu(:,k,t)) ) ));
                L_vv1(g,k,m,t) = (BETA/NU)^2 * (-(mu(:,k,t) .* lambda(:,g,t+1))' *((m==g) + (m==k) - 2*mu(:,m,t)));
                L_vv2(g,k,m,t) = (BETA/NU)^2 * ((k==g) * lambda(:,g,t+1)'*((m==g)-mu(:,m,t)));
                L_vv3(g,k,m,t) = -(BETA/NU)^2 * (lambda(:,g,t+1)'*((k==g)-mu(:,k,t)))*(lambda(:,g,t+1)'*((m==g)-mu(:,m,t)));
                L_vv(g,k,m,t)  = L_vv1(g,k,m,t) + L_vv2(g,k,m,t) + L_vv3(g,k,m,t);
            end
        end
    end
end

for t=1:TIME
    mu_aux = mu(:,:,t);
    for m=1:R*J
        for k=1:R*J
            for n=1:R*J
                v_vv(n,k,m,t) = (BETA^2 / NU) * mu_aux(n,m) * ((k==m) - mu_aux(n,k)) ;
            end
        end
    end
end

v_ww(:,1:TIME)  = (1-CRRA_PAR) .* rw_bar(:,1:TIME).^(1-CRRA_PAR);
v_PP(:,1:TIME)  = (1-CRRA_PAR) .* rw_bar(:,1:TIME).^(1-CRRA_PAR);
v_wP(:,1:TIME)  = -(1-CRRA_PAR) .* rw_bar(:,1:TIME).^(1-CRRA_PAR);

%% SO terms for Temporary eqm 
for t = 1:TIME
    for j=1:J
        for n=1:N
            for k=1:N        
                for m=1:N
                    %diff w.r.t w and w
                    P_ww(j,n,k,m,t) = ((m==k) *  GAMMA(j,k) * GAMMA(j,k) *(-THETA(j)) * pi(n+(j-1)*N,k,t) -  GAMMA(j,k) * GAMMA(j,m) * (-THETA(j)) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t));
                   %diff w.r.t p and p
                    P_pp(j,n,k,m,t) = ((m==k) * (-THETA(j)) * (1-GAMMA(j,k)) * (1-GAMMA(j,k)) * pi(n+(j-1)*N,k,t) - (-THETA(j)) * (1-GAMMA(j,k)) * (1-GAMMA(j,m)) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t));
                    %diff w.r.t T and T
                    P_TT(j,n,k,m,t) = ((m==k) *(-1/THETA(j)) * pi(n+(j-1)*N,k,t) - (-1/THETA(j)) * pi(n+(j-1)*N,k,t)* pi(n+(j-1)*N,m,t));
                    %diff w.r.t w and T
                    P_wT(j,n,k,m,t) = ((m==k) * GAMMA(j,k) * pi(n+(j-1)*N,k,t) - GAMMA(j,k) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t));
                    %diff w.r.t p and T
                    P_pT(j,n,k,m,t) = ((m==k) * (1-GAMMA(j,k)) * pi(n+(j-1)*N,k,t) - (1-GAMMA(j,k)) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t) );
                    %diff w.r.t w and P
                    P_wp(j,n,k,m,t) = ((m==k) * GAMMA(j,k) * (1-GAMMA(j,m))* (-THETA(j)) * pi(n+(j-1)*N,m,t) - (GAMMA(j,k)) * (1-GAMMA(j,m)) * (-THETA(j)) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t));
                end
            end
        end
    end
end  

% Step 4d. Solve for total expenditure
for t = 1:TIME
    for j=1:J
        for i=1:N
            for n=1:N
                for m=1:N
                    %diff w.r.t. X and X
                    X_XX(j,i,n,m,t)   = ((n==m) * (1-GAMMA(j,i)) * varrho(m+(j-1)*N,i,t) - ((1-GAMMA(j,i))^2) * varrho(n+(j-1)*N,i,t)* varrho(m+(j-1)*N,i,t));
                    %diff w.r.t. pi and pi
                    X_pipi(j,i,n,m,t) = ((n==m) * (1-GAMMA(j,i)) * varrho(m+(j-1)*N,i,t) - ((1-GAMMA(j,i))^2) * varrho(n+(j-1)*N,i,t)* varrho(m+(j-1)*N,i,t));
                    %diff w.r.t. pi and X
                    X_piX(j,i,n,m,t)  = ((n==m) * (1-GAMMA(j,i)) * varrho(m+(j-1)*N,i,t) - ((1-GAMMA(j,i))^2) * varrho(n+(j-1)*N,i,t)* varrho(m+(j-1)*N,i,t));
                end
                for s=1:J
                    X_piw(j,i,n,s,t) = -(1-GAMMA(j,i)) * ALPHAS(j,i) * varrho(n+(j-1)*N,i,t) * zeta(i+(s-1)*N,j,t);
                    X_piL(j,i,n,s,t) = -(1-GAMMA(j,i)) * ALPHAS(j,i) * varrho(n+(j-1)*N,i,t) * zeta(i+(s-1)*N,j,t);
                    X_Xw(j,i,n,s,t)  = -(1-GAMMA(j,i)) * ALPHAS(j,i) * varrho(n+(j-1)*N,i,t) * zeta(i+(s-1)*N,j,t);
                    X_XL(j,i,n,s,t)  = -(1-GAMMA(j,i)) * ALPHAS(j,i) * varrho(n+(j-1)*N,i,t) * zeta(i+(s-1)*N,j,t);
                end
            end
            for k=1:J
                for s=1:J
                    %diff w.r.t. w and w
                    X_ww(j,i,k,s,t) = ((k==s) * ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) - (ALPHAS(j,i)^2) * zeta(i+(k-1)*N,j,t) * zeta(i+(s-1)*N,j,t));
                    %diff w.r.t. L and L
                    X_LL(j,i,k,s,t) = ((k==s) * ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) - (ALPHAS(j,i)^2) * zeta(i+(k-1)*N,j,t) * zeta(i+(s-1)*N,j,t));
                    %diff w.r.t. w and L
                    X_wL(j,i,k,s,t) = ((k==s) * ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) - (ALPHAS(j,i)^2) * zeta(i+(k-1)*N,j,t) * zeta(i+(s-1)*N,j,t));
                end
            end
        end
    end
end  
% Step 4e. Solve for an updated w_hat using the labor market clearing

for t = 1:TIME
    for n=1:N
        for i=1:N
            for j=1:J
                for m=1:N
                    w_pipi(j,i,n,m,t) =  ((m==n) * GAMMA(j,i) * chi(n+(j-1)*N,i,t) - (GAMMA(j,i)^2) * chi(n+(j-1)*N,i,t) * chi(m+(j-1)*N,i,t));
                    w_XX(j,i,n,m,t)   =  ((m==n) * GAMMA(j,i) * chi(n+(j-1)*N,i,t) - (GAMMA(j,i)^2) * chi(n+(j-1)*N,i,t) * chi(m+(j-1)*N,i,t));
                    w_piX(j,i,n,m,t)  =  ((m==n) * GAMMA(j,i) * chi(n+(j-1)*N,i,t) - (GAMMA(j,i)^2) * chi(n+(j-1)*N,i,t) * chi(m+(j-1)*N,i,t));
                end
            end
        end
    end
end


SO_dyn  = v2struct(L_LL, L_vL, L_vv, v_vv, v_ww, v_PP, v_wP);
SO_temp = v2struct(P_ww, P_pp, P_TT, P_wT, P_pT, P_wp, X_XX, X_pipi, X_piX, X_piw, X_piL, X_Xw, X_XL, X_ww, X_LL, X_wL, w_pipi, w_XX, w_piX);
end