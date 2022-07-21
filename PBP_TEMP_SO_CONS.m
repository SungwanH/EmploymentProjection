function [w_hat, p_hat, P_hat, pi_hat, X_hat] = PBP_TEMP_SO_CONS(params, t1, T_hat, kappa_hat, L_hat_R, approx, eqm_FO)

%%Temporary problem (=Trade equilibrium) with second order approximation%%
% For second order terms, use constant terms for variables
%% Roll down parameters
v2struct(params.envr);
v2struct(params.tech);
v2struct(params.modl);

%% Roll down approximation points
varrho  =   approx.varrho;
zeta    =   approx.zeta;
chi     =   approx.chi;
pi      =   approx.pi;

%w_hat_FO = eqm_FO.w_hat;
%p_hat_FO = eqm_FO.p_hat;
%X_hat_FO = eqm_FO.X_hat;
%L_hat_FO = eqm_FO.L_hat;
w_hat_FO = eqm_FO.w;
p_hat_FO = eqm_FO.p;
X_hat_FO = eqm_FO.X;
L_hat_FO = NaN(J,N,TIME);
for t=1:TIME
    L_hat_FO(:,:,t) = reshape(eqm_FO.L(:,t),J,N);
end
L_hat_FO(1:J,(R+1):N,1:TIME) = zeros(J,N-R,TIME);
%L_hat_R only includes US states
L_hat = NaN(J,N,TIME);
for t=1:TIME
    L_hat(:,1:R,t)=reshape(L_hat_R(:,t),J,R);
end
L_hat(1:J,(R+1):N,1:TIME) = zeros(J,N-R,TIME); %There is no migration in other countries

ITER_TEMP = 1;
wmax=1;
% Step 4a. Initial guess for the w_tilde_hat
%w_hat = zeros(J,N,TIME);
w_hat = zeros(J,N,TIME);
w_update = zeros(J,N,TIME);
p_hat = zeros(J,N,TIME);
p_old = zeros(J,N,TIME);
pi_hat = zeros(N*J,N,TIME);
X_hat = zeros(J,N,TIME);
P_w = zeros(J,N,N); P_p = zeros(J,N,N); P_T = zeros(J,N,N); 
P_ww = zeros(J,N,N,N); P_pp = zeros(J,N,N,N); P_TT = zeros(J,N,N,N); 
P_wT = zeros(J,N,N,N); P_pT = zeros(J,N,N,N); P_wp = zeros(J,N,N,N); 
%P_ww1 = zeros(J,N,N); P_pp1 = zeros(J,N,N); P_TT1 = zeros(J,N,N); 
%P_wT1 = zeros(J,N,N); P_pT1 = zeros(J,N,N); P_wp1 = zeros(J,N,N); 
%P_ww2 = zeros(J,N,N,N); P_pp2 = zeros(J,N,N,N); P_TT2 = zeros(J,N,N,N); 
%P_wT2 = zeros(J,N,N,N); P_pT2 = zeros(J,N,N,N); P_wp2 = zeros(J,N,N,N); 

X_w = zeros(J,N,J); X_L = zeros(J,N,J);  X_pi = zeros(J,N,N); X_X = zeros(J,N,N);
X_pipi = zeros(J,N,N); X_XX = zeros(J,N,N,N); X_piX = zeros(J,N,N,N);
X_ww = zeros(J,N,J,J); X_LL=zeros(J,N,J,J); X_wL=zeros(J,N,J,J); 
X_piw = zeros(J,N,N,J); X_piL = zeros(J,N,N,J); X_Xw = zeros(J,N,N,J);X_XL = zeros(J,N,N,J);

w_pi = zeros(J,N,N); w_X = zeros(J,N,N); w_L = zeros(J,N); 
w_pipi = zeros(J,N,N,N); w_XX = zeros(J,N,N,N); w_piX = zeros(J,N,N,N);


while (ITER_TEMP <= MAXIT) && (wmax > TOLTEMP)
    
    % Step 4b,c. solve for p_hat and pi_hat
   for t = t1:TIME
        for j=1:J
            for n=1:N
                for k=1:N
                    %diff w.r.t. w
                    P_w(j,n,k) = GAMMA(j,k) * pi(n+(j-1)*N,k,t) * w_hat(j,k,t);
                    %diff w.r.t. p
                    P_p(j,n,k) = (1-GAMMA(j,k)) * pi(n+(j-1)*N,k,t) * p_old(j,k,t);
                    %diff w.r.t. T
                    P_T(j,n,k) = -(1/THETA(j)) * pi(n+(j-1)*N,k,t) * T_hat(j,k,t);                   
                    for m=1:N
                        %diff w.r.t w and w
                        P_ww(j,n,k,m) = ((m==k) *  GAMMA(j,k) * GAMMA(j,k) *(-THETA(j)) * pi(n+(j-1)*N,k,t) -  GAMMA(j,k) * GAMMA(j,m) * (-THETA(j)) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t)) * w_hat_FO(j,k,t) * w_hat_FO(j,m,t);
                       %diff w.r.t p and p
                        P_pp(j,n,k,m) = ((m==k) * (-THETA(j)) * (1-GAMMA(j,k)) * (1-GAMMA(j,k)) * pi(n+(j-1)*N,k,t) - (-THETA(j)) * (1-GAMMA(j,k)) * (1-GAMMA(j,m)) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t)) * p_hat_FO(j,k,t) * p_hat_FO(j,m,t);
                        %diff w.r.t T and T
                        P_TT(j,n,k,m) = ((m==k) *(-1/THETA(j)) * pi(n+(j-1)*N,k,t) - (-1/THETA(j)) * pi(n+(j-1)*N,k,t)* pi(n+(j-1)*N,m,t)) * T_hat(j,k,t) * T_hat(j,m,t);
                        %diff w.r.t w and T
                        P_wT(j,n,k,m) = ((m==k) * GAMMA(j,k) * pi(n+(j-1)*N,k,t) - GAMMA(j,k) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t)) * w_hat_FO(j,k,t) * T_hat(j,m,t);
                        %diff w.r.t p and T
                        P_pT(j,n,k,m) = ((m==k) * (1-GAMMA(j,k)) * pi(n+(j-1)*N,k,t) - (1-GAMMA(j,k)) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t) ) * p_hat_FO(j,k,t) * T_hat(j,m,t);
                        %diff w.r.t w and P
                        P_wp(j,n,k,m) = ((m==k) * GAMMA(j,k) * (1-GAMMA(j,m))* (-THETA(j)) * pi(n+(j-1)*N,m,t) - (GAMMA(j,k)) * (1-GAMMA(j,m)) * (-THETA(j)) * pi(n+(j-1)*N,k,t) * pi(n+(j-1)*N,m,t)) * w_hat_FO(j,k,t) * p_hat_FO(j,m,t);
                    end
                end
            end
        end
    p_hat(:,:,t) = sum(P_w+P_p+P_T,3)  + (1/2) * sum(sum(P_ww+P_pp+P_TT,4),3) + sum(sum(P_wT+P_pT+P_wp,4),3);
%        p_hat(:,:,t) = sum(P_w+P_p+P_T,3); % first-order approximation
    end  
    %price index
    for t=t1:TIME
        P_hat(:,:,t) = sum((ALPHAS.*p_hat(:,:,t)),1);
    end
    % Step 4c. solve for pi_hat
    for t = t1:TIME
        for n=1:N
            for j=1:J
                pi_hat(n+(j-1)*N,:,t) = -THETA(j)*(GAMMA(j,:) .* w_hat(j,:,t) + (1-GAMMA(j,:)).*p_hat(j,:,t) - p_hat(j,n,t) + kappa_hat(n+(j-1)*N,:,t)) + T_hat(j,:,t);
             end
         end
    end        

    % Step 4d. Solve for total expenditure
    for t = t1:TIME
        for j=1:J
            for i=1:N
                for n=1:N
                    %diff w.r.t. pi
                    X_pi(j,i,n)   = (1-GAMMA(j,i)) * varrho(n+(j-1)*N,i,t)* pi_hat(n+(j-1)*N,i,t);
                    %diff w.r.t. X
                    X_X(j,i,n)    = (1-GAMMA(j,i)) * varrho(n+(j-1)*N,i,t) * X_hat(j,n,t);
                    for m=1:N
                        %diff w.r.t. X and X
                        X_XX(j,i,n,m)   = ((n==m) * (1-GAMMA(j,i)) * varrho(m+(j-1)*N,i,t) - ((1-GAMMA(j,i))^2) * varrho(n+(j-1)*N,i,t)* varrho(m+(j-1)*N,i,t)) * X_hat_FO(j,n,t) * X_hat_FO(j,m,t);
                        %diff w.r.t. pi and pi
                        X_pipi(j,i,n,m) = ((n==m) * (1-GAMMA(j,i)) * varrho(m+(j-1)*N,i,t) - ((1-GAMMA(j,i))^2) * varrho(n+(j-1)*N,i,t)* varrho(m+(j-1)*N,i,t)) * pi_hat(n+(j-1)*N,i,t) * pi_hat(m+(j-1)*N,i,t);
                        %diff w.r.t. pi and X
                        X_piX(j,i,n,m)  = ((n==m) * (1-GAMMA(j,i)) * varrho(m+(j-1)*N,i,t) - ((1-GAMMA(j,i))^2) * varrho(n+(j-1)*N,i,t)* varrho(m+(j-1)*N,i,t)) * pi_hat(n+(j-1)*N,i,t) * X_hat_FO(j,m,t);
                    end
                    for s=1:J
                        X_piw(j,i,n,s) = -(1-GAMMA(j,i)) * ALPHAS(j,i) * varrho(n+(j-1)*N,i,t) * zeta(i+(s-1)*N,j,t) * pi_hat(n+(j-1)*N,i,t) * w_hat_FO(s,i,t);
                        X_piL(j,i,n,s) = -(1-GAMMA(j,i)) * ALPHAS(j,i) * varrho(n+(j-1)*N,i,t) * zeta(i+(s-1)*N,j,t) * pi_hat(n+(j-1)*N,i,t) * L_hat_FO(s,i,t);
                        X_Xw(j,i,n,s)  = -(1-GAMMA(j,i)) * ALPHAS(j,i) * varrho(n+(j-1)*N,i,t) * zeta(i+(s-1)*N,j,t) * X_hat_FO(j,n,t) * w_hat_FO(s,i,t);
                        X_XL(j,i,n,s)  = -(1-GAMMA(j,i)) * ALPHAS(j,i) * varrho(n+(j-1)*N,i,t) * zeta(i+(s-1)*N,j,t) * X_hat_FO(j,n,t) * L_hat_FO(s,i,t);
                    end
                end
                for k=1:J
                    %diff w.r.t. w
                    X_w(j,i,k)  = ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) * w_hat(k,i,t);
                    %diff w.r.t. L
                    X_L(j,i,k)  = ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) * L_hat(k,i,t);
                    for s=1:J
                        %diff w.r.t. w and w
                        X_ww(j,i,k,s) = ((k==s) * ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) - (ALPHAS(j,i)^2) * zeta(i+(k-1)*N,j,t) * zeta(i+(s-1)*N,j,t)) * w_hat_FO(k,i,t) * w_hat_FO(s,i,t);
                        %diff w.r.t. L and L
                        X_LL(j,i,k,s) = ((k==s) * ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) - (ALPHAS(j,i)^2) * zeta(i+(k-1)*N,j,t) * zeta(i+(s-1)*N,j,t)) * L_hat_FO(k,i,t) * L_hat_FO(s,i,t);
                        %diff w.r.t. w and L
                        X_wL(j,i,k,s) = ((k==s) * ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) - (ALPHAS(j,i)^2) * zeta(i+(k-1)*N,j,t) * zeta(i+(s-1)*N,j,t)) * w_hat_FO(k,i,t) * L_hat_FO(s,i,t);
                    end
                end
            end
        end
        X_hat(:,:,t) = sum(X_pi+X_X,3) +sum(X_w+X_L,3) + (1/2) * (sum(sum(X_XX+X_pipi,4),3) + sum(sum(X_ww+X_LL,4),3)) + sum(sum(X_piX,4),3) + sum(sum(X_piw+X_piL+X_Xw+X_XL,4),3) +  sum(sum(X_wL,4),3);
%         X_hat(:,:,t) = sum(X_pi+X_X,3) + sum(X_w+X_L,3); %first-order approx
    end  
    % Step 4e. Solve for an updated w_hat using the labor market clearing

    for t = t1:TIME
        for n=1:N
            for i=1:N
                for j=1:J
                    w_pi(j,i,n)   = GAMMA(j,i) * chi(n+(j-1)*N,i,t) * pi_hat(n+(j-1)*N,i,t) ;
                    w_X(j,i,n)    = GAMMA(j,i) * chi(n+(j-1)*N,i,t) * X_hat(j,n,t) ;
                    w_L(j,i)      = -1 * L_hat(j,i,t);
                    for m=1:N
                        w_pipi(j,i,n,m) = ((m==n) * GAMMA(j,i) * chi(n+(j-1)*N,i,t) - (GAMMA(j,i)^2) * chi(n+(j-1)*N,i,t) * chi(m+(j-1)*N,i,t)) * pi_hat(n+(j-1)*N,i,t) * pi_hat(m+(j-1)*N,i,t) ;
                        w_XX(j,i,n,m)   =  ((m==n) * GAMMA(j,i) * chi(n+(j-1)*N,i,t) - (GAMMA(j,i)^2) * chi(n+(j-1)*N,i,t) * chi(m+(j-1)*N,i,t)) * X_hat_FO(j,n,t) * X_hat_FO(j,m,t);
                        w_piX(j,i,n,m)  =  ((m==n) * GAMMA(j,i) * chi(n+(j-1)*N,i,t) - (GAMMA(j,i)^2) * chi(n+(j-1)*N,i,t) * chi(m+(j-1)*N,i,t)) * pi_hat(n+(j-1)*N,i,t) * X_hat_FO(j,m,t);
                    end
                end
            end
        end
        w_update(:,:,t) = w_L + sum(w_pi+w_X,3) + (1/2) * sum(sum(w_pipi+w_XX,4),3) + +sum(sum(w_piX,4),3);
%        w_update(:,:,t) = w_L+ sum(w_pi+w_X,3);
        w_update(:,:,t) = w_update(:,:,t)-w_update(1,1,t); %normalize the first wage to be a constant across periods; should not change anything
    end
    for t=t1:TIME
        checkw(t,1)=max(max(abs(w_hat(:,:,t)-w_update(:,:,t))));
    end
    [wmax loc]=max(checkw);
    wmax;
    %pmax=max(max(abs(p_hat(:,:,t1)-p_old(:,:,t1))));
    %pmax
    
    if ITER_TEMP >30000 || sum(isnan(w_update(:)))>0
        checkw(t1:TIME)
        t1
        disp('Inner loop err')
        ITER_TEMP
        stop
    end
    
    % Step 4e. update the guess
    w_hat = (1-UPDT_W)*w_hat+UPDT_W*w_update;
    %p_old = 0.9*p_old + 0.1*p_hat;
    p_old = p_hat;
    ITER_TEMP=ITER_TEMP+1;
end
end