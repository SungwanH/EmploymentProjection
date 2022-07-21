function [w_hat, p_hat, P_hat, pi_hat, X_hat] = PBP_TEMP_SO_SIM(params, t1, T_hat, kappa_hat, L_hat_R, approx, SO_temp)

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

v2struct(SO_temp);
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
X_w = zeros(J,N,J); X_L = zeros(J,N,J);  X_pi = zeros(J,N,N); X_X = zeros(J,N,N);
w_pi = zeros(J,N,N); w_X = zeros(J,N,N); w_L = zeros(J,N); 

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
                end
            end
        end
    p_hat(:,:,t) = sum(P_w+P_p+P_T,3) + SO_p(:,:,t);
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
                end
                for k=1:J
                    %diff w.r.t. w
                    X_w(j,i,k)  = ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) * w_hat(k,i,t);
                    %diff w.r.t. L
                    X_L(j,i,k)  = ALPHAS(j,i) * zeta(i+(k-1)*N,j,t) * L_hat(k,i,t);
                end
            end
        end
        X_hat(:,:,t) = sum(X_pi+X_X,3) +sum(X_w+X_L,3) + SO_X(:,:,t);
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
                end
            end
        end
        w_update(:,:,t) = w_L + sum(w_pi+w_X,3) + SO_w(:,:,t);
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