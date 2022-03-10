function [w_hat, p_hat, P_hat, pi_hat, X_hat] = PBP_TEMP(params, t1, T_hat, kappa_hat, W, L_hat_R, approx)

%%%Temporary problem (=Trade equilibrium)%%%
%% Roll down parameters
v2struct(params.envr);
v2struct(params.tech);
v2struct(params.modl);
%{
N      =   params.N;
J      =   params.J;
THETA  =   params.THETA;
TIME    =	params.TIME;
UPDT_W  =	params.UPDT_W;
TOLTEMP =	params.TOLTEMP;
MAXIT   =	params.MAXIT;
%}
%% Roll down approximation points
varrho  =   approx.varrho;
zeta    =   approx.zeta;
chi     =   approx.chi;
pi      =   approx.pi;

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
w_hat = W;
w_update = zeros(J,N,TIME);
p_hat = zeros(J,N,TIME);
pi_hat = zeros(N*J,N,TIME);
X_hat = zeros(J,N,TIME);
while (ITER_TEMP <= MAXIT) && (wmax > TOLTEMP)
    
    % Step 4b,c. solve for p_hat and pi_hat

    for t = t1:TIME
        % Step 4b. solve for p_hat
        % for given j, 
        % A*p_hat = RHS
        % the first two rows of A are:
        % ( 1-(1-B_1j)pi_1j1   -(1-B_2j)pi_1j2    -(1-B_3j)pi_1j3 ...)
        % ( -(1-B_1j)pi_2j1   1-(1-B_2j)pi_2j2    -(1-B_3j)pi_2j3 ...)
        RHS = zeros(J,N);
    
        pi_aux =  pi(:,:,t); %Note that pi and kappa are organized by sector; i.e., two consecutive rows are same sectors in different countries
        w_temp = w_hat(:,:,t);
        T_temp = T_hat(:,:,t);
        kappa_temp = kappa_hat(:,:,t);
        for j=1:J
            for n=1:N
                RHS(j,n) = pi_aux(n+(j-1)*N,:)*(GAMMA(j,:).*w_temp(j,:) + kappa_temp(n+(j-1)*N,:) - (1/THETA(j)) .* T_temp(j,:))'; 
            end
        end
        A = zeros(N,N);
        
        for j=1:J
            for n=1:N
                for i=1:N
                    A(n,i) = -(1-GAMMA(j,i))*pi_aux(n+(j-1)*N, i);
                end
                A(n,n) = 1-(1-GAMMA(j,n))*pi_aux(n+(j-1)*N, n);
            end
            p_temp(j,:) = (A\(RHS(j,:)'))';
        end
        p_hat(:,:,t) = p_temp;
    end

    %price index
    for t=t1:TIME
        P_hat(:,:,t) = sum((ALPHAS.*p_hat(:,:,t)),1);
    end
    % Step 4c. solve for pi_hat
    for t = t1:TIME
        for n=1:N
            for j=1:J
%               pi_hat(n+(j-1)*N,ii,t) = -THETA(j)*(B(j,ii) * w_hat(j,ii,t) - (1-B(j,ii)*p_hat(j,ii,t)) + TAU(n+(j-1)*N,ii)) + T_hat(j,ii,t);
                pi_hat(n+(j-1)*N,:,t) = -THETA(j)*(GAMMA(j,:) .* w_hat(j,:,t) + (1-GAMMA(j,:)).*p_hat(j,:,t) - p_hat(j,n,t) - kappa_hat(n+(j-1)*N,:,t)) + T_hat(j,:,t);
             end
         end
    end        

    % Step 4d. Solve for total expenditure
    for t = t1:TIME
        RHS = NaN(J,N);
        RHS1=NaN(J,N,N);
        RHS2=NaN(J,N,J);
        
        varrho_aux = varrho(:,:,t);
        pi_temp = pi_hat(:,:,t);
        zeta_aux = zeta(:,:,t);
        w_temp = w_hat(:,:,t);
        L_temp = L_hat(:,:,t);
        for n=1:N
            for j=1:J
                RHS1(j,:,n) = (1-GAMMA(j,:)).*varrho_aux(n+(j-1)*N,:) .* pi_temp(n+(j-1)*N,:);
            end
        end
        for ii=1:N
            for j=1:J
                for k=1:J
                    RHS2(j,ii,k) = ALPHAS(j,ii)*(zeta_aux(ii+(k-1)*N,j) *(w_temp(k,ii)+L_temp(k,ii))); 
                end
            end
        end
        RHS = sum(RHS1,3) +sum(RHS2,3);
        A = zeros(N,N);
        for j=1:J
            for ii=1:N
                for n=1:N
                    A(ii,n) = -(1-GAMMA(j,ii))*varrho_aux(n+(j-1)*N, ii);
                end
                A(ii,ii) = 1-(1-GAMMA(j,ii))*varrho_aux(ii+(j-1)*N, ii);
            end
            X_temp(j,:) = (A\(RHS(j,:)'))';
        end

        X_hat(:,:,t) = X_temp;        
    end
    % Step 4e. Solve for an updated w_hat using the labor market clearing
    for t = t1:TIME
        RHS_temp=NaN(J,N,N);
        
        chi_aux = chi(:,:,t);
        X_temp = X_hat(:,:,t);
        L_temp = L_hat(:,:,t);
        pi_temp = pi_hat(:,:,t);
        for ii=1:N
            for j=1:J
                for n=1:N
                    RHS_temp(j,ii,n) = chi_aux(n+(j-1)*N, ii) * (pi_temp(n+(j-1)*N, ii) + X_temp(j,n));
                end
            end
        end
        w_update(:,:,t)= GAMMA.*sum(RHS_temp,3) -L_temp;
        w_update(:,:,t)=w_update(:,:,t)-w_update(1,1,t); %normalize the first wage to be a constant across periods; should not change anything
    end

    for t=t1:TIME
        checkw(t,1)=max(max(abs(w_hat(:,:,t)-w_update(:,:,t))));
%    checkw(t,1) = norm((w_update(:,t) - w_hat(:,t))./w_update(:,t),1);
    end
    [wmax loc]=max(checkw);
    wmax;
    if ITER_TEMP >3000 || sum(isnan(w_update(:)))>0
        checkw(t1:TIME)
        t1
        disp('Inner loop err')
        ITER_TEMP
        stop
    end
    
    % Step 4e. update the guess
    w_hat = (1-UPDT_W)*w_hat+UPDT_W*w_update;
    ITER_TEMP=ITER_TEMP+1;
end
end