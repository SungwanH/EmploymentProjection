function [eqm] = PBP_DYN_HETERO(params, t1, t2, E_T_hat, kappa_hat, L_own, L_other, approx, mat_pbp, type)
% Period by period equilibrium
% At every period t1, the perceived path of the shock in fundamental
% (E_T_hat) is updated
% For given L(t1) and the approximation matrices (mu, lambda, chi, varrho, zeta, pi),
% we can solve the expected path of (deviation from PF of) mu, v, L 
% and take those of t1 period outcome.
%%%%%%%%%Algorithm%%%%%%%%%%%%%
%% Roll down parameters

v2struct(params.envr);
v2struct(params.tech);
v2struct(params.modl);
v2struct(params.hetero);

%% Roll down approximation points
mu      =   approx.mu; % note that mu is organized by country, i.e., two consecutive observations are different sectors within the same country
lambda  =   approx.lambda;

T_hat = E_T_hat(:,:,:,t2); %belief at t2 period
%T_hat(:,:,t1,t2) = zeros(J,N);

v_own_hat = zeros(R*J,TIME);
%v_own_hat(:,t1:TIME) = V(:,t1:TIME); % initial value
if type == 1
    L_own_share = share_A; % When this code runs for type A
else
    L_own_share = 1 - share_A; % When this code runs for type B
end
L_hat = NaN(R*J,TIME);
L_own_hat = NaN(R*J,TIME);
L_other_hat = NaN(R*J,TIME);
p_hat = NaN(R*J,TIME);
pi_hat = NaN(R*J,R,TIME);
rw_hat = NaN(J,N,TIME);
mu_own_hat = NaN(R*J,R*J,TIME);
mu_other_hat = NaN(R*J,R*J,TIME);

%%%%%%%%%Algorithm%%%%%%%%%%%%%

%%%Dynamic problem%%%

VMAX =1; ITER_DYN =1;
%fprintf('%s=', 'VMAX')
while (ITER_DYN <= MAXIT) && (VMAX > TOLDYN)
    % Step 2. given the path of value(v_hat), solve for the path of mu_hat

    for t =t1:TIME-1
        mu_aux = mu(:,:,t);
        mu_own_hat(:,:,t) = (BETA/NU) * ones(R*J,1)*v_own_hat(:,t+1)' - (BETA/NU) * mu_aux*v_own_hat(:,t+1);
    end
    mu_own_hat(:,:,TIME) = mu_own_hat(:,:,TIME-1);
    % Step 3. given mu_hat, solve for the path of labor(L_hat)
    % L(nj)_t+1 = lambda(ik,nj)_t+1 * (mu(ik,nj)_t + L(ik)_t)
    L_own_hat(:,t1) = L_own; % initial value of own type
    L_other_hat(:,t1) = zeros(R*J,1);
    L_hat(:,t1) = zeros(R*J,1);
    for t = t1:TIME-1
        lambda_aux = lambda(:,:,t+1);
        if t==t1
            for i=1:R*J
                L_own_hat(i,t+1) = lambda_aux(:,i)'*mu_own_hat(:,i,t) + lambda_aux(:,i)' * L_own_hat(:,t);
                L_other_hat(i,t+1) = L_other(i); % other type's t1+1 period labor
                L_hat(i,t+1) = L_own_share(i) * L_own_hat(i,t+1) + (1-L_own_share(i)) * L_other_hat(i,t+1);
            end
        else
            for i = 1:R*J
                L_own_hat(i,t+1) = lambda_aux(:,i)'*mu_own_hat(:,i,t) + lambda_aux(:,i)' * L_own_hat(:,t);
                L_other_hat(i,t+1) = lambda_aux(:,i)'*mu_own_hat(:,i,t) + lambda_aux(:,i)' * L_other_hat(:,t);
                L_hat(i,t+1) = L_own_share(i) * L_own_hat(i,t+1) + (1-L_own_share(i)) * L_other_hat(i,t+1);
            end
        end 
    end
    
    % Step 4. Inner loop (Temporary problem: solves w, P, pi)
%    [w_hat_iter, p_hat_iter, P_hat_iter, pi_hat_iter, X_hat_iter] = PBP_TEMP(params, t1, T_hat, kappa_hat, W, L_hat, approx); %iterative method
%    W= w_hat; % will be used as next period's initial value
    [w_hat, p_hat, P_hat, pi_hat, X_hat] = PBP_TEMP_MAT(params, t1, T_hat, kappa_hat, L_hat, approx, mat_pbp); % matrix inversion

    for t=t1:TIME
        rw_hat(:,:,t) = w_hat(:,:,t) - ones(J,1)*P_hat(1,:,t); %real wage
    end
    rw_hat_RJ = reshape(rw_hat(:,1:R,:),[R*J,TIME]); 
    % Step 5. solve for a new path of v_hat
    v_own_hat_update = NaN(R*J,TIME);
    %v_hat_SS(:,1) = (eye(R*J) - BETA * mu(:,:,TIME))\(w_hat_RJ(1:R*J,TIME) - P_hat(1:R*J,TIME));
    v_own_hat_SS = (eye(R*J) - BETA * mu(:,:,TIME))\(rw_hat_RJ(:,TIME));
    v_own_hat_update(:,TIME) = v_own_hat_SS;

    for t=TIME-1:-1:t1 
        mu_aux = mu(:,:,t);
        v_own_hat_update(:,t) = rw_hat_RJ(:,t) + BETA * mu_aux(:,:)*v_own_hat_update(:,t+1);
    end

    for t=t1:TIME
        checkV(t,1)=max(abs(v_own_hat_update(:,t)-v_own_hat(:,t)));
    end
    VMAX=max(checkV);
     %fprintf('%3.2e, ', VMAX);
        if ITER_DYN >1000 || sum(sum(isnan(v_own_hat_update(:,t1:end))))>0
            checkV(t1:TIME)
            t1
            stop
            disp('Outer loop err')
        end

    v_own_hat(:,t1:TIME)=(1-UPDT_V)*v_own_hat(:,t1:TIME)+(UPDT_V)*v_own_hat_update(:,t1:TIME);


    ITER_DYN=ITER_DYN+1;
end

eqm.w       = w_hat;
eqm.p       = p_hat;
eqm.P       = P_hat;
eqm.pi      = pi_hat;
eqm.mu      = mu_own_hat;
eqm.L       = L_hat; 
eqm.L_own   = L_own_hat; 
eqm.L_other = L_other_hat; 
eqm.v       = v_own_hat;
eqm.X       = X_hat;
end