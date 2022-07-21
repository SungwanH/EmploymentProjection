function [eqm] = PBP_DYN_SO(params, t1, t2, E_T_hat, kappa_hat, L, V, eqm, approx)
% Period by period equilibrium with second order approximation using
% iterative method (other SO functions use first-order approx. values to
% compute the second order terms)
% At every period t1, the perceived path of the shock in fundamental
% (E_T_hat) is updated
% For given L(t1) and the approximation matrices (mu, lambda, chi, varrho, zeta, pi),
% we can solve the expected path of (deviation from SS of) mu, v, L 
% and take those of t1 period outcome.
%%%%%%%%%Algorithm%%%%%%%%%%%%%
%% Roll down parameters

v2struct(params.envr);
v2struct(params.tech);
v2struct(params.modl);

%% Roll down approximation points
mu      =   approx.mu; % note that mu is organized by country, i.e., two consecutive observations are different sectors within the same country
lambda  =   approx.lambda;
w_bar   =   eqm.w_lev;
P_bar   =   eqm.P_bar;
rw_bar  =   reshape(w_bar(:,1:R,:)./P_bar(:,1:R,:),[R*J,TIME]); 
T_hat = E_T_hat(:,:,:,t2); %belief at t2 period
%T_hat(:,:,t1,t2) = zeros(J,N);

v_hat = NaN(R*J,TIME);
v_hat(:,t1:TIME) = V(:,t1:TIME); % initial value

L_hat = NaN(R*J,TIME);
p_hat = NaN(R*J,TIME);
pi_hat = NaN(R*J,R,TIME);
rw_hat = NaN(J,N,TIME);
mu_hat = NaN(R*J,R*J,TIME);

%%%%%%%%%Algorithm%%%%%%%%%%%%%

%%%Dynamic problem%%%

VMAX =1; ITER_DYN =1;
fprintf('%s=', 'VMAX')
while (ITER_DYN <= MAXIT) && (VMAX > TOLDYN)
    % Step 2. given the path of value(v_hat), solve for the path of mu_hat
    %{
    for t =t1:TIME-1
        mu_aux = mu(:,:,t);
        % first order
        mu_hat(:,:,t) = (BETA/NU) * ones(R*J,1)*v_hat(:,t+1)' - (BETA/NU) * mu_aux*v_hat(:,t+1);
        % second order
        %
        for i=1:R*J
            for g=1:R*J
                for k=1:R*J
                        mu_v(i,g,k) = -(BETA/NU) * ((k==g) - mu_aux(i,k)) * v_hat(k,t+1); 
                    for n=1:R*J
                        mu_vv(i,g,k,n) = -(BETA/NU)^2 * ( (k==n)- mu_aux(i,n) ) * v_hat(k,t+1) * v_hat(n,t+1); 
%                        mu_vv(i,g,k,n) = (BETA/NU)^2 * ( (k==g)*((g==n) *mu_aux(i,g) - mu_aux(i,g)*mu_aux(i,n)) ...
%                            - ((n==g)*mu_aux(i,g) - mu_aux(i,g)*mu_aux(i,n))*mu_aux(i,k) ...
%                            - ((n==k)*mu_aux(i,k) - mu_aux(i,k)*mu_aux(i,n))*mu_aux(i,g) ) * v_hat(k,t+1) * v_hat(n,t+1); 
                    end
                end
            end
        end
        mu_hat(:,:,t) = sum(mu_v,3) + (1/2) * sum(sum(mu_vv,3),4);    
    end
    %}
    mu_hat(:,:,TIME) = mu_hat(:,:,TIME-1);
    %}
    % Step 3. given mu_hat, solve for the path of labor(L_hat)
    L_hat(:,t1) = L;
   
    % derive parts of second order approximation for L_hat 
    for t = t1:TIME-1
        for k=1:R*J
            for g=1:R*J
                L_L(g,k,t) = lambda(k,g,t+1) * L_hat(k,t);
                L_v(g,k,t) = ( (BETA/NU) * (lambda(:,g,t+1)'*((g==k)-mu(:,k,t))) ) * v_hat(k,t+1);
                for m=1:R*J
                    L_LL(g,k,m,t) = ((k==m)*lambda(k,g,t+1) - ( lambda(k,g,t+1)*lambda(m,g,t+1) )) * L_hat(k,t) * L_hat(m,t);
                    L_vL(g,k,m,t) = ((BETA/NU) * (lambda(m,g,t+1)*((k==g)-mu(m,k,t))) - (BETA/NU) * ( lambda(m,g,t+1)*( lambda(:,g,t+1)'*((k==g)-mu(:,k,t)) ) )) * v_hat(k,t+1) * L_hat(m,t);
                    L_vv1(g,k,m,t) = (BETA/NU)^2 * (-(mu(:,k,t) .* lambda(:,g,t+1))' *((m==g) + (m==k) - 2*mu(:,m,t)));
                    L_vv2(g,k,m,t) = (BETA/NU)^2 * ((k==g) * lambda(:,g,t+1)'*((m==g)-mu(:,m,t)));
                    L_vv3(g,k,m,t) = -(BETA/NU)^2 * (lambda(:,g,t+1)'*((k==g)-mu(:,k,t)))*(lambda(:,g,t+1)'*((m==g)-mu(:,m,t)));
                    L_vv(g,k,m,t) = (L_vv1(g,k,m,t) + L_vv2(g,k,m,t) + L_vv3(g,k,m,t)) * v_hat(k,t+1) * v_hat(m,t+1);
 %                    L_vv4(g,m,t) = (lambda(:,g,t+1)'*((m==g)-mu(:,m,t)));
%                    L_LL(g,k,m) = - ( mu(k,g,t) * mu(m,g,t) / (L_bar(g,t+1))^2 ) * L_bar(k,t) * L_hat(k,t) * L_bar(m,t) * L_hat(m,t);
%                    L_vL1(g,k,m) = (BETA/NU) * (lambda(m,g,t+1)*((k==g)-mu(m,k,t)));
%                    for i=1:R*J
                       %L_vL2(g,m,i,t) = - (BETA/NU) * (lambda(m,g,t+1)*(lambda(i,g,t+1).*((m==g)-mu(i,m,t)));
                       %L_vv1(g,k,m,i,t) = (BETA/NU)^2 * (-mu(i,k,t) * lambda(i,g,t+1) *((m==g) + (m==k) - 2*mu(i,m,t)));
                       %L_vv2(g,k,m,i,t) = (BETA/NU)^2 * ((k==g) * lambda(i,g,t+1)*((m==g)-mu(i,m,t)));
                       %L_vv3(g,k,i,t) = -(BETA/NU)^2 * (lambda(i,g,t+1)*((k==g)-mu(i,k,t)));
                       %L_vv4(g,m,i,t) = (lambda(i,g,t+1)*((m==g)-mu(i,m,t)));
                       %L_vv(g,k,m,i) = (BETA/NU)^2 *( (-L_bar(i,t) * mu(i,k,t) * mu(i,g,t) *((m==g) + (m==k) - 2*mu(i,m,t)) + (k==g) * L_bar(i,t) * mu(i,g,t)*((m==g)-mu(i,m,t)))/L_bar(g,t+1)...
                       %               - (L_bar(i,t)*mu(i,g,t).*((k==g)-mu(i,k,t))*(L_bar(i,t)*mu(i,g,t)*((m==g)-mu(i,m,t))))/(L_bar(g,t+1)^2) ) * v_hat(k,t+1) * v_hat(m,t+1);
                       %L_vL(g,k,m,i) = ( (BETA/NU) * (mu(m,g,t)*((k==g)-mu(m,k,t)))/L_bar(g,t+1) - (BETA/NU) * (mu(m,g,t)*L_bar(i,t)*mu(i,g,t).*((m==g)-mu(i,m,t)))/(L_bar(g,t+1)^2) )...
                       %             * v_hat(k,t+1) * L_hat(m,t) * L_bar(m,t);
%                    end
%                    L_vv(g,k,m,t) = (sum(L_vv1(g,k,m,:,t)) + sum(L_vv2(g,k,:,t)) + sum(L_vv3(g,k,:,t))*sum(L_vv4(g,m,:,t))) * v_hat(k,t+1) * v_hat(m,t+1);
%                    L_vL(g,k,m,t) = (L_vL1(g,k,m,t) + sum(L_vL2(g,m,:,t))) * v_hat(k,t+1) * L_hat(m,t);
                   %L_vL(g,k,m,t) = (L_vL1(g,k,m,t) + L_vL2(g,k,m,t)) * v_hat(k,t+1) * L_hat(m,t);
                end
                
            end
        end
%        L_hat(:,t+1) = sum(sum(L_v(:,:,t),3),2) + sum(L_L(:,:,t),2) ;% first-order
% Fixed version:
        L_hat(:,t+1) = sum(L_v(:,:,t),2) + sum(L_L(:,:,t),2) + (1/2) * sum(sum(L_LL(:,:,:,t),3),2) + (1/2) * sum(sum(L_vv(:,:,:,t),3),2) + sum(sum(L_vL(:,:,:,t),3),2);
    end
    %}
   
    %first order version
    %{
    for t = t1:TIME-1
        lambda_aux = lambda(:,:,t+1);
        for i = 1:R*J
    %    L_hat(:,t+1) = sum(lambda_aux.*mu_hat(:,:,t),1)' + lambda_aux' * L_hat(:,t);
            L_hat(i,t+1) = lambda_aux(:,i)'*mu_hat(:,i,t) + lambda_aux(:,i)' * L_hat(:,t);
        end
    end
    %}
    %L_hat(:,t1+1) - L_hat_SO(:,t1+1)
    

    % Step 4. Inner loop (Temporary problem: solves w, P, pi)
    %[w_hat_iter, p_hat_iter, P_hat_iter, pi_hat_iter, X_hat_iter] = PBP_TEMP(params, t1, T_hat, kappa_hat, W, L_hat, approx); %iterative method
%    [w_hat, p_hat, P_hat, pi_hat, X_hat] = PBP_TEMP(params, t1, T_hat, kappa_hat, W, L_hat, approx); %iterative method
%    W= w_hat; % will be used as next period's initial value
%    [w_hat_SO, p_hat_SO, P_hat_SO, pi_hat_SO, X_hat_SO] = PBP_TEMP_SO(params, t1, T_hat, kappa_hat, W, L_hat, approx); 

    [w_hat, p_hat, P_hat, pi_hat, X_hat] = PBP_TEMP_SO(params, t1, T_hat, kappa_hat, L_hat, approx); % matrix inversion

    for t=t1:TIME
        rw_hat(:,:,t) = w_hat(:,:,t) - ones(J,1)*P_hat(1,:,t); %real wage
    end

    % Step 5. solve for a new path of v_hat
    v_hat_update = NaN(R*J,TIME);
    % first / second order terms of CRRA of v (due to CRRA)    
    rw_hat_RJ = reshape(rw_hat(:,1:R,:),[R*J,TIME]); 
    w_hat_RJ = reshape(w_hat(:,1:R,:),[R*J,TIME]);
    P_hat_RJ = reshape(repmat(P_hat(1,1:R,:),[J,1,1]),[R*J,TIME]);
    v_w(:,t1:TIME) = rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* w_hat_RJ(:,t1:TIME);
    v_P(:,t1:TIME) = -rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* P_hat_RJ(:,t1:TIME);
    v_ww(:,t1:TIME)  = (1-CRRA_PAR) .* rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* w_hat_RJ(:,t1:TIME).^2;
    v_PP(:,t1:TIME)  = (1-CRRA_PAR) .* rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* P_hat_RJ(:,t1:TIME).^2;
    v_wP(:,t1:TIME)  = -(1-CRRA_PAR) .* rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* w_hat_RJ(:,t1:TIME) .* P_hat_RJ(:,t1:TIME);
    
    v_hat_SS = (eye(R*J) - BETA * mu(:,:,TIME))\(v_w(:,TIME) + v_P(:,TIME)); % first-order
    v_hat_update(:,TIME) = v_hat_SS; % consider second-order terms (fixed point?)

     for t=TIME-1:-1:t1 
        mu_aux = mu(:,:,t);
        v_v  = BETA * mu_aux(:,:) * v_hat_update(:,t+1); %first order
        for m=1:R*J
            for k=1:R*J
                for n=1:R*J
                    v_vv(n,k,m) = (BETA^2 / NU) * mu_aux(n,m) * ((k==m) - mu_aux(n,k)) * v_hat_update(m,t+1) * v_hat_update(k,t+1);
                end
            end
        end
        v_hat_update(:,t) = v_w(:,t) + v_P(:,t) + v_v + (1/2) * v_ww(:,t) + (1/2) * v_PP(:,t)+ v_wP(:,t) + (1/2) * sum(sum(v_vv,3),2); % second-order
%        v_hat_update(:,t) = rw_hat_RJ(:,t) + BETA * mu_aux(:,:)*v_hat_update(:,t+1);  %first order
    end
    for t=t1:TIME
        checkV(t,1)=max(abs(v_hat_update(:,t)-v_hat(:,t)));
    %    checkV(t,1) = norm((v_hat_update(:,t) - v_hat(:,t))./(v_hat_update(:,t)),1);
    end
    VMAX=max(checkV);
%     fprintf('%3.2e, ', VMAX);
        if ITER_DYN >1000 || sum(sum(isnan(v_hat_update(:,t1:end))))>0
            checkV(t1:TIME)
            t1
            stop
            disp('Outer loop err')
        end

    v_hat(:,t1:TIME)=(1-UPDT_V)*v_hat(:,t1:TIME)+(UPDT_V)*v_hat_update(:,t1:TIME);


    ITER_DYN=ITER_DYN+1;
end

eqm.w = w_hat;
eqm.p = p_hat;
eqm.P = P_hat;
eqm.pi = pi_hat;
eqm.mu = mu_hat;
eqm.L = L_hat; 
eqm.v = v_hat;
eqm.X = X_hat;
end