function [eqm] = PBP_DYN_SO_CONS(params, t1, t2, E_T_hat, kappa_hat, L, V, eqm, approx, eqm_FO)
% Period by period equilibrium with second order approximation using
% constant value for the variables in second-order terms
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
w_bar   =   eqm.w_lev;
P_bar   =   eqm.P_bar;
rw_bar  =   reshape(w_bar(:,1:R,:)./P_bar(:,1:R,:),[R*J,TIME]); 
mu      =   approx.mu; % note that mu is organized by country, i.e., two consecutive observations are different sectors within the same country
lambda  =   approx.lambda;

%w_hat_FO = cons_FO.w_hat;
%p_hat_FO = cons_FO.p_hat;
%v_hat_FO = reshape(eqm_FO.v_hat(:,:,t2),J*R,TIME);
%L_hat_FO = reshape(eqm_FO.L_hat(:,:,:,t2),J*R,TIME);
v_hat_FO = eqm_FO.v;
L_hat_FO = eqm_FO.L;
w_hat_FO = eqm_FO.w;
P_hat_FO = eqm_FO.P;
T_hat = E_T_hat(:,:,:,t2); %belief at t2 period

v_hat = NaN(R*J,TIME);
v_hat(:,t1:TIME) = V(:,t1:TIME); % initial value

L_hat  = NaN(R*J,TIME);
p_hat  = NaN(R*J,TIME);
pi_hat = NaN(R*J,R,TIME);
rw_hat = NaN(J,N,TIME);
mu_hat = NaN(R*J,R*J,TIME);

%%%%%%%%%Algorithm%%%%%%%%%%%%%

%%%Dynamic problem%%%

VMAX =1; ITER_DYN =1;
fprintf('%s=', 'VMAX')
while (ITER_DYN <= MAXIT) && (VMAX > TOLDYN)
    % Step 3. given mu_hat, solve for the path of labor(L_hat)
    L_hat(:,t1) = L;
   
    % derive parts of second order approximation for L_hat 
    for t = t1:TIME-1
        for k=1:R*J
            for g=1:R*J
                L_L(g,k,t) = lambda(k,g,t+1) * L_hat(k,t);
                L_v(g,k,t) = ( (BETA/NU) * (lambda(:,g,t+1)'*((g==k)-mu(:,k,t))) ) * v_hat(k,t+1);

                for m=1:R*J
                    L_LL(g,k,m,t)  = ((k==m)*lambda(k,g,t+1) - ( lambda(k,g,t+1)*lambda(m,g,t+1) )) * L_hat_FO(k,t) * L_hat_FO(m,t);
                    L_vL(g,k,m,t)  = ((BETA/NU) * (lambda(m,g,t+1)*((k==g)-mu(m,k,t))) - (BETA/NU) * ( lambda(m,g,t+1)*( lambda(:,g,t+1)'*((k==g)-mu(:,k,t)) ) )) * v_hat_FO(k,t+1) * L_hat_FO(m,t);
                    L_vv1(g,k,m,t) = (BETA/NU)^2 * (-(mu(:,k,t) .* lambda(:,g,t+1))' *((m==g) + (m==k) - 2*mu(:,m,t)));
                    L_vv2(g,k,m,t) = (BETA/NU)^2 * ((k==g) * lambda(:,g,t+1)'*((m==g)-mu(:,m,t)));
                    L_vv3(g,k,m,t) = -(BETA/NU)^2 * (lambda(:,g,t+1)'*((k==g)-mu(:,k,t)))*(lambda(:,g,t+1)'*((m==g)-mu(:,m,t)));
                    L_vv(g,k,m,t)  = (L_vv1(g,k,m,t) + L_vv2(g,k,m,t) + L_vv3(g,k,m,t)) * v_hat_FO(k,t+1) * v_hat_FO(m,t+1);
                end
                
            end
        end
        L_hat(:,t+1) = sum(L_v(:,:,t),2) + sum(L_L(:,:,t),2) + (1/2) * sum(sum(L_LL(:,:,:,t),3),2) + (1/2) * sum(sum(L_vv(:,:,:,t),3),2) + sum(sum(L_vL(:,:,:,t),3),2);
    end
   
    % Step 4. Inner loop (Temporary problem: solves w, P, pi)
    [w_hat, p_hat, P_hat, pi_hat, X_hat] = PBP_TEMP_SO_CONS(params, t1, T_hat, kappa_hat, L_hat, approx, eqm_FO); % matrix inversion

    for t=t1:TIME
        rw_hat(:,:,t) = w_hat(:,:,t) - ones(J,1)*P_hat(1,:,t); %real wage
    end
    rw_hat_RJ = reshape(rw_hat(:,1:R,:),[R*J,TIME]); 
    w_hat_RJ = reshape(w_hat(:,1:R,:),[R*J,TIME]);
    P_hat_RJ = reshape(repmat(P_hat(1,1:R,:),[J,1,1]),[R*J,TIME]);
    w_hat_RJ_FO = reshape(w_hat_FO(:,1:R,:),[R*J,TIME]);
    P_hat_RJ_FO = reshape(repmat(P_hat_FO(1,1:R,:),[J,1,1]),[R*J,TIME]);
    % Step 5. solve for a new path of v_hat
    v_hat_update = NaN(R*J,TIME);
    % first / second order terms of CRRA of v (due to CRRA)
    v_w(:,t1:TIME) = rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* w_hat_RJ(:,t1:TIME);
    v_P(:,t1:TIME) = -rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* P_hat_RJ(:,t1:TIME);
    v_ww(:,t1:TIME)  = (1-CRRA_PAR) .* rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* w_hat_RJ_FO(:,t1:TIME).^2;
    v_PP(:,t1:TIME)  = (1-CRRA_PAR) .* rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* P_hat_RJ_FO(:,t1:TIME).^2;
    v_wP(:,t1:TIME)  = -(1-CRRA_PAR) .* rw_bar(:,t1:TIME).^(1-CRRA_PAR) .* w_hat_RJ_FO(:,t1:TIME) .* P_hat_RJ_FO(:,t1:TIME);
 
    %v_hat_SS = (eye(R*J) - BETA * mu(:,:,TIME))\(rw_hat_RJ(:,TIME)); % first-order
    v_hat_SS = (eye(R*J) - BETA * mu(:,:,TIME))\(v_w(:,TIME) + v_P(:,TIME)); % first-order
    v_hat_update(:,TIME) = v_hat_SS; % consider second-order terms (fixed point?)


    for t=TIME-1:-1:t1 
        mu_aux = mu(:,:,t);
        v_v  = BETA * mu_aux(:,:) * v_hat_update(:,t+1); %first order
        for m=1:R*J
            for k=1:R*J
                for n=1:R*J
                    v_vv(n,k,m) = (BETA^2 / NU) * mu_aux(n,m) * ((k==m) - mu_aux(n,k)) * v_hat_FO(m,t+1) * v_hat_FO(k,t+1);
                end
            end
        end
%        v_hat_update(:,t) = rw_hat_RJ(:,t) + v_v + (1/2) * sum(sum(v_vv,3),2); % second-order
        v_hat_update(:,t) = v_w(:,t) + v_P(:,t) + v_v + (1/2) * sum(sum(v_vv,3),2) + (1/2) * v_ww(:,t) + (1/2) * v_PP(:,t)+ v_wP(:,t); % second-order
    end
    for t=t1:TIME
        checkV(t,1)=max(abs(v_hat_update(:,t)-v_hat(:,t)));
    end
    VMAX=max(checkV);
     fprintf('%3.2e, ', VMAX);
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