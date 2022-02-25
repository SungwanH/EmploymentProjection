function [eqm] = PBP_DYN(params, t1, t2, E_T_hat, kappa_hat, L, V, approx)
% Period by period equilibrium
% At every period t1, the perceived path of the shock in fundamental
% (E_T_hat) is updated
% For given L(t1) and the approximation matrices (mu, lambda, s),
% we can solve the expected path of (deviation from SS of) mu, v, L 
% and take those of t1 period outcome.
%%%%%%%%%Algorithm%%%%%%%%%%%%%
%% Roll down parameters

v2struct(params.envr);
v2struct(params.tech);
v2struct(params.modl);

%{
N       =   params.N;
J       =   params.J;
R       =   params.R;
THETA   =   params.THETA;
TIME    =	params.TIME;
BETA    =	params.BETA;
NU      =   params.NU;
UPDT_V  =	params.UPDT_V;
TOLDYN  =	params.TOLDYN;
MAXIT   =   params.MAXIT;
%}
%% Roll down approximation points
mu      =   approx.mu;
lambda  =   approx.lambda;

%t1: starting period 

T_hat = E_T_hat(:,:,:,t2); %belief at t2 period

v_hat = zeros(R*J,TIME);
v_hat(:,t1:TIME) = V(:,t1:TIME); % initial value
W = zeros(J,N,TIME);
%w_hat(:,t1:TIME) = W(:,t1:TIME);
%W = zeros(J,N,TIME);
%for t=1:TIME
%    for j=1:J
%        for n=1:N
%            W(j,n,t) = W_NJ(n+(j-1)*N,t);
%        end
%    end
%end

L_hat = zeros(R*J,TIME);
p_hat = zeros(R*J,TIME);
pi_hat = zeros(R*J,R,TIME);
rw_hat = zeros(J,N,TIME);
%%%%%%%%%Algorithm%%%%%%%%%%%%%

%%%Dynamic problem%%%

VMAX =1; ITER_DYN =1;

while (ITER_DYN <= MAXIT) && (VMAX > TOLDYN)
% Step 2. given the path of value(v_hat), solve for the path of mu_hat
mu_hat = zeros(R*J,R*J,TIME);

for t =t1:TIME-1
    mu_aux = mu(:,:,t);
    mu_hat(:,:,t) = (BETA/NU) * v_hat(:,t+1)' - (BETA/NU) * mu_aux*v_hat(:,t+1);
end

% Step 3. given mu_hat, solve for the path of labor(L_hat)
L_hat(:,t1) = L;
for t = t1:TIME-1
        lambda_aux = lambda(:,:,t);
    for i = 1:R*J
%        L_hat(i,t+1) = lambda_aux(:,i)'*(mu_hat(:,i,t) + L_aux(:,1,t));
%        L_hat(:,t+1) = lambda_aux'*mu_hat(:,:,t) + lambda_aux' * L_hat(:,t);
% should it be mu(t+1) or mu(t)?
        L_hat(i,t+1) = lambda_aux(:,i)'*mu_hat(:,i,t) + lambda_aux(:,i)' * L_hat(:,t);
    end
end

% Step 4. Inner loop (solves w, P, pi)
%%%Temporary problem (=Trade equilibrium)%%%
[w_hat, p_hat, P_hat, pi_hat, X_hat] = PBP_TEMP(params, t1, T_hat, kappa_hat, W, L_hat, approx);
W= w_hat;
for t=t1:TIME
    rw_hat(:,:,t) = w_hat(:,:,t) - ones(J,1)*P_hat(:,:,t); %real wage
end
rw_hat_RJ = reshape(rw_hat(:,1:R,:),[R*J,TIME]);
% Step 5. solve for a new path of v_hat
v_hat_update = zeros(R*J,TIME);
%v_hat_SS(:,1) = (eye(R*J) - BETA * mu(:,:,TIME))\(w_hat_RJ(1:R*J,TIME) - P_hat(1:R*J,TIME));
v_hat_SS(:,1) = (eye(R*J) - BETA * mu(:,:,TIME))\(rw_hat_RJ(:,TIME));
v_hat_update(:,TIME) = v_hat_SS;

for t=TIME-1:-1:t1 
    mu_aux = mu(1:R*J,1:R*J,t);
    % here p_hat should be p_hat (R) not R*J 
%    v_hat_update(:,t) = w_hat_RJ(:,t) - p_hat_RJ(:,t) + BETA * mu_aux(:,:)*v_hat_update(:,t+1);
    v_hat_update(:,t) = rw_hat_RJ(:,t) + BETA * mu_aux(:,:)*v_hat_update(:,t+1);
end

for t=t1:TIME
    checkV(t,1)=max(abs(v_hat_update(:,t)-v_hat(:,t)));
%    checkV(t,1) = norm((v_hat_update(:,t) - v_hat(:,t))./(v_hat_update(:,t)),1);
end
VMAX=max(checkV)
    if ITER_DYN >1000 || sum(isnan(v_hat_update(:)))>0
        checkV(t1:TIME)
        t1
        stop
        disp('Outer loop err')
    end

v_hat=(1-UPDT_V)*v_hat+(UPDT_V)*v_hat_update;


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