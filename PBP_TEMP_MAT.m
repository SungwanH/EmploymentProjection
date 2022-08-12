function [w_hat, p_hat, P_hat, pi_hat, X_hat] = PBP_TEMP_MAT(params, t1, T_hat, kappa_hat, L_hat_R, approx, mat_pbp)
% This function solves for the temporary equilibrium using matrix inversion
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

%% Roll down matrices
v2struct(mat_pbp);
%

% Reformulate the order of shocks
% L_hat_R only includes US states and the order is j+(i-1)*J
% We reformulate it to have an order of i+(j-1)*N
L_hat_T     = zeros(J*N,TIME);
T_hat_T     = zeros(J*N,TIME);
kappa_hat_T = zeros(J*N,TIME);
L_hat = zeros(J,N,TIME);
for t=t1:TIME
    L_hat(1:J,1:R,t)=reshape(L_hat_R(:,t),J,R);
end
L_hat(1:J,(R+1):N,1:TIME) = zeros(J,N-R,TIME); %There is no migration in other countries

for t=t1:TIME
    for i=1:N
        for j=1:J
            L_hat_T(i+(j-1)*N,t) = L_hat(j,i,t);
            T_hat_T(i+(j-1)*N,t) = T_hat(j,i,t);
            for n=1:N
               kappa_hat_T(i+(n-1)*N+(j-1)*N*N,t) = kappa_hat(i+(j-1)*N,n,t);
            end
        end
    end
end

w_hat_NJ = zeros(N*J,TIME);
for t=t1:TIME
    temp_right=((T_tilde(:,:,t)-eye(N*J)) * L_hat_T(:,t) + Q_tilde(:,:,t) * T_hat_T(:,t) + F_tilde(:,:,t) * kappa_hat_T(:,t));
    temp_right_truncation=temp_right(2:end);
    temp_left=(eye(N*J) - M_tilde(:,:,t) - T_tilde(:,:,t) );
    temp_left_truncation=temp_left(2:end, 2:end);
    w_hat_NJ(2:N*J,t)=temp_left_truncation\temp_right_truncation;
%    w_hat_NJ(:,t) = (eye(N*J) - M_tilde(:,:,t) - T_tilde(:,:,t) ) \ ((T_tilde(:,:,t)-eye(N*J)) * L_hat_T(:,t) + Q_tilde(:,:,t) * T_hat_T(:,t) + F_tilde(:,:,t) * kappa_hat_T(:,t));
%    w_hat_NJ(:,t) = pinv(eye(N*J) - M_tilde(:,:,t) - T_tilde(:,:,t) ) * ((T_tilde(:,:,t)-eye(N*J)) * L_hat_T(:,t) + Q_tilde(:,:,t) * T_hat_T(:,t) + F_tilde(:,:,t) * kappa_hat_T(:,t));
end
%w_hat_NJ(:,:)-w_hat_NJ(1,:)
%reshape w_hat_NJ to J by N matrix
w_hat = zeros(J,N,TIME);
for t=t1:TIME
    for n=1:N
        for j=1:J
            w_hat(j,n,t) = w_hat_NJ(n+(j-1)*N,t) - w_hat_NJ(1,t); % normalize
        end
    end
end
%Recovering p
psum_temp = zeros(N*J,TIME,N);
for t=t1:TIME
    for o=1:N
        for j=1:J
            for i=1:N
%                psum_temp(o,j,i,t)= pi(o+(j-1)*N,i,t)*(GAMMA(j,i)*w_hat_NJ(i+(j-1)*N,t)+ kappa_hat(o+(j-1)*N,i,t) - (1/THETA(j)) * T_hat(j,i,t));
                psum_temp(o+(j-1)*N,t,i)= pi(o+(j-1)*N,i,t)*(GAMMA(j,i)*w_hat(j,i,t)+ kappa_hat(o+(j-1)*N,i,t) - (1/THETA(j)) * T_hat(j,i,t));
            end
        end
    end
end
psum=sum(psum_temp,3);
p_hat = zeros(J,N,TIME);
for t=t1:TIME
    psum_temp = psum(:,t);
    DELTA_temp = DELTA(:,:,t);
    for j=1:J
        for n=1:N
%            for o=1:N
%                p_hat_temp(j,n,t,o) = DELTA(n+(j-1)*N,o+(j-1)*N,t) * (psum(o,j,t) );
                p_hat(j,n,t) = DELTA_temp(n+(j-1)*N,:) * (psum_temp );
%            end
        end
    end
end

%price index
P_hat = zeros(1,N,TIME);
for t=t1:TIME
    P_hat(:,:,t) = sum((ALPHAS.*p_hat(:,:,t)),1);
end

% pi_hat
pi_hat = zeros(N*J,N,TIME);
for t = t1:TIME
    for n=1:N
        for j=1:J
%               pi_hat(n+(j-1)*N,ii,t) = -THETA(j)*(B(j,ii) * w_hat(j,ii,t) - (1-B(j,ii)*p_hat(j,ii,t)) + TAU(n+(j-1)*N,ii)) + T_hat(j,ii,t);
            pi_hat(n+(j-1)*N,:,t) = -THETA(j)*(GAMMA(j,:) .* w_hat(j,:,t) + (1-GAMMA(j,:)).*p_hat(j,:,t) - p_hat(j,n,t) + kappa_hat(n+(j-1)*N,:,t)) + T_hat(j,:,t);
         end
     end
end        
X_temp_G = zeros(J,N,N,N,TIME);
X_temp_H = zeros(J,N,N,J,TIME);
X_hat = zeros(J,N,TIME);
for t=t1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                for l=1:N
                    X_temp_G(j,i,o,l,t)=G(i,j,o,l,t)*pi_hat(l+(j-1)*N,o,t);
                end
                for k=1:J
                    X_temp_H(j,i,o,k,t)=H(i,j,o,k,t)*(w_hat(k,o,t)+L_hat(k,o,t));
                end
             end
         end
     end
end        
X_hat = sum(sum(X_temp_G,3),4) + sum(sum(X_temp_H,3),4);
%{
w_hat(:,:,1,1)
p_hat(:,:,1,1)
P_hat(:,:,1,1)
pi_hat(:,:,1,1)
X_hat(:,:,1,1)
stop
%}
end