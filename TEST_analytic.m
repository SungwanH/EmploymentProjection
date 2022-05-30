%This is a test code for computing the temporary equilibrium analyticially

clear 
close all
clc;
%data_4_sector=load('DATA/BASE_FOURSECTOR.mat', 'VALjn00', 'Din00','mu0','L0');
%% Load parameter values and approximation points
load('DATA_onesector/DGP.mat', 'eqm_dgp','approx_dgp');
N=10;
J=1;
E_T_hat = eqm_dgp.E_T_hat;
pi=approx_dgp.pi;
chi=approx_dgp.chi;
varrho=approx_dgp.varrho;
params=PARAMS(0.5);
v2struct(params.envr);
THETA  = params.modl.THETA;
GAMMA  = params.modl.GAMMA;
ALPHAS = params.modl.ALPHAS;
clear approx_dgp VALjn00 Din00 mu0 L0
ALPHAS=1/J*ones(J,N);
THETA=4*ones(J,1);
GAMMA=0.5*ones(J,N);
%pi=ones(N*J,N,TIME) / N ;
BETA = 0.8;
NU = 1.3;
TIME=1;

%% Shocks
T_hat = zeros(J,N,TIME);
%T_hat(1,1,1:TIME) = 0.1;
T_hat = E_T_hat(1:J,1:N,2,1);
L_hat = ones(N*J,TIME) *0.1;
L_hat = zeros(N*J,TIME);
Kappa_hat = zeros(N*J,N*J,TIME);

%% 
BTHETA = zeros(N*J,N*J,TIME);
% pi(n*j,i): country n's expenditure share on good j from country i
% sum(pi(1,:))=1
tic
for t=1:TIME
    for i=1:N
        for n=1:N
            for j=1:J
                BTHETA(n+(j-1)*N,i+(j-1)*N,t)= pi(n+(j-1)*N,i,t)*(1-GAMMA(j,i));
            end
        end
    end
end
toc
%Leontief inverse
DELTA = zeros(N*J,N*J,TIME);
for t=1:TIME
    DELTA(:,:,t) = inv(eye(N*J)-BTHETA(:,:,t));
end

%C = -theta^j Gamma^ij
for i=1:N
    for j=1:J
        C(i,j) = -THETA(j) * GAMMA(j,i);
    end
end
%D(ijnom) = -theta^j{(1-Gamma^ij)*DELTA(ijoj)-DETLA(njoj))*pi(ojmj)GAMMA(mj)}
%E(nji) = T(ij) -theta(j){kappa(nji) + sum_o sum_m (1-Gamma^ij)*DELTA(ijoj)-DETLA(njoj))*pi(ojmj)kappa(ojmj)-1/theta(j)T(mj)} 
tic
for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                for o=1:N
                    for m=1:N
                        D_temp(n,j,i,m,o,t) = -THETA(j) * ((1-GAMMA(j,i))*DELTA(i+(j-1)*N,o+(j-1)*N,t)-DELTA(n+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t) * GAMMA(j,m);
                        E_temp(n,j,i,o,m,t) =  - THETA(j) * Kappa_hat(n+(j-1)*N,i) + T_hat(j,i) ...
                         - THETA(j) * ((1-GAMMA(j,i)) * DELTA(i+(j-1)*N,o+(j-1)*N,t) - DELTA(n+(j-1)*N,o+(j-1)*N,t))* pi(o+(j-1)*N,m,t) * (Kappa_hat(o+(j-1)*N,m,t)-1/THETA(j)*T_hat(j,m,t)) ;
                    end
                end
            end
        end
    end
end
toc
D = sum(D_temp,5);
E = sum(sum(E_temp,5),4);


for i=1:N
    for j=1:J
        for o=1:N
            U(i+(j-1)*N,o+(j-1)*N,t) = (1-GAMMA(j,i))*varrho(o+(j-1)*N,i,t); 
        end
    end
end
BGAMMA = zeros(N*J,N*J,TIME);
for t=1:TIME
    BGAMMA(:,:,t) = inv(eye(N*J)-U(:,:,t));
end
%G(ijol) = BGAMMA(ijo)(1-Gamma(oj))varrho(ljo)
for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                for l=1:N
                    G(i,j,o,l,t) = BGAMMA(i+(j-1)*N,o+(j-1)*N,t) * (1-GAMMA(j,o))* varrho(l+(j-1)*N,o,t);
                end
            end
        end
    end
end
%H(ijok) = BGAMMA(ijo)(ALPHA(oj)xi(okj))
for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                for k=1:J
                    H(i,j,o,k,t) = BGAMMA(i+(j-1)*N,o+(j-1)*N,t) * (ALPHAS(j,o))* chi(o+(k-1)*N,j,t);
                end
            end
        end
    end
end

%M(ijmj) = GAMMA(ij) sum_n chi(nji)(D(njim)+sum_o sum_l G(njol)D(ljom))
for t=1:TIME
    for n=1:N
        for j=1:J
            for o=1:N
                for l=1:N
                    for m=1:N
                        GD_temp(n,j,m,t,o,l) = G(n,j,o,l,t)*D(l,j,o,m,t);
                    end
                end
            end
        end
    end
end
GD_sum = sum(sum(GD_temp,6),5);
for t=1:TIME
    for i=1:N
        for j=1:J
            for m=1:N
                for n=1:N
                    M_temp(i+(j-1)*N,m+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t) * (D(n,j,i,m,t) + GD_sum(n,j,m,t));
                end
            end
        end
    end
end
M = sum(M_temp,4);

%Q(ijij) = GAMMA(ij)sum_n chi(nji) C(ij) sum_o sum_l G(njol)C(oj)
for t=1:TIME
    for i=1:N
        for j=1:J
            for m=1:N
                for n=1:N
%                    Q_temp(i+(j-1)*N,i+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t) * (C(i,j,t) + GC_sum(n,j,t));
                    Q_temp(i+(j-1)*N,i+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t) * C(i,j,t);
                end
            end
        end
    end
end
Q = sum(Q_temp,4);

%R(ijoj) = GAMMA(ij)sum_n chi(nji) sum_l G(njol)C(oj)
for t=1:TIME
    for n=1:N
        for j=1:J
            for o=1:N
                for l=1:N
                    GC_temp(n,j,o,t,l) = G(n,j,o,l,t)*C(o,j,t);
                end
            end
        end
    end
end
GC_sum = sum(GC_temp,5);
for t=1:TIME
    for i=1:N
        for j=1:J
            for m=1:N
                for n=1:N
                    R_temp(i+(j-1)*N,o+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t) * GC_sum(n,j,o,t);
                end
            end
        end
    end
end
R = sum(R_temp,4);

%T(ijok) = GAMMA(ij)sum_n chi(nji) (sum_i sum_k H(njok)) sum_i H(njok)
T_temp = zeros(N*J,N*J,TIME,N);
for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                for k=1:J
                    for n=1:N
                        T_temp(i+(j-1)*N,o+(k-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t)*H(n,j,o,k,t);
                    end
                end
            end
        end
    end
end
T= sum(T_temp,4);
for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                for o=1:N
                    for l=1:N
                       E_STAR_temp(i+(j-1)*N,t,n,o,l) = GAMMA(j,i)*chi(n+(j-1)*N,i,t) * (E(n,j,i,t) + G(n,j,o,l,t)*E(l,j,o,t));
                    end
                end
            end
        end
    end
end
E_STAR = sum(sum(sum(E_STAR_temp,5),4),3);
% Final equation:
for t=1:TIME
    w_hat(:,t) = (eye(N*J) - M(:,:,t) - Q(:,:,t)- R(:,:,t) - T(:,:,t) ) \ ( (T-eye(N*J))*L_hat(:,t) + E_STAR(:,t));
end
w_hat