%This is a test code for computing the temporary equilibrium analyticially

clear 
close all
clc;
%data_4_sector=load('DATA/BASE_FOURSECTOR.mat', 'VALjn00', 'Din00','mu0','L0');
%{
load('DATA_onesector/DGP.mat', 'eqm_dgp','approx_dgp');
%% Load parameter values and approximation points
N=10;
J=1;
E_T_hat = eqm_dgp.E_T_hat;
pi=approx_dgp.pi;
chi=approx_dgp.chi;
zeta=approx_dgp.zeta;
varrho=approx_dgp.varrho;
params=PARAMS(0.5);
v2struct(params.envr);
THETA  = params.modl.THETA;
GAMMA  = params.modl.GAMMA;
ALPHAS = params.modl.ALPHAS;
clear approx_dgp VALjn00 Din00 mu0 L0
%}
rng(12093248)
TIME=1;
N=5;
J=2;
pi = rand(N*J,N,TIME);
varrho = rand(N*J,N,TIME);
chi = rand(N*J,N,TIME);
zeta = ones(N*J,J,TIME) * 1/J;
for t=1:TIME
    for i=1:N
        pi(i,:,t) = pi(i,:,t)./sum(pi(i,:,t)); % should sum up to 1
        varrho(i,:,t) = varrho(i,:,t)./sum(varrho(i,:,t));% should sum up to 1
    end
end
ALPHAS=1/J*ones(J,N);
THETA=4*ones(J,1);
%GAMMA=0.5*ones(J,N);
GAMMA=rand(J,N);

%% Shocks
T_hat = zeros(J,N,TIME);
T_hat = rand(J,N,TIME);
%T_hat(1,1,1:TIME) = 0.1;
%T_hat = E_T_hat(1:J,1:N,2,1);
%L_hat = ones(N*J,TIME) *0.1;
L_hat = zeros(N*J,TIME);
kappa_hat = zeros(N*J,N,TIME);
%kappa_hat = rand(N*J,N,TIME);
w_hat = rand(J,N,TIME)*0.1; %For test


BTHETA = zeros(N*J,N*J,TIME);
% pi(n*j,i): country n's expenditure share on good j from country i
% sum(pi(1,:))=1

for t=1:TIME
    for i=1:N
        for n=1:N
            for j=1:J
                BTHETA(n+(j-1)*N,i+(j-1)*N,t)= pi(n+(j-1)*N,i,t)*(1-GAMMA(j,i));
            end
        end
    end
end

%Leontief inverse
DELTA = zeros(N*J,N*J,TIME);
for t=1:TIME
    DELTA(:,:,t) = inv(eye(N*J)-BTHETA(:,:,t));
end

%% Test: obtain p from the new code
for t=1:TIME
    for i=1:N
        for o=1:N
            for j=1:J
                Psum_temp(o,j,i,t)= pi(o+(j-1)*N,i,t)*(GAMMA(j,i)*w_hat(j,i,t)+ kappa_hat(o+(j-1)*N,i,t) - (1/THETA(j)) * T_hat(j,i,t));
            end
        end
    end
end
Psum=sum(Psum_temp,3);
for t=1:TIME
    for j=1:J
        for n=1:N
            for o=1:N
                p_hat_temp(j,n,o,t) = DELTA(n+(j-1)*N,o+(j-1)*N,t) * (Psum(o,j,t) );
            end
        end
    end
end
p_hat_nj = sum(p_hat_temp,3);
%% TEST Compare with the old code
for t = 1:TIME
        % Step 4b. solve for p_hat
        % for given j, 
        % A*p_hat = RHS
        % the first two rows of A are:
        % ( 1-(1-B_1j)pi_1j1   -(1-B_2j)pi_1j2    -(1-B_3j)pi_1j3 ...)
        % ( -(1-B_1j)pi_2j1   1-(1-B_2j)pi_2j2    -(1-B_3j)pi_2j3 ...)
        RHS = zeros(J,N);
    
        pi_aux =  pi(1:N*J,1:N,t); %Note that pi and kappa are organized by sector; i.e., two consecutive rows are same sectors in different countries
        w_temp = w_hat(1:J,1:N,t);
        T_temp = T_hat(1:J,1:N,t);
        kappa_temp = kappa_hat(1:N*J,1:N,t);
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

%T_hat_temp =zeros(J,N,TIME);
w_hat_temp = zeros(J,N,TIME);
% Pi_hat
for t = 1:TIME
    for n=1:N
        for j=1:J
            for i=1:N
%               pi_hat(n+(j-1)*N,ii,t) = -THETA(j)*(B(j,ii) * w_hat(j,ii,t) - (1-B(j,ii)*p_hat(j,ii,t)) + TAU(n+(j-1)*N,ii)) + T_hat(j,ii,t);
                pi_hat(n+(j-1)*N,i,t) = -THETA(j)*(GAMMA(j,i) .* w_hat(j,i,t) + (1-GAMMA(j,i)).*p_hat(j,i,t) - p_hat(j,n,t) - kappa_hat(n+(j-1)*N,i,t)) + T_hat(j,i,t);
             end
         end
     end
end        

p_hat_nj % p from the new code
p_hat % p from the old code
%% 
%C = -THETA^j Gamma^ij
C= zeros(N,J);
for i=1:N
    for j=1:J
        C(i,j) = -THETA(j) * GAMMA(j,i);
    end
end
%D(ijnom) =
%-THETA^j{(1-Gamma^ij)*DELTA(ijoj)-DETLA(njoj))*pi(ojmj)GAMMA(mj)}
%E(nji) = T(ij) -THETA(j){kappa(nji) + sum_o sum_m (1-Gamma^ij)*DELTA(ijoj)-DETLA(njoj))*pi(ojmj)kappa(ojmj)-1/THETA(j)T(mj)} 
D_temp = zeros(N,J,N,N,N,TIME);
E_temp = zeros(N,J,N,N,N,TIME);

for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                for o=1:N
                    for m=1:N
                        D_temp(n,j,i,m,o,t) = -THETA(j) * ((1-GAMMA(j,i))*DELTA(i+(j-1)*N,o+(j-1)*N,t)-DELTA(n+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t) * GAMMA(j,m);
                        D_temp_ij(j,i,m,o,t) = -THETA(j) * ((1-GAMMA(j,i))*DELTA(i+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t) * GAMMA(j,m);
                        D_temp_nj(j,n,m,o,t) = -THETA(j) * (-DELTA(n+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t) * GAMMA(j,m);
%                        E_temp(n,j,i,m,o,t) =  - THETA(j) * kappa_hat(n+(j-1)*N,i) + T_hat(j,i) ...
%                         - THETA(j) * ((1-GAMMA(j,i))* DELTA(i+(j-1)*N,o+(j-1)*N,t) - DELTA(n+(j-1)*N,o+(j-1)*N,t))* pi(o+(j-1)*N,m,t) * (kappa_hat(o+(j-1)*N,m,t) - 1/THETA(j) * T_hat(j,m,t)) ;
                        E_temp(n,j,i,m,o,t) = ((1-GAMMA(j,i))* DELTA(i+(j-1)*N,o+(j-1)*N,t) - DELTA(n+(j-1)*N,o+(j-1)*N,t))* pi(o+(j-1)*N,m,t);
                        E_temp_ij(j,i,m,o,t) = ((1-GAMMA(j,i))* DELTA(i+(j-1)*N,o+(j-1)*N,t))* pi(o+(j-1)*N,m,t);
                        E_temp_nj(j,n,m,o,t) = (- DELTA(n+(j-1)*N,o+(j-1)*N,t))* pi(o+(j-1)*N,m,t);
                    end
                end
            end
        end
    end
end
   
D = sum(D_temp,5);%
D_ij = sum(D_temp_ij,4);
D_nj = sum(D_temp_nj,4);
E_ij = sum(E_temp_ij,4);
E_nj = sum(E_temp_nj,4);
E = sum(E_temp,5);
F = E_temp;
%% Test to check whether using D_A and D_B can correctly recover P
D_temp_A = zeros(J,N,N,N,TIME);
D_temp_B = zeros(J,N,N,N,TIME);
for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                for m=1:N
                    D_temp_A(j,i,m,o,t) = (DELTA(i+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t) * GAMMA(j,m);
                    D_temp_B(j,i,m,o,t) = (DELTA(i+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t);
                end
            end
        end
    end
end
%                    
D_A = sum(D_temp_A,4);
D_B = sum(D_temp_B,4);
p_recov_temp= zeros(J,N,N,TIME);
for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                for m=1:N
                    for o=1:N
                        p_recov_temp(j,i,m,t) = D_A(j,i,m,t)*w_hat(j,m,t) + D_B(j,i,m,t)* (- (1/THETA(j))*T_hat(j,m,t)); %Assume kappa is zero
                    end
                end
             end
        end
    end
end
p_recov=sum(p_recov_temp,3);
%E = sum(sum(E_temp,5),4);

%% Test (new code) Pi
%T_hat_temp =zeros(J,N,TIME);
for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                for m=1:N
                        Pi_temp(n+(j-1)*N,i,m,t) = C(i,j) * w_hat(j,i,t) + D_ij(j,i,m,t)*w_hat(j,m,t)+ D_nj(j,n,m,t)*w_hat(j,m,t) + T_hat(j,i,t) - THETA(j)* kappa_hat(n+(j-1)*N,i,t) + E_ij(j,i,m,t)*T_hat(j,m,t) + E_nj(j,n,m,t)*T_hat(j,m,t);%-F(n,j,i,m,o,t)*THETA(j)*kappa_hat(o+(j-1)*N,m,t) ;
%                    for o=1:N
%                        Pi_temp(n+(j-1)*N,i,m,o,t) = C(i,j) * w_hat(j,i,t) + D(n,j,i,m,t)*w_hat(j,m,t)...
%                            + T_hat(j,i,t) - THETA(j)* kappa_hat(n+(j-1)*N,i,t) + E(n,j,i,m,t)*T_hat(j,m,t)-F(n,j,i,m,o,t)*THETA(j)*kappa_hat(o+(j-1)*N,m,t) ;
%                    end
                end
             end
        end
    end
end
%pi_nj= sum(sum(Pi_temp,3),4);
pi_nj= (sum(Pi_temp,3));

pi_nj %pi from the new code
pi_hat %pi from the old code

%% Parts in X (expenditure)
U = zeros(N*J,N*J,TIME);
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
G = zeros(N,J,N,N,TIME);
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
H = zeros(N,J,N,J,TIME);
for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                for k=1:J
                    H(i,j,o,k,t) = BGAMMA(i+(j-1)*N,o+(j-1)*N,t) * (ALPHAS(j,o))* zeta(o+(k-1)*N,j,t);
                end
            end
        end
    end
end

%M(ijmj) = GAMMA(ij) sum_n chi(nji)(D(njim)+sum_o sum_l G(njol)D(ljom))
GD_temp = zeros(N,J,N,TIME,N,N);
M_temp = zeros(N*J,N*J,TIME,N);
for t=1:TIME
    for n=1:N
        for j=1:J
            for m=1:N
                for o=1:N
                    for l=1:N
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
Q_temp = zeros(N*J,N*J,TIME,N);
for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
%               Q_temp(i+(j-1)*N,i+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t) * (C(i,j,t) + GC_sum(n,j,t));
                Q_temp(i+(j-1)*N,i+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t) * C(i,j,t);
            end
        end
    end
end
Q = sum(Q_temp,4);

%R(ijoj) = GAMMA(ij)sum_n chi(nji) sum_l G(njol)C(oj)
GC_temp = zeros(N,J,N,TIME,N);
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
R_temp = zeros(N*J,N*J,TIME,N);
for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
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
E_STAR_temp = zeros(N*J,TIME,N,N,N);
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
    w_hat_nj(:,t) = (eye(N*J) - M(:,:,t) - Q(:,:,t)- R(:,:,t) - T(:,:,t) ) \ ( (T-eye(N*J))*L_hat(:,t) + E_STAR(:,t));
end
w_hat
w_hat_nj