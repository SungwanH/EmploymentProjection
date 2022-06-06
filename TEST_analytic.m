%This is a test code for computing the temporary equilibrium analyticially

clear 
close all
clc;
load('DATA/DGP.mat', 'eqm_dgp','approx_dgp');
load('DATA/NLPF_DD.mat', 'eqm_nlpf_dd','approx_nlpf_dd');
%eqm_nlpf_dd, approx_nlpf_dd)
%data_4_sector=load('DATA/BASE_FOURSECTOR.mat', 'VALjn00', 'Din00','mu0','L0');
%{
load('DATA_onesector/DGP.mat', 'eqm_dgp','approx_dgp');
%% Load parameter values and approximation points
N=10;
J=1;
%}
%T_hat = eqm_dd.E_T_hat;
pi=approx_nlpf_dd.pi;
chi=approx_nlpf_dd.chi;
zeta=approx_nlpf_dd.zeta;
varrho=approx_nlpf_dd.varrho;

T_hat = eqm_dgp.E_T_hat;
pi=approx_dgp.pi;
chi=approx_dgp.chi;
zeta=approx_dgp.zeta;
varrho=approx_dgp.varrho;
%}
params=PARAMS_TEST(0.5);
v2struct(params.envr);
THETA  = params.modl.THETA;
GAMMA  = params.modl.GAMMA;
ALPHAS = params.modl.ALPHAS;
clear approx_dgp VALjn00 Din00 mu0 L0
%}
rng(20220605)
TIME=2;
N=5;
J=2;
%{
pi = rand(N*J,N,TIME);
varrho = rand(N*J,N,TIME);
chi = rand(N*J,N,TIME);
zeta = ones(N*J,J,TIME) * 1/J;
for t=1:TIME
    for i=1:N*J
        pi(i,:,t) = pi(i,:,t)./sum(pi(i,:,t)); % should sum up to 1
        varrho(i,:,t) = varrho(i,:,t)./sum(varrho(i,:,t));% should sum up to 1
    end
end
ALPHAS=1/J*ones(J,N);
THETA=4*ones(J,1);
%GAMMA=0.5*ones(J,N);
GAMMA=rand(J,N);
%}
%pi = rand(N*J,N,TIME);
%% Shocks
%T_hat = zeros(J,N,TIME);
T_hat = rand(J,N,TIME)*0.1;
%T_hat(1,1,1:TIME) = 0.1;
%T_hat = E_T_hat(1:J,1:N,2,1);
L_hat = rand(J,N,TIME)*0.1;
%L_hat = zeros(J,N,TIME);
%kappa_hat = zeros(N*J,N,TIME);
kappa_hat = rand(N*J,N,TIME)*0.1;
w_hat_test = rand(J,N,TIME)*0.1; %For test


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
                Psum_temp(o,j,i,t)= pi(o+(j-1)*N,i,t)*(GAMMA(j,i)*w_hat_test(j,i,t)+ kappa_hat(o+(j-1)*N,i,t) - (1/THETA(j)) * T_hat(j,i,t));
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
        w_temp = w_hat_test(1:J,1:N,t);
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

% Pi_hat
for t = 1:TIME
    for n=1:N
        for j=1:J
            for i=1:N
%               pi_hat(n+(j-1)*N,ii,t) = -THETA(j)*(B(j,ii) * w_hat(j,ii,t) - (1-B(j,ii)*p_hat(j,ii,t)) + TAU(n+(j-1)*N,ii)) + T_hat(j,ii,t);
                pi_hat(n+(j-1)*N,i,t) = -THETA(j)*(GAMMA(j,i) .* w_hat_test(j,i,t) + (1-GAMMA(j,i)).*p_hat(j,i,t) - p_hat(j,n,t) + kappa_hat(n+(j-1)*N,i,t)) + T_hat(j,i,t);
             end
         end
     end
end        
p_hat_nj % p from the new code
p_hat % p from the old code

% Step 4d. Solve for total expenditure
for t = 1:TIME
    RHS = NaN(J,N);
    RHS1=NaN(J,N,N);
    RHS2=NaN(J,N,J);
    
    varrho_aux = varrho(:,:,t);
    pi_temp = pi_hat(:,:,t);
    zeta_aux = zeta(:,:,t);
    w_temp = w_hat_test(:,:,t);
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

for t = 1:TIME
    RHS_temp=NaN(J,N,N);
    
    chi_aux = chi(:,:,t);
    X_temp = X_hat(:,:,t);
    L_temp = L_hat(:,:,t);
    pi_temp = pi_hat(:,:,t);
    for ii=1:N
        for j=1:J
            for n=1:N
                RHS_temp_old(j,ii,n) = chi_aux(n+(j-1)*N, ii) * (pi_temp(n+(j-1)*N, ii) + X_temp(j,n));
            end
        end
    end
end

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
F = zeros(N,J,N,N,N,TIME);
for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                for o=1:N
                    for m=1:N
                        D_temp(n,j,i,m,o,t) = -THETA(j) * ((1-GAMMA(j,i))*DELTA(i+(j-1)*N,o+(j-1)*N,t)-DELTA(n+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t) * GAMMA(j,m);
                        E_temp(n,j,i,m,o,t) = ((1-GAMMA(j,i))* DELTA(i+(j-1)*N,o+(j-1)*N,t) - DELTA(n+(j-1)*N,o+(j-1)*N,t))* pi(o+(j-1)*N,m,t);
                        F(n,j,i,m,o,t) = THETA(j) * ((1-GAMMA(j,i))* DELTA(i+(j-1)*N,o+(j-1)*N,t) - DELTA(n+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t);
%                        D_temp_ij(j,i,m,o,t) = -THETA(j) * ((1-GAMMA(j,i))*DELTA(i+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t) * GAMMA(j,m);
%                        D_temp_nj(j,n,m,o,t) = -THETA(j) * (-DELTA(n+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t) * GAMMA(j,m);
%                        E_temp_ij(j,i,m,o,t) = ((1-GAMMA(j,i))* DELTA(i+(j-1)*N,o+(j-1)*N,t))* pi(o+(j-1)*N,m,t);
%                        E_temp_nj(j,n,m,o,t) = (- DELTA(n+(j-1)*N,o+(j-1)*N,t))* pi(o+(j-1)*N,m,t);
                    end
                end
            end
        end
    end
end
   
D = sum(D_temp,5);%
E = sum(E_temp,5);
%D_ij = sum(D_temp_ij,4);
%D_nj = sum(D_temp_nj,4);
%E_ij = sum(E_temp_ij,4);
%E_nj = sum(E_temp_nj,4);
%F = E_temp

%% Test (new code) Pi
%T_hat_temp =zeros(J,N,TIME);

Pi_temp=zeros(N*J,N,TIME);
for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                for m=1:N
                    for o=1:N
                       FK_temp(n,j,i,t,m,o) = F(n,j,i,m,o,t) * kappa_hat(o+(j-1)*N,m,t);  
                    end
                end
            end
        end
    end
end
FK = sum(sum(FK_temp,6),5);             


for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                temptemp=0;
                for m=1:N
                    temptemp= temptemp+D(n,j,i,m,t)*w_hat_test(j,m,t) + E(n,j,i,m,t)*T_hat(j,m,t);
                end
                Pi_temp(n+(j-1)*N,i,t) = C(i,j) * w_hat_test(j,i,t) +  T_hat(j,i,t) - THETA(j)* kappa_hat(n+(j-1)*N,i,t) + temptemp - FK(n,j,i,t);  
            end
        end
    end
end
pi_nj=Pi_temp;
pi_nj %pi from the new code
pi_hat %pi from the old code

%% Parts in X (expenditure)
U = zeros(N*J,N*J,TIME);
for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                U(i+(j-1)*N,o+(j-1)*N,t) = (1-GAMMA(j,i))*varrho(o+(j-1)*N,i,t); 
            end
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
%H(ijok) = BGAMMA(ijo)(ALPHA(oj)chi(okj))
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
GC_temp = zeros(N,J,N,N,TIME);
for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                for l=1:N
                    temptemp=0;
                    for m=1:N
                        temptt=0;
                        for h=1:N
                            temptt = temptt - G(i,j,o,l,t)*F(l,j,o,m,h,t)*kappa_hat(h+(j-1)*N,m,t);
                        end
                        temptemp = temptemp + G(i,j,o,l,t)*(D(l,j,o,m,t)*w_hat_test(j,m,t) + E(l,j,o,m,t) *T_hat(j,m,t)) + temptt;
                    end
                    Gsum_temp(j,i,o,l,t) = G(i,j,o,l,t)*(C(o,j)*w_hat_test(j,o,t) + T_hat(j,o,t) - THETA(j)*kappa_hat(l+(j-1)*N,o,t)) + temptemp;
                end
            end
        end
    end
end
Gsum = sum(sum(Gsum_temp,3),4);

for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                for k=1:J
                    Hsum_temp(j,i,o,k,t) = H(i,j,o,k,t) * (w_hat_test(k,o,t)+L_hat(k,o,t));
                end
            end
        end
    end
end
Hsum = sum(sum(Hsum_temp,4),3);

X_hat_ij = Gsum + Hsum;

%X_hat_ij %New code
X_hat    %old code

for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                for l=1:N
                    X_temp_G(j,i,o,l,t)=G(i,j,o,l,t)*pi_nj(l+(j-1)*N,o,t);
                end
                for k=1:J
                    X_temp_H(j,i,o,k,t)=H(i,j,o,k,t)*(w_hat_test(k,o,t)+L_hat(k,o,t));
                end
             end
         end
     end
end        
X_hat_test = sum(sum(X_temp_G,3),4) + sum(sum(X_temp_H,3),4)
%% Recovering labor market clearing condition

%Q(ijij) = GAMMA(ij)sum_n chi(nji) C(ij) sum_o sum_l G(njol)C(oj)
M1_temp = zeros(N*J,N*J,TIME,N);
for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                M1_temp(i+(j-1)*N,i+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t) * C(i,j);
            end
        end
    end
end
M1 = sum(M1_temp,4);
%M(ijmj) = GAMMA(ij) sum_n chi(nji)(D(njim)+sum_o sum_l G(njol)D(ljom))
GD_temp = zeros(N,J,N,TIME,N,N);
M2_temp = zeros(N*J,N*J,TIME,N);
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
    for n=1:N
        for j=1:J
            for m=1:N
                for i=1:N
                    M2_temp(i+(j-1)*N,m+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t) * (D(n,j,i,m,t) + GD_sum(n,j,m,t));
                end
            end
        end
    end
end
M2 = sum(M2_temp,4);


%R(ijoj) = GAMMA(ij)sum_n chi(nji) sum_l G(njol)C(oj)
GC_temp = zeros(N,J,N,TIME,N);
for t=1:TIME
    for n=1:N
        for j=1:J
            for o=1:N
                for l=1:N
                    GC_temp(n,j,o,t,l) = G(n,j,o,l,t)*C(o,j);
                end
            end
        end
    end
end
GC_sum = sum(GC_temp,5);
M3_temp = zeros(N*J,N*J,TIME,N);
for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                for n=1:N
                    M3_temp(i+(j-1)*N,o+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t) * GC_sum(n,j,o,t);
                end
            end
        end
    end
end
M3 = sum(M3_temp,4);

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

Q1_temp = zeros(N*J,N*J,TIME,N);
F1 = zeros(N*J, N*N*J,TIME);
for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                Q1_temp(i+(j-1)*N,i+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t);
                F1(i+(j-1)*N,n+(i-1)*N+(j-1)*N*N,t)  = -THETA(j)*GAMMA(j,i)*chi(n+(j-1)*N,i,t);
            end
        end
    end
end
Q1= sum(Q1_temp,4);

Q2_temp = zeros(N*J,N*J,TIME,N);
F2_temp = zeros(N*J,N*N*J,TIME,N);
for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                for m=1:N
                    Q2_temp(i+(j-1)*N,m+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t)*E(n,j,i,m,t);
                    for h=1:N
                        F2_temp(i+(j-1)*N,h+(m-1)*N+(j-1)*N*N,t,n) = -GAMMA(j,i)*chi(n+(j-1)*N,i,t)*F(n,j,i,m,h,t);
                    end
                end
            end
        end
    end
end
Q2= sum(Q2_temp,4);
F2= sum(F2_temp,4);


Q3_temp = zeros(N*J,N*J,TIME,N,N);
F3_temp = zeros(N*J,N*N*J,TIME,N);
Q4_temp = zeros(N*J,N*J,TIME,N,N,N);
F4_temp = zeros(N*J,N*J,TIME,N,N,N);
for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                for o=1:N
                    for l=1:N
                        Q3_temp(i+(j-1)*N,o+(j-1)*N,t,n,l) = GAMMA(j,i)*chi(n+(j-1)*N,i,t)*G(n,j,o,l,t);
                        F3_temp(i+(j-1)*N,l+(o-1)*N+(j-1)*N*N,t,n) = -THETA(j)*GAMMA(j,i)*chi(n+(j-1)*N,i,t)*G(n,j,o,l,t);
                        for m=1:N
                            Q4_temp(i+(j-1)*N,m+(j-1)*N,t,n,o,l) = GAMMA(j,i)*chi(n+(j-1)*N,i,t)*G(n,j,o,l,t)*E(l,j,o,m,t);
                            for h=1:N
                                F4_temp(i+(j-1)*N,h+(m-1)*N+(j-1)*N*N,t,n,o,l) = -GAMMA(j,i)*chi(n+(j-1)*N,i,t)*G(n,j,o,l,t)*F(l,j,o,m,h,t);
                            end
                        end
                    end
                end
            end
        end
    end
end
Q3= sum(sum(Q3_temp,4),5);
F3= sum(F3_temp,4);
Q4= sum(sum(sum(Q4_temp,4),5),6);
F4= sum(sum(sum(F4_temp,4),5),6);

% Final equation:
M_tilde = M1+M2+M3;
Q_tilde = Q1+Q2+Q3+Q4;
F_tilde = F1+F2+F3+F4;
% Reformulate the order of shocks
for t=1:TIME
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
for t=1:TIME
    M_tilde_temp = M_tilde(:,:,t);
    T_tilde_temp = T(:,:,t);
    Q_tilde_temp = Q_tilde(:,:,t);
    F_tilde_temp = F_tilde(:,:,t);
    L_hat_T_temp = L_hat_T(:,t);
    T_hat_T_temp = T_hat_T(:,t);
    kappa_hat_T_temp = kappa_hat_T(:,t);
    w_hat_nj(:,t) = (eye(N*J) - M_tilde_temp - T_tilde_temp ) \ ((T_tilde_temp-eye(N*J)) * L_hat_T_temp + Q_tilde_temp * T_hat_T_temp + F_tilde_temp * kappa_hat_T_temp);
%        w_hat_nj(:,t) = (eye(N*J) - M_tilde(:,:,t) - T(:,:,t) ) \ ((T(:,:,t)-eye(N*J)) * L_hat_T(:,t) + Q_tilde(:,:,t) * T_hat_T(:,t) + F_tilde(:,:,t) * kappa_hat_T(:,t));
end
w_hat_nj
%% check: compare RHS of Labor mkt clearing condition from old vs new code
% old code; pi_hat and X_hat is coming from the old code
for t=1:TIME
    for i=1:N
        for n=1:N
             for j=1:J
                RHS_old_temp(i+(j-1)*N,t,n) = GAMMA(j,i) * chi(n+(j-1)*N,i,t) * (pi_hat(n+(j-1)*N,i,t)+X_hat(j,n,t)); 
            end
        end
    end
end
%for t=1:TIME
%    RHS_old(:,t) = sum(RHS_temp(:,t,:),3);
%end
RHS_old = sum(RHS_old_temp,3);
RHS_old
% new code
for t=1:TIME
    for i=1:N
        for j=1:J
            w_hat_test_T(i+(j-1)*N,t) = w_hat_test(j,i,t);
        end
    end
end
for t=1:TIME
    RHS_new(:,t) = M_tilde(:,:,t) * w_hat_test_T(:,t) + T(:,:,t) * (w_hat_test_T(:,t) + L_hat_T(:,t)) + Q_tilde(:,:,t) * T_hat_T(:,t) + F_tilde(:,:,t) * kappa_hat_T(:,t);
end
RHS_new

%% Find w_hat_iter from iterative method (old code)
ITER_TEMP = 1;
wmax=1;
% Step 4a. Initial guess for the w_tilde_hat
%w_hat_iter = zeros(J,N,TIME);
w_hat_iter = zeros(J,N,TIME);
w_update   = zeros(J,N,TIME);
p_hat      = zeros(J,N,TIME);
pi_hat     = zeros(N*J,N,TIME);
X_hat      = zeros(J,N,TIME);
UPDT_W     = 0.1; %update speed for wage loop (lower value->conservative)
TOLTEMP    = 1E-15;  % tolerance rate for linear temporary equilibrium
MAXIT      = 1E+8; %maximum number of iterations
while (ITER_TEMP <= MAXIT) && (wmax > TOLTEMP)

    for t = 1:TIME
        % Step 4b. solve for p_hat
        RHS = zeros(J,N);
    
        pi_aux =  pi(:,:,t); %Note that pi and kappa are organized by sector; i.e., two consecutive rows are same sectors in different countries
        w_temp = w_hat_iter(:,:,t);
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
%toc
    %price index
    for t=1:TIME
        P_hat(:,:,t) = sum((ALPHAS.*p_hat(:,:,t)),1);
    end
    % Step 4c. solve for pi_hat
    for t = 1:TIME
        for n=1:N
            for j=1:J
                pi_hat(n+(j-1)*N,:,t) = -THETA(j)*(GAMMA(j,:) .* w_hat_iter(j,:,t) + (1-GAMMA(j,:)).*p_hat(j,:,t) - p_hat(j,n,t) + kappa_hat(n+(j-1)*N,:,t)) + T_hat(j,:,t);
             end
         end
    end        

    % Step 4d. Solve for total expenditure
    for t = 1:TIME
        RHS = NaN(J,N);
        RHS1=NaN(J,N,N);
        RHS2=NaN(J,N,J);
        
        varrho_aux = varrho(:,:,t);
        pi_temp = pi_hat(:,:,t);
        zeta_aux = zeta(:,:,t);
        w_temp = w_hat_iter(:,:,t);
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
    % Step 4e. Solve for an updated w_hat_iter using the labor market clearing
    for t = 1:TIME
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

    for t=1:TIME
        checkw(t,1)=max(max(abs(w_hat_iter(:,:,t)-w_update(:,:,t))));
    end
    [wmax loc]=max(checkw);
    wmax;
    if ITER_TEMP >10000 || sum(isnan(w_update(:)))>0
        checkw(1:TIME)
        1
        disp('Inner loop err')
        ITER_TEMP
        stop
    end
    
    % Step 4e. update the guess
    w_hat_iter = (1-UPDT_W) * w_hat_iter + UPDT_W * w_update;
    ITER_TEMP = ITER_TEMP + 1;
end

for t=1:TIME
    for i=1:N
        for j=1:J
            w_hat_iter_old(i+(j-1)*N,t) = w_hat_iter(j,i,t);
        end
    end
end

%Comparison 
w_hat_nj(:,1) - w_hat_nj(1,1) % wage from new code using inversion
%w_hat_nj(:,2) - w_hat_nj(1,2)
w_hat_iter_old  % wage from old code using iteration

%% Iteration with new code
w_hat_iter_new = zeros(N*J,TIME);
w_update_new = zeros(N*J,TIME);
ITER_TEMP_NEW = 0;
wmax_new=1;
while (ITER_TEMP_NEW <= MAXIT) && (wmax_new > TOLTEMP)
for t=1:TIME
    w_update_new(:,t) = M_tilde(:,:,t) * w_hat_iter_new(:,t) + T(:,:,t) * (w_hat_iter_new(:,t) + L_hat_T(:,t)) + Q_tilde(:,:,t) * T_hat_T(:,t) + F_tilde(:,:,t) * kappa_hat_T(:,t) - L_hat_T(:,t);
    w_update_new(:,t)=w_update_new(:,t)-w_update_new(1,t); %normalize the first wage to be a constant across periods; should not change anything
end

    for t=1:TIME
        checkw(t,1)=max(abs(w_hat_iter_new(:,t)-w_update_new(:,t)));
    end
    [wmax_new loc]=max(checkw);
    wmax_new;
    if ITER_TEMP_NEW >10000 || sum(isnan(w_update_new(:)))>0
      checkw(1:TIME)
      1 
      disp('Inner loop err')
      ITER_TEMP_NEW
      stop
    end
w_hat_iter_new = (1-UPDT_W) * w_hat_iter_new + UPDT_W * w_update_new;
ITER_TEMP_NEW = ITER_TEMP_NEW+1;
end
w_hat_iter_new % wage from new code using iteration