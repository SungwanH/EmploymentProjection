function mat_pbp = MAT(params, approx)
% This function provides matrices which are used in deriving wage
% using matrix inversion in temporary equilibrium
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

%This is used in deriving p

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
D_temp = zeros(N,J,N,N,TIME,N);
E_temp = zeros(N,J,N,N,TIME,N);
F = zeros(N,J,N,N,N,TIME);
for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                for o=1:N
                    for m=1:N
                        D_temp(n,j,i,m,t,o) = -THETA(j) * ((1-GAMMA(j,i))*DELTA(i+(j-1)*N,o+(j-1)*N,t)-DELTA(n+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t) * GAMMA(j,m);
                        E_temp(n,j,i,m,t,o) = ((1-GAMMA(j,i))* DELTA(i+(j-1)*N,o+(j-1)*N,t) - DELTA(n+(j-1)*N,o+(j-1)*N,t))* pi(o+(j-1)*N,m,t);
                        F(n,j,i,m,o,t) = THETA(j) * ((1-GAMMA(j,i))* DELTA(i+(j-1)*N,o+(j-1)*N,t) - DELTA(n+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t);
                    end
                end
            end
        end
    end
end
D = sum(D_temp,6);
E = sum(E_temp,6);
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


%% Recovering labor market clearing condition

%M1(ijij) = GAMMA(ij)sum_n chi(nji) C(ij) sum_o sum_l G(njol)C(oj)
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
%M2(ijmj) = GAMMA(ij) sum_n chi(nji)(D(njim)+sum_o sum_l G(njol)D(ljom))
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

%M3(ijoj) = GAMMA(ij)sum_n chi(nji) sum_l G(njol)C(oj)
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
                                F4_temp(i+(j-1)*N,h+(m-1)*N+(j-1)*N*N,t,n,o,l) = -GAMMA(j,i)*chi(n+(j-1)*N,i,t)*G(n,j,o,l,t)*F(l,j,o,m,h);
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

mat_pbp = v2struct(M_tilde, Q_tilde, F_tilde, T, DELTA, G, H);
end