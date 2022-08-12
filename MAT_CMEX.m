function mat_pbp = MAT_CMEX(params, approx)
% This function provides matrices which are used in deriving wage
% using matrix inversion in temporary equilibrium

%% Roll down parameters
v2struct(params.envr);
v2struct(params.tech);
v2struct(params.modl);

%% Set the number of cores
N_THREAD=6;



%% Roll down approximation points
varrho  =   approx.varrho;
zeta    =   approx.zeta;
chi     =   approx.chi;
pi      =   approx.pi;


%% The following block is to be replaced by STEP1
BTHETA = zeros(N*J,N*J,TIME);
U = zeros(N*J,N*J,TIME);

% pi(n*j,i): country n's expenditure share on good j from country i
% sum(pi(1,:))=1

%{
for t=1:TIME
    for i=1:N
        for n=1:N
            for j=1:J
                BTHETA(n+(j-1)*N,i+(j-1)*N,t)= pi(n+(j-1)*N,i,t)*(1-GAMMA(j,i));
            end
        end
    end
end

for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                U(i+(j-1)*N,o+(j-1)*N,t) = (1-GAMMA(j,i))*varrho(o+(j-1)*N,i,t); 
            end
        end
    end
end
%}

MAT_CPP_S1

%% Back to Matlab for Leontief inverse
DELTA = zeros(N*J,N*J,TIME);

for t=1:TIME
    DELTA(:,:,t) = inv(eye(N*J)-BTHETA(:,:,t));
end

BGAMMA = zeros(N*J,N*J,TIME);
for t=1:TIME
    BGAMMA(:,:,t) = inv(eye(N*J)-U(:,:,t));
end





%% Stacking in C, Step II
C= zeros(N,J);
G = zeros(N,J,N,N,TIME);

%D_temp = zeros(N,J,N,N,TIME,N);
%E_temp = zeros(N,J,N,N,TIME,N);
D=zeros(N,J,N,N,TIME);
E=zeros(N,J,N,N,TIME);
F = zeros(N,J,N,N,N,TIME);

%GD_temp = zeros(N,J,N,TIME,N,N);
GD_sum=zeros(N,J,N,TIME);

%T_tilde_temp = zeros(N*J,N*J,TIME,N);
T_tilde = zeros(N*J,N*J,TIME);


%GC_temp = zeros(N,J,N,TIME,N);
GC_sum = zeros(N,J,N,TIME);

H = zeros(N,J,N,J,TIME);


%M1_temp = zeros(N*J,N*J,TIME,N);
%M2_temp = zeros(N*J,N*J,TIME,N);
%M3_temp = zeros(N*J,N*J,TIME,N);
M1 = zeros(N*J,N*J,TIME);
M2 = zeros(N*J,N*J,TIME);
M3 = zeros(N*J,N*J,TIME);


%Q1_temp = zeros(N*J,N*J,TIME,N);
%Q2_temp = zeros(N*J,N*J,TIME,N);
%Q3_temp = zeros(N*J,N*J,TIME,N,N);
%Q4_temp = zeros(N*J,N*J,TIME,N,N,N);

Q1=  zeros(N*J,N*J,TIME);
Q2= zeros(N*J,N*J,TIME);
Q3=zeros(N*J,N*J,TIME);
Q4= zeros(N*J,N*J,TIME);

%F2_temp = zeros(N*J,N*N*J,TIME,N);
%F3_temp = zeros(N*J,N*N*J,TIME,N);
%F4_temp = zeros(N*J,N*J*N,TIME,N,N,N);


F = zeros(N,J,N,N,N,TIME);
F1 = zeros(N*J, N*N*J,TIME);
F2= zeros(N*J,N*N*J,TIME); % F2: N*J by N*N*J by TIME
F3=zeros(N*J,N*N*J,TIME);
F4=zeros(N*J,N*J*N,TIME);




%{
for i=1:N
    for j=1:J
        C(i,j) = -THETA(j) * GAMMA(j,i);
    end
end

for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                for m=1:N
                    for o=1:N
                        D_temp(n,j,i,m,t,o) = -THETA(j) * ((1-GAMMA(j,i))*DELTA(i+(j-1)*N,o+(j-1)*N,t)-DELTA(n+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t) * GAMMA(j,m);
                        E_temp(n,j,i,m,t,o) = ((1-GAMMA(j,i))* DELTA(i+(j-1)*N,o+(j-1)*N,t) - DELTA(n+(j-1)*N,o+(j-1)*N,t))* pi(o+(j-1)*N,m,t);
                        F(n,j,i,m,o,t) = THETA(j) * ((1-GAMMA(j,i))* DELTA(i+(j-1)*N,o+(j-1)*N,t) - DELTA(n+(j-1)*N,o+(j-1)*N,t)) * pi(o+(j-1)*N,m,t);
                        D(n,j,i,m,t)=D(n,j,i,m,t)+D_temp(n,j,i,m,t,o); 
                        E(n,j,i,m,t)=E(n,j,i,m,t)+E_temp(n,j,i,m,t,o); 
                    end
                end
            end
        end
    end
end

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



for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                M1_temp(i+(j-1)*N,i+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t) * C(i,j);
            end
        end
    end
end





for t=1:TIME
    for n=1:N
        for j=1:J
            for m=1:N
                for o=1:N
                    for l=1:N
                        GD_temp(n,j,m,t,o,l) = G(n,j,o,l,t)*D(l,j,o,m,t);
                        GD_sum(n, j, m, t) = GD_sum(n, j, m, t) + GD_temp(n, j, m, t, o, l);                        
                    end
                end
            end
        end
    end
end
%GD_sum = sum(sum(GD_temp,6),5);


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


%M3(ijoj) = GAMMA(ij)sum_n chi(nji) sum_l G(njol)C(oj)
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
%GC_sum = sum(GC_temp,5);

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


%T_tilde(ijok) = GAMMA(ij)sum_n chi(nji) (sum_i sum_k H(njok)) sum_i H(njok)

for t=1:TIME
    for i=1:N
        for j=1:J
            for o=1:N
                for k=1:J
                    for n=1:N
                        T_tilde_temp(i+(j-1)*N,o+(k-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t)*H(n,j,o,k,t);
                    end
                end
            end
        end
    end
end



for t=1:TIME
    for i=1:N
        for j=1:J
            for n=1:N
                Q1_temp(i+(j-1)*N,i+(j-1)*N,t,n) = GAMMA(j,i)*chi(n+(j-1)*N,i,t);
                F1(i+(j-1)*N,n+(i-1)*N+(j-1)*N*N,t)  = -THETA(j)*GAMMA(j,i)*chi(n+(j-1)*N,i,t); % F1: N*J by N*N*J by TIME
            end
        end
    end
end


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


%F2= sum(F2_temp,4); % F2: N*J by N*N*J by TIME
%M1 = sum(M1_temp,4); %M1: N*J by N*J by TIME 
%M2 = sum(M2_temp,4); %M2: N*J by N*J by TIME  
%M3 = sum(M3_temp,4); %M3: N*J by N*J by TIME 
%Q1= sum(Q1_temp,4); % Q1: N*J by N*J by TIME
%Q2= sum(Q2_temp,4); % Q2: N*J by N*J by TIME
%T_tilde= sum(T_tilde_temp,4); %T_Tilde: N*J by N*J by TIME 

%Q3= sum(sum(Q3_temp,4),5);        % Q3: N*J by N*J by TIME
%F3= sum(F3_temp,4);               % F3: N*J by N*N*J by TIME 

%Q4= sum(sum(sum(Q4_temp,4),5),6); % Q4: N*J by N*J by TIME
%F4= sum(sum(sum(F4_temp,4),5),6); % F4: N*J by N*N*J by TIME 
%}


MAT_CPP_S2




%% EXIT LOOP HERE
% Final equation:
M_tilde = M1+M2+M3;
Q_tilde = Q1+Q2+Q3+Q4;
F_tilde = F1+F2+F3+F4;

mat_pbp = v2struct(M_tilde, Q_tilde, F_tilde, T_tilde, DELTA, G, H);
end