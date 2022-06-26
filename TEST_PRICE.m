%This is a test code for comparing first-order approximation / second-order
%approximation with the nonlinear outcome.

clear 
close all
clc;

rng(20220620)
N=5; % Since I assume kappa =1 for all regions, increasing number of regions does not give variances
J=50;

ALPHAS=1/J*ones(J,N);
THETA=5*ones(J,1);  
GAMMA=max(rand(J,N),0.1);
%GAMMA=ones(J,N)*1/2;
%% Shocks
T = rand(J,N); % baseline productivity
%T_prime = rand(J,N); %counterfactual productivity
T_prime = T*1.5;
T_prime(1:J-5,1) = T(1:J-5,1)*0.9;
T_prime(1:J-5,2) = T(1:J-5,2)*0.5;
T_hat = log(T_prime) - log(T); %log difference
w = rand(J,N);
w_prime = w; % Assume no wage response
w_hat = log(w_prime) - log(w); % = zero
%w_hat = zeros(J,N);
kappa=rand(N*J,N)+1;
kappa=ones(N*J,N);

%% Level value: BASELINE
ITER_LEV = 1;
pmax_lev = 1;
p_update = zeros(J,N);
MAXIT = 1e+8;
p_old = ones(J,N);
TOLTEMP = 1e-14;
while (ITER_LEV <= MAXIT) && (pmax_lev > TOLTEMP)    
    for n=1:N
        for j=1:J
            for i=1:N
                temp(j,n,i) = (w(j,i)^GAMMA(j,i) * p_old(j,i)^(1 - GAMMA(j,i))*kappa(n+(j-1)*N,i))^(-THETA(j)) * T(j,i);
            end
        end
    end
    p_update = (sum(temp,3)).^(-1/THETA(j));
    pmax_lev=max(max(abs(p_old(:,:)-p_update(:,:))));
    pmax_lev

    p_old = p_update;
    ITER_LEV = ITER_LEV+1;
    
end

for i=1:N
    for j=1:J
        for n=1:N
            pi(n+(j-1)*N,i) = ( (w(j,i)^GAMMA(j,i) * p_old(j,i)^(1 - GAMMA(j,i))*kappa(n+(j-1)*N,i))/p_old(j,n))^(-THETA(j)) * T(j,i);
        end
    end
end

%% Level value: Counterfactual
ITER_LEV_prime = 1;
pmax_lev_prime = 1;
p_update_prime = zeros(J,N);
p_old_prime = ones(J,N);
while (ITER_LEV_prime <= MAXIT) && (pmax_lev_prime > TOLTEMP)
    
    
    for i=1:N
        for j=1:J
            for n=1:N
                temp(j,n,i) = ((w(j,i)^GAMMA(j,i) * p_old_prime(j,i)^(1 - GAMMA(j,i))*kappa(n+(j-1)*N,i))^(-THETA(j)) * T_prime(j,i));
            end
        end
    end
    p_update_prime = sum(temp,3).^(-1/THETA(j));

    pmax_lev_prime=max(max(abs(p_old_prime(:,:)-p_update_prime(:,:))));

    p_old_prime = p_update_prime;
    ITER_LEV_prime = ITER_LEV_prime+1;

end
% compute pi_prime
for i=1:N
    for j=1:J
        for n=1:N
            pi_prime(n+(j-1)*N,i) = ( (w(j,i)^GAMMA(j,i) * p_old_prime(j,i)^(1 - GAMMA(j,i))*kappa(n+(j-1)*N,i))/p_old_prime(j,n))^(-THETA(j)) * T_prime(j,i);
        end
    end
end
%save the outcome
p_lev = p_old; % baseline
p_prime_lev = p_old_prime; %counterfactual
p_hat_lev = p_prime_lev./p_lev; %hat

%% First-order (FO)
ITER_FO = 1;
pmax_FO = 1;
p_hat_update = zeros(J,N);
p_hat_old = zeros(J,N);
while (ITER_FO <= MAXIT) && (pmax_FO > TOLTEMP)
    for i=1:N
        for j=1:J
            for n=1:N
                p_hat_temp(j,n,i) = pi(n+(j-1)*N,i) * ( GAMMA(j,i) * w_hat(j,i) + (1 - GAMMA(j,i) ) * p_hat_old(j,i) - (1/THETA(j)) * T_hat(j,i));
            end
        end
    end
    p_hat_update = sum(p_hat_temp,3);

    pmax_FO=max(max(abs(p_hat_old(:,:)-p_hat_update(:,:))));
    pmax_FO

    p_hat_old = p_hat_update;
    ITER_FO = ITER_FO+1;

end
%save the outcome
p_hat_FO = exp(p_hat_old);
p_FO = p_lev .* exp(p_hat_old);


%% Second-order(SO)
ITER_SO = 1;
pmax_SO = 1;
p_hat_old_SO=p_hat_lev;%p_hat_old;%zeros(J,N);
P_ww = zeros(J,N,N,N); P_pp = zeros(J,N,N,N); P_TT = zeros(J,N,N,N);
P_wT = zeros(J,N,N,N); P_pT = zeros(J,N,N,N); P_wp = zeros(J,N,N,N);
while (ITER_SO <= MAXIT) && (pmax_SO > TOLTEMP)
    for k=1:N
        for n=1:N
            for j=1:J
                %diff w.r.t. w
                P_w(j,n,k) = GAMMA(j,k) * pi(n+(j-1)*N,k) * w_hat(j,k);
                %diff w.r.t. p
                P_p(j,n,k) = (1-GAMMA(j,k)) * pi(n+(j-1)*N,k) * p_hat_old_SO(j,k);
                %diff w.r.t. T
                P_T(j,n,k) = -(1/THETA(j))* pi(n+(j-1)*N,k) * T_hat(j,k);
                %P_pp1(j,n,k)= -THETA(j) * (1 -GAMMA(j,k))^2 * pi(n+(j-1)*N,k) * p_hat_old_SO(j,k) * p_hat_old_SO(j,k);
                %P_pp2(j,n,k)= THETA(j) * (1-GAMMA(j,k)) * pi(k+(j-1)*N,k) * p_hat_old_SO(j,k) * p_hat_old_SO(j,k);
                for m=1:N
                    %diff w.r.t w and w
                    P_ww(j,n,k,m) = ((m==k) * GAMMA(j,k) * (-THETA(j)) * GAMMA(j,k) * pi(n+(j-1)*N,k) - (-THETA(j)) * pi(n+(j-1)*N,k) * pi(n+(j-1)*N,m) * GAMMA(j,k) * GAMMA(j,m)) * w_hat(j,k) * w_hat(j,m);
                   %diff w.r.t p and p
                    P_pp(j,n,k,m) = ((m==k) * (-THETA(j)) * (1-GAMMA(j,k)) * (1-GAMMA(j,k)) * pi(n+(j-1)*N,k) - (-THETA(j)) * (1-GAMMA(j,k)) * (1-GAMMA(j,m)) * pi(n+(j-1)*N,k) * pi(n+(j-1)*N,m)) * p_hat_old_SO(j,k) * p_hat_old_SO(j,m);
                    %diff w.r.t T and T
                    P_TT(j,n,k,m) = ((m==k) * (-1/THETA(j)) * pi(n+(j-1)*N,k) - (-1/THETA(j)) * pi(n+(j-1)*N,k) * pi(n+(j-1)*N,m)) * T_hat(j,k) * T_hat(j,m);
                    %diff w.r.t w and T
                    P_wT(j,n,k,m) = ((m==k) * GAMMA(j,k) * pi(n+(j-1)*N,k) - pi(n+(j-1)*N,k) * pi(n+(j-1)*N,m) * GAMMA(j,k))* w_hat(j,k) * T_hat(j,m);
                    %diff w.r.t p and T
                    P_pT(j,n,k,m) = ((m==k) * (1-GAMMA(j,k)) * pi(n+(j-1)*N,k) - pi(n+(j-1)*N,k) * pi(n+(j-1)*N,m) * (1-GAMMA(j,k))) * p_hat_old_SO(j,k) * T_hat(j,m);
                   %diff w.r.t w and P
                    P_wp(j,n,k,m) = ((m==k) * (1-GAMMA(j,k)) * (-THETA(j)) * GAMMA(j,k) * pi(n+(j-1)*N,k) - (-THETA(j)) * (GAMMA(j,k)) * (1-GAMMA(j,m)) * pi(n+(j-1)*N,k) * pi(n+(j-1)*N,m)) * w_hat(j,k) * p_hat_old_SO(j,m);
                end
            end
        end
    end
    p_hat_new_SO(:,:) = sum(P_w+P_p+P_T,3)  + (1/2) * sum(sum(P_ww+P_pp+P_TT+P_wT+P_pT+P_wp,4),3);
    %p_hat_new_SO(:,:) = sum(P_w+P_p+P_T,3); % first-order approximation
  
    pmax_SO=max(max(abs(p_hat_new_SO(:,:)-p_hat_old_SO(:,:))));
    pmax_SO
    if ITER_SO >50000 || sum(isnan(p_hat_new_SO(:)))>0
        disp('Second Order err')
        ITER_SO
        stop
    end
    p_hat_old_SO = 0.9*p_hat_old_SO + 0.1*p_hat_new_SO;
%    p_hat_old_SO = p_hat_new_SO ;
    ITER_SO = ITER_SO+1;
end
p_hat_SO = exp(p_hat_old_SO);
p_SO = p_lev .* exp(p_hat_old_SO);

%% Plot the results
price = [p_hat_lev(:) p_hat_FO(:) p_hat_SO(:)];
%price = [p_prime_lev(:) p_FO(:)];
disp('########################################')
disp('Approximation error for price')
sum(abs(price(:,1)-price(:,2))) % error of the first order
sum(abs(price(:,1)-price(:,3))) % error of the second order
x = linspace(min(price(:)),max(price(:)));
y = x;
figure
hold on
plot(x,y,'black')
scatter(price(:,1),price(:,2),'red')
scatter(price(:,1),price(:,3),'blue','d')
xlabel('Model Outcome (price)') 
ylabel('Approximation (price)') 
legend('45 degree line', 'First Order', 'Second Order', 'location', 'best')
title('Level / FO / SO Comparison')
saveas(gcf,'FIGURES/Priceloop.png')
(abs(price(:,1)-price(:,2))) - (abs(price(:,1)-price(:,3)))