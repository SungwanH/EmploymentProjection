function [eqm_nlpf_HAT,approx_nlpf_HAT] = NLPF_HAT(params, starting_point, hat_fundamentals,initial_guess)
% Non-linear baseline perefect foresight equilibrium (time-difference)
% 50 regions + 37 other countries
%% Roll down parameters
v2struct(params.envr);
v2struct(params.tech);
v2struct(params.modl);

%% Roll down the initial guess for the solution
v2struct(initial_guess);

%% Roll down time differences in fundamentals; only T_HAT in this vintage
v2struct(hat_fundamentals); %this overwrites TIME in params.envr

%%%%%%%%%Algorithm%%%%%%%%%%%%%
%%%Dynamic problem%%%
Ymax = 1; ITER_DYN = 1;
while (ITER_DYN <= MAXIT) && (Ymax > TOL_NL)
    % Rolling out the starting point for the solution. Note that this
    % initial value CANNOT be moved to outside the while loop
    v2struct(starting_point);

    %Solving for the path of migration flows mu
   
    Vaux=reshape(v_td,1,J*R,TIME);
    Vaux2=repmat(Vaux,R*J,1,1); % Vaux2, first dimension is repeation, second dimension is the value of each cell in a period, ordered by country, i.e., two consecutive entries are obs from different sectors for the same country; note that this needs be consistent with wage and value function update below
   
    %Computing mu1
    num       = mu0.*(Vaux2(:,:,2).^BETA);
    den       = sum(num,2);
    den       = den*ones(1,R*(J));
    mu00      = num./den;
    mu        = zeros(R*(J),R*(J),TIME);
    mu(:,:,1) = mu00;  %this is mu1 (what we call mu0 in the paper)
    
    %Solving for mu(t+1)
    for t=1:TIME-2
        num=mu(:,:,t).*(Vaux2(:,:,t+2).^BETA);
        den=sum(num,2);
        den=den*ones(1,R*(J));
        mu(:,:,t+1)=num./den;  %this is mu(t+1)
    end
       mu(:,:,TIME)=mu(:,:,TIME-1);
       
    %Solving for the path of employment
    L00         = reshape(L0,J,R);
    L00_aux2    = reshape(L00,R*(J),1); %two consecutive obs are sectors from the same country
    L00_aux4    = repmat(L00_aux2,1,R*(J)); % each row is a constant; two consecutive obs in a column are different sectors from the same country
    L00_aux5    = mu00.*L00_aux4;
    L1          = sum(L00_aux5)';
    L1          = reshape(L1,J,R);
    Ldyn        = zeros(J,R,TIME);
    Ldyn(:,:,2) = L1;
    Ldyn(:,:,1) = reshape(L0,(J),R);
    for t=2:TIME-1
        aux  = reshape(Ldyn(:,:,t),R*(J),1);
        aux3 = repmat(aux,1,R*(J));
        aux4 = mu(:,:,t).*aux3;
        aux5 = sum(aux4)';
        Ldyn(:,:,t+1)=reshape(aux5,J,R);   %this is the path for employment
    end
    %Ldyn(:,:,TIME)=0;

    %%%%Temporary Equilibrium%%%%%  
    Ltemp=Ldyn;   %path of employment 

    realwages      = ones(J,N,TIME);    %real wages. This matrix will store equilibrium real wages from the temporary equilibrium at each time t
    wf00           = ones(J,N,TIME);    %wages
    pf00           = ones(J,N,TIME);    %prices
    Pf00           = ones(N,TIME);    %price index
    pi             = zeros(N*J,N,TIME); %trade share
    pi(:,:,1)      = Din0;
    VALjn00        = zeros(J,N,TIME);   %labor income
    VALjn00(:,:,1) = VALjn0;
    X              = zeros(J,N,TIME);   %expenditure
    X(:,:,1)       = X0;
    
    w_guess   = ones(J,N); %initial guess for wage
    p_guess   = ones(J,N); %initial guess for good prices
    kappa_hat = ones(J*N,N); % relative change in trade cost

    %static sub-problem at each time t
    for t=1:TIME-1
        if mod(t,10)==0
            fprintf('%d, ', t);
        end
        %Shocks
        T_temp           = T_HAT(:,:,t+1); 
        Ljn_hat          = ones(J,N); %change in employment in the US
        Ljn_hat(:,1:R)   = Ltemp(:,:,t+1)./Ltemp(:,:,t);
        
        [wf0, pf0, Pf0, pi_temp, X_temp, VALjn] = NLPF_TEMP_HAT(params, VALjn0, Din0, kappa_hat, T_temp, Ljn_hat, w_guess, p_guess);
        wf00(:,:,t+1)    = wf0;
        pf00(:,:,t+1)    = pf0;
        Pf00(:,t+1)    = Pf0; 
        pi(:,:,t+1)      = pi_temp;
        X(:,:,t+1)       = X_temp;
        VALjn00(:,:,t+1) = VALjn;
%     w_guess          = wf0;
%     p_guess          = pf0;

        
        %updating the initial conditions
        VALjn0=VALjn;
        Din0=pi_temp;
     
        %storing equilibrium real wages
        realwages(:,:,t+1) = wf0(:,:)./(ones(J,1)*Pf0);
    end
    
    %%%%Solving for the new path of values in time differences%%%        
    realwages_us=realwages(:,1:R,:); 
    realwages_us_aux2=reshape(realwages_us,R*J,TIME); % after this step, realwages_us_aux2 is a 200X200; in the first dimension, two consecutive observations are different sectors of the same country
    %realwages_us_nuu=realwages_us_aux2.^(1/NU);
    realwages_us_nuu=exp((1/(1-CRRA_PAR)) * realwages_us_aux2.^(1-CRRA_PAR)).^(1/NU); % This term is due to CRRA utility
    
    % we now update value function backwards; Y standard for the changes in
    % value function from this iteration
    Y=NaN(R*(J),TIME);    
    Y(:,TIME)=realwages_us_nuu(:,TIME);    %fix this part
    for tt=TIME-1:-1:2
        temp0=ones(J*R,1)*(Y(:,tt+1).^BETA)';
        temp=sum(mu(:,:,tt).*temp0,2); 
        Y(:,tt)=realwages_us_nuu(:,tt).*temp;
    end    
    
    %Excess function
    checkY = zeros(TIME-1,1); 
    checkY = max(abs(log(Y(:,2:end))-log(v_td(:,2:end))));
    Ymax = max(checkY)    
    v_td(:,2:end) = UPDT_V_NL * Y(:,2:end) + (1-UPDT_V_NL) * v_td(:,2:end);
    ITER_DYN = ITER_DYN+1;
end

% approximation points
varrho = zeros(N*J,N,TIME); % period t of these variables are based on period t levels.
zeta   = zeros(N*J,J,TIME);
chi    = zeros(N*J,N,TIME);
lambda = zeros(R*J,R*J,TIME);

for t=1:TIME
    for n=1:N
        for j=1:J
            for ii=1:N
                varrho(n+(j-1)*N,ii,t) = pi(n+(j-1)*N,ii,t) * X(j,n,t) / X(j,ii,t);
                chi(n+(j-1)*N,ii,t)    = pi(n+(j-1)*N,ii,t) * X(j,n,t) / VALjn00(j,ii,t);
            end
        end
    end
end
for t=1:TIME
    for k=1:J
        for ii=1:N
            for j=1:J
                zeta(ii+(k-1)*N,j,t) = VALjn00(k,ii,t)./ X(j,ii,t);
            end
        end
    end
end
for t=1:TIME
    for k=1:J
        for ii=1:R
            for j=1:J
                for n=1:R
                    if t==1
                        lambda(k+(ii-1)*J,j+(n-1)*J,1) = mu(k+(ii-1)*J,j+(n-1)*J,1)*Ldyn(k,ii,1)/Ldyn(j,n,1);
                    else
                        lambda(k+(ii-1)*J,j+(n-1)*J,t) = mu(k+(ii-1)*J,j+(n-1)*J,t-1)*Ldyn(k,ii,t-1)/Ldyn(j,n,t);
                    end
                end
            end
        end
    end
end

%normalize
%for t=1:TIME
%    for i=1:R*J
%        lambda(:,i,t) =    lambda(:,i,t)./sum(sum(lambda(:,i,t)));
%    end
%end


eqm_nlpf_HAT=v2struct(v_td, Ldyn, realwages, wf00, pf00, VALjn00, X);
approx_nlpf_HAT= v2struct(mu, pi, varrho, chi, zeta, lambda);

%{
if SS==1
    eqm_nlpf_HAT_SS = v2struct(v_td, Ldyn, realwages, wf00, pf00, VALjn00, X);
    approx_nlpf_HAT_SS = v2struct(mu, pi, varrho, chi, zeta, lambda);
    eqm_nlpf_HAT = eqm_nlpf_HAT_SS;
    approx_nlpf_HAT = approx_nlpf_HAT_SS;
    save('DATA/NLPF_HAT_SS.mat', 'eqm_nlpf_HAT_SS','approx_nlpf_HAT_SS'); 

elseif SS==0
    eqm_nlpf_HAT = v2struct(v_td, Ldyn, realwages, wf00, pf00, VALjn00, X);
    approx_nlpf_HAT = v2struct(mu, pi, varrho, chi, zeta, lambda);
    save('DATA/NLPF_HAT.mat', 'eqm_nlpf_HAT','approx_nlpf_HAT'); 

else
    eqm_nlpf_HAT_belief = v2struct(v_td, Ldyn, realwages, wf00, pf00, VALjn00, X);
    approx_nlpf_HAT_belief = v2struct(mu, pi, varrho, chi, zeta, lambda);
    save('DATA/NLPF_HAT_BELIEF.mat', 'eqm_nlpf_HAT_belief','approx_nlpf_HAT_belief'); 
end
%}


end