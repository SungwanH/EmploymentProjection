function [eqm_nlpf_CROSS,approx_nlpf_CROSS] = NLPF_DD_NEW(params, starting_point, base_point, hat_fundamentals, initial_guess)
% Non-linear double-difference (time & cross difference) perefect foresight equilibrium
% 50 regions + 37 other countries
% Timing assumption: t=1 corresponds to the first time index in the code
%% Roll down parameters
v2struct(params.envr);
v2struct(params.tech);
v2struct(params.modl);

%% Roll down the initial guess for the solution
v2struct(initial_guess);

%% Roll down time differences in fundamentals; only T_HAT in this vintage
v2struct(hat_fundamentals); %this overwrites TIME in params.envr
v2struct(base_point);
%%%%%%%%%Algorithm%%%%%%%%%%%%%
%%%Dynamic problem%%%
Ymax = 1; ITER_DYN = 1;
while (ITER_DYN <= MAXIT) && (Ymax > TOL_NL)
    % Rolling out the starting point for the solution. Note that this
    % initial value CANNOT be moved to outside the while loop
    v2struct(starting_point);
    v2struct(base_point);

    %Solving for the path of migration flows mu
   
    Vaux=reshape(v_cd,1,J*R,TIME);
    Vaux2=repmat(Vaux,R*J,1,1); % Vaux2, first dimension is repeation, second dimension is the value of each cell in a period, ordered by country, i.e., two consecutive entries are obs from different sectors for the same country; note that this needs be consistent with wage and value function update below
    
    %Computing mu0
    %% 
    mu        = zeros(R*(J),R*(J),TIME);
    %Computing mu1
%    vartheta  = mu0.*mu_base_dot(:,:,1).*(Vaux2(:,:,1).^BETA);
%    den       = sum(vartheta,2);
%    den       = den*ones(1,R*(J));
%    mu00      = vartheta./den;

    vartheta  = mu0.*mu_base_dot(:,:,1).*(Vaux2(:,:,1).^BETA); %coming from CDP
    num       = vartheta.*(Vaux2(:,:,2).^BETA);
    den       = sum(num,2);
    den       = den*ones(1,R*(J));
    mu(:,:,1) = num./den; % this is mu'(1)

    %Solving for mu'(t+1)
    for t=1:TIME-2
        num         = mu(:,:,t).*mu_base_dot(:,:,t+1).*(Vaux2(:,:,t+2).^BETA);
        den         = sum(num,2);
        den         = den*ones(1,R*(J));
        mu(:,:,t+1) = num./den;  %this is mu(t+1)
    end
    mu(:,:,TIME)=mu(:,:,TIME-1);
    %mu(1,:,1:10)
    %Solving for the path of employment
    L00         = reshape(L0,J,R);
    L00_aux2    = reshape(L00,R*(J),1); %two consecutive obs are sectors from the same country
    L00_aux4    = repmat(L00_aux2,1,R*(J)); % each row is a constant; two consecutive obs in a column are different sectors from the same country
%    L00_aux5    = mu_base(:,:,1).*L00_aux4;
    L00_aux5    = mu0.*L00_aux4;
    L1          = sum(L00_aux5)';
    L1          = reshape(L1,J,R);
    Ldyn        = zeros(J,R,TIME);
    Ldyn(:,:,1) = L1;
%    Ldyn(:,:,1) = reshape(L0,(J),R);
    for t=1:TIME-1
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
    % pi, VALjn00, X should come from next period from pre-period(observed)
    pi(:,:,1)      = Din0; %%%%
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
        T_temp           = T_HAT(:,:,t+1); %double difference
        Ljn_hat          = ones(J,N); %change in employment in the US
        Ljn_hat(:,1:R)   = (Ltemp(:,:,t+1)./Ltemp(:,:,t))./Ldyn_base_dot(:,:,t+1); % double difference
        VALjn_dot        = VALjn_base_dot(:,:,t+1); %Time difference of labor income in baseline economy
        Din_dot          = Din_base_dot(:,:,t+1);

        [wf0, pf0, Pf0, pi_temp, X_temp, VALjn] = NLPF_DD_TEMP(params, VALjn_dot, Din_dot, VALjn0, Din0, kappa_hat, T_temp, Ljn_hat, w_guess, p_guess);
        wf00(:,:,t+1)    = wf0;
        pf00(:,:,t+1)    = pf0;
        Pf00(:,t+1)      = Pf0; 
        pi(:,:,t+1)      = pi_temp;
        X(:,:,t+1)       = X_temp;
        VALjn00(:,:,t+1) = VALjn;
        w_guess          = wf0;
        p_guess          = pf0;
        %updating the initial conditions
        VALjn0=VALjn;
        Din0=pi_temp;
     
        %storing equilibrium real wages
        realwages(:,:,t+1) = wf0(:,:)./(ones(J,1)*Pf0);
    end

    %%%%Solving for the new path of values in cross differences%%%        
    realwages_us=realwages(:,1:R,:); 
    realwages_us_aux2=reshape(realwages_us,R*J,TIME); % after this step, realwages_us_aux2 is a 200X200; in the first dimension, two consecutive observations are different sectors of the same country
%    realwages_us_nuu=realwages_us_aux2.^(1/NU);
    rw_lev(:,:,1) = realwages_us(:,:,1);
    rw_lev_base = ones(J,R,TIME); % This holds only when the initial period starts from steady state
    for t=1:TIME-1
        rw_lev(:,:,t+1) = realwages_us(:,:,t+1) .* rw_lev_base(:,1:R,t+1) .* rw_lev(:,1:R,t) ./ rw_lev_base(:,1:R,t);
    end
    rw_lev_RJ = reshape(rw_lev,R*J,TIME);
    realwages_us_nuu=exp((1/(1-CRRA_PAR)) * rw_lev_RJ.^(1-CRRA_PAR)).^(1/NU); % This term is due to CRRA utility
    
    % we now update value function backwards; Y standard for the changes in
    % value function from this iteration
    %Y=NaN(R*(J),TIME);    
    Y=ones(R*(J),TIME);    
    %Y(:,TIME)=realwages_us_nuu(:,TIME);    
    for tt=TIME-1:-1:2
        temp0=ones(J*R,1)*(Y(:,tt+1).^BETA)';
        if tt==1
%            temp=sum(mu0.*mu_base_dot(:,:,1).*temp0,2); 
            temp=sum(vartheta.*temp0,2); 
        else
            temp=sum(mu(:,:,tt-1).*mu_base_dot(:,:,tt).*temp0,2); 
        end
        Y(:,tt)=realwages_us_nuu(:,tt).*temp;
    end 
    %for period 1 (index: 1)
    temp0=ones(J*R,1)*(Y(:,2).^BETA)';
    temp=sum(vartheta.*temp0,2); 
    Y(:,1)=realwages_us_nuu(:,1).*temp;

    %Excess function
    checkY = zeros(TIME,1); 
    for t=1:TIME
        checkY(t,1) = max(abs(log(Y(:,t))-log(v_cd(:,t))));
    end
    Ymax = max(checkY)    
    v_cd = UPDT_V_NL * Y + (1-UPDT_V_NL) * v_cd;
    ITER_DYN = ITER_DYN+1;
end

%% Recover level wage and price
w_lev(:,:,1) = ones(J,N);
p_lev(:,:,1) = ones(J,N);
for t=1:TIME-1
    w_td_cf(:,:,t+1) = wf00(:,:,t+1).* wf00_base(:,:,t+1);
    p_td_cf(:,:,t+1) = pf00(:,:,t+1).* pf00_base(:,:,t+1);
    w_lev(:,:,t+1)   = w_lev(:,:,t) .* w_td_cf(:,:,t+1);
    p_lev(:,:,t+1)   = p_lev(:,:,t) .* p_td_cf(:,:,t+1);
end
for t=1:TIME
    P_lev(1,:,t) = prod((p_lev(:,:,t)./ALPHAS).^ALPHAS,1);
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
                        lambda(k+(ii-1)*J,j+(n-1)*J,1) = mu(k+(ii-1)*J,j+(n-1)*J,1)*L0(k,ii,1)/Ldyn(j,n,1);
                    else
                        lambda(k+(ii-1)*J,j+(n-1)*J,t) = mu(k+(ii-1)*J,j+(n-1)*J,t-1)*Ldyn(k,ii,t-1)/Ldyn(j,n,t);
                    end
                end
            end
        end
    end
end
%{
for t=2:TIME
    v_cd(:,t-1)        = v_cd(:,t);
    Ldyn(:,:,t-1)      = Ldyn(:,:,t);
    realwages(:,:,t-1) = realwages(:,:,t);
    wf00(:,:,t-1)      = wf00(:,:,t);
    pf00(:,:,t-1)      = pf00(:,:,t);
    Pf00(:,t-1)        = pf00(:,t);
    w_lev(:,:,t-1)     = w_lev(:,:,t);
    p_lev(:,:,t-1)     = p_lev(:,:,t);
    VALjn00(:,:,t-1)   = VALjn00(:,:,t);
    mu(:,:,t-1)        = mu(:,:,t);
    pi(:,:,t-1)        = pi(:,:,t);
    varrho(:,:,t-1)    = varrho(:,:,t);
    chi(:,:,t-1)       = chi(:,:,t);
    zeta(:,:,t-1)      = zeta(:,:,t);
    lambda(:,:,t-1)    = lambda(:,:,t);
end
%}
eqm_nlpf_CROSS=v2struct(v_cd, Ldyn, realwages, wf00, pf00, Pf00, w_lev, p_lev, P_lev, VALjn00, X);
approx_nlpf_CROSS= v2struct(mu, pi, varrho, chi, zeta, lambda);

end