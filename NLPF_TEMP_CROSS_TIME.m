function [wf0, pf0, Pf0, Dinp, X, VALjn] = NLPF_TEMP_CROSS_TIME(params, VALjn_dot, Din_dot, VALjn0, Din, kappa_hat, T_hat, Ljn_hat, w_guess, p_guess)
% This function solves temporary equilibrium given Labor.
% Updates price and wage simultaneously but at a different speed
%% Roll down parameters
v2struct(params.tech)
v2struct(params.envr)
v2struct(params.modl)

%wf0 = ones(J,N);
%pf0 = ones(J,N);
wf0 = w_guess;
pf0 = p_guess;
X = zeros(J,N);
Dinp = zeros(N*J,N);
ITER_TEMP = 0;
wfmax = 1;
pfmax = 1;

while (ITER_TEMP <= MAXIT) && ((wfmax > TOL_NL_TEMP) || (pfmax > TOL_NL_TEMP))  %|| (Xmax > TOL_NL))
    
    % step  1. get an updated p
    lw = log(wf0);
    lp = log(pf0);
    % calculating log cost
    lc = GAMMA.*lw + (1-GAMMA).*lp;
    c  = exp(lc);

    for j=1:J
        idx = 1+(j-1)*N:1:N*j;
        LT(idx,1) = ones(N,1)*(1/THETA(j));
    end 
    Din_k = Din.*(kappa_hat.^(-1./(LT*ones(1,N))));      

    %Calculate change in prices 
    for j=1:1:J
        for n=1:1:N
%        phat(j,n)=Din_k(n+(j-1)*N,:)*((T_hat(j,:).^(THETA(j)*GAMMA(j,:))).*(c(j,:).^(-THETA(j))))';
        phat(j,n) = Din_dot(n+(j-1)*N,:) .* Din_k(n+(j-1)*N,:)*(T_hat(j,:).*(c(j,:).^(-THETA(j))))';
        phat(j,n) = phat(j,n)^(-1/THETA(j));
        end              
    end

    % step 2 calculate for Dinp (bilateral trade shares)
    for n=1:1:N
        cp(:,n)    = c(:,n).^( - THETA );
        phatp(:,n) = phat(:,n).^(-THETA);
    end
    
    for n=1:1:N
        idx=n:N:length(Din)-(N-n);
        DD(idx,:) = Din_dot(idx,:) .* Din_k(idx,:).*(cp.*T_hat);   
        %DD(isnan(DD)) = 0;
    end
    for n=1:1:N
        idx = n:N:length(Din)-(N-n);
        Dinp(idx,:) = DD(idx,:)./(phatp(:,n)*ones(1,N)); 
    end
    %Dinp(1:20,1:20)
    % step 3: Solve for total expenditure (by taking inverse)

    % for given j, 
    % A*X = RHS
    % the first two rows of A are:
    % ( 1-(1-B_1j)pi_1j1   -(1-B_1j)pi_2j1    -(1-B_1j)pi_3j1 ...)
    % ( -(1-B_2j)pi_1j2   1-(1-B_2j)pi_2j2    -(1-B_2j)pi_3j2 ...)
    
    RHS = zeros(J,N);
    for j=1:J
        RHS(j,:) = ALPHAS(j,:).*sum((wf0.*Ljn_hat.*VALjn0.*VALjn_dot),1); 
    end
    

    for j=1:J
        idx=(1+(j-1)*N):(j*N);
        PI=Dinp(idx,:);
        X(j,:)=((eye(N)-(1-GAMMA(j,:)')*ones(1,N).*PI')\(RHS(j,:)'))';
    end
    

    %}    
    % Alternative formulation below to check the stacking: results checks
    % out
    %{
    RHS = zeros(J,N);
    for j=1:J
        RHS(j,:) = ALPHAS(j,:).*sum((wf0.*Ljn_hat.*VALjn0.*VALjn_dot),1); 
    end
        A = zeros(N,N);
    for j=1:J
        for i=1:N
            for n=1:N
                A(i,n) = -(1-GAMMA(j,i))*Dinp(n+(j-1)*N, i);
            end
            A(i,i) = 1-(1-GAMMA(j,i))*Dinp(i+(j-1)*N, i);
        end
        X(j,:) = (A\(RHS(j,:)'))';
    end
    %}

    for n=1:N
        for j=1:J
            for i=1:N
                piX(j,i,n) = Dinp(n+(j-1)*N, i) * X(j,n);
            end
        end
    end
    sumpiX = sum(piX,3);
 
    % Update wage using labor market clearing
    w_new = (GAMMA.*sumpiX)./ (Ljn_hat .*VALjn0 .*VALjn_dot);
    
    %update p and w     
%    pfdev2    = norm((phat - pf0)./(pf0),1);
    pfdev    = log(phat) - log(pf0);
    pf0      = phat;
    pfmax    = max(max(pfdev));
    
%    wfdev    = abs(w_new - wf0); % Checking tolerance
    wfdev    = abs(log(w_new) - log(wf0));
    wf0      = w_new*UPDT_W_NL + wf0*(1-UPDT_W_NL);
    wfmax    = max(max(wfdev));
    ITER_TEMP       = ITER_TEMP + 1;
    
% normalize nominal wage change in 1 to 0.    
     wf0=wf0./wf0(1);
     pf0=pf0./wf0(1);
end

% Price index
Pf0 = prod(pf0.^(ALPHAS),1);
% Update labor income in level
VALjn = VALjn0 .* (wf0 .* Ljn_hat).*VALjn_dot;
%disp('Number of iteration')
%disp(ITER_TEMP)
%sum(L(:))

end