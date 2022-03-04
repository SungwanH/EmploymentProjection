function [eqm_nlpf_HAT,approx_nlpf_HAT] = NLPF_HAT(params, initial, SS)
% Non-linear baseline perefect foresight equilibrium
% 50 regions + 37 other countries
%% Roll down parameters
v2struct(params.envr);
v2struct(params.tech);
v2struct(params.modl);

%% Roll down initial values
v2struct(initial);

%%%%%%%%%Algorithm%%%%%%%%%%%%%
%%%Dynamic problem%%%
Ymax = 1; ITER_DYN = 1;
while (ITER_DYN <= MAXIT) && (Ymax > TOL_NL)
    %Solving for the path of migration flows mu
    for t=1:TIME
        V(:,:,t) = reshape(v_td(:,t),J,R);
    end
    
    for t=1:TIME
        Vaux0(:,:,t) = reshape(V(:,:,t),(J)*R,1);
    end
    Vaux = permute(Vaux0,[2,1,3]);
    
    Vaux2 = zeros((J)*R,R*(J),TIME);
    for t=1:TIME
        Vaux2(:,:,t) = repmat(Vaux(:,:,t),R*(J),1);
    end
    
    %Computing mu1
    num       = mu0.*(Vaux2(:,:,2).^BETA);
    den       = sum(num')';
    den       = den*ones(1,R*(J));
    mu00      = num./den;
    mu        = zeros(R*(J),R*(J),TIME);
    mu(:,:,1) = mu00;  %this is mu1 (what we call mu0 in the paper)
    
    %Solving for mu(t+1)
    for t=1:TIME-2
        num=mu(:,:,t).*(Vaux2(:,:,t+2).^BETA);
        den=sum(num')';
        den=den*ones(1,R*(J));
        mu(:,:,t+1)=num./den;  %this is mu(t+1)
    end
   
    %Solving for the path of employment
    L00         = reshape(L0,J,R);
    L00_aux2    = reshape(L00,R*(J),1);
    L00_aux4    = repmat(L00_aux2,1,R*(J));
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
    % run the model once to match the initial data with our setting
    w_guess   = ones(J,N); %initial guess for wage
    p_guess   = ones(J,N); %initial guess for good prices
    kappa_hat = ones(J*N,N); % relative change in trade cost
    
    load('DATA/BASE_FOURSECTOR.mat', 'VALjn00', 'Din00')
    VALjn00   = VALjn00*8; % change quarterly data to bi-annual
    Ljn_hat00 = ones(J,N);
    T_hat00   = ones(J,N);
    [~, ~, ~, Din00_matched, X_matched, VALjn00_matched] = NLPF_TEMP_HAT(params, VALjn00, Din00, kappa_hat, T_hat00, Ljn_hat00, w_guess, p_guess);
    % update the initial values                                                
    VALjn0    = VALjn00_matched; %Labor compensation (w*L)
    Din0      = Din00_matched; %bilateral trade shares
    X0        = X_matched; %Total expenditure 

    %%%%Temporary Equilibrium%%%%%  
    Ltemp=Ldyn;   %path of employment 

    realwages      = ones(J,N,TIME);    %real wages. This matrix will store equilibrium real wages from the temporary equilibrium at each time t
    wf00           = ones(J,N,TIME);    %wages
    pf00           = ones(J,N,TIME);    %prices
    Pf00           = ones(J,N,TIME);    %price index
    pi             = zeros(N*J,N,TIME); %trade share
    pi(:,:,1)      = Din0;
    VALjn00        = zeros(J,N,TIME);   %labor income
    VALjn00(:,:,1) = VALjn0;
    X              = zeros(J,N,TIME);   %expenditure
    X(:,:,1)       = X0;
    %static sub-problem at each time t
    for t=1:TIME-2
        disp(t);
        %Shocks
        T_temp           = T_HAT(:,:,t+1); 
        Ljn_hat          = ones(J,N); %change in employment in the US
        Ljn_hat(:,1:R)   = Ltemp(:,:,t+1)./Ltemp(:,:,t);
        
        [wf0, pf0, Pf0, pi_temp, X_temp, VALjn] = NLPF_TEMP_HAT(params, VALjn0, Din0, kappa_hat, T_temp, Ljn_hat, w_guess, p_guess);
        wf00(:,:,t+1)    = wf0;
        pf00(:,:,t+1)    = pf0;
        Pf00(1,:,t+1)    = Pf0; 
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
    
    %%%%Solving for the new path of values in time differences%%%
    realwagesaux = zeros(J,N,TIME);
    for t = 1:TIME
        realwagesaux(1,:,t) = 1;
        realwagesaux(:,:,t) = realwages(:,:,t);
    end
    
    rwage=ones(R,J,TIME);
    
    for t = 1:TIME
        rwage(:,:,t) = permute(realwagesaux(:,1:R,t),[2,1,3]);
    end
    
    rw     = rwage.^(1/NU);
    rw     = permute(rw,[2,1,3]);
    rw_aux = zeros(R*(J),1,TIME);
    for t=1:TIME
        rw_aux(:,:,t) = reshape(rw(:,:,t),R*(J),1);
    end
    
    rw_aux2 = zeros(R*(J),R*(J),TIME);
    for t=1:TIME
        rw_aux2(:,:,t) = repmat(rw_aux(:,:,t),1,R*(J));
    end
    
    rwagenu=zeros(R*(J),R*(J),TIME);
    
    for t=1:TIME-1
        rwagenu(:,:,t+1) = mu(:,:,t).*rw_aux2(:,:,t+1);
    end
    
    num=zeros(size(rwagenu));
    for t=1:TIME-1
        num(:,:,t) = rwagenu(:,:,t).*(Vaux2(:,:,t+1).^BETA);
    end
    
    Y = zeros(R*(J),1,TIME);
    for t=1:TIME
        Y(:,:,t) = sum(num(:,:,t)')';  %this is the new path of Ys
    end
    Y(:,:,TIME)=1;
    
        
    Ynew = zeros((J)*R, TIME);
    for t=1:TIME
        Ynew(:,t) = Y(:,:,t);
    end
    Ynew(:,TIME) = 1;
    v_td(:,TIME) = 1;
    
    %Excess function
    checkY = zeros(TIME,1); 
    for t=2:TIME
%        checkY(t,1)=max(abs(Ynew(:,t)-v_td(:,t)));
        checkY(t,1) = max(abs(log(Ynew(:,t))-log(v_td(:,t))));
    end
    Ymax = max(checkY)
    
    v_td = UPDT_V_NL * Ynew + (1-UPDT_V_NL) * v_td;
    ITER_DYN = ITER_DYN+1;
end

% approximation points
varrho = zeros(N*J,N,TIME-2);
zeta   = zeros(N*J,J,TIME-2);
chi    = zeros(N*J,N,TIME-2);
lambda = zeros(R*J,R*J,TIME);

for t=1:TIME-2
    for n=1:N
        for j=1:J
            for ii=1:N
                varrho(n+(j-1)*N,ii,t) = pi(n+(j-1)*N,ii,t) * X(j,n,t) / X(j,ii,t);
                chi(n+(j-1)*N,ii,t)    = pi(n+(j-1)*N,ii,t) * X(j,n,t) / VALjn00(j,ii,t);
            end
        end
    end
end
for t=1:TIME-2
    for k=1:J
        for ii=1:N
            for j=1:J
                zeta(ii+(k-1)*N,j,t) = VALjn00(k,ii,t)./ X(j,ii,t);
            end
        end
    end
end
for t=1:TIME-2
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
for t=TIME-1:TIME
varrho(:,:,t)  = varrho(:,:,TIME-2);
pi(:,:,t)      = pi(:,:,TIME-2);
chi(:,:,t)     = chi(:,:,TIME-2);
zeta(:,:,t)    = zeta(:,:,TIME-2);
lambda(:,:,t)  = lambda(:,:,TIME-2);
mu(:,:,t)      = mu(:,:,TIME-2);
wf00(:,:,t)    = wf00(:,:,TIME-2);
pf00(:,:,t)    = pf00(:,:,TIME-2);
VALjn00(:,:,t) = VALjn00(:,:,TIME-2);
X(:,:,t)       = X(:,:,TIME-2);
end
%normalize
%for t=1:TIME
%    for i=1:R*J
%        lambda(:,i,t) =    lambda(:,i,t)./sum(sum(lambda(:,i,t)));
%    end
%end


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
eqm_nlpf_HAT = eqm_nlpf_HAT_belief;
approx_nlpf_HAT = approx_nlpf_HAT_belief;
save('DATA/NLPF_HAT_BELIEF.mat', 'eqm_nlpf_HAT_belief','approx_nlpf_HAT_belief'); 
end



end