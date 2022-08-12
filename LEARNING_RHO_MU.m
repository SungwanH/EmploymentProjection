function [MU_HAT,RHO_HAT]= LEARNING_RHO_MU(params, W)
% This function generates estimates of MU and RHO using OLS
%% Roll down parameters
v2struct(params.prod);
v2struct(params.envr);
v2struct(params.tech);

RHO_OLS = zeros(J,TIME);
MU_OLS  = zeros(J,TIME);
RHO_OLS(:,:) = RHO;
MU_OLS(:,:) = MU;

%t_bef = 10; %length of period before model begins
if W==1

    RHO_HAT = RHO * ones(J,TIME);
    MU_HAT  = MU * ones(J,TIME);

else
    gap_before = reshape(log(T(:,US,1))-log(T_PREV),J,t_bef); % the twenty periods before 
    gap= reshape(log(T(:,US,:)) - log(T(:,CHINA,:)),J,TIME);
    
            
    if ESTM_BOTH==1
        
        for t=1:ENDT % The pre-period realizations are always used in estimation
            for j=1:J
                Y = [gap_before(j,1:t_bef)';gap(j,1:t)'];
                X = [ones(t_bef+t-1,1),[gap_before(j,1:t_bef,1)';gap(j,1:t-1)']];
                B = X\Y;
                MU_OLS(j,t)  = B(1)/(1-B(2));
                RHO_OLS(j,t) = B(2);
            end
        end
        
    elseif ESTM_BOTH==0
    
        RHO_OLS(:,:) = RHO;
        for t=1:ENDT
            for j=1:J
                TEMP = [gap_before(j,1:t_bef)';gap(j,1:t)'];
                Y= TEMP(2:end)-RHO*TEMP(1:end-1);
    %            B = X\Y;
                %    MU_OLS(t) = B(1)/(1-RHO);
                MU_OLS(j,t) =mean(Y)/(1-RHO);
            end
        end
        
        
    end
    RHO_HAT = W * RHO + (1-W) .* RHO_OLS(:,1:TIME);
    MU_HAT  = W * MU  + (1-W) .* MU_OLS(:,1:TIME);
    RHO_HAT(:,ENDT+1:TIME) = RHO;
    MU_HAT(:,ENDT+1:TIME)  = MU;
end
end
