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

t_bef = 20; %length of period before model begins
for t=1:t_bef
    for j=1:J
        gap_temp(j,t) = log(T(j,US,1))-log(T_PREV(j,t));
    end
end
%gap(j,1) = log(T(US,1)) - log(T(CHINA,1));
for t=1:TIME
    for j=1:j
        gap(j,t) = log(T(j,US,t)) - log(T(j,CHINA,t));
    end
end
if ESTM_BOTH==1
    for t=1:2
        for j=1:J
            Y = [gap_temp(j,2:t_bef)';gap(j,1:t)'];
            X = [ones(t_bef+t-1,1),[gap_temp(j,1:t_bef,1)';gap(j,1:t-1)']];
            B = X\Y;
            MU_OLS(j,t)  = B(1)/(1-B(2));
            RHO_OLS(j,t) = B(2);
        end
    end
    % To make the belief smooth, I put initial values as follows:
    MU_OLS(:,1) =0.35;
    RHO_OLS(:,1) =0.99;
    MU_OLS(:,2) =0.3;
    RHO_OLS(:,2) =0.98;
    
    for t=1:ENDT-2
        for j=1:J
        %    Y = [gap_temp(2:t_bef);gap(1:t+2,1)];
        %    X = [ones(t+1+t_bef,1),[gap_temp(1:t_bef);gap(1:t+1,1)]];
            Y = [gap(j,2:t+2)'];
            X = [ones(t+1,1),[gap(j,1:t+1)']];
            B = X\Y;
            MU_OLS(j,t+2) = B(1)/(1-B(2));
            RHO_OLS(j,t+2) = B(2);
        end
    end
else
    RHO_OLS(:,:) = RHO;
    for t=1:2
        for j=1:J
            Y = [gap_temp(j,2)';gap(j,1:t)'];
            X = [ones(2+t-1,1)];
            B = X\Y;
            %    MU_OLS(t) = B(1)/(1-RHO);
            MU_OLS(j,t) = B(1);
        end
    end
    
    %MU_OLS(1) =0.35;
    %RHO_OLS(1) =0.99;
    %MU_OLS(2) =0.3;
    %RHO_OLS(2) =0.98;
    
    for t=1:ENDT-2
        for j=1:J
            Y = [gap_temp(j,2:t_bef)';gap(j,1:t+2)'];
            X = [ones(t+1+t_bef,1)];
            %Y = [gap(2:t+2,1)];
            %X = [ones(t+1,1),[gap(1:t+1,1)]];
            B = X\Y;
            MU_OLS(j,t+2) = B(1);
        end
    end
end
RHO_HAT = W * RHO + (1-W) .* RHO_OLS(:,1:TIME);
MU_HAT  = W * MU  + (1-W) .* MU_OLS(:,1:TIME);
RHO_HAT(:,ENDT+1:TIME) = RHO;
MU_HAT(:,ENDT+1:TIME)  = MU;

end
