function T_belief = BELIEF(params, W)
% This function generates the belief (at each period) on the productivity
% path in each region
%% Roll down parameters
v2struct(params.prod);
v2struct(params.envr);

[MU_HAT,RHO_HAT]= LEARNING_RHO_MU(params, W);
% use this RHO_HAT to generate belief future productivity/fundamental path
%Initialize belief:
for t=1:TIME
    T_belief(:,:,:,t) = T;
end
RHO_HAT(:,ENDT+1:TIME) = RHO; % perfect foresight after t=ENDT
MU_HAT(:,ENDT+1:TIME) = MU; % perfect foresight after t=ENDT
for j=1:J
    for t2=1:ENDT
        summ = 0;
        for t1=t2+1:TIME
            summ = summ + (1-RHO_HAT(j,t2))*MU_HAT(j,t2)*(RHO_HAT(j,t2))^(t1-t2-1);
            temp(j,t1,t2) = summ;
        end
    end
end
for j=1:J
    for t2=1:ENDT
        for t1=t2+1:TIME 
         T_belief(j,CHINA,t1,t2)= exp(log(T(j,US,t1)) - temp(j,t1,t2) - ((RHO_HAT(j,t2))^(t1-t2)) * (log(T(j,US,t2))-log(T(j,CHINA,t2))));
        end
    end
end
end
