function T = PRODUCTIVITY_DGP(params, flag)
% This function generates the objective DGP for the productivities in each
% region
%% Roll down parameters
v2struct(params.prod);
v2struct(params.envr);

% Initialize productivity
T = T_BASE(:,:,1:TIME);
%draw error term for productivity
EPS = PRODUCTIVITY_DRAW(params);

 %now this is for AR(1) for productvity gap (in log)
% ln(z_it+1^US) - ln(z_it+1^CHN) = RHO * (ln(z_it^US) - ln(z_it^CHN)) + EPS(it+1)
gap = zeros(J,1,TIME);
gap(:,1,1) = log(T(:,US,1)) - log(T(:,CHINA,1));
for t=1:TIME-1
    for k=1:J
        gap(k,1,t+1) = (1-RHO) * MU + RHO * gap(k,t) + EPS(k,t+1);
    end
end

for t=1:TIME
    T(:,CHINA,t) = exp(log(T(:,US,t)) - gap(:,1,t));
end
if nargin>1 %now this is for AR(1) for growth rate    
    for t=1:TIME-1
        for k=1:J
            GR(k,1,t+1) = (1-RHO) * MU + RHO * GR(k,t) + EPS(k,t+1);
        end
    end
end
%Here I force to converge at ENDT+1
%T(CHINA,ENDT+1:TIME) = T(CHINA,TIME);
end
