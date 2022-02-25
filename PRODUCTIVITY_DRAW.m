function EPS = PRODUCTIVITY_DRAW(params)
% This function generates the objective DGP for the productivities in each
% region
%% Roll down parameters
%N      =   params.N;
%J      =   params.J;
%TIME   =   params.TIME;
%EPST   =   params.EPST;
%SIGMA  =   params.SIGMA; % Error term in objective DGP

v2struct(params.envr);
v2struct(params.modl);
v2struct(params.prod);

rng(123457)
%productivity gap draw
EPS = zeros(J,TIME);

for t=1:EPST
    EPS(:,t) = normrnd(0,SIGMA,[J,1]);
end

end
