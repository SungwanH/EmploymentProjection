function [eqm] = FIXED_POINT_HETERO(params, E_A_T_hat, E_B_T_hat, kappa_A_hat, kappa_B_hat, L_A, L_B, V, approx, mat_pbp)
% This function finds fixed point of Labor share of A and B in each period.
% At every period t1, the perceived path of the shock in fundamental
% (E_T_hat) is updated
% For given labor of other type in t1+1,
% we obtain the labor allocation of own type, and this is an input for
% solving for other type's labor allocation. Iterate this to find a fixed
% point.
% We solve each type's  expected path of (deviation from PF of) mu, v, L. 
%%%%%%%%%Algorithm%%%%%%%%%%%%%
%% Roll down parameters

v2struct(params.envr);
v2struct(params.tech);
v2struct(params.modl);


% Solve for A given L_B(t+1)

E_A_T_hat = params.belief.E_T_hat;
E_B_T_hat = params.belief.E_T_hat;
kappa_A_hat = zeros(N*J,N,TIME);
kappa_B_hat = zeros(N*J,N,TIME);

ITER_FP = 1;
Lmax=1;
L_A = zeros(N*J,1);
L_B = zeros(N*J,1);
L_A_2 = zeros(N*J,1);
L_B_2 = zeros(N*J,1);
while (ITER_FP <= MAXIT) && (Lmax > TOLFP)
    for t1=1:1
    %for t1=1:ENDT+1
        [eqm_A] = PBP_DYN_HETERO(params, t1, t1, E_A_T_hat, kappa_A_hat, L_A, L_B_2, V, approx, mat_pbp);
        L_A_new = eqm_A.L_A;
        L_A_2_new = eqm_A.L_A(:,2);
        [eqm_B] = PBP_DYN_HETERO(params, t1, t1, E_B_T_hat, kappa_B_hat, L_B, L_A_2, V, approx, mat_pbp);
        L_B_new = eqm_B.L_B;
        L_B_2_new = eqm_B.L_B(:,2);
    
        
    %    for t=t1:TIME
            checkL(t,1)=max(abs(L_A_2_new(:,t)-L_A_2(:,t)));
            checkL(t,2)=max(abs(L_B_2_new(:,t)-L_B_2(:,t)));
    %    end
        [Lmax loc]=max(max(checkL));
        Lmax
        if ITER_FP >3000 || sum(isnan(w_update(:)))>0
            checkL(t1:TIME)
            t1
            disp('Fixed Point loop err')
            ITER_FP
            stop
        end
        
        %L_hat = (1-UPDT_L)*L_A_temp+UPDT_L*L_A_update;
        L_A_2 = L_A_2_new;
        L_B_2 = L_B_2_new;
        ITER_FP=ITER_FP+1;
    end
end

end