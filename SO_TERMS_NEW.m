function [SO_temp, SO_dyn] = SO_TERMS_NEW(params, t1, eqm_sim, SO_COEF_temp, SO_COEF_dyn)
% This function computes second order terms using the cofficients of
% second order terms (input (SO_COEF_temp, SO_COEF_dyn)) and second order
% moments (calling SO_MOMENTS later)
%% Roll down parameters
v2struct(params.envr);
v2struct(params.tech);
v2struct(params.modl);

%% Roll down approximation points


[mom_dyn, mom_temp] = SO_MOMENTS(params, t1, eqm_sim);
v2struct(mom_dyn);
v2struct(mom_temp);
v2struct(SO_COEF_temp);
v2struct(SO_COEF_dyn);
%% SO terms for Dynamic eqm 
SO_L = zeros(R*J,TIME); SO_v = zeros(R*J,TIME);
for t = t1:TIME-1
    for m=1:R*J
        for k=1:R*J
            for g=1:R*J
                SO_L_LL(g,k,m) = L_LL(g,k,m,t) * ll(k,m,t);
                SO_L_vL(g,k,m) = L_vL(g,k,m,t) * vl(k,m,t);
                SO_L_vv(g,k,m) = L_vv(g,k,m,t) * vv(k,m,t);
            end
        end
    end
    SO_L(:,t+1) = (1/2) * sum(sum(SO_L_LL,3),2) + (1/2) * sum(sum(SO_L_vv,3),2) + sum(sum(SO_L_vL,3),2);
end

for t=t1:TIME
    for m=1:R*J
        for k=1:R*J
            for n=1:R*J
               SO_v_vv(n,k,m) =  v_vv(n,k,m,t) * vv(m,k,t);
            end
        end
    end
    SO_v(:,t) = (1/2) * sum(sum(SO_v_vv,3),2); % This is R*J by 1 
end

SO_v_ww(:,t1:TIME) = v_ww(:,t1:TIME) .* ww_dyn(:,t1:TIME);
SO_v_PP(:,t1:TIME) = v_PP(:,t1:TIME) .* PP_dyn(:,t1:TIME);
SO_v_wP(:,t1:TIME) = v_wP(:,t1:TIME) .* wP_dyn(:,t1:TIME);
SO_rw = (1/2) * SO_v_ww + (1/2) * SO_v_PP + SO_v_wP; % This is J by R 

%% SO terms for Temporary eqm 
SO_p = zeros(J,N,TIME); SO_X = zeros(J,N,TIME); SO_w = zeros(J,N,TIME);
for t = t1:TIME
    for j=1:J
        for n=1:N
            for k=1:N        
                for m=1:N
                    %diff w.r.t w and w
                    SO_P_ww(j,n,k,m) = P_ww(j,n,k,m,t) * ww_price(j,k,m,t);
                   %diff w.r.t p and p
                    SO_P_pp(j,n,k,m) = P_pp(j,n,k,m,t) * pp(j,k,m,t);
                    %diff w.r.t T and T
                    SO_P_TT(j,n,k,m) = P_TT(j,n,k,m,t) * TT(j,k,m,t);
                    %diff w.r.t w and T
                    SO_P_wT(j,n,k,m) = P_wT(j,n,k,m,t) * wT(j,k,m,t);
                    %diff w.r.t p and T
                    SO_P_pT(j,n,k,m) = P_pT(j,n,k,m,t) * pT(j,k,m,t);
                    %diff w.r.t w and P
                    SO_P_wp(j,n,k,m) = P_wp(j,n,k,m,t) * wp(j,k,m,t);
                end
            end
        end
    end
SO_p(:,:,t) = (1/2) * sum(sum(SO_P_ww+SO_P_pp+SO_P_TT,4),3) + sum(sum(SO_P_wT+SO_P_pT+SO_P_wp,4),3);
end  

% Step 4d. Solve for total expenditure
for t = t1:TIME
    for j=1:J
        for i=1:N
            for n=1:N
                for m=1:N
                    %diff w.r.t. X and X
                    SO_X_XX(j,i,n,m) = X_XX(j,i,n,m,t) * XX(j,n,m,t);
                    %diff w.r.t. pi and pi
                    SO_X_pipi(j,i,n,m) = X_pipi(j,i,n,m,t) * pipi(n+(j-1)*N,m+(i-1)*N,t);
                    %diff w.r.t. pi and X
                    SO_X_piX(j,i,n,m) = X_piX(j,i,n,m,t) * piX(n+(j-1)*N,m+(i-1)*N,t);
                end
                for s=1:J
                    SO_X_piw(j,i,n,s) = X_piw(j,i,n,s,t) * piw(n+(j-1)*N,i+(s-1)*N,t);
                    SO_X_piL(j,i,n,s) = X_piL(j,i,n,s,t) * pil_ca(n+(j-1)*N,i+(s-1)*N,t);
                    SO_X_Xw(j,i,n,s) = X_Xw(j,i,n,s,t) * Xw(j,n,s,i,t);
                    SO_X_XL(j,i,n,s) = X_XL(j,i,n,s,t) * Xl(j,n,s,i,t);
                end
            end
            for k=1:J
                for s=1:J
                    %diff w.r.t. w and w
                    SO_X_ww(j,i,k,s) = X_ww(j,i,k,s,t) * ww_ca(k,i,s,t);
                    %diff w.r.t. L and L
                    SO_X_LL(j,i,k,s) = X_LL(j,i,k,s,t) * ll_ca(k,i,s,t);
                    %diff w.r.t. w and L
                    SO_X_wL(j,i,k,s) = X_wL(j,i,k,s,t) * wl_ca(k,i,s,t);
                end
            end
        end
    end
    SO_X(:,:,t) = (1/2) * (sum(sum(SO_X_XX+SO_X_pipi,4),3) + sum(sum(SO_X_ww+SO_X_LL,4),3)) + sum(sum(SO_X_piX,4),3) + sum(sum(SO_X_piw+SO_X_piL+SO_X_Xw+SO_X_XL,4),3) +  sum(sum(SO_X_wL,4),3);
end  
% Step 4e. Solve for an updated w_hat using the labor market clearing

for t = t1:TIME
    for n=1:N
        for i=1:N
            for j=1:J
                for m=1:N
                    SO_w_pipi(j,i,n,m) = w_pipi(j,i,n,m,t) * pipi(n+(j-1)*N,m+(i-1)*N,t);
                    SO_w_XX(j,i,n,m) = w_XX(j,i,n,m,t) * XX(j,n,m,t);
                    SO_w_piX(j,i,n,m) = w_piX(j,i,n,m,t) * piX(n+(j-1)*N,m+(i-1)*N,t);
                end
            end
        end
    end
    SO_w(:,:,t) = (1/2) * sum(sum(SO_w_pipi+SO_w_XX,4),3) + +sum(sum(SO_w_piX,4),3);
end

SO_dyn  = v2struct(SO_L, SO_v, SO_rw);
SO_temp = v2struct(SO_p, SO_X, SO_w);
end