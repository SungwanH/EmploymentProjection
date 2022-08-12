function [mom_dyn, mom_temp] = SO_MOMENTS(params, t1, eqm_sim)
% This code gives Second order moments from the first-order approximations


%% Roll down parameters
v2struct(params.envr)
v2struct(params.belief)


%% Roll down Baseline outcomes
v2struct(eqm_sim)
%v2struct(approx_sim)

% Remove all the first simulation results
%{
L_sim = L_sim(:,:,:,2:NUM_SIM);
v_sim = v_sim(:,:,2:NUM_SIM);
w_sim = w_sim(:,:,:,2:NUM_SIM);
P_sim = P_sim(:,:,:,2:NUM_SIM);
p_sim = p_sim(:,:,:,2:NUM_SIM);
X_sim = X_sim(:,:,:,2:NUM_SIM);
pi_sim = pi_sim(:,:,:,2:NUM_SIM);
E_T_hat_sim = E_T_hat_sim(:,:,:,:,2:NUM_SIM,:);
NUM_SIM = NUM_SIM-1;
%}
%% Dynamic part
% ww, wp, pp are derived in the temporary part
%L = reshape(L_sim,J*R,TIME,NUM_SIM);
for s=1:NUM_SIM
    for t=1:TIME
        L(:,t,s) = reshape(L_sim(:,:,t,s),[J*R,1,1,1]);
    end
end

for s=1:NUM_SIM
    for t=1:TIME-1
        v_temp = v_sim(:,t+1,s);
        L_temp = L(:,t,s);
        vv_temp(:,:,t,s) = v_temp * v_temp';
        ll_temp(:,:,t,s) = L_temp * L_temp';
        vl_temp(:,:,t,s) = v_temp * L_temp';
        
    end
end
vv_temp(:,:,TIME,:) = vv_temp(:,:,TIME-1,:);
ll_temp(:,:,TIME,:) = ll_temp(:,:,TIME-1,:);
vl_temp(:,:,TIME,:) = vl_temp(:,:,TIME-1,:);

vv = sum(vv_temp,4)/NUM_SIM;
ll = sum(ll_temp,4)/NUM_SIM;
vl = sum(vl_temp,4)/NUM_SIM;

for s=1:NUM_SIM
    for t=1:TIME
        %for jj = 1:J                
        w_dyn_temp = reshape(w_sim(:,:,t,s),[R*J,1]);
        P_temp = reshape(repmat(P_sim(1,:,t,s),[J,1,1]),[R*J,1]);
        ww_dyn_temp(:,t,s) = w_dyn_temp .* w_dyn_temp;
        PP_dyn_temp(:,t,s) = P_temp .* P_temp;
        wP_dyn_temp(:,t,s) = w_dyn_temp .* P_temp;
        %end
    end
end

ww_dyn = sum(ww_dyn_temp,3)/NUM_SIM;
PP_dyn = sum(PP_dyn_temp,3)/NUM_SIM;
wP_dyn = sum(wP_dyn_temp,3)/NUM_SIM;
%% Temporary part
% Price equation
for s=1:NUM_SIM
    for t=1:TIME
        if t<ENDT+1
            t2 = t;
        else
            t2 = ENDT+1;
        end
        for jj = 1:J                
        w_temp = w_sim(jj,:,t,s);
        p_temp = p_sim(jj,:,t,s);
        T_temp = E_T_hat_sim(jj,:,t,t2,s,t1); %Check this part; when the economy move to t1, we compute the belief deviation at t in which the belief is made at t1
        ww_temp(jj,:,:,t,s) = w_temp' * w_temp;
        pp_temp(jj,:,:,t,s) = p_temp' * p_temp;
        TT_temp(jj,:,:,t,s) = T_temp' * T_temp;
        wT_temp(jj,:,:,t,s) = w_temp' * T_temp;
        pT_temp(jj,:,:,t,s) = p_temp' * T_temp;
        wp_temp(jj,:,:,t,s) = w_temp' * p_temp;
        end
    end
end
ww_price = sum(ww_temp,5)/NUM_SIM;
pp = sum(pp_temp,5)/NUM_SIM;
TT = sum(TT_temp,5)/NUM_SIM;
wT = sum(wT_temp,5)/NUM_SIM;
pT = sum(pT_temp,5)/NUM_SIM;
wp = sum(wp_temp,5)/NUM_SIM;

% Current account 
ww_temp   = NaN(J,N,J,TIME,NUM_SIM);
ll_temp   = NaN(J,N,J,TIME,NUM_SIM);
wl_temp   = NaN(J,N,J,TIME,NUM_SIM);
Xw_temp   = NaN(J,N,J,N,TIME,NUM_SIM);
Xl_temp   = NaN(J,N,J,N,TIME,NUM_SIM);
XX_temp   = NaN(J,N,N,TIME,NUM_SIM);
pipi_temp = NaN(N*J,N*N,TIME,NUM_SIM);
piX_temp  = NaN(N*J,N*N,TIME,NUM_SIM);
piw_temp  = NaN(N*J,N*J,TIME,NUM_SIM);
pil_temp  = NaN(N*J,N*J,TIME,NUM_SIM);
for jj = 1:J                
    for nn=1:N
        for ss=1:J
                ww_temp(jj,nn,ss,:,:) = w_sim(jj,nn,:,:) .* w_sim(ss,nn,:,:);
                ll_temp(jj,nn,ss,:,:) = L_sim(jj,nn,:,:) .* L_sim(ss,nn,:,:);
                wl_temp(jj,nn,ss,:,:) = w_sim(jj,nn,:,:) .* L_sim(ss,nn,:,:);
            for ii=1:N
                Xw_temp(jj,nn,ss,ii,:,:) = X_sim(jj,nn,:,:) .* w_sim(ss,ii,:,:);
                Xl_temp(jj,nn,ss,ii,:,:) = X_sim(jj,nn,:,:) .* L_sim(ss,ii,:,:);
            end
        end
        for kk=1:N
            XX_temp(jj,nn,kk,:,:) = X_sim(jj,nn,:,:) .* X_sim(jj,kk,:,:);
        end
        for ii=1:N
            for mm=1:N
                pipi_temp(nn+(jj-1)*N,mm+(ii-1)*N,:,:) = pi_sim(nn+(jj-1)*N,ii,:,:) .* pi_sim(mm+(jj-1)*N,ii,:,:);
                piX_temp(nn+(jj-1)*N,mm+(ii-1)*N,:,:) = pi_sim(nn+(jj-1)*N,ii,:,:) .* X_sim(jj,mm,:,:);
            end
            for ss=1:J
                piw_temp(nn+(jj-1)*N,ii+(ss-1)*N,:,:) = pi_sim(nn+(jj-1)*N,ii,:,:) .* w_sim(ss,ii,:,:);
                pil_temp(nn+(jj-1)*N,ii+(ss-1)*N,:,:) = pi_sim(nn+(jj-1)*N,ii,:,:) .* L_sim(ss,ii,:,:);
            end
        end
        
    end
end

ww_ca = sum(ww_temp,5)/NUM_SIM;
ll_ca = sum(pp_temp,5)/NUM_SIM;
wl_ca = sum(TT_temp,5)/NUM_SIM;
Xw = sum(Xw_temp,6)/NUM_SIM;
Xl = sum(Xl_temp,6)/NUM_SIM;
XX = sum(XX_temp,5)/NUM_SIM;
pipi = sum(pipi_temp,4)/NUM_SIM;
piX = sum(piX_temp,4)/NUM_SIM;
piw = sum(piw_temp,4)/NUM_SIM;
pil_ca = sum(pil_temp,4)/NUM_SIM;

% Mkt clearing condition
%pipi, XX, piX are derived in current account
piw_temp_mc  = NaN(N*J,N*N,TIME,NUM_SIM);
for jj = 1:J                
    for nn=1:N
        for mm=1:N
            piw_temp_mc(nn+(jj-1)*N,ii+(mm-1)*N,:,:) = pi_sim(nn+(jj-1)*N,ii,:,:) .* w_sim(jj,mm,:,:);
            pil_temp_mc(nn+(jj-1)*N,ii+(mm-1)*N,:,:) = pi_sim(nn+(jj-1)*N,ii,:,:) .* L_sim(jj,mm,:,:);
        end
    end
end
piw_mc = sum(piw_temp,4)/NUM_SIM;
pil_mc = sum(pil_temp,4)/NUM_SIM;
mom_dyn    = v2struct(vv, ll, vl, ww_dyn, PP_dyn, wP_dyn);
mom_temp   = v2struct(ww_price, pp, wp, TT, wT, pT, ww_ca, ll_ca, wl_ca, Xw, Xl, XX, pipi, piX, piw, pil_ca, piw_mc, pil_mc);

%toc

