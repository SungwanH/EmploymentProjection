function FIGURES(params, eqm_nlpf_HAT_SS, eqm_nlpf_HAT, eqm_nlpf_dd, eqm_nlpf_dd_belief, eqm_dgp, eqm_recur)

v2struct(params.envr);
v2struct(params.prod);
%% Belief 
T_belief = BELIEF(params, W_TRUE);



figure
hold on
title('Objective Productivity and Belief (sector 1)')
plot(1:TIME, permute(T(1,CHINA,1:TIME),[2,3,1]))
plot(1:TIME, permute(T_belief(1,CHINA,1:TIME,1),[2,3,4,1]),'--')
plot(3:TIME, permute(T_belief(1,CHINA,3:TIME,3),[2,3,4,1]),'--')
plot(5:TIME, permute(T_belief(1,CHINA,5:TIME,5),[2,3,4,1]),'--')
plot(10:TIME, permute(T_belief(1,CHINA,10:TIME,10),[2,3,4,1]),'--')
legend('Productivity','Belief at t=1','Belief at t=5','Belief at t=10','location','best')
saveas(gcf,'figures/prod_and_belief_sec1.png')
%{
figure
hold on
title('Productivity Belief across Sector')
plot(2:TIME, permute(T_belief(1,CHINA,2:TIME,2),[2,3,4,1]))
plot(2:TIME, permute(T_belief(2,CHINA,2:TIME,2),[2,3,4,1]))
plot(2:TIME, permute(T_belief(3,CHINA,2:TIME,2),[2,3,4,1]))
plot(2:TIME, permute(T_belief(4,CHINA,2:TIME,2),[2,3,4,1]))
saveas(gcf,'figures/belief_sectors.png')

figure
hold on
title('Objective Productivity and Belief (sector 2)')
plot(1:TIME, permute(T(2,CHINA,1:TIME),[2,3,1]))
plot(1:TIME, permute(T_belief(2,CHINA,1:TIME,1),[2,3,4,1]),'--')
plot(3:TIME, permute(T_belief(2,CHINA,3:TIME,3),[2,3,4,1]),'--')
plot(5:TIME, permute(T_belief(2,CHINA,5:TIME,5),[2,3,4,1]),'--')
plot(7:TIME, permute(T_belief(2,CHINA,7:TIME,7),[2,3,4,1]),'--')
plot(10:TIME, permute(T_belief(2,CHINA,10:TIME,10),[2,3,4,1]),'--')
legend('Productivity','Belief at t=1','Belief at t=5','Belief at t=10','location','best')
saveas(gcf,'figures/prod_and_belief_sec2.png')
%}
%% nonlinear perfect foresight
v_td = eqm_nlpf_HAT.v_td;
v_td_SS = eqm_nlpf_HAT_SS.v_td;
Ldynamic = permute(sum(eqm_nlpf_HAT.Ldyn,1),[2,3,1]);

US_manu_share_in_world=reshape(sum(eqm_nlpf_HAT.VALjn00(1,1:R,:),2)./sum(eqm_nlpf_HAT.VALjn00(1,:,:),2),TIME,1);
%LdynamicManu= reshape(sum(eqm_nlpf_HAT.Ldyn(1,:,:),2),TIME,1);
%LdynamicServ= reshape(sum(eqm_nlpf_HAT.Ldyn(4,:,:),2),TIME,1);

Ldynamic_SS = permute(sum(eqm_nlpf_HAT_SS.Ldyn,1),[2,3,1]);
for t=1:TIME-1
    Ltrend(:,t)=Ldynamic(:,t+1)-Ldynamic(:,t);
end
for t=1:TIME_SS-1
    Ltrend_SS(:,t)=Ldynamic_SS(:,t+1)-Ldynamic_SS(:,t);
end
realwages_fig = permute(sum(eqm_nlpf_HAT.realwages,1),[2,3,1]);
realwages_fig_SS = permute(sum(eqm_nlpf_HAT_SS.realwages,1),[2,3,1]);

figure
hold on
title('Perfect Foresight (Deriving SS): Labor (level)')
plot(1:TIME_SS-1,Ldynamic_SS(1:R,1:TIME_SS-1))
saveas(gcf,'figures/labor_level_SS.png')

figure
hold on
title('Perfect Foresight (Deriving SS):Labor converegence(Ldyn(t+1)-Ldyn(t))')
plot(1:TIME_SS-1,Ltrend_SS(:,1:TIME_SS-1))
saveas(gcf,'figures/labor_convergence_level_SS.png')

figure
hold on
title('Perfect Foresight: Value (time difference)')
plot(2:TIME,v_td(:,2:TIME))
saveas(gcf,'figures/Value_level.png')

figure
hold on
title('Perfect Foresight: Labor (level)')
plot(1:TIME-1,Ldynamic(1:R,1:TIME-1))
saveas(gcf,'figures/labor_level.png')


figure
hold on
title('Perfect Foresight: Labor converegence (Ldyn(t+1)-Ldyn(t))')
plot(1:TIME-1,Ltrend(1:R,1:TIME-1))
saveas(gcf,'figures/labor_convergence_level.png')

figure
hold on
title('Perfect Foresight: Real wages (time difference)')
 plot(1:TIME,realwages_fig(1:R,1:TIME))
saveas(gcf,'figures/realwage_level.png')


figure
hold on
title('Perfect Foresight: Labor Alabama (level)')
plot(1:TIME-1,Ldynamic(1,1:TIME-1))
saveas(gcf,'figures/labor_level_ALA.png')

figure
hold on
title('Perfect Foresight: Labor California (level)')
plot(1:TIME-1,Ldynamic(5,1:TIME-1))
saveas(gcf,'figures/labor_level_CAL.png')

figure
hold on
title('Perfect Foresight: Labor converegence California')
plot(1:TIME-1,Ltrend(5,1:TIME-1))
saveas(gcf,'figures/labor_convergence_level_CAL.png')

figure
hold on
title('Perfect Foresight: Real wages California')
 plot(1:TIME,realwages_fig(5,1:TIME))
saveas(gcf,'figures/realwage_level_CAL.png')
%{
figure
hold on
title('Perfect Foresight: US Manufacture Share in the World')
 plot(1:TIME,US_manu_share_in_world(1:TIME))
saveas(gcf,'figures/PF_US_manu_share.png')

figure
hold on
title('Perfect Foresight: US Manufacture Employment')
 plot(1:TIME,LdynamicManu(1:TIME))
saveas(gcf,'figures/PF_US_manu_emp.png')

figure
hold on
title('Perfect Foresight: US Service Employment')
 plot(1:TIME,LdynamicServ(1:TIME))
saveas(gcf,'figures/PF_US_serv_emp.png')
%}

%% DATA path (linearized)
%Ldyn, L_hat, v_hat, w_hat, p_hat, P_hat, pi_hat, mu_hat, X_hat, E_T_hat, L_dgp, w_dgp, p_dgp, P_dgp, pi_dgp, mu_dgp, X_dgp, VALjn00);
Ldyn_dgp = eqm_dgp.Ldyn;
Lhat_dgp = eqm_dgp.L_hat;
Ldynamic_dgp = permute(sum(Ldyn_dgp,1),[2,3,1]);
L_belief_dgp = eqm_dgp.L_belief_dgp;
L_belief_agg_dgp = permute(sum(L_belief_dgp,1),[2,3,4,1]);

figure
hold on
title('DATA: Labor in California (sector1)')
plot(1:TIME-1,permute(eqm_nlpf_HAT.Ldyn(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(Ldyn_dgp(1,5,1:TIME-1),[1,3,2]))
plot(1:TIME-1,permute(L_belief_dgp(1,5,1:TIME-1,1),[2,3,4,1]),'--')
plot(5:TIME-1,permute(L_belief_dgp(1,5,5:TIME-1,5),[2,3,4,1]),'--')
%plot(10:TIME-1,permute(L_belief_dgp(1,5,10:TIME-1,10),[2,3,4,1]),'--')
saveas(gcf,'figures/dgp_labor_CAL_sec1.png')


figure
hold on
title('DATA: Labor in Alabama (sector-aggregated)')
plot(1:TIME-1,Ldynamic(1,1:TIME-1))
plot(1:TIME-1,Ldynamic_dgp(1,1:TIME-1),'--')
plot(1:TIME-1,L_belief_agg_dgp(1,1:TIME-1,1),':')
plot(2:TIME-1,L_belief_agg_dgp(1,2:TIME-1,2),':')
plot(3:TIME-1,L_belief_agg_dgp(1,3:TIME-1,3),':')
plot(4:TIME-1,L_belief_agg_dgp(1,4:TIME-1,4),':')
plot(10:TIME-1,L_belief_agg_dgp(1,10:TIME-1,10),':')
legend('Nonlin PF','DATA','Belief','location','best')
saveas(gcf,'figures/dgp_labor_agg_ALA.png')
figure
hold on
title('DATA: Labor in Colorado (sector-aggregated)')
plot(1:TIME-1,Ldynamic(6,1:TIME-1))
plot(1:TIME-1,Ldynamic_dgp(6,1:TIME-1),'--')
plot(1:TIME-1,L_belief_agg_dgp(6,1:TIME-1,1),':')
plot(2:TIME-1,L_belief_agg_dgp(6,2:TIME-1,2),':')
plot(3:TIME-1,L_belief_agg_dgp(6,3:TIME-1,3),':')
plot(4:TIME-1,L_belief_agg_dgp(6,4:TIME-1,4),':')
plot(10:TIME-1,L_belief_agg_dgp(6,10:TIME-1,10),':')
legend('Nonlin PF','DATA','Belief','location','best')
saveas(gcf,'figures/dgp_labor_agg_COL.png')

figure
hold on
title('DATA: Labor in California (sector-aggregated)')
plot(1:TIME-1,Ldynamic(5,1:TIME-1))
plot(1:TIME-1,Ldynamic_dgp(5,1:TIME-1),'--')
plot(1:TIME-1,L_belief_agg_dgp(5,1:TIME-1,1),':')
plot(2:TIME-1,L_belief_agg_dgp(5,2:TIME-1,2),':')
plot(3:TIME-1,L_belief_agg_dgp(5,3:TIME-1,3),':')
plot(4:TIME-1,L_belief_agg_dgp(5,4:TIME-1,4),':')
plot(10:TIME-1,L_belief_agg_dgp(5,10:TIME-1,10),':')
legend('Nonlin PF','DATA','Belief','location','best')
saveas(gcf,'figures/dgp_labor_agg_CAL.png')
%{
figure
hold on
title('DATA: Labor in Texas (sector-aggregated)')
plot(1:TIME-1,Ldynamic(43,1:TIME-1))
plot(1:TIME-1,Ldynamic_dgp(43,1:TIME-1),'--')
plot(1:TIME-1,L_belief_agg_dgp(43,1:TIME-1,1),':')
plot(2:TIME-1,L_belief_agg_dgp(43,2:TIME-1,2),':')
plot(3:TIME-1,L_belief_agg_dgp(43,3:TIME-1,3),':')
plot(4:TIME-1,L_belief_agg_dgp(43,4:TIME-1,4),':')
plot(10:TIME-1,L_belief_agg_dgp(43,10:TIME-1,10),':')
legend('Nonlin PF','DATA','Belief','location','best')
saveas(gcf,'figures/dgp_labor_agg_TEX.png')
%}
%% Recursive method
%eqm_recur = v2struct(L_pf_lev, L_belief_lev, wf_pf, pf_pf, L_pf, w_pf, p_pf, pi_pf, mu_pf,...
%     wf_belief, pf_belief, L_belief, w_belief, p_belief, pi_belief, mu_belief);
L_pf_recur = eqm_recur.L_pf_lev;
L_belief_recur = eqm_recur.L_belief_lev;
L_pf_agg_recur = permute(sum(L_pf_recur,1),[2,3,4,1]);
L_belief_agg_recur = permute(sum(L_belief_recur,1),[2,3,4,1]);

figure
hold on
title('RECURSIVE: pct. dev. Labor in California (sector 1)')
plot(1:TIME-1,permute(eqm_recur.L_pf(1,5,1:TIME-1,1),[3,2,4,1]),':')
plot(1:TIME-1,permute(eqm_recur.L_belief(1,5,1:TIME-1,1),[3,2,4,1]))
plot(2:TIME-1,permute(eqm_recur.L_pf(1,5,2:TIME-1,2),[3,2,4,1]),':')
plot(2:TIME-1,permute(eqm_recur.L_belief(1,5,2:TIME-1,2),[3,2,4,1]))
plot(3:TIME-1,permute(eqm_recur.L_pf(1,5,3:TIME-1,3),[3,2,4,1]),':')
plot(3:TIME-1,permute(eqm_recur.L_belief(1,5,3:TIME-1,3),[3,2,4,1]))
plot(5:TIME-1,permute(eqm_recur.L_pf(1,5,5:TIME-1,5),[3,2,4,1]),':')
plot(5:TIME-1,permute(eqm_recur.L_belief(1,5,5:TIME-1,5),[3,2,4,1]))
plot(7:TIME-1,permute(eqm_recur.L_pf(1,5,7:TIME-1,7),[3,2,4,1]),':')
plot(7:TIME-1,permute(eqm_recur.L_belief(1,5,7:TIME-1,7),[3,2,4,1]))
plot(10:TIME-1,permute(eqm_recur.L_pf(1,5,10:TIME-1,10),[3,2,4,1]),':')
plot(10:TIME-1,permute(eqm_recur.L_belief(1,5,10:TIME-1,10),[3,2,4,1]))
%plot(10:TIME-1,L_belief_agg_dgp(5,10:TIME-1,10),'--')
legend('PF (RECOVERED)','Belief (RECOVERED)','location','best')
saveas(gcf,'figures/recur_labor_sec1_CAL_pct.png')

%{
figure
hold on
title('RECURSIVE: pct. dev. Labor in TEXAS (sector 1)')
plot(1:TIME-1,permute(eqm_recur.L_pf(1,43,1:TIME-1,1),[3,2,4,1]),':')
plot(1:TIME-1,permute(eqm_recur.L_belief(1,43,1:TIME-1,1),[3,2,4,1]))
plot(3:TIME-1,permute(eqm_recur.L_pf(1,43,3:TIME-1,3),[3,2,4,1]),':')
plot(3:TIME-1,permute(eqm_recur.L_belief(1,43,3:TIME-1,3),[3,2,4,1]))
plot(5:TIME-1,permute(eqm_recur.L_pf(1,43,5:TIME-1,5),[3,2,4,1]),':')
plot(5:TIME-1,permute(eqm_recur.L_belief(1,43,5:TIME-1,5),[3,2,4,1]))
plot(7:TIME-1,permute(eqm_recur.L_pf(1,43,7:TIME-1,7),[3,2,4,1]),':')
plot(7:TIME-1,permute(eqm_recur.L_belief(1,43,7:TIME-1,7),[3,2,4,1]))
plot(10:TIME-1,permute(eqm_recur.L_pf(1,43,10:TIME-1,10),[3,2,4,1]),':')
plot(10:TIME-1,permute(eqm_recur.L_belief(1,43,10:TIME-1,10),[3,2,4,1]))
%plot(10:TIME-1,L_belief_agg_dgp(5,10:TIME-1,10),'--')
legend('PF (RECOVERED)','Belief (RECOVERED)','location','best')
saveas(gcf,'figures/recur_labor_sec1_TEX_pct.png')
%}

figure
hold on
title('RECURSIVE: Labor in California (sector aggregated)')
plot(1:TIME-1,Ldynamic(5,1:TIME-1,1))
plot(1:TIME-1,Ldynamic_dgp(5,1:TIME-1))
plot(1:TIME-1,L_belief_agg_recur(5,1:TIME-1,1),':')
plot(1:TIME-1,L_belief_agg_dgp(5,1:TIME-1,1))
plot(1:TIME-1,L_pf_agg_recur(5,1:TIME-1,1),'--')
%plot(2:TIME-1,L_pf_agg_recur(43,2:TIME-1,2),'--')
%plot(3:TIME-1,L_pf_agg_recur(43,3:TIME-1,3),'--')
%plot(5:TIME-1,L_pf_agg_recur(5,5:TIME-1,5),'--')
%plot(5:TIME-1,L_belief_agg_recur(5,5:TIME-1,5),':')
plot(10:TIME-1,L_pf_agg_recur(5,10:TIME-1,10),'--')
%plot(1:TIME-1,L_belief_agg_recur(43,1:TIME-1,1),':')
%plot(2:TIME-1,L_belief_agg_recur(43,2:TIME-1,2),':')
%plot(3:TIME-1,L_belief_agg_recur(43,3:TIME-1,3),':')
plot(10:TIME-1,L_belief_agg_recur(5,10:TIME-1,10),':')
%plot(10:TIME-1,L_belief_agg_dgp(5,10:TIME-1,10))
%plot(10:TIME-1,L_belief_agg_dgp(5,10:TIME-1,10),'--')
legend('Nonlin PF','DATA','Belief (RECOVERED)','Belief (DGP)','PF (RECOVERED)','location','best')
saveas(gcf,'figures/recur_labor_agg_CAL.png')


%{
figure
hold on
title('RECURSIVE: Labor in TEXAS (sector aggregated)')
plot(1:TIME-1,Ldynamic(43,1:TIME-1,1))
plot(1:TIME-1,Ldynamic_dgp(43,1:TIME-1))
plot(1:TIME-1,L_belief_agg_recur(43,1:TIME-1,1),':')
plot(1:TIME-1,L_belief_agg_dgp(43,1:TIME-1,1))
plot(1:TIME-1,L_pf_agg_recur(43,1:TIME-1,1),'--')
%plot(2:TIME-1,L_pf_agg_recur(43,2:TIME-1,2),'--')
%plot(3:TIME-1,L_pf_agg_recur(43,3:TIME-1,3),'--')
%plot(5:TIME-1,L_pf_agg_recur(43,5:TIME-1,5),'--')
%plot(5:TIME-1,L_belief_agg_recur(43,5:TIME-1,5),':')
plot(10:TIME-1,L_pf_agg_recur(43,10:TIME-1,10),'--')
%plot(1:TIME-1,L_belief_agg_recur(43,1:TIME-1,1),':')
%plot(2:TIME-1,L_belief_agg_recur(43,2:TIME-1,2),':')
%plot(3:TIME-1,L_belief_agg_recur(43,3:TIME-1,3),':')
plot(10:TIME-1,L_belief_agg_recur(43,10:TIME-1,10),':')
%plot(10:TIME-1,L_belief_agg_dgp(43,10:TIME-1,10))
%plot(10:TIME-1,L_belief_agg_dgp(5,10:TIME-1,10),'--')
legend('Nonlin PF','DATA','Belief (RECOVERED)','Belief (DGP)','PF (RECOVERED)','location','best')
saveas(gcf,'figures/recur_labor_agg_TEX.png')
%}

figure
hold on
title('RECURSIVE: Labor in California (sector 1)')
plot(1:TIME-1,permute(eqm_nlpf_dd.Ldyn(1,5,1:TIME-1,1),[3,2,4,1]))
plot(1:TIME-1,permute(Ldyn_dgp(1,5,1:TIME-1),[3,2,1]))
plot(1:TIME-1,permute(L_belief_recur(1,5,1:TIME-1,1),[3,2,4,1]),':')
plot(1:TIME-1,permute(L_belief_dgp(1,5,1:TIME-1,1),[3,2,4,1]))
plot(1:TIME-1,permute(L_pf_recur(1,5,1:TIME-1,1),[3,2,4,1]),'--')
%plot(2:TIME-1,L_pf_agg_recur(43,2:TIME-1,2),'--')
%plot(3:TIME-1,L_pf_agg_recur(43,3:TIME-1,3),'--')
%plot(5:TIME-1,permute(L_pf_recur(1,5,5:TIME-1,5),[3,2,4,1]),'--')
%plot(5:TIME-1,permute(L_belief_recur(1,5,5:TIME-1,5),[3,2,4,1]),':')
plot(10:TIME-1,permute(L_pf_recur(1,5,10:TIME-1,10),[3,2,4,1]),'--')
%plot(1:TIME-1,L_belief_agg_recur(43,1:TIME-1,1),':')
%plot(2:TIME-1,L_belief_agg_recur(43,2:TIME-1,2),':')
%plot(3:TIME-1,L_belief_agg_recur(43,3:TIME-1,3),':')
plot(10:TIME-1,permute(L_belief_recur(1,5,10:TIME-1,10),[3,2,4,1]),':')
plot(10:TIME-1,permute(L_belief_dgp(1,5,10:TIME-1,10),[3,2,4,1]))
%plot(10:TIME-1,L_belief_agg_dgp(5,10:TIME-1,10),'--')
legend('Nonlin PF','DATA','Belief (RECOVERED)','Belief (DGP)','PF (RECOVERED)','location','best')
saveas(gcf,'figures/recur_labor_sec1_CAL.png')


%{
figure
hold on
title('RECURSIVE: Labor in California (sector 2)')
plot(1:TIME-1,permute(eqm_nlpf_HAT.Ldyn(2,5,1:TIME-1,1),[3,2,4,1]))
plot(1:TIME-1,permute(Ldyn_dgp(2,5,1:TIME-1),[3,2,1]))
plot(1:TIME-1,permute(L_belief_recur(2,5,1:TIME-1,1),[3,2,4,1]),':')
plot(1:TIME-1,permute(L_belief_dgp(2,5,1:TIME-1,1),[3,2,4,1]))
plot(1:TIME-1,permute(L_pf_recur(2,5,1:TIME-1,1),[3,2,4,1]),'--')
%plot(2:TIME-1,L_pf_agg_recur(43,2:TIME-1,2),'--')
%plot(3:TIME-1,L_pf_agg_recur(43,3:TIME-1,3),'--')
%plot(5:TIME-1,permute(L_pf_recur(2,5,5:TIME-1,5),[3,2,4,1]),'--')
plot(10:TIME-1,permute(L_pf_recur(2,5,10:TIME-1,10),[3,2,4,1]),'--')
%plot(1:TIME-1,L_belief_agg_recur(43,1:TIME-1,1),':')
%plot(2:TIME-1,L_belief_agg_recur(43,2:TIME-1,2),':')
%plot(3:TIME-1,L_belief_agg_recur(43,3:TIME-1,3),':')
%plot(5:TIME-1,permute(L_belief_recur(2,5,5:TIME-1,5),[3,2,4,1]),':')
plot(10:TIME-1,permute(L_belief_recur(2,5,10:TIME-1,10),[3,2,4,1]),':')
%plot(10:TIME-1,permute(L_belief_dgp(2,5,10:TIME-1,10),[3,2,4,1]))
%plot(10:TIME-1,L_belief_agg_dgp(5,10:TIME-1,10),'--')
legend('Nonlin PF','DATA','Belief (RECOVERED)','Belief (DGP)','PF (RECOVERED)','location','best')
saveas(gcf,'figures/recur_labor_sec2_CAL.png')



figure
hold on
title('RECURSIVE: Labor in Alabama (sector 3)')
plot(1:TIME-1,permute(eqm_nlpf_HAT.Ldyn(3,1,1:TIME-1,1),[3,2,4,1]))
plot(1:TIME-1,permute(Ldyn_dgp(3,1,1:TIME-1),[3,2,1]))
plot(1:TIME-1,permute(L_belief_recur(3,1,1:TIME-1,1),[3,2,4,1]),':')
plot(1:TIME-1,permute(L_belief_dgp(3,1,1:TIME-1,1),[3,2,4,1]))
plot(1:TIME-1,permute(L_pf_recur(3,1,1:TIME-1,1),[3,2,4,1]),'--')
%plot(2:TIME-1,L_pf_agg_recur(43,2:TIME-1,2),'--')
%plot(3:TIME-1,L_pf_agg_recur(43,3:TIME-1,3),'--')
%plot(5:TIME-1,permute(L_pf_recur(3,1,5:TIME-1,5),[3,2,4,1]),'--')
plot(10:TIME-1,permute(L_pf_recur(3,1,10:TIME-1,10),[3,2,4,1]),'--')
%plot(1:TIME-1,L_belief_agg_recur(43,1:TIME-1,1),':')
%plot(2:TIME-1,L_belief_agg_recur(43,2:TIME-1,2),':')
%plot(3:TIME-1,L_belief_agg_recur(43,3:TIME-1,3),':')
%plot(5:TIME-1,permute(L_belief_recur(3,1,5:TIME-1,5),[3,2,4,1]),':')
plot(10:TIME-1,permute(L_belief_recur(3,1,10:TIME-1,10),[3,2,4,1]),':')
plot(10:TIME-1,permute(L_belief_dgp(3,1,10:TIME-1,10),[3,2,4,1]))
%plot(10:TIME-1,L_belief_agg_dgp(5,10:TIME-1,10),'--')
legend('Nonlin PF','DATA','Belief (RECOVERED)','Belief (DGP)','PF (RECOVERED)','location','best')
saveas(gcf,'figures/recur_labor_sec3_ALA.png')

figure
hold on
title('RECURSIVE: Labor in TEXAS (sector 1)')
plot(1:TIME-1,permute(eqm_nlpf_HAT.Ldyn(1,43,1:TIME-1,1),[3,2,4,1]))
plot(1:TIME-1,permute(Ldyn_dgp(1,43,1:TIME-1),[3,2,1]))
plot(1:TIME-1,permute(L_belief_recur(1,43,1:TIME-1,1),[3,2,4,1]),':')
plot(1:TIME-1,permute(L_belief_dgp(1,43,1:TIME-1,1),[3,2,4,1]))
plot(1:TIME-1,permute(L_pf_recur(1,43,1:TIME-1,1),[3,2,4,1]),'--')
plot(10:TIME-1,permute(L_pf_recur(1,43,10:TIME-1,10),[3,2,4,1]),'--')
plot(10:TIME-1,permute(L_belief_recur(1,43,10:TIME-1,10),[3,2,4,1]),':')
legend('Nonlin PF','DATA','Belief (RECOVERED)','Belief (DGP)','PF (RECOVERED)','location','best')
saveas(gcf,'figures/recur_labor_sec1_TEX.png')
%}
end