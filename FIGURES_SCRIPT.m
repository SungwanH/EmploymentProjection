%This script makes figures

%% Arrnage data within the same structure 
% Aggregation by region: different DGP assumptions
eqm_dgp.L_region=reshape(sum(eqm_dgp.Ldyn,1),N,TIME);
eqm_dgp_w_01.L_region=reshape(sum(eqm_dgp_w_01.Ldyn,1),N,TIME);
eqm_dgp_w_09.L_region=reshape(sum(eqm_dgp_w_09.Ldyn,1),N,TIME);
eqm_dgp_w_re.L_region=reshape(sum(eqm_dgp_w_re.Ldyn,1),N,TIME);
eqm_nlpf_dd.L_region=reshape(sum(eqm_nlpf_dd.Ldyn,1),N,TIME);

% Assume the true DGP is W=0.5; collect belief and the corresponding
% counterfactual for PF
eqm_dgp.L_belief_region=reshape(sum(eqm_dgp.L_belief_dgp,1),N,TIME,ENDT+1);

eqm_dgp_w_09.L_belief_region=reshape(sum(eqm_dgp_w_09.L_belief_dgp,1),N,TIME,ENDT+1);
eqm_dgp_w_01.L_belief_region=reshape(sum(eqm_dgp_w_01.L_belief_dgp,1),N,TIME,ENDT+1);

eqm_recur.L_belief_region=reshape(sum(eqm_recur.L_belief_lev,1),N,TIME,ENDT);
eqm_recur.L_pf_region=reshape(sum(eqm_recur.L_pf_lev,1),N,TIME,ENDT);

eqm_recur_w_01.L_belief_region=reshape(sum(eqm_recur_w_01.L_belief_lev,1),N,TIME,ENDT);
eqm_recur_w_01.L_pf_region=reshape(sum(eqm_recur_w_01.L_pf_lev,1),N,TIME,ENDT);

eqm_recur_w_09.L_belief_region=reshape(sum(eqm_recur_w_09.L_belief_lev,1),N,TIME,ENDT);
eqm_recur_w_09.L_pf_region=reshape(sum(eqm_recur_w_09.L_pf_lev,1),N,TIME,ENDT);

eqm_recur_w_re.L_belief_region=reshape(sum(eqm_recur_w_re.L_belief_lev,1),N,TIME,ENDT);
eqm_recur_w_re.L_pf_region=reshape(sum(eqm_recur_w_re.L_pf_lev,1),N,TIME,ENDT);

%% Compare the evoluation of employment under different assumption on the DGP
Figure1=figure;
plot(1:TIME,eqm_dgp.L_region(1,1:TIME),'.','linewidth',2)
hold on
plot(1:TIME,eqm_dgp_w_01.L_region(1,1:TIME),':','linewidth',2)
plot(1:TIME,eqm_dgp_w_09.L_region(1,1:TIME),'--','linewidth',2)
plot(1:TIME,eqm_dgp_w_re.L_region(1,1:TIME),'Marker','+','linewidth',2)
plot(1:TIME,eqm_nlpf_dd.L_region(1,1:TIME),'linewidth',2)
xline(30);
legend('W=0,5','W=0.1','W=0.9','W=1','PF','Shocks stop','location','southeast');
print(Figure1,'figures/response_under_different_learning.png','-dpng','-r600');


%% Assume the true DGP W=0.5, comparing the recovered belief and the actual belief
Figure2=figure;
h(1)=plot(1:TIME,eqm_dgp.L_region(1,1:TIME),'linewidth',2);
hold on
h(2)=plot(1:TIME,eqm_dgp.L_belief_region(1,1:TIME,1),'.','linewidth',0.1,'Marker','+');
h(2)=plot(10:TIME,eqm_dgp.L_belief_region(1,10:TIME,10),'.','linewidth',0.1,'Marker','+');
h(2)=plot(15:TIME,eqm_dgp.L_belief_region(1,15:TIME,15),'.','linewidth',0.1,'Marker','+');
h(3)=plot(1:TIME,eqm_recur.L_belief_region(1,1:TIME,1),'--','linewidth',2);
h(3)=plot(10:TIME,eqm_recur.L_belief_region(1,10:TIME,10),'--','linewidth',2);
h(3)=plot(15:TIME,eqm_recur.L_belief_region(1,15:TIME,15),'--','linewidth',2);
legend(h([1,2 ,3]),{'Actual evolution','Belief held by agent', 'Belief from recursive algorithm'},'location','southeast');
print(Figure2,'figures/recover_belief.png','-dpng','-r600');
h(4)=plot(1:TIME,eqm_recur.L_pf_region(1,1:TIME,1),':','linewidth',2);
h(4)=plot(10:TIME,eqm_recur.L_pf_region(1,10:TIME,10),':','linewidth',2);
h(4)=plot(15:TIME,eqm_recur.L_pf_region(1,15:TIME,15),':','linewidth',2);
legend(h([1,2 ,3,4]),{'Actual evolution (Data)','Belief held by agent', 'Belief recovered (W=0.5)','Counterfactual PF'},'location','southeast');
print(Figure2,'figures/recover_belief_with_pf.png','-dpng','-r600');

Figure4=figure;
h(1)=plot(1:TIME,eqm_dgp.L_region(1,1:TIME),'linewidth',2);
hold on
h(2)=plot(1:TIME,eqm_dgp.L_belief_region(1,1:TIME,1),'.','linewidth',0.1,'Marker','+');
h(3)=plot(1:TIME,eqm_recur.L_belief_region(1,1:TIME,1),'--','linewidth',2);
h(4)=plot(1:TIME,eqm_recur_w_01.L_belief_region(1,1:TIME,1),'.','linewidth',0.1,'Marker','+');
h(5)=plot(1:TIME,eqm_recur_w_09.L_belief_region(1,1:TIME,1),'.','linewidth',0.1,'Marker','+');
legend(h([1,2 ,3,4,5]),{'Actual evolution (Data)','Belief held by agent','Belief recovered (W=0.5)', 'Belief recovered (W=0.1)','Belief recovered (W=0.9)'},'location','southeast');
print(Figure4,'figures/identify_belief_T1.png','-dpng','-r600');


Figure4=figure;
h(1)=plot(1:TIME,eqm_dgp.L_region(1,1:TIME),'linewidth',2);
hold on
h(2)=plot(10:TIME,eqm_dgp.L_belief_region(1,10:TIME,10),':', 'linewidth',2);
h(3)=plot(10:TIME,eqm_recur.L_belief_region(1,10:TIME,10),'--','linewidth',2);
h(4)=plot(10:TIME,eqm_recur_w_01.L_belief_region(1,10:TIME,10),'.','linewidth',0.1,'Marker','+');
h(5)=plot(10:TIME,eqm_recur_w_09.L_belief_region(1,10:TIME,10),'.','linewidth',0.1,'Marker','o');
legend(h([1,2 ,3,4,5]),{'Actual evolution (Data)','Belief held by agent','Belief recovered (W=0.5)', 'Belief recovered (W=0.1)','Belief recovered (W=0.9)'},'location','southeast');
print(Figure4,'figures/identify_belief_T10.png','-dpng','-r600');


Figure5=figure;
h(1)=plot(1:TIME,eqm_dgp.L_region(1,1:TIME),'linewidth',2);
hold on
h(2)=plot(15:TIME,eqm_dgp.L_belief_region(1,15:TIME,15),':', 'linewidth',2);
h(3)=plot(15:TIME,eqm_recur.L_belief_region(1,15:TIME,15),'--','linewidth',2);
h(4)=plot(15:TIME,eqm_recur_w_01.L_belief_region(1,15:TIME,15),'.','linewidth',0.1,'Marker','+');
h(5)=plot(15:TIME,eqm_recur_w_09.L_belief_region(1,15:TIME,15),'.','linewidth',0.1,'Marker','o');
legend(h([1,2 ,3,4,5]),{'Actual evolution (Data)','Belief held by agent','Belief recovered (W=0.5)', 'Belief recovered (W=0.1)','Belief recovered (W=0.9)'},'location','southeast');
print(Figure5,'figures/identify_belief_T15.png','-dpng','-r600');


%% Assume that the True W=0.5; econometricians making different assumption about agent information set end up recovering different impacts of a PF China shock
Figure6=figure;
h(1)=plot(1:TIME,eqm_dgp.L_region(1,1:TIME),'linewidth',2);
hold on
h(2)=plot(1:TIME,eqm_recur.L_pf_region(1,1:TIME,1),'--','linewidth',2);
h(3)=plot(1:TIME,eqm_recur_w_01.L_pf_region(1,1:TIME,1),':','linewidth',2);
h(4)=plot(1:TIME,eqm_recur_w_re.L_pf_region(1,1:TIME,1),'.','linewidth',0.1,'Marker','o');
legend(h([1,2 ,3,4]),{'Data (The effect of PF shock assuming DGP is PF)','The effect of PF shock assuming W=0.5','The effect of PF shock assuming W=0.1', 'The effect of PF shock assuming W=1'},'location','southeast');
print(Figure6,'figures/identify_the_effect_of_PF_shock.png','-dpng','-r600');






%{
Ldynamic = permute(sum(eqm_nlpf_HAT.Ldyn,1),[2,3,1]);
Ldynamic_sec = permute(eqm_nlpf_HAT.Ldyn(:,CHINA,:),[1,3,2]);

Ldynamic_sec_dd = permute(eqm_nlpf_dd.Ldyn(1,:,:),[1,3,2]);
Ldynamic_sec_dd_belief = permute(eqm_nlpf_dd_belief.Ldyn(:,CHINA,:),[1,3,2]);
Ldynamic_dgp = permute(eqm_dgp.Ldyn,[2,3,1]);
LdynamicManu= reshape(sum(eqm_nlpf_HAT.Ldyn(1,:,:),2),TIME,1);
LdynamicManu_cross= reshape(sum(eqm_nlpf_dd.Ldyn(1,:,:),2),TIME,1);
L_belief_dgp = eqm_dgp.L_belief_dgp;
L_belief_agg_dgp = permute(sum(L_belief_dgp,1),[2,3,4,1]);

Ldyn_dgp = eqm_dgp.Ldyn;
Lhat_dgp = eqm_dgp.L_hat;
L_belief_dgp = eqm_dgp.L_belief_dgp;

L_pf_recur = eqm_recur.L_pf_lev;
L_belief_recur = eqm_recur.L_belief_lev;
L_belief_recur_w_01 = eqm_recur_w_01.L_belief_lev;
L_belief_recur_w_09 = eqm_recur_w_09.L_belief_lev;
L_pf_agg_recur = permute(sum(L_pf_recur,1),[2,3,4,1]);
L_belief_agg_recur = permute(sum(L_belief_recur,1),[2,3,4,1]);
L_belief_agg_recur_w_01 = permute(sum(L_belief_recur_w_01,1),[2,3,4,1]);
L_belief_agg_recur_w_09 = permute(sum(L_belief_recur_w_09,1),[2,3,4,1]);



%% Plot the three different belief
figure
plot(2:TIME,L_belief_agg_dgp(1,2:TIME,1))
hold on
plot(2:TIME,L_belief_agg_dgp(1,2:TIME,3))
plot(2:TIME,L_belief_agg_dgp(1,2:TIME,5))
plot(2:TIME,L_belief_agg_dgp(1,2:TIME,7))
plot(2:TIME,L_belief_agg_dgp(1,2:TIME,10))

plot(2:TIME,L_pf_agg_recur(1,2:TIME,1),'--')
plot(2:TIME,L_pf_agg_recur(1,2:TIME,3),'--')
plot(2:TIME,L_pf_agg_recur(1,2:TIME,5),'--')
plot(2:TIME,L_pf_agg_recur(1,2:TIME,7),'--')
plot(2:TIME,L_pf_agg_recur(1,2:TIME,10),'--')


plot(1:TIME,Ldynamic_sec_dd(1,1:TIME,1),'LineWidth',2)



plot(2:TIME,L_belief_agg_recur(1,2:TIME,1),':')
hold on
plot(2:TIME,L_belief_agg_dgp(1,2:TIME,1))
plot(2:TIME,L_belief_agg_recur_w_09(1,2:TIME,1))
plot(2:TIME,L_belief_agg_recur_w_01(1,2:TIME,1))



%%

figure
hold on
title('DATA: Labor in region 1')
plot(1:TIME,Ldynamic_sec_dd(1,1:TIME))
plot(1:TIME,Ldynamic_dgp(1,1:TIME),'--')
plot(2:TIME,L_belief_agg_dgp(1,2:TIME,2),':')
plot(3:TIME,L_belief_agg_dgp(1,3:TIME,3),':')
plot(4:TIME,L_belief_agg_dgp(1,4:TIME,4),':')
plot(4:TIME,L_belief_agg_dgp(1,4:TIME,10),':')

plot(4:TIME,L_belief_agg_dgp(1,4:TIME,20),':')
%plot(10:TIME,L_belief_agg_dgp(1,10:TIME,10),':')
legend('Nonlin PF','DATA','Belief','location','best')
saveas(gcf,'figures/dgp_labor_region1.png')

figure
hold on
title('CHINA Employment')
plot(1:TIME,Ldynamic_sec(1,1:TIME))
plot(1:TIME,Ldynamic_sec_dd(1,1:TIME),'--')
plot(1:TIME,Ldynamic_sec_dd_belief(1,1:TIME),'--')
plot(1:TIME,L_belief_agg_dgp(1,1:TIME,1),':')
%plot(1:TIME,L_belief_agg_dgp(5,1:TIME,1),':')
legend('Baseline(Nonlinear TimeDiff)','Nonlinear Cross-Time(Objective T)','Nonlinear Cross-Time(Belief)','Linear Belief','location','best')
saveas(gcf,'figures/NLPF_All_Region1.png')
%}

%{
figure
hold on
title('RECURSIVE: Labor in region 1')
plot(1:TIME,Ldynamic_sec_dd(1,1:TIME,1))
plot(1:TIME,Ldynamic_dgp(1,1:TIME))
plot(2:TIME,L_belief_agg_recur(1,2:TIME,1),':')
plot(2:TIME,L_belief_agg_dgp(1,2:TIME,1))
plot(2:TIME,L_pf_agg_recur(1,2:TIME,1),'--')
%plot(10:TIME,L_belief_agg_recur(1,10:TIME,10),':')
%plot(10:TIME,L_belief_agg_dgp(1,10:TIME,10))
%plot(10:TIME,L_pf_agg_recur(1,10:TIME,10),'--')
legend('Nonlin PF','DATA','Belief (RECOVERED)','Belief (DGP)','PF (RECOVERED)','location','best')
figure
hold on
title('RECURSIVE: Labor in region 5')
plot(1:TIME,Ldynamic_sec_dd(1,1:TIME,5))
plot(1:TIME,Ldynamic_dgp(5,1:TIME))
plot(1:TIME,L_belief_agg_recur(5,1:TIME,1),':')
plot(1:TIME,L_belief_agg_dgp(5,1:TIME,1))
plot(1:TIME,L_pf_agg_recur(5,1:TIME,1),'--')
legend('Nonlin PF','DATA','Belief (RECOVERED)','Belief (DGP)','PF (RECOVERED)','location','best')

%FIGURES(params, eqm_nlpf_HAT_SS, eqm_nlpf_HAT, eqm_nlpf_dd, eqm_nlpf_dd_belief, eqm_dgp, eqm_recur)
%}

%{
logdiff = -(log(Ldynamic_sec(1,:)) - log(Ldynamic_sec_dd)); 
logdiff_new = -(log(Ldynamic_sec(1,:)) - log(Ldynamic_sec_cross_new)); 
logdiff_belief = -(log(Ldynamic_sec(1,:)) - log(L_belief_agg_dgp(:,:,2))); 

figure
hold on
title('change in log percentage CHINA Employment - NLPF and Linear Approx')
%plot(1:TIME,Ldynamic_sec(1,1:TIME))
plot(1:TIME,logdiff(1,1:TIME))
plot(1:TIME,logdiff_new(1,1:TIME),'--')
plot(1:TIME,logdiff_belief(1,1:TIME),':')
%plot(1:TIME,L_belief_agg_dgp(5,1:TIME,1),':')
legend('Nonlinear Time-Cross(Objective T)','Nonlinear Time-Cross(New code)','Linear Approximation','location','best')
saveas(gcf,'figures/Comparison_Approx_logdiff.png')
saveas(gcf,'figures/Comparison_Approx_logdiff.fig')

figure
hold on
title('CHINA Employment - NLPF and Linear Approx')
%plot(1:TIME,Ldynamic_sec(1,1:TIME))
plot(1:TIME,Ldynamic_sec_dd(1,1:TIME))
plot(1:TIME,Ldynamic_sec_cross_new(1,1:TIME),'--')
plot(1:TIME,L_belief_agg_dgp(1,1:TIME,1),'--')
%plot(1:TIME,L_belief_agg_dgp(5,1:TIME,1),':')
legend('Nonlinear Time-Cross(Objective T)','Nonlinear Time-Cross(New code)','Linear Approximation','location','best')
saveas(gcf,'figures/Comparison_Approx_Sym.png')
saveas(gcf,'figures/Comparison_Approx_Sym.fig')


figure
hold on
title('CHINA Employment -No Surprise Adj')
%plot(1:TIME,Ldynamic_sec(1,1:TIME))
plot(1:TIME,Ldynamic_sec_td(1,1:TIME))
plot(1:TIME,Ldynamic_sec_dd(1,1:TIME))
plot(1:TIME,Ldynamic_sec_cross_new(1,1:TIME),'--')
%plot(1:TIME,L_belief_agg_dgp(1,1:TIME,2),'--')
%plot(1:TIME,L_belief_agg_dgp(5,1:TIME,1),':')
legend('Nonlinear TimeDiff','Nonlinear Time-Cross(Base Code)','Nonlinear Time-Cross(New code)','location','best')
saveas(gcf,'figures/Comparison_No_Surprise.png')
saveas(gcf,'figures/Comparison_No_Surprise.fig')


%}

