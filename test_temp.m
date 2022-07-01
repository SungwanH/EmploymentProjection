%% This test whether linearization or log linearization is a better fit for the LOM for labor
rng(22020612)
N=20;
l_0=exp(rand(N,1)*2);
l_0=5*l_0/sum(l_0);

v_1=rand(N,1);

BETA=0.95; % DISCOUNT RATE
RHO=1; %
rand_temp=rand(N,N)+1;
KAPPA=exp(rand_temp-diag(diag(rand_temp))); %MIGRATINO COST

params=v2struct(l_0,v_1,BETA,RHO,rand_temp,KAPPA,N);


%% baseline l_1;
[l_1_base,D,E]=test_temp_nonlinear(params);

%% counterfactual v_1 set to another random variable; l_0 set to another
% random variable
v_1_counter=v_1;%+0.1*rand(N,1) ;%  set this equal to v_1 to shut down the change in v_1; alternatively, set this to 0.1*rand(N,1) to test the level approximation
l_0_counter=exp(rand(N,1));
l_0_counter=5*l_0_counter/sum(l_0_counter);
params_counter=params;
params_counter.l_0=l_0_counter;
params_counter.v_1=v_1_counter;
l_1_counter=test_temp_nonlinear(params_counter);
pc_change_nonlinear=log(l_1_counter)-log(l_1_base);

%% log linearization under the counterfactual:
% NOTE THAT I DID NOT CODE UP THE V PART, FOR LOG LINEARIZATYION, SO THIS
% CODE IS NOT GOOD FOR TESTING V SHOCK; SEE LINE 20
pc_change_log=sum(E.*log(l_0_counter./l_0),1)';  
l_1_counter_log=exp(pc_change_log).*l_1_base;


pc_change_log2=sum(E.*(l_0_counter./l_0-1),1)';  
l_1_counter_log2=(pc_change_log2+1).*l_1_base;


%% linearization under the counterfactual;
l_1_counter_level_c1=sum(D.*(l_0_counter-l_0),1)';
l_1_counter_level_c2=zeros(N,1);
for jj=1:N
    for ii=1:N
        for kk=1:N
        l_1_counter_level_c2(jj,1)=l_1_counter_level_c2(jj,1)+l_0(ii,1)*BETA/RHO*D(ii,jj)*((kk==jj)-D(ii,kk))*(v_1_counter(kk)-v_1(kk));
        end
    end
end

l_1_counter_level=l_1_counter_level_c1+l_1_counter_level_c2+l_1_base;
pc_change_level=log(l_1_counter_level)-log(l_1_base);

% percent difference between log, level linearization, and nonlinear
figure
scatter(pc_change_nonlinear,pc_change_log)
hold on
scatter(pc_change_nonlinear,pc_change_log2)
plot((min(pc_change_nonlinear):0.001:max(pc_change_nonlinear)),(min(pc_change_nonlinear):0.001:max(pc_change_nonlinear)))
legend('percentage impact of the shock: log approximation', 'percentage error: percentage approximation')
hold off


% log error as a share of log changes
error_pc=(pc_change_log-pc_change_nonlinear)./(abs(pc_change_nonlinear))
error_pc2=(pc_change_log2-pc_change_nonlinear)./(abs(pc_change_nonlinear))


% level prediction
figure
scatter(l_1_counter,l_1_counter_log)
hold on
scatter(l_1_counter,l_1_counter_log2)
plot((min(l_1_counter):0.001:max(l_1_counter)),(min(l_1_counter):0.001:max(l_1_counter)))
legend('log deviation as input and output', 'percentage deviation as input and output')



