function [l_1,D,E] = test_temp_nonlinear(params)
v2struct(params)
D=(exp(BETA*v_1'*ones(N,1))./KAPPA)./(sum(exp(BETA*v_1'*ones(N,1))./KAPPA,2)); % verify that sum(D,2)=1

l_1=sum(l_0.*D,1)';

E=l_0.*D./l_1';

end