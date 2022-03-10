function data = FAKE_DATA(params)
% GET_DATA returns data for a test problem
v2struct(params.envr);
VALjn00=ones(J,N);

% 20% import share
damp = 0.2 / (N-1);
pi_ni_0 = damp*ones(N,N);
pi_ni_0(eye(N,'logical')) = 1 - damp*(N-1);

Din00=pi_ni_0;
for jj=2:J
    Din00=[Din00;pi_ni_0];
end

L0=0.1*ones(J,R);

damp = 0.05 / (N*J-1);
mu0 = damp*ones(N*J,N*J);
mu0(eye(N*J,'logical')) = 1 - damp*(N*J-1);



data=v2struct(VALjn00,Din00,mu0,L0)
end