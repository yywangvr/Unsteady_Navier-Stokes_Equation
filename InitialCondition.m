function [ v_0,rho_0 ] = InitialCondition( X )
% Set initial condition for the coupled nonlinear dynamic problem
len=size(X,1);

%initial velocity is 0
v_0=zeros(len*2,1);

%initial density is 1-0.5*x
rho_0=ones(len,1)-0.5*X(:,1);


end

