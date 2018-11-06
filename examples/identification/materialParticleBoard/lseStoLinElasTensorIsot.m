function [lambda,err,exitflag,output] = lseStoLinElasTensorIsot(C_data)
% function [lambda,err,exitflag,output] = lseStoLinElasTensorIsot(C_data)
% Least-squares estimation for stochastic linear elastic tensor with
% isotropic symmetry
% C_data(:,1): data for random coordinate C1
% C_data(:,2): data for random coordinate C2
% lambda: parameters (la1,la2,la)
% lambda(1) = la1 > 0
% lambda(2) = la2 > 0
% lambda(3) = la < 1/5

% empirical estimates
mC_data = mean(C_data,1);
% vC_data = var(C_data,0,1);
% % vC_data = size(C_data,1)/(size(C_data,1)-1)*moment(C_data,2,1);
% sC_data = sqrt(norm(vC_data));
% dC_data = sC_data/norm(mC_data);
phiC_data = log(96*C_data(:,1).*C_data(:,2).^5);
nuC_data = mean(phiC_data,1);

% initial guess
la = 0; % la < 1/5
la1 = (1-la)/mC_data(1); % la1 > 0
la2 = (1-5*la)/mC_data(2); % la2 > 0

lambda0 = [la1 la2 la];
lb = [0 0 -Inf];
ub = [Inf Inf 1/5];

% display = 'off';
% display = 'iter';
% display = 'iter-detailed';
display = 'final';
% display = 'final-detailed';

tolX = 1e-12; % tolerance on the parameter value
tolFun = 1e-12; % tolerance on the function value
% maxFunEvals = 3e3; % maximum number of function evaluations

% options  = optimoptions('lsqnonlin','Display',display);
% options  = optimoptions('lsqnonlin','Display',display,'TolX',tolX,'TolFun',tolFun);
options  = optimoptions('lsqnonlin','Display',display,'StepTolerance',tolX,'FunctionTolerance',tolFun);

fun = @(lambda) funoptimlseIsot(lambda,mC_data,nuC_data);
[lambda,err,~,exitflag,output] = lsqnonlin(fun,lambda0,lb,ub,options);

end
