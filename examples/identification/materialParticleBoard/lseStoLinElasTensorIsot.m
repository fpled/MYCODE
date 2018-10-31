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
mC = mean(C_data,1);
vC = var(C_data,0,1);
% vC = size(C_data,1)/(size(C_data,1)-1)*moment(C_data,2,1);
sC = sqrt(norm(vC));
dC = sC/norm(mC);

% initial guess
la = 0; % la < 1/5
la1 = (1-la)/mC(1); % la1 > 0
la2 = (1-5*la)/mC(2); % la2 > 0

lambda0 = [la1 la2 la];
lb = [0 0 -Inf];
ub = [Inf Inf 1/5];

display = 'off';
% display = 'iter';
% display = 'iter-detailed';
% display = 'final';
% display = 'final-detailed';

tolX = 1e-6; % tolerance on the parameter value
tolFun = 1e-6; % tolerance on the function value
tolOptimality = 1e-6; % tolerance on the first-order optimality
tolConstraint = 1e-6; % tolerance on the constraint violation
maxFunEvals = 3e3; % maximum number of function evaluations

options  = optimoptions('fmincon','Display',display);
% options  = optimoptions('fmincon','Display',display,'TolX',tolX,'TolFun',tolFun,'MaxFunEvals',maxFunEvals,'OptimalityTolerance',tolOptimality,'ConstraintTolerance',tolConstraint);
% options  = optimoptions('fmincon','Display',display,'StepTolerance',tolX,'FunctionTolerance',tolFun,'MaxFunctionEvaluations',maxFunEvals,'OptimalityTolerance',tolOptimality,'ConstraintTolerance',tolConstraint);

fun = @(lambda) funoptimlseIsot(lambda,mC,dC);
[lambda,err,exitflag,output] = fmincon(fun,lambda0,[],[],[],[],lb,ub,[],options);

end
