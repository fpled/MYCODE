function [lambda,err,exitflag,output] = lseStoLinElasTensorIsotTrans(C_data,MCMCalg)
% function [lambda,err,exitflag,output] = lseElasTensorIsot(C_data,MCMCalg)
% Least-squares estimation for stochastic linear elastic tensor with
% transversely isotropic symmetry
% C_data(:,1): data for random coordinate C1
% C_data(:,2): data for random coordinate C2
% C_data(:,3): data for random coordinate C3
% C_data(:,4): data for random coordinate C4
% C_data(:,5): data for random coordinate C5
% MCMCalg: algorithm for Markov-Chain Monte Carlo (MCMC) method
% MCMCalg = 'MH', 'BUM', 'CUM' or 'SS' ('SS' by default)
% lambda: parameters (la1,la2,la3,la4,la5,la)
% lambda(1) = la1 > 0
% lambda(2) = la2 > 0
% lambda(3) = la3 in R such that 2*sqrt(la1*la2)-la3 > 0
% lambda(4) = la4 > 0
% lambda(5) = la5 > 0
% lambda(6) = la < 1/2

if nargin==1 || isempty(MCMCalg)
    MCMCalg = 'SS';
end

% empirical estimates
mC = mean(C_data,1);
vC = var(C_data,0,1);
% vC = length(C_data)/(length(C_data)-1)*moment(C_data,2,1);
sC = sqrt(norm(vC));
dC = sC/norm(mC);

% initial guess
la = -100; % la < 1/2
la1 = -(mC(2)*la)/(mC(1)*mC(2)-mC(3)^2); % la1 > 0
la2 = -(mC(1)*la)/(mC(1)*mC(2)-mC(3)^2); % la2 > 0
la3 = (2*mC(3)*la)/(mC(1)*mC(2)-mC(3)^2); % la3 in R such that 2*sqrt(la1*la2)-la3 > 0
a = 1-2*la; % a > 0
la4 = a/mC(4); % la4 > 0
la5 = a/mC(5); % la5 > 0

lambda0 = [la1 la2 la3 la4 la5 la];
% f0 = funoptimlseIsotTrans(lambda0,C_data,mC,dC,MCMCalg);

lb = [0 0 -Inf 0 0 -Inf];
ub = [Inf Inf Inf Inf Inf 1/2];
nonlcon = @funconIsotTrans;

% display = 'off';
% display = 'iter';
display = 'iter-detailed';
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

fun = @(lambda) funoptimlseIsotTrans(lambda,C_data,mC,dC,MCMCalg);
[lambda,err,exitflag,output] = fmincon(fun,lambda0,[],[],[],[],lb,ub,nonlcon,options);

end
