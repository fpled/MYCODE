function [lambda,resnorm,residual,exitflag,output] = lseStoLinElasIsot(C_data,lambda0,varargin)
% function [lambda,resnorm,residual,exitflag,output] = lseStoLinElasIsot(C_data,lambda0,varargin)
% Least-squares estimation for stochastic linear elastic tensor with isotropic symmetry
% C_data: data set for random vector C=(C1,C2)
% C_data(:,i): data for random coordinate Ci
% lambda0: initial parameter vector [la1, la2, la] for full parameterization
%                                   [la]           for reduced parameterization
% lambda: optimal parameter vector [la1, la2, la] for full parameterization
%                                  [la]           for reduced parameterization
% For full parameterization:
% lambda(1) = la1 > 0
% lambda(2) = la2 > 0
% lambda(3) = la < 1/5
% For reduced parameterization:
% lambda(1) = la < 1/5
% resnorm: squared 2-norm of the residual at the solution lambda: sum(fun(lambda).^2) 
% with fun = @(lambda) funlsqnonlinElasIsot(lambda,mC_data,nuC_data)
% residual: residual (objective function value) at the solution lambda: residual = fun(lambda)
% exitflag: reason the solver (lsqnonlin) stopped, exit condition of lsqnonlin
% output: information about the optimization process

% Empirical estimates
% N_data = size(C_data,1);
mC_data = mean(C_data,1);
% vC_data = var(C_data,0,1);
% % vC_data = N_data/(N_data-1)*moment(C_data,2,1);
% stdC_data = std(C_data,1);
% cvC_data = stdC_data./mC_data;
% sC_data = sqrt(norm(vC_data));
% mCnorm_data = norm(mC_data);
% dC_data = sC_data/mCnorm_data;
phiC_data = log(96*C_data(:,1).*C_data(:,2).^5);
% phiC_data = log(96) + log(C_data(:,1)) + 5*log(C_data(:,2));
nuC_data = mean(phiC_data,1);

% Parameter bounds
useRedParam = isscalar(lambda0); % reduced parameterization
if useRedParam
    % reduced parameterization
    lb = -Inf; % lower bounds
    ub = 1/5;  % upper bounds
else
    % full parameterization
    lb = [0 0 -Inf];    % lower bounds
    ub = [Inf Inf 1/5]; % upper bounds
end

display = getcharin('display',varargin,'off');
% display = 'off';
% display = 'iter';
% display = 'iter-detailed';
% display = 'final'; % default for lsqnonlin
% display = 'final-detailed';

algo = 'trust-region-reflective'; % default for lsqnonlin
% algo = 'levenberg-marquardt';
% algo = 'interior-point';

tolX = eps; % tolerance on the parameter value (default: 1e-6 for lsqnonlin)
tolFun = eps; % tolerance on the function value (default: 1e-6 for lsqnonlin)
tolOpt = eps; % tolerance on the first-order optimality (default: 1e-6 for lsqnonlin)
maxIters = Inf; % maximum number of iterations
maxFunEvals = Inf; % maximum number of function evaluations

% options = optimoptions('lsqnonlin','Display',display,'Algorithm',algo);
% options = optimoptions('lsqnonlin','Display',display,'Algorithm',algo,...
%     'TolX',tolX,'TolFun',tolFun,...
%     'MaxIter',maxIters,'MaxFunEvals',maxFunEvals);
options = optimoptions('lsqnonlin','Display',display,'Algorithm',algo,...
    'StepTolerance',tolX,'FunctionTolerance',tolFun,'OptimalityTolerance',tolOpt,...
    'MaxIterations',maxIters,'MaxFunctionEvaluations',maxFunEvals);

fun = @(lambda) funlsqnonlinElasIsot(lambda,mC_data,nuC_data);
[lambda,resnorm,residual,exitflag,output] = lsqnonlin(fun,lambda0,lb,ub,options);

end
