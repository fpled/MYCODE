function [lambda,loglfval,loglfval0] = mleStoLinElasIsot(C_data,lambda0,varargin)
% function [lambda,loglfval,loglfval0] = mleStoLinElasIsot(C_data,lambda0,varargin)
% Maximum likelihood estimation for stochastic linear elastic tensor with isotropic symmetry
% C_data: data set for random vector C=(C1,C2)
% C_data(:,i): data for random coordinate Ci
% lambda0: initial parameter vector [la1, la2, la]
% lambda: optimal parameter vector [la1, la2, la]
% lambda(1) = la1 > 0
% lambda(2) = la2 > 0
% lambda(3) = la < 1/5
% loglfval: log-likelihood function evaluated at lambda: loglfval = loglfval = -nloglf(lambda,C_data)
% loglfval0: log-likelihood function evaluated at lambda0: loglfval0 = -nloglf(lambda0,C_data)

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
% display = 'off'; % default for gamfit and mle
% display = 'iter';
% display = 'final';

tolX = 1e-8; % tolerance on the parameter value (1e-8 by default for gamfit and 1e-6 for mle)
tolFun = 1e-8; % tolerance on the function value (1e-8 by default for gamfit and 1e-6 for mle)
tolBnd = 1e-6; % tolerance on the parameter bound (1e-6 by default for gamfit and mle)
maxIters = Inf; % maximum number of iterations
maxFunEvals = Inf; % maximum number of function evaluations
% optimFun = 'fminsearch'; % optimization function, 'fminsearch' (default) or 'fmincon'

options = statset(statset('mlecustom'),'Display',display,...
    'TolX',tolX,'TolFun',tolFun,'TolBnd',tolBnd,...
    'MaxIter',maxIters,'MaxFunEvals',maxFunEvals);

% n_data = size(C_data,1);
% C1_data = C_data(:,1);
% C2_data = C_data(:,2);
data = C_data(:); % data = [C1_data; C2_data];

% nloglf = @(lambda,data,cens,freq) -n_data*( (1-lambda(3))*log(lambda(1)) - gammaln(1-lambda(3)) )...
%     - n_data*( (1-5*lambda(3))*log(lambda(2)) - gammaln(1-5*lambda(3)) )...
%     + lambda(3)*( sum(log(data(1:n_data))) + 5*sum(log(data(n_data+1:end))) )...
%     + lambda(1)*sum(data(1:n_data))...
%     + lambda(2)*sum(data(n_data+1:end));
nloglf = @(lambda,data,cens,freq) nloglfElasIsot(lambda,data); % negative log-likelihood function
lambda = mle(data,'nloglf',nloglf,'Start',lambda0,...
    'LowerBound',lb,'UpperBound',ub,'Options',options);%,'OptimFun',optimFun);

if nargout>=2
    loglfval = -nloglf(lambda,data);
end
if nargout>=3
    loglfval0 = -nloglf(lambda0,data);
end

end

