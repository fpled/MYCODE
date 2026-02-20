function [lambda,fval,exitflag,output,varargout] = lseStoLinElasIsotTrans(C_data,lambda0,N,MCMCalgo,s,varargin)
% function [lambda,fval,exitflag,output,varargout] = lseStoLinElasIsotTrans(C_data,lambda0,N,MCMCalgo,s,varargin)
% Least-squares estimation for stochastic linear elastic tensor with transversely isotropic symmetry
% C_data: data set for random vector C=(C1,C2,C3,C4,C5)
% C_data(:,i): data for random coordinate Ci
% lambda0: initial parameter vector [la1, la2, la3, la4, la5, la] for full parameterization
%                                   [la1, la2, la3, la]           for reduced parameterization
% N: number of samples
% MCMCalgo: algorithm for Markov Chain Monte Carlo (MCMC) method
% MCMCalgo = 'IMH', 'RWMH' or 'SS' ('RWMH' by default)
% s: random number generator settings, specified as a structure with fields Type, Seed, and State (current settings by default)
% lambda: optimal parameter vector [la1, la2, la3, la4, la5, la] for full parameterization
%                                  [la1, la2, la3, la]           for reduced parameterization
% lambda(1) = la1 > 0
% lambda(2) = la2 > 0
% lambda(3) = la3 in R such that 2*sqrt(la1*la2)-la3 > 0
% For full parameterization:
% lambda(4) = la4 > 0
% lambda(5) = la5 > 0
% lambda(6) = la < 1/2, la < 0 for MH sampling using a trivariate normal proposal pdf with positive-definite covariance matrix Sigma
% For reduced parameterization:
% lambda(4) = la < 1/2, la < 0 for MH sampling using a trivariate normal proposal pdf with positive-definite covariance matrix Sigma
% fval: objective function value (residual) at the solution lambda: fval = fun(lambda)
% with fun = @(lambda) funoptimlseElasIsotTrans(lambda,C_data,mC_data,nuC_data,N,MCMCalgo,s)
% exitflag: reason the solver (fmincon) stopped, exit condition of fmincon
% output: information about the optimization process

if nargin<4 || isempty(MCMCalgo)
    MCMCalgo = 'RWMH'; % Random Walk Metropolis-Hastings algorithm
end
if nargin<5 || isempty(s)
    s = rng; % get the current random number generator settings
end

% Empirical estimates
% N_data = size(C_data,1);
mC_data = mean(C_data,1);
% vC_data = var(C_data,0,1);
% % vC_data = N_data/(N_data-1)*moment(C_data,2,1);
% stdC_data = std(C_data,0,1);
% cvC_data = stdC_data./mC_data;
% sC_data = sqrt(norm(vC_data));
% mCnorm_data = norm(mC_data);
% dC_data = sC_data/mCnorm_data;
phiC_data = log((C_data(:,1).*C_data(:,2)-C_data(:,3).^2).*(C_data(:,4).^2).*(C_data(:,5).^2));
% phiC_data = log(C_data(:,1).*C_data(:,2)-C_data(:,3).^2) + 2*log(C_data(:,4)) + 2*log(C_data(:,5));
nuC_data = mean(phiC_data,1);

% Parameter bounds
useRedParam = (length(lambda0)==4); % reduced parameterization
if useRedParam
    % reduced parameterization
    lb = [0 0 -Inf -Inf];   % lower bounds
    % ub = [Inf Inf Inf 1/2]; % upper bounds
    ub = [Inf Inf Inf 0];   % upper bounds
else
    % full parameterization
    lb = [0 0 -Inf 0 0 -Inf];       % lower bounds
    % ub = [Inf Inf Inf Inf Inf 1/2]; % upper bounds
    ub = [Inf Inf Inf Inf Inf 0];   % upper bounds
end

nonlcon = @funconIsotTrans; % nonlinear constraints

display = getcharin('display',varargin,'off');
% display = 'off';
% display = 'iter';
% display = 'iter-detailed';
% display = 'notify';
% display = 'notify-detailed';
% display = 'final'; % default for fmincon
% display = 'final-detailed';

algo = 'interior-point'; % default for fmincon
% algo = 'trust-region-reflective'; % does not work, as it requires to provide the gradient in the objective function, it allows only bounds or only linear equality constraints (but not both), and it does not accept nonlinear constraints.
% algo = 'sqp';
% algo = 'sqp-legacy';
% algo = 'active-set';

tolX = eps; % tolerance on the parameter value (1e-10 by default for fmincon and 1e-6 for patternsearch)
tolFun = eps; % tolerance on the function value (1e-6 by default for fmincon and patternsearch)
tolOpt = eps; % tolerance on the first-order optimality (1e-6 by default for fmincon)
% tolCon = eps; % tolerance on the constraint violation (1e-6 by default for fmincon and patternsearch)
tolCon = 0; % tolerance on the constraint violation (1e-6 by default for fmincon and patternsearch)
maxIters = Inf; % maximum number of iterations
maxFunEvals = Inf; % maximum number of function evaluations
finDiffType = 'central'; % finite differences, 'forward' or 'central' ('forward' by default)
useParallel = getcharin('useParallel',varargin,false); % parallel computing to estimates gradients of objective function (false by default)
plotFcn = getcharin('plotFcn',varargin,[]); % plot functions called at each iteration ([] by default)
% plotFcn = {'optimplot',...
%     'optimplotx','optimplotfunccount','optimplotfval','optimplotfvalconstr',...
%     'optimplotconstrviolation','optimplotstepsize','optimplotfirstorderopt'}; % plot functions called at each iteration

% options = optimoptions('fmincon','Display',display,'Algorithm',algo,'SpecifyConstraintGradient',true,'UseParallel',useParallel);
% options = optimoptions('fmincon','Display',display,'Algorithm',algo,...
%     'TolX',tolX,'TolFun',tolFun,'TolCon',tolCon,...
%     'MaxIter',maxIters,'MaxFunEvals',maxFunEvals,...
%     'SpecifyConstraintGradient',true,'UseParallel',useParallel);
options = optimoptions('fmincon','Display',display,'Algorithm',algo,...
    'StepTolerance',tolX,'FunctionTolerance',tolFun,'OptimalityTolerance',tolOpt,'ConstraintTolerance',tolCon,...
    'MaxIterations',maxIters,'MaxFunctionEvaluations',maxFunEvals,...
    'SpecifyConstraintGradient',true,'FiniteDifferenceType',finDiffType,'UseParallel',useParallel,'PlotFcn',plotFcn);
if strcmpi(algo,'active-set')
    options.TolConSQP = 0;
end

% Check first derivative (gradient) function against finite-difference
% approximation for the nonlinear constraint function funconIsotTrans
% [valid,err] = checkGradients(@funconIsotTrans,lambda0,'ISConstraint',true,Display="on");

if useParallel
    myparallel('start');
end

% Objective function
fun = @(lambda) funoptimlseElasIsotTrans(lambda,C_data,mC_data,nuC_data,N,MCMCalgo,s);

globalSolver = any(ischarin({'PatternSearch','GlobalSearch','MultiStart'},varargin));

if ~globalSolver
    % Single Local Solution (Single Local Minimum) via fmincon (Optimization Toolbox) solver
    [lambda,fval,exitflag,output] = fmincon(fun,lambda0,[],[],[],[],lb,ub,nonlcon,options);
    
    % problem = createOptimProblem('fmincon','objective',fun, ...
    %     'x0',lambda0, 'lb',lb,'ub',ub, ...
    %     'nonlcon',nonlcon,'options',options);
    % [lambda,fval,exitflag,output] = fmincon(problem);
    
    % Check if the purported (candidate) solution lambda is a local solution with patternsearch
    % % displayp = 'off';
    % displayp = 'iter';
    % % displayp = 'diagnose';
    % % displayp = 'final'; % default for patternsearch
    % 
    % algop = 'classic'; % default for patternsearch
    % % algop = 'nups';
    % % algop = 'nups-gps';
    % % algop = 'nups-mads';
    % 
    % meshTol = 1e-6; % tolerance on the mesh size (1e-6 by default for patternsearch)
    % plotFcnp = {'psplotbestf','psplotfuncount','psplotmeshsize','psplotbestx','psplotmaxconstr'}; % plot functions for patternsearch
    % 
    % optionsp = optimoptions('patternsearch','Display',displayp,'Algorithm',algop,...
    %     'StepTolerance',tolX,'FunctionTolerance',tolFun,'ConstraintTolerance',tolCon,'MeshTolerance',meshTol,...
    %     'MaxIterations',maxIters,'MaxFunctionEvaluations',maxFunEvals,...
    %     'UseParallel',useParallel,'UseCompletePoll',true,'UseCompleteSearch',true,'UseVectorized',false,'PlotFcn',plotFcnp);
    % 
    % % Set the candidate solution lambda as the start point
    % [lambdap,fvalp,exitflagp,outputp] = patternsearch(fun,lambda,[],[],[],[],lb,ub,nonlcon,optionsp);
    % 
    % % problem.solver = 'patternsearch';
    % % problem.x0 = lambda;
    % % problem.options = optionsp;
    % % [lambdap,fvalp,exitflagp,outputp] = patternsearch(problem);
elseif ischarin('PatternSearch',varargin)
    % Single Local Solution via patternsearch (Global Optimization Toolbox) solver
    
    % Properties (Global Options) for patternsearch
    % displayp = 'off';
    % displayp = 'iter';
    displayp = 'diagnose';
    % displayp = 'final'; % default for patternsearch
    
    algop = 'classic'; % default for patternsearch
    % algop = 'nups';
    % algop = 'nups-gps';
    % algop = 'nups-mads';
    
    meshTol = 1e-6; % tolerance on the mesh size (1e-6 by default for patternsearch)
    plotFcnp = {'psplotbestf','psplotfuncount','psplotmeshsize','psplotbestx','psplotmaxconstr'}; % plot functions for patternsearch
    
    optionsp = optimoptions('patternsearch','Display',displayp,'Algorithm',algop,...
        'StepTolerance',tolX,'FunctionTolerance',tolFun,'ConstraintTolerance',tolCon,'MeshTolerance',meshTol,...
        'MaxIterations',maxIters,'MaxFunctionEvaluations',maxFunEvals,...
        'UseParallel',useParallel,'UseCompletePoll',true,'UseCompleteSearch',true,'UseVectorized',false,'PlotFcn',plotFcnp);
    
    [lambda,fval,exitflag,output] = patternsearch(fun,lambda0,[],[],[],[],lb,ub,nonlcon,optionsp);
    
    % problem = createOptimProblem('patternsearch','objective',fun, ...
    %     'x0',lambda0, 'lb',lb,'ub',ub, ...
    %     'nonlcon',nonlcon,'options',optionsp);
    % [lambda,fval,exitflag,output] = patternsearch(problem);
else
    % Single Global Solution (Single Global Minimum) or Multiple Local Solutions (Multiple Local Minima)
    % via GlobalSearch or MultipleStart (Global Optimization Toolbox) solver
    
    % Create Problem Structure
    problem = createOptimProblem('fmincon','objective',fun, ...
        'x0',lambda0, 'lb',lb,'ub',ub, ...
        'nonlcon',nonlcon,'options',options);
    
    % Validate the problem structure by running the solver on the structure
    % [lambda,fval,exitflag,output] = fmincon(problem);
    
    % Properties (Global Options) for GlobalSearch and MultiStart
    % displayg = 'off';
    displayg = 'iter';
    % displayg = 'final'; % default for GlobalSearch and MultiStart
    
    tolXg = 1e-12; % tolerance on the parameter value (1e-6 by default for GlobalSearch and MultiStart)
    tolFung = 1e-12; % tolerance on the function value (1e-6 by default for GlobalSearch and MultiStart)
    maxTime = Inf; % maximum time (Inf by default for GlobalSearch and MultiStart)
    startPointsToRun = 'bounds-ineqs'; % start points to run, 'all', 'bounds' or 'bounds-ineqs' ('all' by default for GlobalSearch)
    numTrialPoints = 1e3; % number of trial points to examine (only for GlobalSearch, 1000 by default)
    useParallelg = getcharin('useParallel',varargin,false); % distribute local solver calls to multiple processors (only for MultiStart, false by default)
    if useParallelg && ischarin('MultiStart',varargin)
        plotFcng = [];
        problem.options.PlotFcn = []; % disable plot functions for fmincon
        problem.options.UseParallel = false; % disable parallel computing for fmincon
    else
        plotFcng = {'gsplotbestf','gsplotfunccount'}; % plot functions for GlobalSearch and MultiStart
    end
    
    % Create Solver Object
    if ischarin('GlobalSearch',varargin)
        % Construct a GlobalSearch object
        gs = GlobalSearch('Display',displayg,...
            'XTolerance',tolXg,'FunctionTolerance',tolFung,...
            'MaxTime',maxTime,'StartPointsToRun',startPointsToRun,...
            'NumTrialPoints',numTrialPoints,'PlotFcn',plotFcng); % GlobalSearch object
        % ms = MultiStart(gs); % MultiStart object based on our GlobalSearch attributes
        % ms.UseParallel = useParallelg;
    elseif ischarin('MultiStart',varargin)
        % Construct a MultiStart object
        ms = MultiStart('Display',displayg,...
            'XTolerance',tolXg,'FunctionTolerance',tolFung,...
            'MaxTime',maxTime,'StartPointsToRun',startPointsToRun,...
            'UseParallel',useParallelg,'PlotFcn',plotFcng); % MultiStart object
    end
    
    % Run multi-start object
    if ischarin('GlobalSearch',varargin)
        [lambda,fval,exitflag,output,solutions] = run(gs,problem);
    elseif ischarin('MultiStart',varargin)
        k = 50; % number of start points
        [lambda,fval,exitflag,output,solutions] = run(ms,problem,k);
    end
end

if useParallel
    myparallel('stop');
end

if nargout>4 && any(ischarin({'GlobalSearch','MultiStart'},varargin))
    varargout{1} = solutions;
end

end
