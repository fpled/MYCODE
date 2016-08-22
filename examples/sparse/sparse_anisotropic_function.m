%% Adaptive sparse polynomial approximation - Anisotropic function %%
%%-----------------------------------------------------------------%%
% [Chkifa Cohen Schwab 2014]

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');

%% Filename and Pathname
filename = 'sparse_anisotropic_function';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Scalar-valued anisotropic function
d = 16; % parametric dimension

% y = x_3*sin(x_4+x_16)
% fun = @(x) x(:,3).*sin(x(:,4)+x(:,16));
fun = vectorize('x3*sin(x4+x16)');

% y = 1/(1 + sum_{j=1}^d (g_j*x_j)) with g_j = 3/(5*j^3)
% fun = @(x) 1./(ones(size(x,1),1)+x*(3./(5.*(1:size(x,2)).^3))');

% y = 1/(1 + (sum_{j=1}^d (g_j*x_j))^2) with g_j = 5/(j^3)
% fun = @(x) 1./(ones(size(x,1),1)+(x*(5./((1:size(x,2)).^3))').^2);

% y = 1/(1 + sum_{j=1}^d (g_j*x_j)) with g_j = 10^(-j)
% fun = @(x) 1./(ones(size(x,1),1)+x*(10.^(-(1:size(x,2))))');

v = UniformRandomVariable(0,1);
rv = RandomVector(v,d);

% [fun,rv] = multivariateFunctionsBenchmarks('anisotropic',d);

fun = MultiVariateFunction(fun,d);
fun.evaluationAtMultiplePoints = true;

%% Adaptive sparse approximation using least-squares
p = 50;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,d);

s = AdaptiveSparseTensorAlgorithm();
% s.nbSamples = 1;
% s.addSamplesFactor = 0.1;
s.tol = 1e-6;
s.tolStagnation = 5e-2;
% s.tolOverfit = 1.1;
% s.bulkParameter = 0.5;
% s.adaptiveSampling = true;
% s.adaptationRule = 'reducedmargin';
s.maxIndex = p;
% s.display = true;
% s.displayIterations = true;

ls = LeastSquaresSolver();
ls.regularization = false;
% ls.regularizationType = 'l1';
ls.errorEstimation = true;
% ls.errorEstimationType = 'leaveout';
% ls.errorEstimationOptions.correction = true;

t = tic;
[f,err,~,y] = s.leastSquares(fun,bases,ls,rv);
time = toc(t);

%% Outputs
fprintf('\n')
fprintf('parametric dimension = %d\n',ndims(f.basis))
% fprintf('parametric dimension = %d\n',numel(rv))
fprintf('basis dimension = %d\n',numel(f.basis))
fprintf('order = [ %s ]\n',num2str(max(f.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(f.basis.indices.array)
fprintf('nb samples = %d\n',size(y,1))
fprintf('CV error = %d\n',norm(err))
fprintf('elapsed time = %f s\n',time)

Ntest = 1000;
[errtest,xtest,fxtest,ytest] = computeTestError(f,fun,Ntest,rv);
fprintf('test error = %d\n',errtest)

%% Display multi-index set
dim = [3 4 16];
plotMultiIndexSet(f,'dim',dim,'legend',false)
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);

dim = [1 2 4];
plotMultiIndexSet(f,'dim',dim,'legend',false)
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);
