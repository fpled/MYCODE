%% Adaptive sparse polynomial approximation - Anisotropic function %%
%%-----------------------------------------------------------------%%
% [Chkifa Cohen Schwab 2014]

% clc
clear all
close all

%% Filename and Pathname
d = 16; % number of random variables
filename = ['sparse_approx_anisotropic_function_nbvar_' num2str(d)];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'RESULTS',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
% set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object

%% Scalar-valued Anisotropic function
% y = x_3*sin(x_4+x_16)
% fun = @(x) (x(:,3).*sin(x(:,4)+x(:,16)));
fun = vectorize('x3*sin(x4+x16)');

% y = 1/(1 + sum_{j=1}^d (g_j*x_j)) with g_j = 3/(5*j^3)
% fun = @(x) (1./(ones(size(x,1),1)+x*(3./(5.*(1:size(x,2)).^3))'));

% y = 1/(1 + (sum_{j=1}^d (g_j*x_j))^2) with g_j = 5/(j^3)
% fun = @(x) (1./(ones(size(x,1),1)+(x*(5./((1:size(x,2)).^3))').^2));

% y = 1/(1 + sum_{j=1}^d (g_j*x_j)) with g_j = 10^(-j)
% fun = @(x) (1./(ones(size(x,1),1)+x*(10.^(-(1:size(x,2))))'));

fun = MultiVariateFunction(fun,d);
fun.evaluationAtMultiplePoints = true;

rv = RandomVector(UniformRandomVariable(0,1),d);

h = PolynomialFunctionalBasis(LegendrePolynomials(),0:15);
H = FunctionalBases(h,d);

%% Adaptive sparse approximation using least-squares
s = AdaptiveSparseTensorAlgorithm();
% s.addSamplesFactor = 0.1;
s.tol = 1e-4;
% s.tolStagnation = 5e-2;
% s.tolOverfit = 1.1;
% s.bulkParameter = 0.5;
% s.adaptiveSampling = true;
% s.adaptationRule = 'reducedmargin';
% s.display = true;
s.maxIndex = 15;

ls = LeastSquaresSolver();
ls.regularization = false;
% ls.regularization = true;
% ls.regularizationOptions = struct('lambda',0);
% ls.regularizationType = 'l1';
% ls.modelSelection = true;
% ls.modelSelectionOptions.stopIfErrorIncrease = false;
% ls.errorEstimation = true;
% ls.errorEstimationType = 'leaveout';
ls.errorEstimationOptions.correction = false;
% ls.basisAdaptation = false;
% ls.basisAdaptationPath = [];
% ls.solver = '\';
% ls.solver = 'qr';
% ls.options = struct();

% rng('default')

t = tic;
[f,err,N] = s.leastSquares(@(x) functionEval(fun,x),H,ls,rv);
time = toc(t);

%% Results
fprintf('nb rand vars = %d\n',ndims(f.basis))% fprintf('nb rand vars = %d\n',numel(rv))
fprintf('dimension = %d\n',numel(f.basis))
fprintf('order = [ %s ]\n',num2str(max(f.basis.indices.array)))
fprintf('multi-index set = \n')
disp(f.basis.indices.array)
fprintf('nb samples = %d\n',N)
fprintf('CV error = %d\n',err)
fprintf('elapsed time = %f s\n',time)

xtest = random(rv,300,1);xtest=[xtest{:}];
ytest = fun.functionEval(xtest);
fxtest = f.functionEval(xtest);
error = norm(ytest-fxtest)/norm(ytest);
fprintf('test error = %d\n',error)

%% Display multi-index set
dim = [3 4 16];
plot(f.basis.indices,'dim',dim,'legend',false)
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);

dim = [1 2 4];
plot(f.basis.indices,'dim',dim,'legend',false)
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);
