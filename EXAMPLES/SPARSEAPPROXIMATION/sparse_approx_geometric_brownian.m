%% Adaptive sparse polynomial approximation - Geometric brownian %%
%%---------------------------------------------------------------%%

% clc
% clear all
close all

%% Filename and Pathname
d = 10; % number of random variables
filename = ['sparse_approx_geometric_brownian_nbvar_' num2str(d)];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'RESULTS',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
% set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object

%% Vector-valued Geometric brownian
fun = @(x) geometric_brownian_kl(x,-1,0.5,1,100);
% q = 2;
% fun = @(x) (1/(2^(size(x,2)))*prod(3*x.^q+1,2));

fun = MultiVariateFunction(fun,d);
fun.evaluationAtMultiplePoints = true;

rv = RandomVector(UniformRandomVariable(0,1),d);

h = PolynomialFunctionalBasis(LegendrePolynomials(),0:15);
H = FunctionalBases(h,d);

%% Adaptive sparse approximation using least-squares
s = AdaptiveSparseTensorAlgorithm();
% s.addSamplesFactor = 0.1;
s.tol = 1e-3;
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
fprintf('CV error = %d\n',norm(err))
fprintf('elapsed time = %f s\n',time)

xtest = random(rv,300,1);xtest=[xtest{:}];
ytest = fun.functionEval(xtest);
fxtest = f.functionEval(xtest);
error = norm(ytest-fxtest)/norm(ytest);
fprintf('test error = %d\n',error)

%% Display random evaluation of Brownian motion
plot_geometric_brownian_kl(ytest(1,:)',fxtest(1,:)');
mysaveas(pathname,'geometric_brownian_kl.fig','fig');
mymatlab2tikz(pathname,'geometric_brownian_kl.tex');

%% Display multi-index set
dim = [1 3 5];
plot(f.basis.indices,'dim',dim,'legend',false)
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);

dim = [1 7 10];
plot(f.basis.indices,'dim',dim,'legend',false)
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);
