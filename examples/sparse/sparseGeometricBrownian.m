%% Adaptive sparse polynomial approximation - Geometric brownian %%
%%---------------------------------------------------------------%%

% clc
clearvars
close all
rng('default');

%% Input data
filename = 'geometricBrownian';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','sparse',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Vector-valued geometric brownian
d = 10; % parametric dimension
c = -1;
s = 0.5;
T = 1;
m = 101; % output size
n = m-1;

% fun = @(x) geometricBrownian(x,c,s,T,n);
% v = NormalRandomVariable(0,1);
% rv = RandomVector(v,d);

[fun,rv] = multivariateFunctionsBenchmarks('geometricbrownian',d,c,s,T,n);

fun = UserDefinedFunction(fun,d,m);
fun.evaluationAtMultiplePoints = true;

%% Adaptive sparse least-squares approximation
p = 50;
bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
bases = FunctionalBases(bases);

s = AdaptiveSparseTensorAlgorithm();
s.tol = 1e-3;
s.tolStagnation = 5e-2;
s.display = true;
s.displayIterations = true;

ls = LinearModelLearningSquareLoss();
ls.errorEstimation = true;
ls.sharedCoefficients = false;

t = tic;
[f,err,~,y] = s.leastSquares(fun,bases,ls,rv);
time = toc(t);

%% Outputs
fprintf('\n')
fprintf('parametric dimension = %d\n',ndims(f.basis))
fprintf('basis dimension = %d\n',cardinal(f.basis))
fprintf('order = [ %s ]\n',num2str(max(f.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(f.basis.indices.array)
fprintf('nb samples = %d\n',size(y,1))
fprintf('CV error = %e\n',norm(err))
fprintf('elapsed time = %f s\n',time)

%% Test
N = 1000;
errL2 = testError(f,fun,N,rv);
fprintf('mean squared error = %d\n',errL2)

%% Display random evaluations
x = random(rv,1);
plotGeometricBrownian(fun(x),f(x));
mysaveas(pathname,'geometric_brownian','fig');
mymatlab2tikz(pathname,'geometric_brownian.tex');

%% Display multi-index set
dim = [1 3 5];
plotMultiIndexSet(f,'dim',dim,'legend',false);
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);

dim = [1 7 10];
plotMultiIndexSet(f,'dim',dim,'legend',false);
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);
