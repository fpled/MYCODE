%% Adaptive sparse polynomial approximation - Anisotropic function %%
%%-----------------------------------------------------------------%%
% [Chkifa, Cohen, Schwab, 2014]

% clc
clearvars
close all
% rng('default');

%% Input data
filename = 'anisotropicFunction';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','sparse',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Scalar-valued anisotropic function
d = 16; % parametric dimension

% fun(x) = x_3*sin(x_4+x_16)
% fun = @(x) x(:,3).*sin(x(:,4)+x(:,16));
% fun = vectorize('x3*sin(x4+x16)');

% fun(x) = 1/(1 + sum_{j=1}^d (g_j*x_j)) with g_j = 3/(5*j^3)
% fun = @(x) 1./(ones(size(x,1),1)+x*(3./(5.*(1:size(x,2)).^3))');

% fun(x) = 1/(1 + (sum_{j=1}^d (g_j*x_j))^2) with g_j = 5/(j^3)
% fun = @(x) 1./(ones(size(x,1),1)+(x*(5./((1:size(x,2)).^3))').^2);

% fun(x) = 1/(1 + sum_{j=1}^d (g_j*x_j)) with g_j = 10^(-j)
% fun = @(x) 1./(ones(size(x,1),1)+x*(10.^(-(1:size(x,2))))');

% v = UniformRandomVariable(0,1);
% rv = RandomVector(v,d);

[fun,rv] = multivariateFunctionsBenchmarks('anisotropic',d);

fun = UserDefinedFunction(fun,d);
fun.evaluationAtMultiplePoints = true;

%% Adaptive sparse least-squares approximation
p = 50;
bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
bases = FunctionalBases(bases);

s = AdaptiveSparseTensorAlgorithm();
s.tol = 1e-6;
s.tolStagnation = 5e-2;
s.display = true;
s.displayIterations = true;

ls = LeastSquaresSolver();
ls.errorEstimation = true;

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
fprintf('CV error = %d\n',err)
fprintf('elapsed time = %f s\n',time)

%% Test
N = 1000;
errL2 = testError(f,fun,N,rv);
fprintf('mean squared error = %d\n',errL2)

%% Display multi-index set
dim = [3 4 16];
plotMultiIndexSet(f,'dim',dim,'legend',false);
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);

dim = [1 2 4];
plotMultiIndexSet(f,'dim',dim,'legend',false);
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);
