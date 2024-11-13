%% Adaptive sparse polynomial approximation - Polynomial function %%
%%----------------------------------------------------------------%%

% clc
clearvars
close all
% rng('default');

%% Scalar-valued polynomial function
d = 2; % parametric dimension
% fun = @(x) ((x(:,1)+x(:,2)).^2+((x(:,1)-x(:,2))./10).^2);
fun = @(x) (x(:,1).^2+(x(:,2)./10).^2);

fun = UserDefinedFunction(fun,d);
fun.evaluationAtMultiplePoints = true;

v = UniformRandomVariable(0,1);
rv = RandomVector(v,d);

%% Adaptive sparse least-squares approximation
p = 50;
bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
bases = FunctionalBases(bases);

s = AdaptiveSparseTensorAlgorithm();
s.tol = 1e-12;
s.tolStagnation = 5e-2;
s.display = true;
s.displayIterations = true;

ls = LinearModelLearningSquareLoss();
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
fprintf('CV error = %d\n',norm(err))
fprintf('elapsed time = %f s\n',time)

%% Test
N = 1000;
errL2 = testError(f,fun,N,rv);
fprintf('mean squared error = %d\n',errL2)

%% Display multi-index set
dim = 1:2;
plotMultiIndexSet(f,'dim',dim,'legend',false);

%% Quantities of interest: mean, variance, Sobol indices
% Numerical approximate values
num.mean = mean(f);
num.var = variance(f);
num.S1 = SensitivityAnalysis.sobolIndices(f,1,d);
num.S2 = SensitivityAnalysis.sobolIndices(f,2,d);
num.S12 = SensitivityAnalysis.sobolIndices(f,[1,2],d);
num.S1T = SensitivityAnalysis.totalSobolIndices(f,1,d);
num.S2T = SensitivityAnalysis.totalSobolIndices(f,2,d);
% num.S1T = num.S1 + num.S12 + num.S12;
% num.S2T = num.S2 + num.S12 + num.S12;
% Comparative table
ff = '%9.5f';
disp('+-------------------+-----------+')
disp('| Quantity \ Value  | Numerical |')
disp('+-------------------+-----------+')
fprintf(['| Mean E            | ' ff ' |\n'],num.mean)
fprintf(['| Variance V        | ' ff ' |\n'],num.var)
fprintf(['| Sobol index S_1   | ' ff ' |\n'],num.S1)
fprintf(['| Sobol index S_2   | ' ff ' |\n'],num.S2)
fprintf(['| Sobol index S_12  | ' ff ' |\n'],num.S12)
fprintf(['| Sobol index S_1^T | ' ff ' |\n'],num.S1T)
fprintf(['| Sobol index S_2^T | ' ff ' |\n'],num.S2T)
disp('+-------------------+-----------+')
