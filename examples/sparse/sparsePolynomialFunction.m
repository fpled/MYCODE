%% Adaptive sparse polynomial approximation - Polynomial function %%
%%----------------------------------------------------------------%%
% [Sudret, 2008, RESS]

% clc
clearvars
close all
rng('default');

%% Input data
filename = 'polynomialFunction';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','sparse',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Scalar-valued polynomial function of total degree q*d = 2*3 = 6
% fun(x) = 1/(2^d) * prod_{j=1}^d (3*(x_j)^q+1)
d = 3; % parametric dimension
q = 2; % partial degree

% fun = @(x) 1/(2^(size(x,2)))*prod(3*x.^q+1,2);
% v = UniformRandomVariable(0,1);
% rv = RandomVector(v,d);

[fun,rv] = multivariateFunctionsBenchmarks('polynomial',d,q);

fun = UserDefinedFunction(fun,d);
fun.evaluationAtMultiplePoints = true;

%% Adaptive sparse least-squares approximation
p = 50;
% basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
% bases = FunctionalBases.duplicate(basis,d);
bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
bases = FunctionalBases(bases);

s = AdaptiveSparseTensorAlgorithm();
s.tol = 1e-12;
s.tolStagnation = 1e-1;
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
fprintf('CV error = %e\n',err)
fprintf('elapsed time = %f s\n',time)

%% Test
N = 1000;
errL2 = testError(f,fun,N,rv);
fprintf('mean squared error = %d\n',errL2)

%% Display multi-index set
dim = 1:3;
plotMultiIndexSet(f,'dim',dim,'legend',false);
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);

%% Quantities of interest: mean, variance, Sobol indices
if q==2
    % Analytical exact values
    anal.mean = 1;
    anal.var = (6/5)^d - 1;
    anal.S1 = 5^(-1)/anal.var;
    anal.S2 = anal.S1;
    anal.S3 = anal.S1;
    anal.S12 = 5^(-2)/anal.var;
    anal.S13 = anal.S12;
    anal.S23 = anal.S12;
    anal.S123 = 5^(-3)/anal.var;
    anal.S1T = anal.S1 + anal.S12 + anal.S13 + anal.S123;
    anal.S2T = anal.S2 + anal.S12 + anal.S23 + anal.S123;
    anal.S3T = anal.S3 + anal.S13 + anal.S23 + anal.S123;
    % Numerical approximate values
    num.mean = mean(f);
    num.var = variance(f);
    num.S1 = SensitivityAnalysis.sobolIndices(f,1,d);
    num.S2 = SensitivityAnalysis.sobolIndices(f,2,d);
    num.S3 = SensitivityAnalysis.sobolIndices(f,3,d);
    num.S12 = SensitivityAnalysis.sobolIndices(f,[1,2],d);
    num.S13 = SensitivityAnalysis.sobolIndices(f,[1,3],d);
    num.S23 = SensitivityAnalysis.sobolIndices(f,[2,3],d);
    num.S123 = SensitivityAnalysis.sobolIndices(f,[1,2,3],d);
    num.S1T = SensitivityAnalysis.totalSobolIndices(f,1,d);
    num.S2T = SensitivityAnalysis.totalSobolIndices(f,2,d);
    num.S3T = SensitivityAnalysis.totalSobolIndices(f,3,d);
%     num.S1T = num.S1 + num.S12 + num.S13 + num.S123;
%     num.S2T = num.S2 + num.S12 + num.S23 + num.S123;
%     num.S3T = num.S3 + num.S13 + num.S23 + num.S123;
    % Comparative table
    fanal = '%10.5f';
    fnum = '%9.5f';
    ferr = '%14.4e';
    fprintf('\n')
    disp('+-------------------+------------+-----------+----------------+')
    disp('| Quantity \ Value  | Analytical | Numerical | Relative error |')
    disp('+-------------------+------------+-----------+----------------+')
    fprintf(['| Mean E            | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.mean,num.mean,abs((anal.mean-num.mean)/anal.mean))
    fprintf(['| Variance V        | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.var,num.var,abs((anal.var-num.var)/anal.var))
    fprintf(['| Sobol index S_1   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S1,num.S1,abs((anal.S1-num.S1)/anal.S1))
    fprintf(['| Sobol index S_2   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S2,num.S2,abs((anal.S2-num.S2)/anal.S2))
    fprintf(['| Sobol index S_3   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S3,num.S3,abs((anal.S3-num.S3)/anal.S3))
    fprintf(['| Sobol index S_12  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S12,num.S12,abs((anal.S12-num.S12)/anal.S12))
    fprintf(['| Sobol index S_13  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S13,num.S13,abs((anal.S13-num.S13)/anal.S13))
    fprintf(['| Sobol index S_23  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S23,num.S23,abs((anal.S23-num.S23)/anal.S23))
    fprintf(['| Sobol index S_123 | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S123,num.S123,abs((anal.S123-num.S123)/anal.S123))
    fprintf(['| Sobol index S_1^T | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S1T,num.S1T,abs((anal.S1T-num.S1T)/anal.S1T))
    fprintf(['| Sobol index S_2^T | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S2T,num.S2T,abs((anal.S2T-num.S2T)/anal.S2T))
    fprintf(['| Sobol index S_3^T | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S3T,num.S3T,abs((anal.S3T-num.S3T)/anal.S3T))
    disp('+-------------------+------------+-----------+----------------+')
end
