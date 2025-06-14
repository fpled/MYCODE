%% Adaptive sparse polynomial approximation - Ishigami function %%
%%--------------------------------------------------------------%%
% [Ishigami, Homma, 1990, ISUMA]
% [Saltelli, Chan, Scott, 2000, Wiley]
% [Sudret, 2008, RESS]
% [Blatman, Sudret, 2010, PEM]
% [Blatman, Sudret, 2010, RESS]
% [Blatman, Sudret, 2011, JCP]

% clc
clearvars
close all
rng('default');

%% Input data
filename = 'ishigamiFunction';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','sparse',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Scalar-valued ishigami function
% fun(x) = sin(x_1) + a*sin(x_2)^2 + b*(x_3)^4*sin(x_1)
d = 3; % parametric dimension
a = 7;
b = 0.1;

% fun = vectorize('sin(x1) + 7*sin(x2)^2 + 0.1*x3^4*sin(x1)');
% fun = @(x) sin(x(:,1)) + a.*sin(x(:,2)).^2 + b.*x(:,3).^4.*sin(x(:,1));
% v = UniformRandomVariable(-pi,pi);

% fun = vectorize('sin(pi*x1) + 7*sin(pi*x2)^2 + 0.1*(pi*x3)^4*sin(pi*x1)');
% fun = @(x) sin(pi*x(:,1)) + a.*sin(pi*x(:,2)).^2 + b.*(pi*x(:,3)).^4.*sin(pi*x(:,1));
% v = UniformRandomVariable(-1,1);

% fun = vectorize('sin(-pi+2*pi*x1) + 7*sin(-pi+2*pi*x2)^2 + 0.1*(-pi+2*pi*x3)^4*sin(-pi+2*pi*x1)');
% fun = @(x) sin(-pi+2*pi*x(:,1)) + a.*sin(-pi+2*pi*x(:,2)).^2 + b.*(-pi+2*pi*x(:,3)).^4.*sin(-pi+2*pi*x(:,1));
% v = UniformRandomVariable(0,1);

% rv = RandomVector(v,d);

[fun,rv] = multivariateFunctionsBenchmarks('ishigami',d,a,b);

fun = UserDefinedFunction(fun,d);
fun.evaluationAtMultiplePoints = true;

%% Adaptive sparse least-squares approximation
p = 50;
bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
bases = FunctionalBases(bases);

s = AdaptiveSparseTensorAlgorithm();
s.tol = 1e-4;
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
% Analytical exact values
anal.mean = a/2;
anal.var = 1/2 + a^2/8 + b*pi^4/5 + b^2*pi^8/18;
anal.S1 = (1/2*(1 + b*pi^4/5)^2)/anal.var;
anal.S2 = (a^2/8)/anal.var;
anal.S3 = 0;
anal.S12 = 0;
anal.S13 = (b^2*pi^8*(8/225))/anal.var;
anal.S23 = 0;
anal.S123 = 0;
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
% num.S1T = num.S1 + num.S12 + num.S13 + num.S123;
% num.S2T = num.S2 + num.S12 + num.S23 + num.S123;
% num.S3T = num.S3 + num.S13 + num.S23 + num.S123;
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
fprintf(['| Sobol index S_3   | ' fanal ' | ' fnum ' |                |\n'],anal.S3,num.S3)
fprintf(['| Sobol index S_12  | ' fanal ' | ' fnum ' |                |\n'],anal.S12,num.S12)
fprintf(['| Sobol index S_13  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S13,num.S13,abs((anal.S13-num.S13)/anal.S13))
fprintf(['| Sobol index S_23  | ' fanal ' | ' fnum ' |                |\n'],anal.S23,num.S23)
fprintf(['| Sobol index S_123 | ' fanal ' | ' fnum ' |                |\n'],anal.S123,num.S123)
fprintf(['| Sobol index S_1^T | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S1T,num.S1T,abs((anal.S1T-num.S1T)/anal.S1T))
fprintf(['| Sobol index S_2^T | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S2T,num.S2T,abs((anal.S2T-num.S2T)/anal.S2T))
fprintf(['| Sobol index S_3^T | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S3T,num.S3T,abs((anal.S3T-num.S3T)/anal.S3T))
disp('+-------------------+------------+-----------+----------------+')
