%% Adaptive sparse polynomial approximation - Polynomial function %%
%%----------------------------------------------------------------%%
% [Sudret 2008]

% clc
clear all
close all

%% Filename and Pathname
d = 3; % number of random variables
q = 2; % power of random variables
filename = ['sparse_approx_polynomial_function_partial_degree_' num2str(q) '_nbvar_' num2str(d)];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'RESULTS',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
% set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object

%% Polynomial function of total degree q*M = 2*3 = 6
% y = 1/(2^M) * prod_{j=1}^{M}(3*(x_j)^q+1)
fun = @(x) (1/(2^(size(x,2)))*prod(3*x.^q+1,2));

fun = MultiVariateFunction(fun,d);
fun.evaluationAtMultiplePoints = true;

rv = RandomVector(UniformRandomVariable(0,1),d);

h = PolynomialFunctionalBasis(LegendrePolynomials(),0:15);
H = FunctionalBases(h,d);

%% Adaptive sparse approximation using least-squares
s = AdaptiveSparseTensorAlgorithm();
% s.addSamplesFactor = 0.1;
s.tol = 1e-12;
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
dim = 1:3;
plot(f.basis.indices,'dim',dim,'legend',false)
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);

%% Quantities of interest : mean, variance, Sobol indices
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
    num.mean = mean(f.data);
    num.var = var(f.data);
    num.S1 = sobol_indices(f.data,1);
    num.S2 = sobol_indices(f.data,2);
    num.S3 = sobol_indices(f.data,3);
    num.S12 = sobol_indices_group(f.data,[1,2]) - num.S1 - num.S2;
    num.S13 = sobol_indices_group(f.data,[1,3]) - num.S1 - num.S3;
    num.S23 = sobol_indices_group(f.data,[2,3]) - num.S2 - num.S3;
    num.S123 = sobol_indices_group(u,[1,2,3]) - num.S1 - num.S2 - num.S3 - num.S12 - num.S13 - num.S23;
    num.S1T = num.S1 + num.S12 + num.S13 + num.S123;
    num.S2T = num.S2 + num.S12 + num.S23 + num.S123;
    num.S3T = num.S3 + num.S13 + num.S23 + num.S123;
    % Comparative table
    fanal = '%10.5f';
    fnum = '%9.5f';
    ferr = '%14.4e';
    disp('+-------------------+------------+-----------+----------------+')
    disp('| Quantity \ Value  | Analytical | Numerical | Relative error |')
    disp('+-------------------+------------+-----------+----------------+')
    fprintf(['| Mean  E           | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.mean,num.mean,abs((anal.mean-num.mean)/anal.mean))
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
