%% Adaptive sparse polynomial approximation - Ishigami function %%
%%--------------------------------------------------------------%%
% [Ishigami & Homma, 1990], [Saltelli, Chan & Scott, 2000], [Sudret 2008], [Blatman & Sudret 2011]

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');

%% Filename and Pathname
filename = 'sparse_ishigami_function';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Scalar-valued ishigami function
% y = sin(x_1) + a*sin(x_2)^2 + b*(x_3)^4*sin(x_1)
d = 3; % parametric dimension
a = 7;
b = 0.1;

fun = vectorize('sin(x1)+7*sin(x2)^2+0.1*x3^4*sin(x1)');
v = UniformRandomVariable(-pi,pi);
V = RVUNIFORM(-pi,pi);

% fun = vectorize('sin(pi*x1) + 7*sin(pi*x2)^2 + 0.1*(pi*x3)^4*sin(pi*x1)');
% v = UniformRandomVariable(-1,1);
% V = RVUNIFORM(-1,1);

% fun = vectorize('sin(-pi+2*pi*x1) + 7*sin(-pi+2*pi*x2)^2 + 0.1*(-pi+2*pi*x3)^4*sin(-pi+2*pi*x1)');
% v = UniformRandomVariable(0,1);
% V = RVUNIFORM(0,1);

rv = RandomVector(v,d);
RV = RANDVARS(repmat({V},1,d));

% [fun,rv] = multivariateFunctionsBenchmarks('ishigami',d,a,b);

fun = MultiVariateFunction(fun,d);
fun.evaluationAtMultiplePoints = true;

[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',2);

%% Adaptive sparse approximation using least-squares
p = 50;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,d);

s = AdaptiveSparseTensorAlgorithm();
% s.nbSamples = 1;
% s.addSamplesFactor = 0.1;
s.tol = 1e-4;
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
[f,err,y] = s.leastSquares(fun,bases,ls,rv);
time = toc(t);

ind = f.basis.indices.array;
switch gettypebase(PC)
    case 1
        ind(:,ndims(f.basis)+1) = sum(ind(:,1:ndims(f.basis)),2);
    case 2
        ind(:,ndims(f.basis)+1) = max(ind(:,1:ndims(f.basis)),[],2);
end
PC = setindices(PC,ind,'update');
u = f.data';
u = PCMATRIX(u,[size(u,1) 1],PC);

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
dim = 1:3;
plotMultiIndexSet(f,'dim',dim,'legend',false)
mysaveas(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end))],'fig');
mymatlab2tikz(pathname,['multi_index_set_dim' sprintf('_%d',dim(1:end)) '.tex']);

%% Quantities of interest : mean, variance, Sobol indices
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
num.mean = mean(u);
num.var = variance(u);
num.S1 = sobol_indices(u,1);
num.S2 = sobol_indices(u,2);
num.S3 = sobol_indices(u,3);
num.S12 = sobol_indices_group(u,[1,2]) - num.S1 - num.S2;
num.S13 = sobol_indices_group(u,[1,3]) - num.S1 - num.S3;
num.S23 = sobol_indices_group(u,[2,3]) - num.S2 - num.S3;
num.S123 = sobol_indices_group(u,[1,2,3]) - num.S1 - num.S2 - num.S3 - num.S12 - num.S13 - num.S23;
num.S1T = num.S1 + num.S12 + num.S13 + num.S123;
num.S2T = num.S2 + num.S12 + num.S23 + num.S123;
num.S3T = num.S3 + num.S13 + num.S23 + num.S123;
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
