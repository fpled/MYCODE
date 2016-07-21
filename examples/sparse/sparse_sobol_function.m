%% Adaptive sparse polynomial approximation - Sobol function %%
%%-----------------------------------------------------------%%
% [Sobol 2003], [Sudret 2008]

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');

%% Filename and Pathname
filename = 'sparse_sobol_function';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Scalar-valued Sobol function
% y = prod_{j=1}^d (|4*x_j-2|+a_j)/(1+a_j)
d = 8; % parametric dimension
a = [1,2,5,10,20,50,100,500];

% fun = @(x) ((abs(4*x(:,1)-2)+a(1))/(1+a(1)).*(abs(4*x(:,2)-2)+a(2))/(1+a(2)).*...
%     (abs(4*x(:,3)-2)+a(3))/(1+a(3)).*(abs(4*x(:,4)-2)+a(4))/(1+a(4)).*...
%     (abs(4*x(:,5)-2)+a(5))/(1+a(5)).*(abs(4*x(:,6)-2)+a(6))/(1+a(6)).*...
%     (abs(4*x(:,7)-2)+a(7))/(1+a(7)).*(abs(4*x(:,8)-2)+a(8))/(1+a(8)));
fun = @(x) (prod((abs(4*x-2)+repmat(a,size(x,1),1))./repmat(1+a,size(x,1),1),2));

v = UniformRandomVariable(0,1);
rv = RandomVector(v,d);

% [fun,rv] = multivariateFunctionsBenchmarks('sobol',d,a);

fun = MultiVariateFunction(fun,d);
fun.evaluationAtMultiplePoints = true;

V = RVUNIFORM(0,1);
RV = RANDVARS(repmat({V},1,d));
[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',2);

%% Adaptive sparse approximation using least-squares
p = 50;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,d);

s = AdaptiveSparseTensorAlgorithm();
% s.nbSamples = 1;
% s.addSamplesFactor = 0.1;
s.tol = 5e-2;
s.tolStagnation = 1e-2;
s.tolOverfit = 1.05;
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
[f,err,N] = s.leastSquares(fun,bases,ls,rv);
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
fprintf('parametric dimension = %d\n',ndims(f.basis))% fprintf('parametric dimension = %d\n',numel(rv))
fprintf('basis dimension = %d\n',numel(f.basis))
fprintf('order = [ %s ]\n',num2str(max(f.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(f.basis.indices.array)
fprintf('nb samples = %d\n',N)
fprintf('optimal nb samples = %d (according to empirical rule without regularization)\n',(ndims(f.basis)-1)*numel(f.basis))
fprintf('CV error = %d\n',err)
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
anal.mean = 1;
anal.var = 1/(1+a(1))^2*(4/3+2*a(1)+a(1)^2)*1/(1+a(2))^2*(4/3+2*a(2)+a(2)^2)*1/(1+a(3))^2*(4/3+2*a(3)+a(3)^2)...
    *1/(1+a(4))^2*(4/3+2*a(4)+a(4)^2)*1/(1+a(5))^2*(4/3+2*a(5)+a(5)^2)*1/(1+a(6))^2*(4/3+2*a(6)+a(6)^2)...
    *1/(1+a(7))^2*(4/3+2*a(7)+a(7)^2)*1/(1+a(8))^2*(4/3+2*a(8)+a(8)^2) - 1;
anal.S1 = 1/anal.var*1/(3*(1+a(1))^2);
anal.S2 = 1/anal.var*1/(3*(1+a(2))^2);
anal.S3 = 1/anal.var*1/(3*(1+a(3))^2);
anal.S4 = 1/anal.var*1/(3*(1+a(4))^2);
anal.S5 = 1/anal.var*1/(3*(1+a(5))^2);
anal.S6 = 1/anal.var*1/(3*(1+a(6))^2);
anal.S7 = 1/anal.var*1/(3*(1+a(7))^2);
anal.S8 = 1/anal.var*1/(3*(1+a(8))^2);
anal.S12 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(2))^2);
anal.S13 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(3))^2);
anal.S14 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(4))^2);
anal.S15 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(5))^2);
anal.S16 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(6))^2);
anal.S17 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(7))^2);
anal.S18 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(8))^2);
% Numerical approximate values
num.mean = mean(u);
num.var = variance(u);
num.S1 = sobol_indices(u,1);
num.S2 = sobol_indices(u,2);
num.S3 = sobol_indices(u,3);
num.S4 = sobol_indices(u,4);
num.S5 = sobol_indices(u,5);
num.S6 = sobol_indices(u,6);
num.S7 = sobol_indices(u,7);
num.S8 = sobol_indices(u,8);
num.S12 = sobol_indices_group(u,[1,2]) - num.S1 - num.S2;
num.S13 = sobol_indices_group(u,[1,3]) - num.S1 - num.S3;
num.S14 = sobol_indices_group(u,[1,4]) - num.S1 - num.S4;
num.S15 = sobol_indices_group(u,[1,5]) - num.S1 - num.S5;
num.S16 = sobol_indices_group(u,[1,6]) - num.S1 - num.S6;
num.S17 = sobol_indices_group(u,[1,7]) - num.S1 - num.S7;
num.S18 = sobol_indices_group(u,[1,8]) - num.S1 - num.S8;
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
fprintf(['| Sobol index S_4   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S4,num.S4,abs((anal.S4-num.S4)/anal.S4))
fprintf(['| Sobol index S_5   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S5,num.S5,abs((anal.S5-num.S5)/anal.S5))
fprintf(['| Sobol index S_6   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S6,num.S6,abs((anal.S6-num.S6)/anal.S6))
fprintf(['| Sobol index S_7   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S7,num.S7,abs((anal.S7-num.S7)/anal.S7))
fprintf(['| Sobol index S_8   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S8,num.S8,abs((anal.S8-num.S8)/anal.S8))
fprintf(['| Sobol index S_12  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S12,num.S12,abs((anal.S12-num.S12)/anal.S12))
fprintf(['| Sobol index S_13  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S13,num.S13,abs((anal.S13-num.S13)/anal.S13))
fprintf(['| Sobol index S_14  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S14,num.S14,abs((anal.S14-num.S14)/anal.S14))
fprintf(['| Sobol index S_15  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S15,num.S15,abs((anal.S15-num.S15)/anal.S15))
fprintf(['| Sobol index S_16  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S16,num.S16,abs((anal.S16-num.S16)/anal.S16))
fprintf(['| Sobol index S_17  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S17,num.S17,abs((anal.S17-num.S17)/anal.S17))
fprintf(['| Sobol index S_18  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S18,num.S18,abs((anal.S18-num.S18)/anal.S18))
disp('+-------------------+------------+-----------+----------------+')

% Only the three or four first variables have a significant influence on the solution variance
