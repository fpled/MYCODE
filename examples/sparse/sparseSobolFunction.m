%% Adaptive sparse polynomial approximation - Sobol function %%
%%-----------------------------------------------------------%%
% [Saltelli, Sobol, RESS, 1995]
% [Sobol, 2003, RESS]
% [Sudret, 2008, RESS]

% clc
clearvars
close all
rng('default');

%% Input data
filename = 'sobolFunction';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','sparse',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Scalar-valued Sobol function
% fun(x) = prod_{j=1}^d (|4*x_j-2|+a_j)/(1+a_j)
d = 8; % parametric dimension
a = [1,2,5,10,20,50,100,500];

% fun = @(x) prod((abs(4*x-2)+repmat(a,size(x,1),1))./repmat(1+a,size(x,1),1),2);
% v = UniformRandomVariable(0,1);
% rv = RandomVector(v,d);

[fun,rv] = multivariateFunctionsBenchmarks('sobol',d,a);

fun = UserDefinedFunction(fun,d);
fun.evaluationAtMultiplePoints = true;

%% Adaptive sparse least-squares approximation
p = 50;
bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
bases = FunctionalBases(bases);

s = AdaptiveSparseTensorAlgorithm();
s.tol = 5e-2;
s.tolStagnation = 1e-2;
s.tolOverfit = 1.05;
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
fprintf('optimal nb samples = %d (according to empirical rule without regularization)\n',(ndims(f.basis)-1)*cardinal(f.basis))
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
num.mean = mean(f);
num.var = variance(f);
num.S1 = SensitivityAnalysis.sobolIndices(f,1,d);
num.S2 = SensitivityAnalysis.sobolIndices(f,2,d);
num.S3 = SensitivityAnalysis.sobolIndices(f,3,d);
num.S4 = SensitivityAnalysis.sobolIndices(f,4,d);
num.S5 = SensitivityAnalysis.sobolIndices(f,5,d);
num.S6 = SensitivityAnalysis.sobolIndices(f,6,d);
num.S7 = SensitivityAnalysis.sobolIndices(f,7,d);
num.S8 = SensitivityAnalysis.sobolIndices(f,8,d);
num.S12 = SensitivityAnalysis.sobolIndices(f,[1,2],d);
num.S13 = SensitivityAnalysis.sobolIndices(f,[1,3],d);
num.S14 = SensitivityAnalysis.sobolIndices(f,[1,4],d);
num.S15 = SensitivityAnalysis.sobolIndices(f,[1,5],d);
num.S16 = SensitivityAnalysis.sobolIndices(f,[1,6],d);
num.S17 = SensitivityAnalysis.sobolIndices(f,[1,7],d);
num.S18 = SensitivityAnalysis.sobolIndices(f,[1,8],d);
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
