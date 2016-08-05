%% Monoscale stochastic linear diffusion circular inclusions isotropic %%
%%---------------------------------------------------------------------%%
% [Beck, Nobile, Tamellini, Tempone 2011,2014], [Chkifa, Cohen, Migliorati, Tempone 2014]

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');
myparallel('start');

%% Input data

filename = 'monoscale_sto_lin_diff_circ_inclusions_iso';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Domains and meshes

D = DOMAIN(2,[0.0,0.0],[1.0,1.0]);

r = 0.13;
B = cell(1,9);
B{1} = CIRCLE(0.2,0.2,r);
B{2} = CIRCLE(0.2,0.5,r);
B{3} = CIRCLE(0.2,0.8,r);
B{4} = CIRCLE(0.5,0.8,r);
B{5} = CIRCLE(0.8,0.8,r);
B{6} = CIRCLE(0.8,0.5,r);
B{7} = CIRCLE(0.8,0.2,r);
B{8} = CIRCLE(0.5,0.2,r);
B{9} = DOMAIN(2,[0.4,0.4],[0.6,0.6]);

cl = 0.02;
problem.S = gmshdomainwithinclusion(D,B,cl,cl,[pathname 'gmsh_circular_inclusions']);

%% Random variables

d = 8; % parametric dimension d = 2, 4, 8
v = UniformRandomVariable(-0.99,-0.2);
rv = RandomVector(v,d);

V = RVUNIFORM(-0.99,-0.2);
RV = RANDVARS(repmat({V},1,d));
[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',1);

%% Materials

% Deterministic subdomains
% Linear diffusion coefficient
K_det = 1;
mat_det = FOUR_ISOT('k',K_det); % uniform value
mat_det = setnumber(mat_det,0);
k = [0 9];
problem.S = setmaterial(problem.S,mat_det,k+1);

% Stochastic subdomains
p = 1;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,d);
vb = basis.basis.randomVariable;
rvb = getRandomVector(bases);
H = FullTensorProductFunctionalBasis(bases);
I = gaussIntegrationRule(vb,2);
I = I.tensorize(d);

mat_sto = MATERIALS();
for i=1:d
    % Linear diffusion coefficient
    % K_sto(xi) = 1 + xi
    fun = @(x) 1 + x(:,i);
    funtr = @(x) fun(transfer(rvb,rv,x));
    fun = MultiVariateFunction(funtr,d);
    fun.evaluationAtMultiplePoints = true;
    
    K_sto = H.projection(fun,I);
    
    mat_sto{i} = FOUR_ISOT('k',K_sto); % uniform value
    mat_sto{i} = setnumber(mat_sto{i},i);
end
switch d
    case 2
        k = [2 4 6 8];
        problem.S = setmaterial(problem.S,mat_sto{1},k+1);
        k = [1 3 5 7];
        problem.S = setmaterial(problem.S,mat_sto{2},k+1);
    case 4
        k = [1 5];
        problem.S = setmaterial(problem.S,mat_sto{1},k+1);
        k = [2 6];
        problem.S = setmaterial(problem.S,mat_sto{2},k+1);
        k = [3 7];
        problem.S = setmaterial(problem.S,mat_sto{3},k+1);
        k = [4 8];
        problem.S = setmaterial(problem.S,mat_sto{4},k+1);
    case 8
        for i=1:d
            problem.S = setmaterial(problem.S,mat_sto{i},i+1);
        end
    otherwise
        error('Wrong number of variables')
end

%% Dirichlet boundary conditions

problem.S = final(problem.S);
problem.S = addcl(problem.S,[]);

%% Stiffness matrices and sollicitation vectors

if israndom(problem.S)
    problem.A = [];
else
    problem.A = calc_rigi(problem.S);
end

% Source term
f = 100;

if israndom(f)
    problem.b = [];
else
    problem.b = bodyload(keepgroupelem(problem.S,10),[],'QN',f);
end

%% Adaptive sparse approximation using least-squares

p = 50;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,d);
rv = getRandomVector(bases);

s = AdaptiveSparseTensorAlgorithm();
% s.nbSamples = 1;
% s.addSamplesFactor = 0.1;
s.tol = 1e-2;
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

fun = @(xi) solveSystem(calcOperator(funEval(problem,xi)));
fun = MultiVariateFunction(fun,d,getnbddlfree(problem.S));
fun.evaluationAtMultiplePoints = false;

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
% fprintf('multi-index set = \n')
% disp(f.basis.indices.array)
fprintf('order = [ %s ]\n',num2str(max(f.basis.indices.array)))
fprintf('nb samples = %d\n',size(y,1))
fprintf('CV error = %d\n',norm(err))
fprintf('elapsed time = %f s\n',time)

Ntest = 100;
[errtest,xtest,fxtest,ytest] = computeTestError(f,fun,Ntest);
fprintf('test error = %d\n',errtest)

%% Save variables

save(fullfile(pathname,'all.mat'));

%% Display domains and meshes

plotDomain(D,B);
mysaveas(pathname,'domain',formats,renderer);
mymatlab2tikz(pathname,'domain.tex');

plotPartition(problem.S,'legend',false);
mysaveas(pathname,'mesh_partition',formats,renderer);

plotModel(problem.S,'legend',false);
mysaveas(pathname,'mesh',formats,renderer);

%% Display multi-index set

for i=1:2:d
    plotMultiIndexSet(f,'dim',[i i+1],'legend',false)
    mysaveas(pathname,['multi_index_set_dim_' num2str(i) '_' num2str(i+1)],'fig');
    mymatlab2tikz(pathname,['multi_index_set_dim_' num2str(i) '_' num2str(i+1) '.tex']);
end

%% Display statistical outputs of solution

% plotStats(problem.S,u);

plotMean(problem.S,u);
mysaveas(pathname,'mean_solution',formats,renderer);

plotVar(problem.S,u);
mysaveas(pathname,'var_solution',formats,renderer);

plotStd(problem.S,u);
mysaveas(pathname,'std_solution',formats,renderer);

for i=1:d
    plotSobolIndices(problem.S,u,i);
    mysaveas(pathname,['sobol_indices_solution_var_' num2str(i)],formats,renderer);
    
    plotSensitivityIndicesMaxVar(problem.S,u,i);
    mysaveas(pathname,['sensitivity_indices_solution_var_' num2str(i)],formats,renderer);
end

%% Display random evaluations of solution

% nbsamples = 3;
% for i=1:nbsamples
%     Stest = randomeval(problem.S,xtest(i,:)');
%     plotSolution(Stest,ytest(i,:)');
%     plotSolution(Stest,fxtest(i,:)');
% end

myparallel('stop');
