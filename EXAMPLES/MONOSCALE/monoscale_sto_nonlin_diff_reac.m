%% Monoscale stochastic nonlinear diffusion reaction %%
%%---------------------------------------------------%%

% clc
% clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');
myparallel('start');

%% Input data

filename = 'monoscale_sto_nonlin_diff_reac';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Domains and meshes

D = DOMAIN(2,[0.0,0.0],[1.0,1.0]);

nbelem = [20,20];
problem.S = build_model(D,'nbelem',nbelem);
% cl = 0.05;
% problem.S = build_model(D,'cl',cl,'filename',[pathname 'gmsh_domain']);

%% Random variables

d = 2; % parametric dimension
v = UniformRandomVariable(0,1);
rv = RandomVector(v,d);

V = RVUNIFORM(0,1);
RV = RANDVARS(repmat({V},1,d));
[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',2);

%% Materials

% Linear diffusion coefficient
% K(xi) = 1 + xi
p = 1;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,d);
rvb = getRandomVector(bases);
H = FullTensorProductFunctionalBasis(bases);
I = gaussIntegrationRule(v,2);
I = I.tensorize(d);

fun = @(x) 1 + x(:,1);
funtr = @(x) fun(transfer(rvb,rv,x));
fun = MultiVariateFunction(funtr,d);
fun.evaluationAtMultiplePoints = true;

K = H.projection(fun,I);
K = PCMATRIX(K.tensor.data,[1 1],PC);
% K = ones(1,1,PC) + X{1};

% Nonlinear reaction parameter
% R(xi) = xi
fun = @(x) x(:,2);
funtr = @(x) fun(transfer(rvb,rv,x));
fun = MultiVariateFunction(funtr,d);
fun.evaluationAtMultiplePoints = true;

R = H.projection(fun,I);
R = PCMATRIX(R.tensor.data,[1 1],PC);
% R = X{2};

mat = FOUR_ISOT('k',K,'r',R); % uniform value
problem.S = setmaterial(problem.S,mat);

%% Dirichlet boundary conditions

problem.S = final(problem.S);
problem.S = addcl(problem.S,[]);

%% Stiffness matrices and sollicitation vectors

if israndom(problem.S)
    problem.A = [];
    problem.Atang = [];
else
    problem.A = @(u) calc_fint(problem.S,u);
    problem.Atang = @(u) calc_rigitang(problem.S,u);
end

% Source term
f = 100;

if israndom(f)
    problem.b = [];
else
    problem.b = bodyload(problem.S,[],'QN',f);
end

%% Adaptive sparse approximation using least-squares

problem.solver = NEWTONSOLVER('type','tangent','increment',true,...
    'maxiter',100,'tol',1e-12,'display',false,'stopini',true);

p = 50;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,d);
rv = getRandomVector(bases);

s = AdaptiveSparseTensorAlgorithm();
% s.nbSamples = 1;
% s.addSamplesFactor = 0.1;
s.tol = 1e-6;
s.tolStagnation = 1e-1;
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

% u0 = solveSystem(calcOperator(funEval(problem,mean(RANDVARS(PC)))));
% fun = @(xi) solveSystem(calcOperator(funEval(problem,xi)),'inittype',u0);
fun = @(xi) solveSystem(calcOperator(funEval(problem,xi)));
fun = MultiVariateFunction(fun,d,getnbddlfree(problem.S));
fun.evaluationAtMultiplePoints = false;

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
% fprintf('multi-index set = \n')
% disp(f.basis.indices.array)
fprintf('order = [ %s ]\n',num2str(max(f.basis.indices.array)))
fprintf('nb samples = %d\n',N)
fprintf('CV error = %d\n',norm(err))
fprintf('elapsed time = %f s\n',time)

Ntest = 100;
[errtest,xtest,fxtest,ytest] = f.computeError(fun,Ntest);
fprintf('test error = %d\n',norm(errtest))

%% Save variables

save(fullfile(pathname,'all.mat'));

%% Display domains and meshes

plotDomain(D);
mysaveas(pathname,'domain',formats,renderer);
mymatlab2tikz(pathname,'domain.tex');

% plotPartition(problem.S,'legend',false);
% mysaveas(pathname,'mesh_partition',formats,renderer);

plotModel(problem.S,'legend',false);
mysaveas(pathname,'mesh',formats,renderer);

%% Display multi-index set

plotMultiIndexSet(f,'legend',false)
mysaveas(pathname,'multi_index_set','fig');
mymatlab2tikz(pathname,'multi_index_set.tex');

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
