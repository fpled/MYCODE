%% Multiscale stochastic nonlinear diffusion reaction aligned inclusions isotropic %%
%%---------------------------------------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');
myparallel('start');

%% Input data

n = 8; % number of patches
% for rho = [0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8]
% for tol = 1:5

filename = ['multiscale_sto_nonlin_diff_reac_' num2str(n) '_align_inclusions_iso'];
% filename = ['multiscale_sto_nonlin_diff_reac_' num2str(n) '_align_inclusions_iso_tol_3_rho_' num2str(rho)];
% filename = ['multiscale_sto_nonlin_diff_reac_' num2str(n) '_align_inclusions_iso_tol_'  num2str(tol) '_rho_aitken'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
% pathname = fullfile('/Users/Op/Documents/Recherche/GeM/Results',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc2'};
renderer = 'OpenGL';

directSolver = true;
iterativeSolver = true;

%% Domains and meshes

% Global
glob = Global();
glob_out = GlobalOutside();

D = DOMAIN(2,[0.0,0.0],[2.0,2*n]);

nbelem = [20,20*n];
glob.S = build_model(D,'nbelem',nbelem);
% cl = 0.25;
% glob.S = build_model(D,'cl',cl,'filename',[pathname 'gmsh_domain']);

% Patches
patches = Patches(n);

D_patch = cell(1,n);
for k=1:n
    D_patch{k} = DOMAIN(2,[0.5,2*k-1.5],[1.5,2*k-0.5]);
end

nbelem_patch = [20,20];
for k=1:n
    patches.patches{k}.S = build_model(D_patch{k},'nbelem',nbelem_patch);
end
% cl_patch = 0.025;
% for k=1:n
%     patches.patches{k}.S = build_model(D_patch{k},'cl',cl_patch,'filename',[pathname 'gmsh_patch_' num2str(k)]);
% end

% Partition of global mesh
glob = partition(glob,patches);

%% Random variables

d = 2*n; % parametric dimension
v = UniformRandomVariable(0,1);
rv = RandomVector(v,d);

%% Materials

% Linear diffusion coefficient
K_out = 1;
K_patch = cell(1,n);
K_in = cell(1,n);
% Nonlinear reaction parameter
R_patch = cell(1,n);

p = 1;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,[],d);
vb = basis.basis.randomVariable;
rvb = getRandomVector(bases);
H = FullTensorProductFunctionalBasis(bases);
I = gaussIntegrationRule(vb,2);
I = I.tensorize(d);

for k=1:n
    patch = patches.patches{k};
    % K_patch(x,xi) = 1 + f(x) * xi
    % K_in(x)       = 1
    % R_patch(x,xi) = f(x) * xi
    % with f(x) = 1 if ||x-c||_Inf < L
    %           = 0 if ||x-c||_Inf >= L
    L = norm(getsize(D_patch{k}),Inf)/4;
    c = getcenter(D_patch{k});
    f = @(x) distance(x,c,Inf)<L;
    
    fun = @(xi) ones(size(xi,1),patch.S.nbnode) + xi(:,2*k-1) * double(squeeze(f(patch.S.node)))';
    funtr = @(xi) fun(transfer(rvb,rv,xi));
    fun = MultiVariateFunction(funtr,d,patch.S.nbnode);
    fun.evaluationAtMultiplePoints = true;
    
    K_patch{k} = H.projection(fun,I);
    
    fun = @(xi) xi(:,2*k) * double(squeeze(f(patch.S.node)))';
    funtr = @(xi) fun(transfer(rvb,rv,xi));
    fun = MultiVariateFunction(funtr,d,patch.S.nbnode);
    fun.evaluationAtMultiplePoints = true;
    
    R_patch{k} = H.projection(fun,I);
    
    K_in{k} = 1;
end

% Complementary subdomain
mat_out = FOUR_ISOT('k',K_out); % uniform value
mat_out = setnumber(mat_out,0);
glob.S = setmaterial(glob.S,mat_out,getnumgroupelemwithparam(glob.S,'partition',0));

% Patches
mat_patch = MATERIALS();
for k=1:n
    % mat_patch{k} = FOUR_ISOT('k',K_patch{k},'r',R_patch{k}); % uniform value
    mat_patch{k} = FOUR_ISOT('k',FENODEFIELD(K_patch{k}),'r',FENODEFIELD(R_patch{k})); % nodal values
    mat_patch{k} = setnumber(mat_patch{k},k);
    patches.patches{k}.S = setmaterial(patches.patches{k}.S,mat_patch{k});
end

% Fictitious patches
mat_in = MATERIALS();
for k=1:n
    mat_in{k} = FOUR_ISOT('k',K_in{k}); % uniform value
    % mat_in{k} = FOUR_ISOT('k',FENODEFIELD(K_in{k})); % nodal values
    mat_in{k} = setnumber(mat_in{k},k);
    glob.S = setmaterial(glob.S,mat_in{k},getnumgroupelemwithparam(glob.S,'partition',k));
end

%% Dirichlet boundary conditions

% Global
glob.S = final(glob.S);
glob.S = addcl(glob.S,[]);
glob.S_out = getfinalmodelpart(glob.S,0);
% S_in = cell(1,n);
% for k=1:n
%     S_in{k} = getfinalmodelpart(glob.S,k);
% end

% Complementary subdomain
glob_out.S_out = glob.S_out;

% Patches
for k=1:n
    patches.patches{k}.S = final(patches.patches{k}.S);
end

% Interfaces
interfaces = Interfaces(patches);

%% Stiffness matrices and sollicitation vectors

% Source term
f = 1;

% Complementary subdomain
glob_out.A_out = calc_rigi(glob_out.S_out);
glob_out.b_out = bodyload(glob_out.S_out,[],'QN',f);

% Global
glob.A = calc_rigi(glob.S);
for k=1:n
    glob.A_in{k} = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
end
glob.b_out = bodyload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),[],'QN',f);

% Patches
for k=1:n
    if ~israndom(patches.patches{k}.S)
        patches.patches{k}.A = @(u) calc_fint(patches.patches{k}.S,u);
        patches.patches{k}.Atang = @(u) calc_rigitang(patches.patches{k}.S,u);
    end
    if ~israndom(f)
        patches.patches{k}.b = bodyload(patches.patches{k}.S,[],'QN',f);
    end
end

%% Mass matrices

for k=1:n
    interfaces.interfaces{k}.M = calc_massgeom(interfaces.interfaces{k}.S);
end

%% Projection operators

glob.P_out = calcProjection(glob);
for k=1:n
    [interfaces.interfaces{k}.P_glob] = calcProjection(interfaces.interfaces{k},glob);
    [interfaces.interfaces{k}.P_glob_out,numnode] = calcProjection(interfaces.interfaces{k},glob_out);
    interfaces.interfaces{k}.P_patch = calcProjection(patches.patches{k},interfaces.interfaces{k});
    % plotProjectionOperator(glob,patches.patches{k},numnode);
end

%% Parameters for global and local problems

% Global problem
glob.increment = true;

% Local problems
for k=1:n
    patches.patches{k}.changeOfVariable = false;
    patches.patches{k}.increment = true;
    patches.patches{k}.initializationType = 'zero';
    patches.patches{k}.solver = NEWTONSOLVER('type','tangent','increment',patches.patches{k}.increment,...
        'maxiter',100,'tol',1e-12,'display',false,'stopini',true);
end

%% Direct solver

p = 50;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,[],d);
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

DS = DirectSolver();
DS.changeOfVariable = false;
DS.display = true;
DS.solver = NEWTONSOLVER('type','tangent','increment',true,...
    'maxiter',100,'tol',1e-12,'display',false,'stopini',true);
DS.initializationType = 'zero';
if directSolver
    [U_ref,w_ref,lambda_ref,output_ref] = DS.solveRandom(glob_out,patches,interfaces,s,bases,ls,rv);
    save(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','output_ref');
else
    load(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','output_ref');
end

%% Outputs

fprintf('\n')
fprintf('spatial dimension = %d for U_ref\n',U_ref.sz)
for k=1:n
    fprintf('                  = %d for w_ref{%u}\n',w_ref{k}.sz,k)
    fprintf('                  = %d for lambda_ref{%u}\n',lambda_ref{k}.sz,k)
end
fprintf('parametric dimension = %d\n',ndims(U_ref.basis))
% fprintf('parametric dimension = %d\n',numel(rv))
fprintf('basis dimension = %d for U_ref\n',numel(U_ref.basis))
for k=1:n
    fprintf('                = %d for w_ref{%u}\n',numel(w_ref{k}.basis),k)
    fprintf('                = %d for lambda_ref{%u}\n',numel(lambda_ref{k}.basis),k)
end
fprintf('order = [ %s ] for U_ref\n',num2str(max(U_ref.basis.indices.array)))
for k=1:n
    fprintf('      = [ %s ] for w_ref{%u}\n',num2str(max(w_ref{k}.basis.indices.array)),k)
    fprintf('      = [ %s ] for lambda_ref{%u}\n',num2str(max(lambda_ref{k}.basis.indices.array)),k)
end
% fprintf('multi-index set for U_ref = \n')
% disp(num2str(U_ref.basis.indices.array))
% for k=1:n
%     fprintf('multi-index set for w_ref{%u} = \n',k)
%     disp(num2str(w_ref{k}.basis.indices.array))
%     fprintf('multi-index set for lambda_ref{%u} = \n',k)
%     disp(num2str(lambda_ref{k}.basis.indices.array))
% end
fprintf('nb samples = %d\n',output_ref.nbSamples)
fprintf('CV error = %d for U_ref\n',norm(output_ref.CVErrorGlobalSolution))
for k=1:n
    fprintf('         = %d for w_ref{%u}\n',norm(output_ref.CVErrorLocalSolution{k}),k)
    fprintf('         = %d for lambda_ref{%u}\n',norm(output_ref.CVErrorLagrangeMultiplier{k}),k)
end
fprintf('elapsed time = %f s\n',output_ref.time)

%% Global-local Iterative solver

s.tol = 1e-3;
% s.tol = 10^(-tol);
s.tolStagnation = 1e-1;
s.display = true;
s.displayIterations = false;

IS = IterativeSolver();
IS.maxIterations = 20;
IS.tolerance = eps;
IS.relaxation = 'Aitken';
IS.updateRelaxationParameter = true;
% IS.relaxation = rho;
% IS.updateRelaxationParameter = false;
IS.errorCriterion = 'reference';
IS.referenceSolution = {U_ref,w_ref,lambda_ref};
IS.display = true;
IS.displayIterations = true;
if iterativeSolver
    [U,w,lambda,output] = IS.solveRandom(glob,patches,interfaces,s,bases,ls,rv);
    save(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
else
    load(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
end

%% Outputs

fprintf('\n')
fprintf('spatial dimension = %d for U\n',U.sz)
for k=1:n
    fprintf('                  = %d for w{%u}\n',w{k}.sz,k)
    fprintf('                  = %d for lambda{%u}\n',lambda{k}.sz,k)
end
fprintf('parametric dimension = %d\n',ndims(U.basis))
% fprintf('parametric dimension = %d\n',numel(rv))
fprintf('basis dimension = %d for U\n',numel(U.basis))
for k=1:n
    fprintf('                = %d for w{%u}\n',numel(w{k}.basis),k)
    fprintf('                = %d for lambda{%u}\n',numel(lambda{k}.basis),k)
end
fprintf('order = [ %s ] for U\n',num2str(max(U.basis.indices.array)))
for k=1:n
    fprintf('      = [ %s ] for w{%u}\n',num2str(max(w{k}.basis.indices.array)),k)
    fprintf('      = [ %s ] for lambda{%u}\n',num2str(max(lambda{k}.basis.indices.array)),k)
end
% fprintf('multi-index set for U = \n')
% disp(num2str(U.basis.indices.array))
% for k=1:n
%     fprintf('multi-index set for w{%u} = \n',k)
%     disp(num2str(w{k}.basis.indices.array))
%     fprintf('multi-index set for lambda{%u} = \n',k)
%     disp(num2str(lambda{k}.basis.indices.array))
% end
fprintf('elapsed time = %f s\n',output.totalTime)

% end

myparallel('stop');
