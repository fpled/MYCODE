%% Multiscale stochastic nonlinear diffusion reaction aligned inclusions anisotropic %%
%%-----------------------------------------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');
myparallel('start');

%% Input data

n = 8; % number of patches
filename = ['multiscale_sto_nonlin_diff_reac_' num2str(n) '_align_inclusions_aniso'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
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

D = DOMAIN(2,[0.0,0.0],[2*n,2.0]);

nbelem = [10*n,10];
glob.S = build_model(D,'nbelem',nbelem);
% cl = 0.25;
% glob.S = build_model(D,'cl',cl,'filename',[pathname 'gmsh_domain']);

% Patches
patches = Patches(n);

D_patch = cell(1,n);
for k=1:n
    D_patch{k} = DOMAIN(2,[2*k-1.5,0.5],[2*k-0.5,1.5]);
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

V = RVUNIFORM(0,1);
RV = RANDVARS(repmat({V},1,d));
[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',1);

%% Materials

% Linear diffusion coefficient
K_out = 1;
K_patch = cell(1,n);
K_in = cell(1,n);
% Nonlinear reaction parameter
R_patch = cell(1,n);

p = 1;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,d);
vb = basis.basis.randomVariable;
rvb = getRandomVector(bases);
H = FullTensorProductFunctionalBasis(bases);
I = gaussIntegrationRule(vb,2);
I = I.tensorize(d);

g = 0.8:-0.1:0.1;
for k=1:n
    patch = patches.patches{k};
    % K_patch(x,xi) = 1 + f(x) * g * xi
    % K_in(x)       = 1
    % R_patch(x,xi) = f(x) * g * xi
    % with f(x) = 1 if ||x-c||_Inf < L
    %           = 0 if ||x-c||_Inf >= L
    L = norm(getsize(D_patch{k}),Inf)/4;
    c = getcenter(D_patch{k});
    f = @(x) distance(x,c,Inf)<L;
    
    fun = @(xi) ones(size(xi,1),patch.S.nbnode) + g(k) * xi(:,2*k-1) * double(squeeze(f(patch.S.node)))';
    funtr = @(xi) fun(transfer(rvb,rv,xi));
    fun = MultiVariateFunction(funtr,d,patch.S.nbnode);
    fun.evaluationAtMultiplePoints = true;
    
    K_patch{k} = H.projection(fun,I);
    
    fun = @(xi) g(k) * xi(:,2*k) * double(squeeze(f(patch.S.node)))';
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

DS = DirectSolver();
DS.changeOfVariable = false;
DS.display = true;
DS.solver = NEWTONSOLVER('type','tangent','increment',true,...
    'maxiter',100,'tol',1e-12,'display',false,'stopini',true);
DS.initializationType = 'zero';
if directSolver
    [fU_ref,fw_ref,flambda_ref,output_ref] = DS.solveRandom(glob_out,patches,interfaces,s,bases,ls,rv);
    save(fullfile(pathname,'reference_solution.mat'),'fU_ref','fw_ref','flambda_ref','output_ref');
else
    load(fullfile(pathname,'reference_solution.mat'),'fU_ref','fw_ref','flambda_ref','output_ref');
end

ind_U_ref = fU_ref.basis.indices.array;
ind_w_ref = cellfun(@(x) x.basis.indices.array,fw_ref,'UniformOutput',false);
ind_lambda_ref = cellfun(@(x) x.basis.indices.array,flambda_ref,'UniformOutput',false);
switch gettypebase(PC)
    case 1
        ind_U_ref(:,ndims(fU_ref.basis)+1) = sum(ind_U_ref(:,1:ndims(fU_ref.basis)),2);
        ind_w_ref_end = cellfun(@(ind,f) sum(ind(:,1:ndims(f.basis)),2),ind_w_ref,fw_ref,'UniformOutput',false);
        ind_lambda_ref_end = cellfun(@(ind,f) sum(ind(:,1:ndims(f.basis)),2),ind_lambda_ref,flambda_ref,'UniformOutput',false);
    case 2
        ind_U_ref(:,ndims(fU_ref.basis)+1) = max(ind_U_ref(:,1:ndims(fU_ref.basis)),[],2);
        ind_w_ref_end = cellfun(@(ind,f) max(ind(:,1:ndims(f.basis)),[],2),ind_w_ref,fw_ref,'UniformOutput',false);
        ind_lambda_ref_end = cellfun(@(ind,f) max(ind(:,1:ndims(f.basis)),[],2),ind_lambda_ref,flambda_ref,'UniformOutput',false);
end
ind_w_ref = cellfun(@(ind,ind_end) [ind,ind_end],ind_w_ref,ind_w_ref_end,'UniformOutput',false);
ind_lambda_ref = cellfun(@(ind,ind_end) [ind,ind_end],ind_lambda_ref,ind_lambda_ref_end,'UniformOutput',false);
PC_U_ref = setindices(PC,ind_U_ref,'update');
PC_w_ref = cellfun(@(x) setindices(PC,x,'update'),ind_w_ref,'UniformOutput',false);
PC_lambda_ref = cellfun(@(x) setindices(PC,x,'update'),ind_lambda_ref,'UniformOutput',false);
U_ref = fU_ref.data';
U_ref = PCMATRIX(U_ref,[size(U_ref,1) 1],PC_U_ref);
w_ref = cellfun(@(x) x.data',fw_ref,'UniformOutput',false);
w_ref = cellfun(@(x,PC) PCMATRIX(x,[size(x,1) 1],PC),w_ref,PC_w_ref,'UniformOutput',false);
lambda_ref = cellfun(@(x) x.data',flambda_ref,'UniformOutput',false);
lambda_ref = cellfun(@(x,PC) PCMATRIX(x,[size(x,1) 1],PC),lambda_ref,PC_lambda_ref,'UniformOutput',false);

%% Outputs

fprintf('\n')
fprintf('spatial dimension = %d for U_ref\n',fU_ref.sz)
for k=1:n
    fprintf('                  = %d for w_ref{%u}\n',fw_ref{k}.sz,k)
    fprintf('                  = %d for lambda_ref{%u}\n',flambda_ref{k}.sz,k)
end
fprintf('parametric dimension = %d\n',ndims(fU_ref.basis))% fprintf('parametric dimension = %d\n',numel(rv))
fprintf('basis dimension = %d for U_ref\n',numel(fU_ref.basis))
for k=1:n
    fprintf('                = %d for w_ref{%u}\n',numel(fw_ref{k}.basis),k)
    fprintf('                = %d for lambda_ref{%u}\n',numel(flambda_ref{k}.basis),k)
end
fprintf('order = [ %s ] for U_ref\n',num2str(max(fU_ref.basis.indices.array)))
for k=1:n
    fprintf('      = [ %s ] for w_ref{%u}\n',num2str(max(fw_ref{k}.basis.indices.array)),k)
    fprintf('      = [ %s ] for lambda_ref{%u}\n',num2str(max(flambda_ref{k}.basis.indices.array)),k)
end
% fprintf('multi-index set for U_ref = \n')
% disp(num2str(fU_ref.basis.indices.array))
% for k=1:n
%     fprintf('multi-index set for w_ref{%u} = \n',k)
%     disp(num2str(fw_ref{k}.basis.indices.array))
%     fprintf('multi-index set for lambda_ref{%u} = \n',k)
%     disp(num2str(flambda_ref{k}.basis.indices.array))
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
s.tolStagnation = 1e-1;
s.display = true;
s.displayIterations = false;

IS = IterativeSolver();
IS.maxIterations = 20;
IS.tolerance = eps;
IS.relaxation = 'Aitken';
IS.updateRelaxationParameter = true;
IS.errorCriterion = 'reference';
IS.referenceSolution = {fU_ref,fw_ref,flambda_ref};
IS.display = true;
IS.displayIterations = true;
if iterativeSolver
    [fU,fw,flambda,output] = IS.solveRandom(glob,patches,interfaces,s,bases,ls,rv);
    save(fullfile(pathname,'solution.mat'),'fU','fw','flambda','output');
else
    load(fullfile(pathname,'solution.mat'),'fU','fw','flambda','output');
end

ind_U = fU.basis.indices.array;
ind_w = cellfun(@(x) x.basis.indices.array,fw,'UniformOutput',false);
ind_lambda = cellfun(@(x) x.basis.indices.array,flambda,'UniformOutput',false);
switch gettypebase(PC)
    case 1
        ind_U(:,ndims(fU.basis)+1) = sum(ind_U(:,1:ndims(fU.basis)),2);
        ind_w_end = cellfun(@(ind,f) sum(ind(:,1:ndims(f.basis)),2),ind_w,fw,'UniformOutput',false);
        ind_lambda_end = cellfun(@(ind,f) sum(ind(:,1:ndims(f.basis)),2),ind_lambda,flambda,'UniformOutput',false);
    case 2
        ind_U(:,ndims(fU.basis)+1) = max(ind_U(:,1:ndims(fU.basis)),[],2);
        ind_w_end = cellfun(@(ind,f) max(ind(:,1:ndims(f.basis)),[],2),ind_w,fw,'UniformOutput',false);
        ind_lambda_end = cellfun(@(ind,f) max(ind(:,1:ndims(f.basis)),[],2),ind_lambda,flambda,'UniformOutput',false);
end
ind_w = cellfun(@(ind,ind_end) [ind,ind_end],ind_w,ind_w_end,'UniformOutput',false);
ind_lambda = cellfun(@(ind,ind_end) [ind,ind_end],ind_lambda,ind_lambda_end,'UniformOutput',false);
PC_U = setindices(PC,ind_U,'update');
PC_w = cellfun(@(x) setindices(PC,x,'update'),ind_w,'UniformOutput',false);
PC_lambda = cellfun(@(x) setindices(PC,x,'update'),ind_lambda,'UniformOutput',false);
U = fU.data';
U = PCMATRIX(U,[size(U,1) 1],PC_U);
w = cellfun(@(x) x.data',fw,'UniformOutput',false);
w = cellfun(@(x,PC) PCMATRIX(x,[size(x,1) 1],PC),w,PC_w,'UniformOutput',false);
lambda = cellfun(@(x) x.data',flambda,'UniformOutput',false);
lambda = cellfun(@(x,PC) PCMATRIX(x,[size(x,1) 1],PC),lambda,PC_lambda,'UniformOutput',false);

%% Outputs

fprintf('\n')
fprintf('spatial dimension = %d for U\n',fU.sz)
for k=1:n
    fprintf('                  = %d for w{%u}\n',fw{k}.sz,k)
    fprintf('                  = %d for lambda{%u}\n',flambda{k}.sz,k)
end
fprintf('parametric dimension = %d\n',ndims(fU.basis))
% fprintf('parametric dimension = %d\n',numel(rv))
fprintf('basis dimension = %d for U\n',numel(fU.basis))
for k=1:n
    fprintf('                = %d for w{%u}\n',numel(fw{k}.basis),k)
    fprintf('                = %d for lambda{%u}\n',numel(flambda{k}.basis),k)
end
fprintf('order = [ %s ] for U\n',num2str(max(fU.basis.indices.array)))
for k=1:n
    fprintf('      = [ %s ] for w{%u}\n',num2str(max(fw{k}.basis.indices.array)),k)
    fprintf('      = [ %s ] for lambda{%u}\n',num2str(max(flambda{k}.basis.indices.array)),k)
end
% fprintf('multi-index set for U = \n')
% disp(num2str(fU.basis.indices.array))
% for k=1:n
%     fprintf('multi-index set for w{%u} = \n',k)
%     disp(num2str(fw{k}.basis.indices.array))
%     fprintf('multi-index set for lambda{%u} = \n',k)
%     disp(num2str(flambda{k}.basis.indices.array))
% end
fprintf('elapsed time = %f s\n',output.totalTime)

%% Display domains and meshes

plotDomain(D,D_patch);
mysaveas(pathname,'domain_global_patches',formats,renderer);
mymatlab2tikz(pathname,'domain_global_patches.tex');

% plotPartition(glob,'legend',false);
% mysaveas(pathname,'mesh_partition',formats,renderer);

plotModel(glob,patches,'legend',false);
mysaveas(pathname,'mesh_global_patches',formats,renderer);

% plotModel(glob);
% plotModel(patches);
% plotModel(interfaces);

%% Display evolution of error indicator, stagnation indicator, CPU time, relaxation parameter w.r.t. number of iterations

plotError(output);
mysaveas(pathname,'error','fig');
mymatlab2tikz(pathname,'error.tex');

plotStagnation(output);
mysaveas(pathname,'stagnation','fig');
mymatlab2tikz(pathname,'stagnation.tex');

plotErrorGlobalSolution(output);
mysaveas(pathname,'error_global_solution','fig');
mymatlab2tikz(pathname,'error_global_solution.tex');

plotStagnationGlobalSolution(output);
mysaveas(pathname,'stagnation_global_solution','fig');
mymatlab2tikz(pathname,'stagnation_global_solution.tex');

plotCPUTime(output,'legend',false);
mysaveas(pathname,'cpu_time','fig');
mymatlab2tikz(pathname,'cpu_time.tex');

plotRelaxationParameter(output,'legend',false);
mysaveas(pathname,'relaxation_parameter','fig');
mymatlab2tikz(pathname,'relaxation_parameter.tex');

plotNbSamples(output);
mysaveas(pathname,'nb_samples','fig');
mymatlab2tikz(pathname,'nb_samples.tex');

plotDimStochasticBasis(output);
mysaveas(pathname,'dim_stochastic_basis','fig');
mymatlab2tikz(pathname,'dim_stochastic_basis.tex');

plotCVError(output);
mysaveas(pathname,'cv_error','fig');
mymatlab2tikz(pathname,'cv_error.tex');

%% Display multi-index set

plotMultiIndexSet(fU,'legend',false)
mysaveas(pathname,'multi_index_set_global_solution','fig');
mymatlab2tikz(pathname,'multi_index_set_global_solution.tex');

for k=1:n
    plotMultiIndexSet(fw{k},'legend',false)
    mysaveas(pathname,['multi_index_set_local_solution_' num2str(k)],'fig');
    mymatlab2tikz(pathname,['multi_index_set_local_solution_' num2str(k) '.tex']);
    
    plotMultiIndexSet(flambda{k},'legend',false)
    mysaveas(pathname,['multi_index_set_Lagrange_multiplier_' num2str(k)],'fig');
    mymatlab2tikz(pathname,['multi_index_set_Lagrange_multiplier_' num2str(k) '.tex']);
end

%% Display statistical outputs

% plotStatsAllSolutions(glob,patches,interfaces,U,w,lambda);

plotMeanGlobalSolution(glob,U);
mysaveas(pathname,'mean_global_solution',formats,renderer);

% plotMeanLocalSolution(patches,w);
% mysaveas(pathname,'mean_local_solution',formats,renderer);

% plotMeanLagrangeMultiplier(interfaces,lambda);
% mysaveas(pathname,'mean_Lagrange_multiplier',formats,renderer);

plotMeanMultiscaleSolution(glob,patches,interfaces,U,w);
mysaveas(pathname,'mean_multiscale_solution',formats,renderer);

plotMeanGlobalLocalSolution(glob,patches,interfaces,U,w);
mysaveas(pathname,'mean_global_local_solution',formats,renderer);

plotMeanGlobalLocalSolution(glob,patches,interfaces,U,w,'view3',true);
mysaveas(pathname,'mean_global_local_solution_surf',formats,renderer);

plotVarGlobalSolution(glob,U);
mysaveas(pathname,'var_global_solution',formats,renderer);

% plotVarLocalSolution(patches,w);
% mysaveas(pathname,'var_local_solution',formats,renderer);

% plotVarLagrangeMultiplier(interfaces,lambda);
% mysaveas(pathname,'var_Lagrange_multiplier',formats,renderer);

plotVarMultiscaleSolution(glob,patches,interfaces,U,w);
mysaveas(pathname,'var_multiscale_solution',formats,renderer);

plotVarGlobalLocalSolution(glob,patches,interfaces,U,w);
mysaveas(pathname,'var_global_local_solution',formats,renderer);

plotVarGlobalLocalSolution(glob,patches,interfaces,U,w,'view3',true);
mysaveas(pathname,'var_global_local_solution_surf',formats,renderer);

plotStdGlobalSolution(glob,U);
mysaveas(pathname,'std_global_solution',formats,renderer);

% plotStdLocalSolution(patches,w);
% mysaveas(pathname,'std_local_solution',formats,renderer);

% plotStdLagrangeMultiplier(interfaces,lambda);
% mysaveas(pathname,'std_Lagrange_multiplier',formats,renderer);

plotStdMultiscaleSolution(glob,patches,interfaces,U,w);
mysaveas(pathname,'std_multiscale_solution',formats,renderer);

plotStdGlobalLocalSolution(glob,patches,interfaces,U,w);
mysaveas(pathname,'std_global_local_solution',formats,renderer);

plotStdGlobalLocalSolution(glob,patches,interfaces,U,w,'view3',true);
mysaveas(pathname,'std_global_local_solution_surf',formats,renderer);

for i=1:d
    plotSobolIndicesMultiscaleSolution(glob,patches,interfaces,U,w,i);
    mysaveas(pathname,['sobol_indices_multiscale_solution_var_' num2str(i)],formats,renderer);
    
    plotSensitivityIndicesMaxVarMultiscaleSolution(glob,patches,interfaces,U,w,i);
    mysaveas(pathname,['sensitivity_indices_multiscale_solution_var_' num2str(i)],formats,renderer);
end

%% Display random evaluations of solutions

% nbsamples = 3;
% for i=1:nbsamples
%     xi = random(rv,1,1);
%     U_xi = fU.functionEval(xi);
%     w_xi = cellfun(@(x) x.functionEval(xi),fw,'UniformOutput',false);
%     lambda_xi = cellfun(@(x) x.functionEval(xi),flambda,'UniformOutput',false);
%     % plotAllSolutions(glob,patches.patchEval(xi),interfaces,U_xi',cellfun(@(x) x',w_xi,'UniformOutput',false),cellfun(@(x) x',lambda_xi,'UniformOutput',false));
%     plotGlobalSolution(glob,U_xi');
%     % plotLocalSolution(patches,cellfun(@(x) x',w_xi,'UniformOutput',false));
%     % plotLagrangeMultiplier(interfaces,cellfun(@(x) x',lambda_xi,'UniformOutput',false));
%     plotMultiscaleSolution(glob,patches.patchEval(xi),interfaces,U_xi',cellfun(@(x) x',w_xi,'UniformOutput',false));
% end

myparallel('stop');
