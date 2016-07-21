%% Multiscale stochastic linear diffusion square inclusions isotropic %%
%%--------------------------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');
myparallel('start');

%% Input data

n = 8; % number of patches
filename = ['multiscale_sto_lin_diff_' num2str(n) '_square_inclusions_iso'];
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

D = DOMAIN(2,[0.0,0.0],[1.0,1.0]);

nbelem = [20,20];
glob.S = build_model(D,'nbelem',nbelem);
% cl = 0.05;
% glob.S = build_model(D,'cl',cl,'filename',[pathname 'gmsh_domain']);

% Patches
patches = Patches(n);

D_patch = cell(1,n);
D_patch{1} = DOMAIN(2,[0.1,0.1],[0.3,0.3]);
D_patch{2} = DOMAIN(2,[0.1,0.4],[0.3,0.6]);
D_patch{3} = DOMAIN(2,[0.1,0.7],[0.3,0.9]);
D_patch{4} = DOMAIN(2,[0.4,0.7],[0.6,0.9]);
D_patch{5} = DOMAIN(2,[0.7,0.7],[0.9,0.9]);
D_patch{6} = DOMAIN(2,[0.7,0.4],[0.9,0.6]);
D_patch{7} = DOMAIN(2,[0.7,0.1],[0.9,0.3]);
D_patch{8} = DOMAIN(2,[0.4,0.1],[0.6,0.3]);

nbelem_patch = [40,40];
for k=1:n
    patches.patches{k}.S = build_model(D_patch{k},'nbelem',nbelem_patch);
end
% cl_patch = 0.005;
% for k=1:n
%     patches.patches{k}.S = build_model(D_patch{k},'cl',cl_patch,'filename',[pathname 'gmsh_patch_' num2str(k)]);
% end

% Partition of global mesh
glob = partition(glob,patches);

%% Random variables

d = n; % parametric dimension
v = UniformRandomVariable(0,1);
rv = RandomVector(v,d);

V = RVUNIFORM(0,1);
RV = RANDVARS(repmat({V},1,d));
[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',2);

%% Materials

% Linear diffusion coefficient
K_out = 1;
K_patch = cell(1,n);
K_in = cell(1,n);

p = 1;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,d);
vb = basis.basis.randomVariable;
rvb = getRandomVector(bases);
H = FullTensorProductFunctionalBasis(bases);
I = gaussIntegrationRule(vb,2);
I = I.tensorize(d);

for k=1:n
    patch = patches.patches{k};
    % K_patch(x,xi) = 1 + f(x) * xi
    % K_in(x)       = 1
    % with f(x) = 1 if ||x-c||_Inf < L
    %           = 0 if ||x-c||_Inf >= L
    L = norm(getsize(D_patch{k}),Inf)/4;
    c = getcenter(D_patch{k});
    f = @(x) distance(x,c,Inf)<L;
    
    fun = @(xi) ones(size(xi,1),patch.S.nbnode) + xi(:,k) * double(squeeze(f(patch.S.node)))';
    funtr = @(xi) fun(transfer(rvb,rv,xi));
    fun = MultiVariateFunction(funtr,d,patch.S.nbnode);
    fun.evaluationAtMultiplePoints = true;
    
    K_patch{k} = H.projection(fun,I);
    K_patch{k} = PCMATRIX(permute(K_patch{k}.tensor.data,[d+1 1:d]),[patch.S.nbnode 1],PC);
    % K_patch{k} = ones(patch.S.nbnode,1,PC) + double(squeeze(f(patch.S.node))) * X{k};
    
    K_in{k} = 1;
end

% Complementary subdomain
mat_out = FOUR_ISOT('k',K_out); % uniform value
mat_out = setnumber(mat_out,0);
glob.S = setmaterial(glob.S,mat_out,getnumgroupelemwithparam(glob.S,'partition',0));

% Patches
mat_patch = MATERIALS();
for k=1:n
    % mat_patch{k} = FOUR_ISOT('k',K_patch{k}); % uniform value
    mat_patch{k} = FOUR_ISOT('k',FENODEFIELD(K_patch{k})); % nodal values
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
f = 100;

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
        patches.patches{k}.A = calc_rigi(patches.patches{k}.S);
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
    patches.patches{k}.ChangeOfVariable = false;
    patches.patches{k}.increment = true;
end

%% Direct solver

p = 50;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,d);
rv = getRandomVector(bases);

s = AdaptiveSparseTensorAlgorithm();
% s.nbSamples = 1;
% s.addSamplesFactor = 0.1;
s.tol = 1e-8;
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
if directSolver
    [fU_ref,fw_ref,flambda_ref,output_ref] = DS.solveRandom(glob_out,patches,interfaces,s,bases,ls,rv);
    save(fullfile(pathname,'reference_solution.mat'),'fU_ref','fw_ref','flambda_ref','output_ref');
else
    load(fullfile(pathname,'reference_solution.mat'),'fU_ref','fw_ref','flambda_ref','output_ref');
end

ind_ref = output_ref.f.basis.indices.array;
switch gettypebase(PC)
    case 1
        ind_ref(:,ndims(output_ref.f.basis)+1) = sum(ind_ref(:,1:ndims(output_ref.f.basis)),2);
    case 2
        ind_ref(:,ndims(output_ref.f.basis)+1) = max(ind_ref(:,1:ndims(output_ref.f.basis)),[],2);
end
PC_ref = setindices(PC,ind_ref,'update');
u_ref = output_ref.f.data';
u_ref = PCMATRIX(u_ref,[size(u_ref,1) 1],PC_ref);
U_ref = fU_ref.data';
U_ref = PCMATRIX(U_ref,[size(U_ref,1) 1],PC_ref);
w_ref = cellfun(@(x) x.data',fw_ref,'UniformOutput',false);
w_ref = cellfun(@(x) PCMATRIX(x,[size(x,1) 1],PC_ref),w_ref,'UniformOutput',false);
lambda_ref = cellfun(@(x) x.data',flambda_ref,'UniformOutput',false);
lambda_ref = cellfun(@(x) PCMATRIX(x,[size(x,1) 1],PC_ref),lambda_ref,'UniformOutput',false);

%% Outputs

fprintf('\n')
fprintf('parametric dimension = %d\n',ndims(output_ref.f.basis))% fprintf('parametric dimension = %d\n',numel(rv))
fprintf('basis dimension = %d\n',numel(output_ref.f.basis))
fprintf('order = [ %s ]\n',num2str(max(output_ref.f.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(output_ref.f.basis.indices.array)
fprintf('nb samples = %d\n',output_ref.N)
fprintf('CV error = %d\n',norm(output_ref.CVError))
fprintf('elapsed time = %f s\n',output_ref.time)

%% Global-local Iterative solver

s.tol = 1e-4;
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

ind_glob = fU.basis.indices.array;
ind_patch = cellfun(@(x) x.basis.indices.array,output.f,'UniformOutput',false);
switch gettypebase(PC)
    case 1
        ind_glob(:,ndims(fU.basis)+1) = sum(ind_glob(:,1:ndims(fU.basis)),2);
        ind_patch_end = cellfun(@(ind,f) sum(ind(:,1:ndims(f.basis)),2),ind_patch,output.f,'UniformOutput',false);
    case 2
        ind_glob(:,ndims(fU.basis)+1) = max(ind_glob(:,1:ndims(fU.basis)),[],2);
        ind_patch_end = cellfun(@(ind,f) max(ind(:,1:ndims(f.basis)),[],2),ind_patch,output.f,'UniformOutput',false);
end
ind_patch = cellfun(@(ind,ind_end) [ind,ind_end],ind_patch,ind_patch_end,'UniformOutput',false);
PC_glob = setindices(PC,ind_glob,'update');
PC_patch = cellfun(@(x) setindices(PC,x,'update'),ind_patch,'UniformOutput',false);
U = fU.data';
U = PCMATRIX(U,[size(U,1) 1],PC_glob);
w = cellfun(@(x) x.data',fw,'UniformOutput',false);
w = cellfun(@(x,PC) PCMATRIX(x,[size(x,1) 1],PC),w,PC_patch,'UniformOutput',false);
lambda = cellfun(@(x) x.data',flambda,'UniformOutput',false);
lambda = cellfun(@(x,PC) PCMATRIX(x,[size(x,1) 1],PC),lambda,PC_patch,'UniformOutput',false);

%% Outputs

fprintf('\n')
fprintf('parametric dimension = %d\n',ndims(fU.basis))% fprintf('parametric dimension = %d\n',numel(rv))
fprintf('Global solution : basis dimension = %d\n',numel(fU.basis))
fprintf('                  order = [ %s ]\n',num2str(max(fU.basis.indices.array)))
% fprintf('                  multi-index set = \n')
% disp([repmat('                  ',numel(fU.basis),1) num2str(fU.basis.indices.array)])

for k=1:n
    fprintf('Local solution #%2.d : basis dimension = %d\n',k,numel(fw{k}.basis))
    % fprintf('                     multi-index set = \n')
    % disp([repmat('                     ',numel(fw{k}.basis),1) num2str(fw{k}.basis.indices.array)])
    fprintf('                     order = [ %s ]\n',num2str(max(fw{k}.basis.indices.array)))
    fprintf('                     nb samples = %d\n',output.N{k})
    fprintf('                     CV error = %d\n',norm(output.CVError{k}))
    fprintf('                     elapsed time = %f s\n',output.time(k))
end

%% Save variables

save(fullfile(pathname,'all.mat'));

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

%% Display multi-index set

for i=1:2:d
    plotMultiIndexSet(fU,'dim',[i i+1],'legend',false)
    mysaveas(pathname,'multi_index_set_global_solution','fig');
    mymatlab2tikz(pathname,'multi_index_set_global_solution.tex');
end

for k=1:n
    for i=1:2:d
        plotMultiIndexSet(fw{k},'dim',[i i+1],'legend',false)
        mysaveas(pathname,['multi_index_set_local_solution_' num2str(k) '_dim_' num2str(i) '_' num2str(i+1)],'fig');
        mymatlab2tikz(pathname,['multi_index_set_local_solution_' num2str(k) '_dim_' num2str(i) '_' num2str(i+1) '.tex']);
        
        plotMultiIndexSet(flambda{k},'dim',[i i+1],'legend',false)
        mysaveas(pathname,['multi_index_set_Lagrange_multiplier_' num2str(k) '_dim_' num2str(i) '_' num2str(i+1)],'fig');
        mymatlab2tikz(pathname,['multi_index_set_Lagrange_multiplier_' num2str(k) '_dim_' num2str(i) '_' num2str(i+1) '.tex']);
    end
end

%% Display statistical outputs

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
