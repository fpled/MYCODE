%% Multiscale stochastic linear diffusion form %%
%%---------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');
myparallel('start');

%% Input data

n = 4; % number of patches n = 1, 2, 4
filename = ['multiscale_sto_lin_diff_form_' num2str(n) '_patches'];
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
glob.S = final(glob.S);
glob.S = addcl(glob.S,[]);

% Patches
patches = Patches(n);

D_patch = cell(1,n);
switch n
    case 1
        D_patch{1} = DOMAIN(2,[0.4,0.4],[0.6,0.6]);
    case 2
        D_patch{1} = DOMAIN(2,[0.1,0.1],[0.3,0.3]);
        D_patch{2} = DOMAIN(2,[0.7,0.7],[0.9,0.9]);
    case 4
        D_patch{1} = DOMAIN(2,[0.1,0.1],[0.3,0.3]);
        D_patch{2} = DOMAIN(2,[0.1,0.7],[0.3,0.9]);
        D_patch{3} = DOMAIN(2,[0.7,0.7],[0.9,0.9]);
        D_patch{4} = DOMAIN(2,[0.7,0.1],[0.9,0.3]);
    otherwise
        error('Wrong number of patches')
end
nbelem_patch = [40,40];
for k=1:n
    patches.patches{k}.S = build_model(D_patch{k},'nbelem',nbelem_patch);
    patches.patches{k}.S = final(patches.patches{k}.S);
end
% cl_patch = 0.005;
% for k=1:n
%     patches.patches{k}.S = build_model(D_patch{k},'cl',cl_patch,'filename',[pathname 'gmsh_patch_' num2str(k)]);
%     patches.patches{k}.S = final(patches.patches{k}.S);
% end

% Partition of global mesh
glob = partition(glob,patches);
glob.S_out = getfinalmodelpart(glob.S,0);
% S_in = cell(1,n);
% for k=1:n
%     S_in{k} = getfinalmodelpart(glob.S,k);
% end

% Complementary subdomain
glob_out.S_out = glob.S_out;

% Interfaces
interfaces = Interfaces(patches);

%% Random variables

d = n; % parametric dimension
v = UniformRandomVariable(0,1);
rv = RandomVector(v,d);

%% Bilinear and linear forms

% Linear diffusion coefficient
K_out = 1;
K_patch = cell(1,n);
K_in = cell(1,n);

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
    % K_patch(x,xi) = 1 + beta_patch * f(x) * xi
    % K_in(x)       = 1 + beta_in * f(x)
    % with f(x) = alpha*exp( -Amp*||x-c||_2^2/L^2 ) if ||x-c||_Inf < L
    %           = 0                                 if ||x-c||_Inf >= L
%     alpha = 10;
%     Amp = 2;
%     L = norm(getsize(D_patch{k}),Inf)/4;
%     c = getcenter(D_patch{k});
%     f = @(x) (distance(x,c,Inf)<L) * alpha * exp(-Amp*distance(x,c,2).^2/L^2);
%
%     beta_patch = 1;
%     beta_in = 0;
%     fun = @(xi) ones(size(xi,1),patch.S.nbnode) + beta_patch * xi(:,k) * double(squeeze(f(patch.S.node)))';
%     funtr = @(xi) fun(transfer(rvb,rv,xi));
%     fun = MultiVariateFunction(funtr,d,patch.S.nbnode);
%     fun.evaluationAtMultiplePoints = true;
%
%     K_patch{k} = H.projection(fun,I);
%     K_in{k} = 1 + beta_in * squeeze(f(glob.S.node));
    
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
    K_in{k} = 1;
end

% Complementary subdomain
% a_out = BILINFORM(1,1,K_out);
a_out = DIFFUSIONFORM(K_out);
a_out = setfree(a_out,1);

% Patches
a_patch = cell(1,n);
for k=1:n
    % a_patch{k} = BILINFORM(1,1,K_patch{k}); % uniform value
    % a_patch{k} = BILINFORM(1,1,K_patch{k},0); % nodal values
    a_patch{k} = DIFFUSIONFORM(K_patch{k});
    a_patch{k} = setfree(a_patch{k},1);
    patches.patches{k}.a = a_patch{k};
end

% Fictitious patches
a_in = cell(1,n);
for k=1:n
    % a_in{k} = BILINFORM(1,1,K_in{k}); % uniform value
    % a_in{k} = BILINFORM(1,1,K_in{k},0); % nodal values
    a_in{k} = DIFFUSIONFORM(K_in{k});
    a_in{k} = setfree(a_in{k},1);
end

% Source term
f = 100;

% Global
l = LINFORM(0,f);
l = setfree(l,1);

% Patches
l_patch = l;

%% Stiffness matrices and sollicitation vectors

% Complementary subdomain
glob_out.A_out = calc_matrix(a_out,glob_out.S_out);
glob_out.b_out = calc_vector(l,glob_out.S_out);

% Global
a_out = setselgroup(a_out,getnumgroupelemwithparam(glob.S,'partition',0));
A_out = calc_matrix(a_out,glob.S);
glob.A = A_out;
for k=1:n
    a_in{k} = setselgroup(a_in{k},getnumgroupelemwithparam(glob.S,'partition',k));
    glob.A_in{k} = calc_matrix(a_in{k},glob.S);
    glob.A = glob.A + glob.A_in{k};
end
l_out = setselgroup(l,getnumgroupelemwithparam(glob.S,'partition',0));
glob.b_out = calc_vector(l_out,glob.S);

% Patches
for k=1:n
    % patches.patches{k}.A = calc_matrix(a_patch{k},patches.patches{k}.S);
    patches.patches{k}.b = calc_vector(l_patch,patches.patches{k}.S);
end

%% Mass matrices

% a = BILINFORM(0,0,1);
% a = setfree(a,1);
for k=1:n
    % interfaces.interfaces{k}.M = calc_matrix(a,interfaces.interfaces{k}.S);
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
end

%% Direct solver

p = 50;
basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
bases = FunctionalBases(basis,[],d);
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

plotMultiIndexSet(U,'legend',false)
mysaveas(pathname,'multi_index_set_global_solution','fig');
mymatlab2tikz(pathname,'multi_index_set_global_solution.tex');

for k=1:n
    plotMultiIndexSet(w{k},'legend',false)
    mysaveas(pathname,['multi_index_set_local_solution_' num2str(k)],'fig');
    mymatlab2tikz(pathname,['multi_index_set_local_solution_' num2str(k) '.tex']);
    
    plotMultiIndexSet(lambda{k},'legend',false)
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

plotVarianceGlobalSolution(glob,U);
mysaveas(pathname,'var_global_solution',formats,renderer);

% plotVarianceLocalSolution(patches,w);
% mysaveas(pathname,'var_local_solution',formats,renderer);

% plotVarianceLagrangeMultiplier(interfaces,lambda);
% mysaveas(pathname,'var_Lagrange_multiplier',formats,renderer);

plotVarianceMultiscaleSolution(glob,patches,interfaces,U,w);
mysaveas(pathname,'var_multiscale_solution',formats,renderer);

plotVarianceGlobalLocalSolution(glob,patches,interfaces,U,w);
mysaveas(pathname,'var_global_local_solution',formats,renderer);

plotVarianceGlobalLocalSolution(glob,patches,interfaces,U,w,'view3',true);
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
    
    plotSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,U,w,i);
    mysaveas(pathname,['sensitivity_indices_multiscale_solution_var_' num2str(i)],formats,renderer);
end

%% Display random evaluations of solutions

% nbsamples = 3;
% for i=1:nbsamples
%     xi = random(rv,1,1);
%     U_xi = U.functionEval(xi);
%     w_xi = cellfun(@(x) x.functionEval(xi),w,'UniformOutput',false);
%     lambda_xi = cellfun(@(x) x.functionEval(xi),lambda,'UniformOutput',false);
%     % plotAllSolutions(glob,patches.patchEval(xi),interfaces,U_xi',cellfun(@(x) x',w_xi,'UniformOutput',false),cellfun(@(x) x',lambda_xi,'UniformOutput',false));
%     plotGlobalSolution(glob,U_xi');
%     % plotLocalSolution(patches,cellfun(@(x) x',w_xi,'UniformOutput',false));
%     % plotLagrangeMultiplier(interfaces,cellfun(@(x) x',lambda_xi,'UniformOutput',false));
%     plotMultiscaleSolution(glob,patches.patchEval(xi),interfaces,U_xi',cellfun(@(x) x',w_xi,'UniformOutput',false));
% end

myparallel('stop');
