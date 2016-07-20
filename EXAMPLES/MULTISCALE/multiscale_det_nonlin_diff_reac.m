%% Multiscale deterministic nonlinear diffusion reaction %%
%%-------------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% myparallel('start');

%% Input data

n = 4; % number of patches n = 1, 2, 4
filename = ['multiscale_det_nonlin_diff_reac_' num2str(n) '_patches'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc2'};
renderer = 'OpenGL';

solve_reference = true;
solve_multiscale = true;

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
switch n
    case 1
        D_patch{1} = DOMAIN(2,[0.4,0.4],[0.6,0.6]);
    case 2
        D_patch{1} = DOMAIN(2,[0.2,0.4],[0.4,0.6]);
        D_patch{2} = DOMAIN(2,[0.6,0.4],[0.8,0.6]);
    case 4
        D_patch{1} = DOMAIN(2,[0.2,0.2],[0.4,0.4]);
        D_patch{2} = DOMAIN(2,[0.2,0.6],[0.4,0.8]);
        D_patch{3} = DOMAIN(2,[0.6,0.6],[0.8,0.8]);
        D_patch{4} = DOMAIN(2,[0.6,0.2],[0.8,0.4]);
    otherwise
        error('Wrong number of patches')
end

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

%% Materials

% Linear diffusion coefficient
K_out = 1;
K_patch = cell(1,n);
K_in = cell(1,n);
% Nonlinear reaction parameter
R_patch = cell(1,n);
for k=1:n
    patch = patches.patches{k};
    % K_patch(x)  = 1 + f(x)
    % K_in(x)     = 1
    % K2_patch(x) = f(x)
    % R_patch(x)  = f(x)   
    % with f(x) = 1 if ||x-c||_Inf < L
    %           = 0 if ||x-c||_Inf >= L
    L = norm(getsize(D_patch{k}),Inf)/4;
    c = getcenter(D_patch{k});
    f = @(x) distance(x,c,Inf)<L;
    K_patch{k} = ones(patch.S.nbnode,1) + double(squeeze(f(patch.S.node)));
    R_patch{k} = double(squeeze(f(patch.S.node)));
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
    patches.patches{k}.A = @(u) calc_fint(patches.patches{k}.S,u);
    patches.patches{k}.Atang = @(u) calc_rigitang(patches.patches{k}.S,u);
    patches.patches{k}.b = bodyload(patches.patches{k}.S,[],'QN',f);
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

DS = DirectSolver();
DS.changeOfVariable = false;
DS.display = true;
DS.solver = NEWTONSOLVER('type','tangent','increment',true,...
    'maxiter',100,'tol',1e-12,'display',false,'stopini',true);
DS.initializationType = 'zero';
if solve_reference
    [U_ref,w_ref,lambda_ref] = DS.solve(glob_out,patches,interfaces);
    save(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref');
else
    load(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref');
end

%% Global-local Iterative solver

IS = IterativeSolver();
IS.maxIterations = 50;
IS.tolerance = eps;
IS.relaxation = 'Aitken';
IS.updateRelaxationParameter = true;
IS.errorCriterion = 'reference';
IS.referenceSolution = {U_ref,w_ref,lambda_ref};
IS.display = true;
IS.displayIterations = true;
if solve_multiscale
    [U,w,lambda,output] = IS.solve(glob,patches,interfaces);
    save(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
else
    load(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
end
fprintf('\n');

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

%% Display evolution of error indicator, stagnation indicator and CPU time w.r.t. number of iterations

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

%% Display reference and multscale solutions

% if exist('U_ref','var') && exist('w_ref','var') && exist('lambda_ref','var')
%     plotAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref);
% end
% 
% plotAllSolutions(glob,patches,interfaces,U,w,lambda);

plotGlobalSolution(glob,U);
mysaveas(pathname,'global_solution',formats,renderer);

plotMultiscaleSolution(glob,patches,interfaces,U,w);
mysaveas(pathname,'multiscale_solution',formats,renderer);

plotGlobalLocalSolution(glob,patches,interfaces,U,w);
mysaveas(pathname,'global_local_solution',formats,renderer);

plotGlobalLocalSolution(glob,patches,interfaces,U,w,'view3',true);
mysaveas(pathname,'global_local_solution_surf',formats,renderer);

% myparallel('stop');
