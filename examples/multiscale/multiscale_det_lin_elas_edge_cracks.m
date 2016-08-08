%% Multiscale deterministic linear elasticity edge cracks %%
%%--------------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% myparallel('start');

%% Input data

n = 1; % number of patches n = 1
loading = 'pull'; % 'pull' or 'shear'
filename = ['multiscale_det_lin_elas_' num2str(n) '_edge_cracks_' loading];
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

L = 16;
w = 7;
D = DOMAIN(2,[0.0,-L/2],[w,L/2]);

nbelem = [w*2,L*2];
glob.S = build_model(D,'nbelem',nbelem);
% cl = 0.25;
% glob.S = build_model(D,'cl',cl,'filename',[pathname 'gmsh_domain']);

% Patches
patches = Patches(n);

a = w/2;
l = 1;
D_patch = cell(1,n);
B_patch = cell(1,n);
P_patch = cell(1,n);
switch n
    case 1
        P_patch{1} = [a,0.0];
        B_patch{1} = LIGNE([0.0,0.0],P_patch{1});
        D_patch{1} = DOMAIN(2,[0.0,-l],P_patch{1}+[l,l]);
    otherwise
        error('Wrong number of patches')
end

cl_patch_D = 0.25;
cl_patch_P = 0.05;
for k=1:n
    patches.patches{k}.S = gmshdomainwithedgecrack(D_patch{k},P_patch{k},cl_patch_D,cl_patch_P,[pathname 'gmsh_patch_' num2str(k) '_edge_crack']);
end

% Partition of global mesh
glob = partition(glob,patches);

%% Materials

% Poisson ratio
NU = 0.3;
% Thickness
DIM3 = 1;
% Density
RHO = 1;
% Young modulus
E_out = 1;
E_patch = cell(1,n);
E_in = cell(1,n);
for k=1:n
    patch = patches.patches{k};
    % E_patch(x) = 1 + f(x)
    % E_in(x)    = 1
    % with f(x) = 1 if ||x-c||_Inf < L
    %           = 0 if ||x-c||_Inf >= L
    % L = norm(getsize(D_patch{k}),Inf)/4;
    % c = getcenter(D_patch{k});
    % f = @(x) distance(x,c,Inf)<L;
    % E_patch{k} = ones(patch.S.nbnode,1) + double(squeeze(f(patch.S.node)));
    E_patch{k} = 1;
    E_in{k} = 1;
end

% Complementary subdomain
mat_out = ELAS_ISOT('E',E_out,'NU',NU,'RHO',RHO,'DIM3',DIM3);
mat_out = setnumber(mat_out,0);
glob.S = setmaterial(glob.S,mat_out,getnumgroupelemwithparam(glob.S,'partition',0));

% Patches
mat_patch = MATERIALS();
for k=1:n
    mat_patch{k} = ELAS_ISOT('E',E_patch{k},'NU',NU,'RHO',RHO,'DIM3',DIM3); % uniform value
    % mat_patch{k} = ELAS_ISOT('E',FENODEFIELD(E_patch{k}),'NU',NU,'RHO',RHO,'DIM3',DIM3); % nodal values
    mat_patch{k} = setnumber(mat_patch{k},k);
    patches.patches{k}.S = setmaterial(patches.patches{k}.S,mat_patch{k});
end

% Fictitious patches
mat_in = MATERIALS();
for k=1:n
    mat_in{k} = ELAS_ISOT('E',E_in{k},'NU',NU,'RHO',RHO,'DIM3',DIM3); % uniform value
    % mat_in{k} = ELAS_ISOT('E',FENODEFIELD(E_in{k}),'NU',NU,'RHO',RHO,'DIM3',DIM3); % nodal values
    mat_in{k} = setnumber(mat_in{k},k);
    glob.S = setmaterial(glob.S,mat_in{k},getnumgroupelemwithparam(glob.S,'partition',k));
end

%% Dirichlet boundary conditions

LU = LIGNE([0.0,L/2],[w,L/2]);
LL = LIGNE([0.0,-L/2],[w,-L/2]);
LM = LIGNE([a,0.0],[w,0.0]);

% Global
glob.S = final(glob.S);
switch loading
    case 'pull'
        glob.S = addcl(glob.S,LM,'UY');
        glob.S = addcl(glob.S,POINT([w,0.0]),'UX');
    case 'shear'
        glob.S = addcl(glob.S,LL);
    otherwise
        error('Wrong loading case')
end
glob.S_out = getfinalmodelpart(glob.S,0);
% S_in = cell(1,n);
% for k=1:n
%     S_in{k} = getfinalmodelpart(glob.S,k);
% end

% Complementary subdomain
glob_out.S_out = glob.S_out;

% Patches
for k=1:n
    patches.patches{k}.S = final(patches.patches{k}.S,'duplicate');
end

% Interfaces
interfaces = Interfaces(patches,glob);

%% Stiffness matrices and sollicitation vectors

% Traction force density
f = 1;

% Complementary subdomain
glob_out.A_out = calc_rigi(glob_out.S_out);
switch loading
    case 'pull'
        glob_out.b_out = surfload(glob_out.S_out,LU,{'FX','FY'},[0;f]);
        glob_out.b_out = glob_out.b_out + surfload(glob_out.S_out,LL,{'FX','FY'},[0;-f]);
    case 'shear'
        glob_out.b_out = surfload(glob_out.S_out,LU,{'FX','FY'},[f;0]);
    otherwise
        error('Wrong loading case')
end

% Global
glob.A = calc_rigi(glob.S);
for k=1:n
    glob.A_in{k} = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
end
switch loading
    case 'pull'
        glob.b_out = surfload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),LU,{'FX','FY'},[0;f]);
        glob.b_out = glob.b_out + surfload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),LL,{'FX','FY'},[0;-f]);
    case 'shear'
        glob.b_out = surfload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),LU,{'FX','FY'},[f;0]);
    otherwise
        error('Wrong loading case')
end

% Patches
for k=1:n
    patches.patches{k}.A = calc_rigi(patches.patches{k}.S);
    patches.patches{k}.b = sparse(getnbddlfree(patches.patches{k}.S),1);
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
end

%% Direct solver

DS = DirectSolver();
DS.changeOfVariable = false;
DS.display = true;
if directSolver
    [U_ref,w_ref,lambda_ref] = DS.solve(glob_out,patches,interfaces);
    save(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref');
else
    load(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref');
end

%% Outputs

fprintf('\n')
fprintf('spatial dimension = %d for U_ref\n',length(U_ref))
for k=1:n
    fprintf('                  = %d for w_ref{%u}\n',length(w_ref{k}),k)
    fprintf('                  = %d for lambda_ref{%u}\n',length(lambda_ref{k}),k)
end

%% Global-local Iterative solver

IS = IterativeSolver();
IS.maxIterations = 100;
IS.tolerance = eps;
IS.relaxation = 'Aitken';
IS.updateRelaxationParameter = true;
IS.errorCriterion = 'reference';
IS.referenceSolution = {U_ref,w_ref,lambda_ref};
IS.display = true;
IS.displayIterations = true;
if iterativeSolver
    [U,w,lambda,output] = IS.solve(glob,patches,interfaces);
    save(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
else
    load(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
end

%% Outputs

fprintf('\n')
fprintf('spatial dimension = %d for U\n',length(U))
for k=1:n
    fprintf('                  = %d for w{%u}\n',length(w{k}),k)
    fprintf('                  = %d for lambda{%u}\n',length(lambda{k}),k)
end

%% Display domains and meshes

plotDomain(D,cellfun(@(x,y) {x,y},D_patch,B_patch,'UniformOutput',false));
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

%% Display solutions

for i=1:2
    % plotAllSolutions(glob,patches,interfaces,U,w,lambda,'displ',i);
    % mysaveas(pathname,['all_solutions_' num2str(i)],formats,renderer);
    
    plotGlobalSolution(glob,U,'displ',i);
    mysaveas(pathname,['global_solution_' num2str(i)],formats,renderer);
    
    % plotLocalSolution(patches,w,'displ',i);
    % mysaveas(pathname,['local_solution_' num2str(i)],formats,renderer);
    %
    % plotLagrangeMultiplier(interfaces,lambda,'displ',i);
    % mysaveas(pathname,['Lagrange_multiplier_' num2str(i)],formats,renderer);

    plotMultiscaleSolution(glob,patches,interfaces,U,w,'displ',i);
    mysaveas(pathname,['multiscale_solution_' num2str(i)],formats,renderer);
    
    plotGlobalLocalSolution(glob,patches,interfaces,U,w,'displ',i);
    mysaveas(pathname,['global_local_solution_' num2str(i)],formats,renderer);
    
    plotGlobalLocalSolution(glob,patches,interfaces,U,w,'displ',i,'view3',true);
    mysaveas(pathname,['global_local_solution_' num2str(i) '_surf'],formats,renderer);
end

% myparallel('stop');
