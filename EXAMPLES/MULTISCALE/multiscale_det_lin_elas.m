%% Multiscale deterministic linear elasticity %%
%%--------------------------------------------%%

% clc
clear all
close all

%% Input data

n = 4; % number of patches n = 1, 2, 4
filename = ['multiscale_det_lin_elas_' num2str(n) '_patches'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE','RESULTS',filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
% set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
formats = {'fig','epsc2'};
renderer = 'OpenGL';

solve_reference = true;
solve_multiscale = true;

% Parallel computing
% myparallel('start');

%% Domains and meshes

% Global
glob = GLOBAL();
glob_out = GLOBALOUT();

D = DOMAIN(2,[0.0,0.0],[1.0,1.0]);

nbelem = [20,20];
glob.S = build_model(D,'nbelem',nbelem);
% cl = 0.05;
% system.S = build_model(D,'cl',cl,'filename',[pathname 'gmsh_domain']);

% Patches
patches = PATCHES(n);

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
    patches.PATCH{k}.S = build_model(D_patch{k},'nbelem',nbelem_patch);
end
% cl_patch = 0.005;
% for k=1:n
%     patches.PATCH{k}.S = build_model(D_patch{k},'cl',cl_patch,'filename',[pathname 'gmsh_patch_' num2str(k)]);
% end

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
    patch = patches.PATCH{k};
    % E_patch(x) = 1 + beta_patch * f(x)
    % E_in(x)    = 1 + beta_in * f(x)
    % with f(x) = alpha*exp( -Amp*||x-c||_2^2/L^2 ) if ||x-c||_Inf < L
    %           = 0                                 if ||x-c||_Inf >= L
    % alpha = 10;
    % Amp = 2;
    % L = norm(getsize(D_patch{k}),Inf)/4;
    % c = getcenter(D_patch{k});
    % f = @(x) (distance(x,c,Inf)<L) * alpha * exp(-Amp*distance(x,c,2).^2/L^2);
    % beta_patch = 1;
    % beta_in = 0;
    % E_patch{k} = 1 + beta_patch * squeeze(f(patch.S.node));
    % E_in{k} = 1 + beta_in * squeeze(f(glob.S.node));
    
    % E_patch(x) = 1 + f(x)
    % E_in(x)    = 1
    % with f(x) = 1 if ||x-c||_Inf < L
    %           = 0 if ||x-c||_Inf >= L
    L = norm(getsize(D_patch{k}),Inf)/4;
    c = getcenter(D_patch{k});
    f = @(x) distance(x,c,Inf)<L;
    E_patch{k} = ones(patch.S.nbnode,1) + double(squeeze(f(patch.S.node)));
    E_in{k} = 1;
end

% Outside
mat_out = ELAS_ISOT('E',E_out,'NU',NU,'RHO',RHO,'DIM3',DIM3);
mat_out = setnumber(mat_out,0);
glob.S = setmaterial(glob.S,mat_out,getnumgroupelemwithparam(glob.S,'partition',0));

% Patches
mat_patch = MATERIALS();
for k=1:n
    % mat_patch{k} = ELAS_ISOT('E',E_patch{k},'NU',NU,'RHO',RHO,'DIM3',DIM3); % uniform value
    mat_patch{k} = ELAS_ISOT('E',FENODEFIELD(E_patch{k}),'NU',NU,'RHO',RHO,'DIM3',DIM3); % nodal values
    mat_patch{k} = setnumber(mat_patch{k},k);
    patches.PATCH{k}.S = setmaterial(patches.PATCH{k}.S,mat_patch{k});
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

% Global
glob.S = final(glob.S);
glob.S = addcl(glob.S,[]);
glob.S_out = get_final_model_part(glob.S,0);
glob_out.S_out = glob.S_out;
% S_in = cell(1,n);
% for k=1:n
%     S_in{k} = get_final_model_part(glob.S,k);
% end

% Patches
for k=1:n
    patches.PATCH{k}.S = final(patches.PATCH{k}.S);
end

% Interfaces
interfaces = INTERFACES(patches);

%% Stiffness matrices and sollicitation vectors

% Traction force density
f = [0;-100];

% Outside
glob_out.A_out = calc_rigi(glob_out.S_out);
glob_out.b_out = bodyload(glob_out.S_out,[],{'FX','FY'},f);

% Global
glob.A = calc_rigi(glob.S);
for k=1:n
    glob.A_in{k} = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
end
glob.b_out = bodyload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),[],{'FX','FY'},f);

% Patches
for k=1:n
    patches.PATCH{k}.A = calc_rigi(patches.PATCH{k}.S);
    patches.PATCH{k}.b = bodyload(patches.PATCH{k}.S,[],{'FX','FY'},f);
end

%% Mass matrices

for k=1:n
    interfaces.INTERFACE{k}.M = calc_massgeom(interfaces.INTERFACE{k}.S);
end

%% Projection operators

glob.P_out = calc_P_free(glob.S,glob.S_out);
for k=1:n
    [interfaces.INTERFACE{k}.P_glob] = calc_projection(interfaces.INTERFACE{k},glob);
    [interfaces.INTERFACE{k}.P_glob_out,numnode] = calc_projection(interfaces.INTERFACE{k},glob_out);
    interfaces.INTERFACE{k}.P_patch = calc_P_free(patches.PATCH{k}.S,interfaces.INTERFACE{k}.S);
    % plot_projection_operator(glob,patches.PATCH{k},numnode);
end

%% Parameters for global and local problems

% Global problem
glob.param = setparam(glob.param,'increment',true);
glob.param = setparam(glob.param,'inittype','zero');

% Local problems
for k=1:n
    patches.PATCH{k}.param = setparam(patches.PATCH{k}.param,'change_of_variable',false);
    patches.PATCH{k}.param = setparam(patches.PATCH{k}.param,'increment',true);
    patches.PATCH{k}.param = setparam(patches.PATCH{k}.param,'inittype','zero');
end

%% Monoscale resolution

R = REFERENCESOLVER('display',true,'change_of_variable',false,'inittype','zero');
if solve_reference
    [U_ref,w_ref,lambda_ref] = solve(R,glob_out,patches,interfaces);
    save(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref');
else
    load(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref');
end

%% Multiscale resolution

I = ITERATIVESOLVER('display',true,'displayiter',true,...
    'maxiter',50,'tol',eps,'rho','Aitken',...
    'errorindicator','reference','reference',{{U_ref,w_ref,lambda_ref}});
if solve_multiscale
    [U,w,lambda,result] = solve(I,glob,patches,interfaces);
    save(fullfile(pathname,'solution.mat'),'U','w','lambda','result');
else
    load(fullfile(pathname,'solution.mat'),'U','w','lambda','result');
end
fprintf('\n');

%% Save variables

save(fullfile(pathname,'all.mat'));

%% Display domains and meshes

plot_domain(D,D_patch);
mysaveas(pathname,'domain_global_patches',formats,renderer);
mymatlab2tikz(pathname,'domain_global_patches.tex');

% plot_partition(glob,'nolegend');
% mysaveas(pathname,'mesh_partition',formats,renderer);

plot_model(glob,patches,'nolegend');
mysaveas(pathname,'mesh_global_patches',formats,renderer);

% plot_model(glob);
% plot_model(patches);
% plot_model(interfaces);

%% Display evolution of error indicator, stagnation indicator and CPU time w.r.t. number of iterations

plot_error_indicator_deterministic(result);
mysaveas(pathname,'error_indicator','fig');
mymatlab2tikz(pathname,'error_indicator.tex');

plot_stagnation_indicator_deterministic(result);
mysaveas(pathname,'stagnation_indicator','fig');
mymatlab2tikz(pathname,'stagnation_indicator.tex');

plot_error_indicator_U_deterministic(result);
mysaveas(pathname,'error_indicator_U','fig');
mymatlab2tikz(pathname,'error_indicator_U.tex');

plot_stagnation_indicator_U_deterministic(result);
mysaveas(pathname,'stagnation_indicator_U','fig');
mymatlab2tikz(pathname,'stagnation_indicator_U.tex');

plot_cpu_time(result,'nolegend');
mysaveas(pathname,'cpu_time','fig');
mymatlab2tikz(pathname,'cpu_time.tex');

plot_relaxation_parameter(result,'nolegend');
mysaveas(pathname,'relaxation_parameter','fig');
mymatlab2tikz(pathname,'relaxation_parameter.tex');

%% Display reference and multscale solutions

% if exist('U_ref','var') && exist('w_ref','var') && exist('lambda_ref','var')
%     for i=1:2
%         plot_sols_ref(glob,patches,interfaces,U_ref,w_ref,lambda_ref,'displ',i);
%     end
% end

% for i=1:2
%     plot_sols(glob,patches,interfaces,U,w,lambda,'displ',i);
% end

for i=1:2
    plot_U(glob,U,'displ',i);
    mysaveas(pathname,['U_' num2str(i)],formats,renderer);
    
    plot_sol(glob,patches,interfaces,U,w,'displ',i);
    mysaveas(pathname,['sol_' num2str(i)],formats,renderer);
    
    plot_U_w(glob,patches,interfaces,U,w,'displ',i);
    mysaveas(pathname,['U_w' num2str(i)],formats,renderer);
    
    plot_U_w(glob,patches,interfaces,U,w,'displ',i,'view3');
    mysaveas(pathname,['U_w_surf' num2str(i)],formats,renderer);
end

% myparallel('stop');
