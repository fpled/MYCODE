%% Multiscale deterministic linear diffusion form %%
%%------------------------------------------------%%

% clc
clear all
close all

%% Input data
n = 4; % number of patches n = 1, 2, 4
filename = ['multiscale_det_lin_diff_form_' num2str(n) '_patches'];
pathname = [getfemobjectoptions('path') 'MYCODE/RESULTS/' filename '/'];
if ~exist(pathname,'dir')
    mkdir(pathname);
end
set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
renderer = 'OpenGL';

% Reference solution - Direct resolution of initial problem based on non-overlapping domain decomposition
solve_reference = true;

% Multscale solution - Reformulated global-local iterative algorithm based on overlapping domain decomposition
solve_multiscale = true;

% Parallel computing
% myparallel('start');

%% Domain and mesh definition

% Global
glob = GLOBAL();
glob_out = GLOBALOUT();

D = DOMAIN(2,[0.0,0.0],[5.0,5.0]);
nbelem = [20,20];
glob.S = build_model(D,'nbelem',nbelem);
% cl = 0.25;
% glob.S = build_model(D,'cl',cl);
glob.S = final(glob.S);
glob.S = addcl(glob.S,[]);

% Patches
patches = PATCHES(n);

D_patch = cell(1,n);
switch n
    case 1
        D_patch{1} = DOMAIN(2,[2.0,2.0],[3.0,3.0]);
    case 2
        D_patch{1} = DOMAIN(2,[1.0,2.0],[2.0,3.0]);
        D_patch{2} = DOMAIN(2,[3.0,2.0],[4.0,3.0]);
    case 4
        D_patch{1} = DOMAIN(2,[1.0,1.0],[2.0,2.0]);
        D_patch{2} = DOMAIN(2,[1.0,3.0],[2.0,4.0]);
        D_patch{3} = DOMAIN(2,[3.0,3.0],[4.0,4.0]);
        D_patch{4} = DOMAIN(2,[3.0,1.0],[4.0,2.0]);
    otherwise
        error('Wrong number of patches')
end
nbelem_patch = [40,40];
% cl_patch = 0.025;
for k=1:n
    patches.PATCH{k}.S = build_model(D_patch{k},'nbelem',nbelem_patch);
    % patches.PATCH{k}.S = build_model(D_patch{k},'cl',cl_patch);
    patches.PATCH{k}.S = final(patches.PATCH{k}.S);
end

% Partition of global mesh glob.S
glob = partition(glob,patches);
glob.S_out = get_final_model_part(glob.S,0);
glob_out.S_out = glob.S_out;
% S_in = cell(1,n);
% for k=1:n
%     S_in{k} = get_final_model_part(glob.S,k);
% end

% Interfaces
interfaces = INTERFACES(patches);

%% Bilinear forms and linear forms associated to initial problem

% Linear diffusion coefficients K_out, K_patch and K_in
K_out = 1;
K_patch = cell(1,n);
K_in = cell(1,n);
for k=1:n
    patch = patches.PATCH{k};
    % K_patch(x) = 1 + beta_patch * f(x)
    % K_in(x)    = 1 + beta_in * f(x)
    % with f(x) = alpha*exp( -Amp*||x-c||_2^2/L^2 ) if ||x-c||_Inf < L
    %           = 0                                 if ||x-c||_Inf >= L
    % alpha = 10;
    % Amp = 2;
    % L = norm(getsize(D_patch{k}),Inf)/4;
    % c = getcenter(D_patch{k});
    % f = @(x) (distance(x,c,Inf)<L) * alpha * exp(-Amp*distance(x,c,2).^2/L^2);
    % 
    % beta_patch = 1;
    % P = POINT(patch.S.node);
    % K_patch{k} = 1 + beta_patch * squeeze(f(P));
    % 
    % beta_in = 0;
    % P = POINT(glob.S.node);
    % K_in{k} = 1 + beta_in * squeeze(f(P));
    
    % K_patch(x) = 1 + f(x)
    % K_in(x)    = 1
    % with f(x) = 1 if ||x-c||_Inf < L
    %           = 0 if ||x-c||_Inf >= L
    L = norm(getsize(D_patch{k}),Inf)/4;
    c = getcenter(D_patch{k});
    f = @(x) distance(x,c,Inf)<L;
    P = POINT(patch.S.node);
    K_patch{k} = ones(getnbnode(patch.S),1) + double(squeeze(f(P)));
    
    K_in{k} = 1;
end

% Bilinear form a_out associated to outside subdomain
% a(u,v) = int( K.grad(u).grad(v) )
% a_out = BILINFORM(1,1,K_out);
a_out = DIFFUSIONFORM(K_out);
a_out = setfree(a_out,1);

% Bilinear form a_patch associated to patch
% a(u,v) = int( K.grad(u).grad(v) )
a_patch = cell(1,n);
for k=1:n
    % a_patch{k} = BILINFORM(1,1,K_patch{k}); % uniform value
    % a_patch{k} = BILINFORM(1,1,K_patch{k},0); % nodal values
    a_patch{k} = DIFFUSIONFORM(K_patch{k});
    a_patch{k} = setfree(a_patch{k},1);
end

% Bilinear form a_in associated to fictitious patch
% a(u,v) = int( K.grad(u).grad(v) )
a_in = cell(1,n);
for k=1:n
    % a_in{k} = BILINFORM(1,1,K_in{k}); % uniform value
    % a_in{k} = BILINFORM(1,1,K_in{k},0); % nodal values
    a_in{k} = DIFFUSIONFORM(K_in{k});
    a_in{k} = setfree(a_in{k},1);
end

% Source term f
f = 1;

% Linear form l associated to source term f
% l(v) = int( f.v )
l = LINFORM(0,f);
l = setfree(l,1);

% Linear form l_patch associated to source term f
% l(v) = int( f.v )
l_patch = l;

%% Stiffness matrices and sollicitation vectors

% Stiffness matrix glob_out.A_out and sollicitation vector glob_out.b_out associated to mesh glob_out.S_out
glob_out.A_out = calc_matrix(a_out,glob_out.S_out);
glob_out.b_out = calc_vector(l,glob_out.S_out);

% Stiffness matrices glob.A and glob.A_in and sollicitation vector glob.b_out associated to mesh glob.S
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

% Stiffness matrix patch.A and sollicitation vector patch.b associated to mesh patch.S
for k=1:n
    patches.PATCH{k}.A = calc_matrix(a_patch{k},patches.PATCH{k}.S);
    patches.PATCH{k}.b = calc_vector(l_patch,patches.PATCH{k}.S);
end

%% Mass matrices

% Bilinear form a
% a(u,v) = int( u.v )
% a = BILINFORM(0,0,1);
% a = setfree(a,1);

% Mass matrix interface.M associated to boundary mesh interface.S
for k=1:n
    % interfaces.INTERFACE{k}.M = calc_matrix(a,interfaces.INTERFACE{k}.S);
    interfaces.INTERFACE{k}.M = calc_massgeom(interfaces.INTERFACE{k}.S);
end

%% Projection operators

% Projection operator glob.P_out from mesh glob.S to mesh glob.S_out
glob.P_out = calc_P_free(glob.S,glob.S_out);

% Projection operator interface.P_glob from mesh glob.S to boundary mesh interface.S
for k=1:n
    [interfaces.INTERFACE{k}.P_glob] = calc_projection(interfaces.INTERFACE{k},glob);
    [interfaces.INTERFACE{k}.P_glob_out,numnode] = calc_projection(interfaces.INTERFACE{k},glob_out);
    interfaces.INTERFACE{k}.P_patch = calc_P(patches.PATCH{k}.S,interfaces.INTERFACE{k}.S);
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

%% Direct resolution of initial problem based on non-overlapping domain decomposition

R = REFERENCESOLVER('display',true,'change_of_variable',false,'inittype','zero');
if solve_reference
    [U_ref,w_ref,lambda_ref] = solve(R,glob_out,patches,interfaces);
    save(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref');
else
    if ~exist(fullfile(pathname,'reference_solution.mat'),'file')
        error(['File reference_solution.mat does not exist in folder ' pathname]);
    else
        load(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref');
    end
end

%% Reformulated global-local iterative algorithm based on overlapping domain decomposition

I = ITERATIVESOLVER('display',true,'displayiter',true,...
    'maxiter',50,'tol',eps,'rho','Aitken',...
    'errorindicator','reference','reference',{{U_ref,w_ref,lambda_ref}});
if solve_multiscale
    [U,w,lambda,result] = solve(I,glob,patches,interfaces);
    save(fullfile(pathname,'solution.mat'),'U','w','lambda','result');
else
    if ~exist(fullfile(pathname,'solution.mat'),'file')
        error(['File solution.mat does not exist in folder ' pathname]);
    else
        load(fullfile(pathname,'solution.mat'),'U','w','lambda','result');
    end
end
fprintf('\n');

%% Save all variables

save(fullfile(pathname,'all.mat'));

%% Display domain, partition and mesh

% Display global domain and patches
plot_domain(D,D_patch);
mysaveas(pathname,'domain_global_patches',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(pathname,'domain_global_patches',renderer);
mymatlab2tikz(pathname,'domain_global_patches.tex');

% Display partition of global mesh glob.S
% plot_partition(glob);
% mysaveas(pathname,'mesh_partition',{'fig','epsc2','pdf'},renderer);
% mysaveaspdf(pathname,'mesh_partition',renderer);

% Display global mesh glob.S_out and local meshes patch.S
plot_model(glob,patches,'nolegend');
mysaveas(pathname,'mesh_global_patches',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(pathname,'mesh_global_patches',renderer);

% Display all parts of global mesh glob.S
% plot_model(glob);

% Display local meshes patch.S
% plot_model(patches);

% Display boundary meshes interface.S
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

%% Display reference solution u_ref=(U_ref,w_ref) and multscale solution u=(U,w) at final iteration

% if exist('U_ref','var') && exist('w_ref','var') && exist('lambda_ref','var')
%     plot_sols_ref(glob,patches,interfaces,U_ref,w_ref,lambda_ref);
% end
% 
% plot_sols(glob,patches,interfaces,U,w,lambda);

plot_U(glob,U);
mysaveas(pathname,'U',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(pathname,'U',renderer);

plot_sol(glob,patches,interfaces,U,w);
mysaveas(pathname,'sol',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(pathname,'sol',renderer);

plot_U_w(glob,patches,interfaces,U,w);
mysaveas(pathname,'U_w',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(pathname,'U_w',renderer);

plot_U_w(glob,patches,interfaces,U,w,'surface');
mysaveas(pathname,'U_w_surf',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(pathname,'U_w_surf',renderer);

% myparallel('stop');
