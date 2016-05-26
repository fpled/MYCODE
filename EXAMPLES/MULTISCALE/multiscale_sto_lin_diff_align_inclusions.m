%% Multiscale stochastic linear diffusion aligned inclusions %%
%%-----------------------------------------------------------%%

% clc
clear all
close all

% Parallel computing
myparallel('start');

%% Input data

n = 10; % number of patches
filename = ['multiscale_sto_lin_diff_' num2str(n) '_align_inclusions'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE','RESULTS',filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
% set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
formats = {'fig','epsc2'};
renderer = 'OpenGL';

solve_reference = true;
solve_multiscale = true;

calc_MC_error_estimate_ref = false;
calc_MC_error_estimate = false;

%% Domains and meshes

% Global
glob = GLOBAL();
glob_out = GLOBALOUT();

D = DOMAIN(2,[0.0,0.0],[2.0,2*n]);

nbelem = [20,20*n];
glob.S = build_model(D,'nbelem',nbelem);
% cl = 0.25;
% system.S = build_model(D,'cl',cl,'filename',[pathname 'gmsh_domain']);

% Patches
patches = PATCHES(n);

D_patch = cell(1,n);
for k=1:n
    D_patch{k} = DOMAIN(2,[0.5,2*k-1.5],[1.5,2*k-0.5]);
end

nbelem_patch = [40,40];
for k=1:n
    patches.PATCH{k}.S = build_model(D_patch{k},'nbelem',nbelem_patch);
end
% cl_patch = 0.025;
% for k=1:n
%     patches.PATCH{k}.S = build_model(D_patch{k},'cl',cl_patch,'filename',[pathname 'gmsh_patch_' num2str(k)]);
% end

% Partition of global mesh
% glob = partition(glob,patches);
glob = partition(glob,D_patch);

%% Random variables

rv = RVUNIFORM(0,1);
RV = RANDVARS(repmat({rv},1,n));
[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',1);

%% Materials

% Linear diffusion coefficient
K_out = 1;
K_patch = cell(1,n);
K_in = cell(1,n);
for k=1:n
    patch = patches.PATCH{k};
    % K_patch(x,xi) = 1 + f(x) * xi
    % K_in(x)       = 1
    % with f(x) = 1 if ||x-c||_Inf < L
    %           = 0 if ||x-c||_Inf >= L
    L = norm(getsize(D_patch{k}),Inf)/4;
    c = getcenter(D_patch{k});
    f = @(x) distance(x,c,Inf)<L;
    K_patch{k} = ones(patch.S.nbnode,1,PC) + double(squeeze(f(patch.S.node))) * X{k};
    K_in{k} = 1;
end

% Outside
mat_out = FOUR_ISOT('k',K_out); % uniform value
mat_out = setnumber(mat_out,0);
glob.S = setmaterial(glob.S,mat_out,getnumgroupelemwithparam(glob.S,'partition',0));

% Patches
mat_patch = MATERIALS();
for k=1:n
    % mat_patch{k} = FOUR_ISOT('k',K_patch{k}); % uniform value
    mat_patch{k} = FOUR_ISOT('k',FENODEFIELD(K_patch{k})); % nodal values
    mat_patch{k} = setnumber(mat_patch{k},k);
    patches.PATCH{k}.S = setmaterial(patches.PATCH{k}.S,mat_patch{k});
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

% Source term
f = 1;

% Outside
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
    if ~israndom(patches.PATCH{k}.S)
        patches.PATCH{k}.A = calc_rigi(patches.PATCH{k}.S);
    end
    if ~israndom(f)
        patches.PATCH{k}.b = bodyload(patches.PATCH{k}.S,[],'QN',f);
    end
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

initPC = POLYCHAOS(RV,0,'typebase',1);

method_ref = METHOD('type','leastsquares','display',true,'displayiter',true,...
    'basis','adaptive','initPC',initPC,'maxcoeff',Inf,...
    'algorithm','RMS','bulkparam',0.5,...
    'sampling','adaptive','initsample',2,'addsample',0.1,'maxsample',Inf,...
    'regul','','cv','leaveout','k',10,...
    'tol',1e-8,'tolstagn',1e-1,'toloverfit',1.1,'correction',false,...
    'decompKL',false,'tolKL',1e-12,'cvKL','leaveout','kKL',10);

R = REFERENCESOLVER('display',true,'change_of_variable',false,'inittype','zero');
if solve_reference
    [U_ref,w_ref,lambda_ref,result_ref] = solve_random(R,glob_out,patches,interfaces,method_ref);
    save(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','result_ref');
else
    load(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','result_ref');
end

%% Multiscale resolution

method = METHOD('type','leastsquares','display',true,'displayiter',false,...
    'basis','adaptive','initPC',initPC,'maxcoeff',Inf,...
    'algorithm','RMS','bulkparam',0.5,...
    'sampling','adaptive','initsample',2,'addsample',0.1,'maxsample',Inf,...
    'regul','','cv','leaveout','k',10,...
    'tol',1e-4,'tolstagn',1e-1,'toloverfit',1.1,'correction',false,...
    'decompKL',false,'tolKL',1e-12,'cvKL','leaveout','kKL',10);

I = ITERATIVESOLVER('display',true,'displayiter',true,...
    'maxiter',20,'tol',eps,'rho','Aitken',...
    'errorindicator','reference','reference',{{U_ref,w_ref,lambda_ref}});
if solve_multiscale
    [U,w,lambda,result] = solve_random(I,glob,patches,interfaces,method);
    save(fullfile(pathname,'solution.mat'),'U','w','lambda','result');
else
    load(fullfile(pathname,'solution.mat'),'U','w','lambda','result');
end
fprintf('\n');

%% Monte Carlo error estimation

nbsamples = 100;
if calc_MC_error_estimate_ref
    [mean_error_U_ref,mean_error_w_ref,mean_error_lambda_ref,...
        var_error_U_ref,var_error_w_ref,var_error_lambda_ref,...
        std_error_U_ref,std_error_w_ref,std_error_lambda_ref]...
        = calc_Monte_Carlo_error_estimates(R,glob_out,patches,interfaces,U_ref,w_ref,lambda_ref,nbsamples);
end
if calc_MC_error_estimate
    [mean_error_U,mean_error_w,mean_error_lambda,...
        var_error_U,var_error_w,var_error_lambda,...
        std_error_U,std_error_w,std_error_lambda]...
        = calc_Monte_Carlo_error_estimates(I,glob,patches,interfaces,U,w,lambda,nbsamples);
end

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

%% Display evolution of error indicator, stagnation indicator, CPU time, sparsity ratio, number of samples, dimension of stochastic space, cross-validation error indicator w.r.t. number of iterations

plot_error_indicator(result);
mysaveas(pathname,'error_indicator','fig');
mymatlab2tikz(pathname,'error_indicator.tex');

plot_stagnation_indicator(result);
mysaveas(pathname,'stagnation_indicator','fig');
mymatlab2tikz(pathname,'stagnation_indicator.tex');

plot_error_indicator_U(result);
mysaveas(pathname,'error_indicator_U','fig');
mymatlab2tikz(pathname,'error_indicator_U.tex');

plot_stagnation_indicator_U(result);
mysaveas(pathname,'stagnation_indicator_U','fig');
mymatlab2tikz(pathname,'stagnation_indicator_U.tex');

plot_cpu_time(result,'nolegend');
mysaveas(pathname,'cpu_time','fig');
mymatlab2tikz(pathname,'cpu_time.tex');

plot_relaxation_parameter(result,'nolegend');
mysaveas(pathname,'relaxation_parameter','fig');
mymatlab2tikz(pathname,'relaxation_parameter.tex');

plot_sparsity_ratio(result);
mysaveas(pathname,'sparsity_ratio','fig');
mymatlab2tikz(pathname,'sparsity_ratio.tex');

plot_nb_samples(result);
mysaveas(pathname,'nb_samples','fig');
mymatlab2tikz(pathname,'nb_samples.tex');

plot_dim_stochastic_space(result);
mysaveas(pathname,'dim_stochastic_space','fig');
mymatlab2tikz(pathname,'dim_stochastic_space.tex');

plot_cv_error_indicator(result);
mysaveas(pathname,'cv_error_indicator','fig');
mymatlab2tikz(pathname,'cv_error_indicator.tex');

%% Display multi-index set

PC_U = getPC(U);
plot_multi_index_set(PC_U,'nolegend')
mysaveas(pathname,'multi_index_set_U','fig');
mymatlab2tikz(pathname,'multi_index_set_U.tex');

for k=1:n
    PC_w = getPC(w{k});
    plot_multi_index_set(PC_w,'nolegend')
    mysaveas(pathname,['multi_index_set_w_' num2str(k)],'fig');
    mymatlab2tikz(pathname,['multi_index_set_w_' num2str(k) '.tex']);
    
    PC_lambda = getPC(lambda{k});
    plot_multi_index_set(PC_lambda,'nolegend')
    mysaveas(pathname,['multi_index_set_lambda_' num2str(k)],'fig');
    mymatlab2tikz(pathname,['multi_index_set_lambda_' num2str(k) '.tex']);
end

%% Display evolution of multi-index set

% if isfield(result_ref,{'PC_seq_U','PC_seq_w','PC_seq_lambda'})
%     video_indices(result_ref.PC_seq_U,'filename','multi_index_set_U_ref','pathname',pathname)
%     for k=1:n
%         video_indices(result_ref.PC_seq_w{k},'filename',['multi_index_set_w_ref_' num2str(k)],'pathname',pathname)
%         video_indices(result_ref.PC_seq_lambda{k},'filename',['multi_index_set_lambda_ref_' num2str(k)],'pathname',pathname)
%     end
% end

% if isfield(result,{'PC_seq_w','PC_seq_lambda'})
%     for k=1:n
%         video_indices(result.PC_seq_w{k}{end},'filename',['multi_index_set_w_' num2str(k)],'pathname',pathname)
%         video_indices(result.PC_seq_lambda{k}{end},'filename',['multi_index_set_lambda_' num2str(k)],'pathname',pathname)
%     end
% end

%% Display evolution of cross-validation error indicator, dimension of stochastic space and number of samples w.r.t. number of iterations

% if isfield(result_ref,{'cv_error_indicator_seq_U','cv_error_indicator_seq_w','cv_error_indicator_seq_lambda','PC_seq_U','PC_seq_w','PC_seq_lambda','N_seq'})
%     plot_adaptive_algorithm(result_ref.cv_error_indicator_seq_U,result_ref.PC_seq_U,result_ref.N_seq);
%     mysaveas(pathname,'adaptive_algorithm_U_ref.fig','fig');
%     mymatlab2tikz(pathname,'adaptive_algorithm_U_ref.tex');
%     for k=1:n
%         plot_adaptive_algorithm(result_ref.cv_error_indicator_seq_w{k},result_ref.PC_seq_w{k},result_ref.N_seq);
%         mysaveas(pathname,['adaptive_algorithm_w_ref_' num2str(k) '.fig'],'fig');
%         mymatlab2tikz(pathname,['adaptive_algorithm_w_ref_' num2str(k) '.tex']);
%         plot_adaptive_algorithm(result_ref.cv_error_indicator_seq_lambda{k},result_ref.PC_seq_lambda{k},result_ref.N_seq);
%         mysaveas(pathname,['adaptive_algorithm_lambda_ref_' num2str(k) '.fig'],'fig');
%         mymatlab2tikz(pathname,['adaptive_algorithm_lambda_ref_' num2str(k) '.tex']);
%     end
% end

% if isfield(result,{'cv_error_indicator_seq_w','cv_error_indicator_seq_lambda','PC_seq_w','PC_seq_lambda','N_seq'})
%     for k=1:n
%         plot_adaptive_algorithm(result.cv_error_indicator_seq_w{k}{end},result.PC_seq_w{k}{end},result.N_seq{k}{end});
%         mysaveas(pathname,['adaptive_algorithm_w_' num2str(k) '.fig'],'fig');
%         mymatlab2tikz(pathname,['adaptive_algorithm_w_' num2str(k) '.tex']);
%         plot_adaptive_algorithm(result.cv_error_indicator_seq_lambda{k}{end},result.PC_seq_lambda{k}{end},result.N_seq{k}{end});
%         mysaveas(pathname,['adaptive_algorithm_lambda_' num2str(k) '.fig'],'fig');
%         mymatlab2tikz(pathname,['adaptive_algorithm_lambda_' num2str(k) '.tex']);
%     end
% end

%% Display statistical outputs of multiscale solution

% plot_stats_sols(glob,patches,interfaces,U,w,lambda);

plot_mean_U(glob,U);
mysaveas(pathname,'mean_U',formats,renderer);

plot_mean_sol(glob,patches,interfaces,U,w);
mysaveas(pathname,'mean_sol',formats,renderer);

plot_mean_U_w(glob,patches,interfaces,U,w);
mysaveas(pathname,'mean_U_w',formats,renderer);

plot_mean_U_w(glob,patches,interfaces,U,w,'view3');
mysaveas(pathname,'mean_U_w_surf',formats,renderer);

plot_var_U(glob,U);
mysaveas(pathname,'var_U',formats,renderer);

plot_var_sol(glob,patches,interfaces,U,w);
mysaveas(pathname,'var_sol',formats,renderer);

plot_var_U_w(glob,patches,interfaces,U,w);
mysaveas(pathname,'var_U_w',formats,renderer);

plot_var_U_w(glob,patches,interfaces,U,w,'view3');
mysaveas(pathname,'var_U_w_surf',formats,renderer);

plot_std_U(glob,U);
mysaveas(pathname,'std_U',formats,renderer);

plot_std_sol(glob,patches,interfaces,U,w);
mysaveas(pathname,'std_sol',formats,renderer);

plot_std_U_w(glob,patches,interfaces,U,w);
mysaveas(pathname,'std_U_w',formats,renderer);

plot_std_U_w(glob,patches,interfaces,U,w,'view3');
mysaveas(pathname,'std_U_w_surf',formats,renderer);

M = getM(PC);
for m=1:M
    plot_sobol_indices_sol(glob,patches,interfaces,U,w,m);
    mysaveas(pathname,['sobol_indices_sol_var_' num2str(m)],formats,renderer);
    
    plot_sensitivity_indices_max_var_sol(glob,patches,interfaces,U,w,m);
    mysaveas(pathname,['sensitivity_indices_sol_var_' num2str(m)],formats,renderer);
end

%% Display random evaluations of reference and multiscale solutions

% nbsamples = 3;
% for s=1:nbsamples
%     xi = random(RANDVARS(PC));
%     
%     if exist('U_ref','var') && exist('w_ref','var') && exist('lambda_ref','var')
%         U_ref_xi = randomeval(U_ref,xi);
%         w_ref_xi = randomeval(w_ref,xi);
%         lambda_ref_xi = randomeval(lambda_ref,xi);
%         % plot_sols_ref(glob,randomeval(patches,xi),interfaces,U_ref_xi,w_ref_xi,lambda_ref_xi);
%         plot_sol_ref(glob,randomeval(patches,xi),interfaces,U_ref_xi,w_ref_xi);
%         % plot_U_ref(glob,U_ref_xi);
%         plot_w_ref(patches,w_ref_xi);
%         % plot_lambda_ref(interfaces,lambda_ref_xi);
%     end
%     
%     U_xi = randomeval(U,xi);
%     w_xi = randomeval(w,xi);
%     lambda_xi = randomeval(lambda,xi);
%     % plot_sols(glob,randomeval(patches,xi),interfaces,U_xi,w_xi,lambda_xi);
%     plot_sol(glob,randomeval(patches,xi),interfaces,U_xi,w_xi);
%     % plot_U(glob,U_xi);
%     plot_w(patches,w_xi);
%     % plot_lambda(interfaces,lambda_xi);
% end

myparallel('stop');
