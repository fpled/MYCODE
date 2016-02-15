%% Multiscale stochastic linear elasticity %%
%%-----------------------------------------%%

% clc
clear all
close all

%% Input data

n = 4; % number of patches n = 1, 2, 4
filename = ['multiscale_sto_lin_elas_' num2str(n) '_patches'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE','RESULTS',filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
renderer = 'OpenGL';

% Reference solution - Direct resolution of initial problem based on non-overlapping domain decomposition
solve_reference = true;
calc_MC_error_estimate_ref = false;

% Multiscale solution - Reformulated global-local iterative algorithm based on overlapping domain decomposition
solve_multiscale = true;
calc_MC_error_estimate = false;

% Parallel computing
myparallel('start');

%% Domain and mesh definition

% Global
glob = GLOBAL();
glob_out = GLOBALOUT();

D = DOMAIN(2,[0.0,0.0],[1.0,1.0]);

nbelem = [20,20];
glob.S = build_model(D,'nbelem',nbelem);
% cl = 0.05;
% system.S = build_model(D,'cl',cl,'filename','gmsh_domain');

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

% Partition of global mesh glob.S
glob = partition(glob,patches);

%% Random variables

rv = RVUNIFORM(0,1);
RV = RANDVARS(repmat({rv},1,n));
[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',1);

%% Sampling-based method: L2 Projection, Least-squares minimization/Regression, Interpolation/Collocation

initPC = POLYCHAOS(RV,0,'typebase',1);

% For direct resolution of initial problem
method_ref = METHOD('type','leastsquares','display',true,'displayiter',true,...
    'basis','adaptive','initPC',initPC,'maxcoeff',Inf,...
    'algorithm','RMS','bulkparam',0.5,...
    'sampling','adaptive','initsample',2,'addsample',0.1,'maxsample',Inf,...
    'regul','','cv','leaveout','k',10,...
    'tol',1e-8,'tolstagn',1e-1,'toloverfit',1.1,'correction',false,...
    'decompKL',false,'tolKL',1e-12,'cvKL','leaveout','kKL',10);

% For reformulated global-local iterative algorithm
method = METHOD('type','leastsquares','display',true,'displayiter',false,...
    'basis','adaptive','initPC',initPC,'maxcoeff',Inf,...
    'algorithm','RMS','bulkparam',0.5,...
    'sampling','adaptive','initsample',2,'addsample',0.1,'maxsample',Inf,...
    'regul','','cv','leaveout','k',10,...
    'tol',1e-4,'tolstagn',1e-1,'toloverfit',1.1,'correction',false,...
    'decompKL',false,'tolKL',1e-12,'cvKL','leaveout','kKL',10);

%% Materials associated to initial problem

% Poisson ratio
NU = 0.3;
% Thickness
DIM3 = 1;
% Density
RHO = 1;
% Young modulus E_out, E_patch and E_in
E_out = 1;
E_patch = cell(1,n);
E_in = cell(1,n);
for k=1:n
    patch = patches.PATCH{k};
    % E_patch(x,xi) = 1 + beta_patch * f(x) * xi
    % E_in(x)       = 1 + beta_in * f(x)
    % with f(x) = alpha*exp( -Amp*||x-c||_2^2/L^2 ) if ||x-c||_Inf < L
    %           = 0                                 if ||x-c||_Inf >= L
    % alpha = 10;
    % Amp = 2;
    % L = norm(getsize(D_patch{k}),Inf)/4;
    % c = getcenter(D_patch{k});
    % f = @(x) (distance(x,c,Inf)<L) * alpha * exp(-Amp*distance(x,c,2).^2/L^2);
    % 
    % beta_patch = 1;
    % E_patch{k} = ones(patch.S.nbnode,1,PC) + beta_patch * double(squeeze(f(patch.S.node))) * X{k};
    % 
    % beta_in = 0;
    % E_in{k} = 1 + beta_in * squeeze(f(glob.S.node));
    
    % E_patch(x,xi) = 1 + f(x) * xi
    % E_in(x)       = 1
    % with f(x) = 1 if ||x-c||_Inf < L
    %           = 0 if ||x-c||_Inf >= L
    L = norm(getsize(D_patch{k}),Inf)/4;
    c = getcenter(D_patch{k});
    f = @(x) distance(x,c,Inf)<L;
    E_patch{k} = ones(patch.S.nbnode,1,PC) + double(squeeze(f(patch.S.node))) * X{k};
    E_in{k} = 1;
end

% Material mat_out associated to outside subdomain
% a(u,v) = int( epsilon(u) : K : epsilon(v) )
mat_out = ELAS_ISOT('E',E_out,'NU',NU,'RHO',RHO,'DIM3',DIM3);
mat_out = setnumber(mat_out,0);
glob.S = setmaterial(glob.S,mat_out,getnumgroupelemwithparam(glob.S,'partition',0));

% Material mat_patch associated to patch
% a(u,v) = int( epsilon(u) : K : epsilon(v) )
for k=1:n
    % mat_patch = ELAS_ISOT('E',E_patch{k},'NU',NU,'RHO',RHO,'DIM3',DIM3); % uniform value
    mat_patch = ELAS_ISOT('E',FENODEFIELD(E_patch{k}),'NU',NU,'RHO',RHO,'DIM3',DIM3); % nodal values
    mat_patch = setnumber(mat_patch,k);
    patches.PATCH{k}.S = setmaterial(patches.PATCH{k}.S,mat_patch);
end

% Material mat_in associated to fictitious patch
% a(u,v) = int( epsilon(u) : K : epsilon(v) )
for k=1:n
    mat_in = ELAS_ISOT('E',E_in{k},'NU',NU,'RHO',RHO,'DIM3',DIM3); % uniform value
    % mat_in = ELAS_ISOT('E',FENODEFIELD(E_in{k}),'NU',NU,'RHO',RHO,'DIM3',DIM3); % nodal values
    mat_in = setnumber(mat_in,k);
    glob.S = setmaterial(glob.S,mat_in,getnumgroupelemwithparam(glob.S,'partition',k));
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

% Body force field f
f = [0;-100];

% Stiffness matrix glob_out.A_out and sollicitation vector glob_out.b_out associated to mesh glob_out.S_out
glob_out.A_out = calc_rigi(glob_out.S_out);
glob_out.b_out = bodyload(glob_out.S_out,[],{'FX','FY'},f);

% Stiffness matrices glob.A and glob.A_in and sollicitation vector glob.b_out associated to mesh glob.S
glob.A = calc_rigi(glob.S);
for k=1:n
    glob.A_in{k} = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
end
glob.b_out = bodyload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),[],{'FX','FY'},f);

% Stiffness matrix patch.A and sollicitation vector patch.b associated to mesh patch.S
for k=1:n
    if ~israndom(patches.PATCH{k}.S)
        patches.PATCH{k}.A = calc_rigi(patches.PATCH{k}.S);
    end
    if ~israndom(f)
        patches.PATCH{k}.b = bodyload(patches.PATCH{k}.S,[],{'FX','FY'},f);
    end
end

%% Mass matrices

% Mass matrix interface.M associated to boundary mesh interface.S
for k=1:n
    interfaces.INTERFACE{k}.M = calc_massgeom(interfaces.INTERFACE{k}.S);
end

%% Projection operators

% Projection operator glob.P_out from mesh glob.S to mesh glob.S_out
glob.P_out = calc_P_free(glob.S,glob.S_out);

% Projection operator interface.P_glob from mesh glob.S to boundary mesh interface.S
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
    [U_ref,w_ref,lambda_ref,result_ref] = solve_random(R,glob_out,patches,interfaces,method_ref);
    save(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','result_ref');
else
    load(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','result_ref');
end

%% Monte Carlo error estimation of reference solution u_ref=(U_ref,w_ref)

if calc_MC_error_estimate_ref && exist('U_ref','var') && exist('w_ref','var') && exist('lambda_ref','var')
    nbsamples = 100;
    [mean_error_U_ref,mean_error_w_ref,mean_error_lambda_ref,...
        var_error_U_ref,var_error_w_ref,var_error_lambda_ref,...
        std_error_U_ref,std_error_w_ref,std_error_lambda_ref]...
        = calc_Monte_Carlo_error_estimates(R,glob_out,patches,interfaces,U_ref,w_ref,lambda_ref,nbsamples);
end

%% Multiscale resolution using global-local iterative algorithm based on overlapping domain decomposition

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

%% Monte Carlo error estimation of multiscale solution u=(U,w) at final iteration

if calc_MC_error_estimate
    nbsamples = 100;
    [mean_error_U,mean_error_w,mean_error_lambda,...
        var_error_U,var_error_w,var_error_lambda,...
        std_error_U,std_error_w,std_error_lambda]...
        = calc_Monte_Carlo_error_estimates(I,glob,patches,interfaces,U,w,lambda,nbsamples);
end

%% Save all variables

save(fullfile(pathname,'all.mat'));

%% Display domain, partition and mesh

% Display global domain and patches
plot_domain(D,D_patch);
mysaveas(pathname,'domain_global_patches',{'fig','epsc2'},renderer);
mymatlab2tikz(pathname,'domain_global_patches.tex');

% Display partition of global mesh glob.S
% plot_partition(glob,'nolegend');
% mysaveas(pathname,'mesh_partition',{'fig','epsc2'},renderer);

% Display global mesh glob.S_out and local meshes patch.S
plot_model(glob,patches,'nolegend');
mysaveas(pathname,'mesh_global_patches',{'fig','epsc2'},renderer);

% Display all parts of global mesh glob.S
% plot_model(glob);

% Display local meshes patch.S
% plot_model(patches);

% Display boundary meshes interface.S
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

%% Display statistical outputs : mean, variance, standard deviation, Sobol and other sensitivity indices

plot_mean_U(glob,U);
mysaveas(pathname,'mean_U',{'fig','epsc2'},renderer);

plot_mean_sol(glob,patches,interfaces,U,w);
mysaveas(pathname,'mean_sol',{'fig','epsc2'},renderer);

plot_var_U(glob,U);
mysaveas(pathname,'var_U',{'fig','epsc2'},renderer);

plot_var_sol(glob,patches,interfaces,U,w);
mysaveas(pathname,'var_sol',{'fig','epsc2'},renderer);

for i=1:2
    % plot_stats_sols(glob,patches,interfaces,U,w,lambda,'displ',i);
    
    plot_mean_U(glob,U,'displ',i);
    mysaveas(pathname,['mean_U_' num2str(i)],{'fig','epsc2'},renderer);
    
    plot_mean_sol(glob,patches,interfaces,U,w,'displ',i);
    mysaveas(pathname,['mean_sol_' num2str(i)],{'fig','epsc2'},renderer);
    
    plot_mean_U_w(glob,patches,interfaces,U,w,'displ',i);
    mysaveas(pathname,['mean_U_w' num2str(i)],{'fig','epsc2'},renderer);
    
    plot_mean_U_w(glob,patches,interfaces,U,w,'displ',i,'view3');
    mysaveas(pathname,['mean_U_w_surf' num2str(i)],{'fig','epsc2'},renderer);
    
    plot_var_U(glob,U,'displ',i);
    mysaveas(pathname,['var_U_' num2str(i)],{'fig','epsc2'},renderer);
    
    plot_var_sol(glob,patches,interfaces,U,w,'displ',i);
    mysaveas(pathname,['var_sol_' num2str(i)],{'fig','epsc2'},renderer);
    
    plot_var_U_w(glob,patches,interfaces,U,w,'displ',i);
    mysaveas(pathname,['var_U_w' num2str(i)],{'fig','epsc2'},renderer);
    
    plot_var_U_w(glob,patches,interfaces,U,w,'displ',i,'view3');
    mysaveas(pathname,['var_U_w_surf' num2str(i)],{'fig','epsc2'},renderer);

    plot_std_U(glob,U,'displ',i);
    mysaveas(pathname,['std_U_' num2str(i)],{'fig','epsc2'},renderer);
    
    plot_std_sol(glob,patches,interfaces,U,w,'displ',i);
    mysaveas(pathname,['std_sol_' num2str(i)],{'fig','epsc2'},renderer);
    
    plot_std_U_w(glob,patches,interfaces,U,w,'displ',i);
    mysaveas(pathname,['std_U_w' num2str(i)],{'fig','epsc2'},renderer);
    
    plot_std_U_w(glob,patches,interfaces,U,w,'displ',i,'view3');
    mysaveas(pathname,['std_U_w_surf' num2str(i)],{'fig','epsc2'},renderer);
    
    M = getM(PC);
    for m=1:M
        plot_sobol_indices_sol(glob,patches,interfaces,U,w,m,'displ',i);
        mysaveas(pathname,['sobol_indices_sol_' num2str(i) '_var_' num2str(m)],{'fig','epsc2'},renderer);
        
        plot_sensitivity_indices_max_var_sol(glob,patches,interfaces,U,w,m,'displ',i);
        mysaveas(pathname,['sensitivity_indices_sol_' num2str(i) '_var_' num2str(m)],{'fig','epsc2'},renderer);
    end
end

%% Display relative error in statistical outputs : mean, variance, standard deviation

% if exist('U_ref','var') && exist('w_ref','var') && exist('lambda_ref','var')
%     for i=1:2
%         % plot_error_stats_sols(glob,patches,interfaces,U,w,lambda,U_ref,w_ref,lambda_ref,'displ',i);
%         
%         plot_error_mean_sol(glob,patches,interfaces,U,w,U_ref,w_ref,'displ',i);
%         mysaveas(pathname,['error_mean_sol_' num2str(i)],{'fig','epsc2'},renderer);
%         
%         plot_error_var_sol(glob,patches,interfaces,U,w,U_ref,w_ref,'displ',i);
%         mysaveas(pathname,['error_var_sol_' num2str(i)],{'fig','epsc2'},renderer);
%         
%         plot_error_std_sol(glob,patches,interfaces,U,w,U_ref,w_ref,'displ',i);
%         mysaveas(pathname,['error_std_sol_' num2str(i)],{'fig','epsc2'},renderer);
%     end
% end

%% Display random evaluations of reference solution u_ref=(U_ref,w_ref) and multiscale solution u=(U,w) at final iteration

% nbsamples = 3;
% for s=1:nbsamples
%     xi = random(RANDVARS(PC));
%     
%     if exist('U_ref','var') && exist('w_ref','var') && exist('lambda_ref','var')
%         U_ref_xi = randomeval(U_ref,xi);
%         w_ref_xi = randomeval(w_ref,xi);
%         lambda_ref_xi = randomeval(lambda_ref,xi);
%         for i=1:2
%             % plot_sols_ref(glob,randomeval(patches,xi),interfaces,U_ref_xi,w_ref_xi,lambda_ref_xi,'displ',i);
%             plot_sol_ref(glob,randomeval(patches,xi),interfaces,U_ref_xi,w_ref_xi,'displ',i);
%             % plot_U_ref(glob,U_ref_xi,'displ',i);
%             plot_w_ref(patches,w_ref_xi,'displ',i);
%             % plot_lambda_ref(interfaces,lambda_ref_xi,'displ',i);
%         end
%     end
%     
%     U_xi = randomeval(U,xi);
%     w_xi = randomeval(w,xi);
%     lambda_xi = randomeval(lambda,xi);
%     for i=1:2
%         % plot_sols(glob,randomeval(patches,xi),interfaces,U_xi,w_xi,lambda_xi,'displ',i);
%         plot_sol(glob,randomeval(patches,xi),interfaces,U_xi,w_xi,'displ',i);
%         % plot_U(glob,U_xi,'displ',i);
%         plot_w(patches,w_xi,'displ',i);
%         % plot_lambda(interfaces,lambda_xi,'displ',i);
%     end
% end

myparallel('stop');
