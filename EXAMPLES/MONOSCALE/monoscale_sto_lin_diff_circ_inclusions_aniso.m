%% Monoscale stochastic linear diffusion circular inclusions anisotropic %%
%%-----------------------------------------------------------------------%%
% [Beck, Nobile, Tamellini, Tempone 2011,2014]

% clc
clear all
close all

%% Input data
M = 4; % number of random variables
FileName = ['monoscale_sto_lin_diff_' num2str(M) '_circ_inclusions_aniso'];
PathName = [getfemobjectoptions('path') 'MYCODE/RESULTS/' FileName '/'];
if ~exist(PathName,'dir')
    dos(['mkdir ' PathName]);
end
set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
renderer = 'opengl';

% Parallel computing
% myparallel('start');

%% Domain and mesh definition

if exist([PathName 'gmsh_circular_inclusions.msh'],'file')
    system.S = gmsh2femobject(2,[PathName 'gmsh_circular_inclusions.msh'],2);
elseif exist([PathName 'gmsh_circular_inclusions.geo'],'file')
    system.S = gmsh2femobject(2,[PathName 'gmsh_circular_inclusions.geo'],2);
else
    D = DOMAIN(2,[0.0,0.0],[1.0,1.0]);
    r = 0.13;
    B = cell(1,9);
    B{1} = CIRCLE(0.2,0.2,r);
    B{2} = CIRCLE(0.2,0.5,r);
    B{3} = CIRCLE(0.2,0.8,r);
    B{4} = CIRCLE(0.5,0.8,r);
    B{5} = CIRCLE(0.8,0.8,r);
    B{6} = CIRCLE(0.8,0.5,r);
    B{7} = CIRCLE(0.8,0.2,r);
    B{8} = CIRCLE(0.5,0.2,r);
    B{9} = DOMAIN(2,[0.4,0.4],[0.6,0.6]);
    cl = 0.02;
    system.S = gmshdomainwithinclusion(D,B,cl,cl,[PathName 'gmsh_circular_inclusions']);
end

%% Random variables

rv = RVUNIFORM(-0.99,0);
RV = RANDVARS(repmat({rv},1,M));
[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',1);

%% Materials

% Linear diffusion coefficients K_square and K_circle
K_square = 1;
K_circle = cell(1,M);
g = [1 0.9 0.75 0.6];
for m=1:M
    % K_circle(xi) = 1 + g * xi
    K_circle{m} = ones(1,1,PC) + g(m) * X{m};
end

% Material mat_square associated to squares
% a(u,v) = int( K.grad(u).grad(v) )
mat_square = FOUR_ISOT('k',K_square); % uniform value
mat_square = setnumber(mat_square,0);
k = [0 2 4 6 8 9];
system.S = setmaterial(system.S,mat_square,k+1);

% Material mat_circle associated to circles
% a(u,v) = int( K.grad(u).grad(v) )
mat_circle = cell(1,M);
for m=1:M
    mat_circle{m} = FOUR_ISOT('k',K_circle{m}); % uniform value
    mat_circle{m} = setnumber(mat_circle{m},m);
end
k=[1 3 5 7];
for m=1:M
    system.S = setmaterial(system.S,mat_circle{m},k(m)+1);
end

%% Finalization and application of Dirichlet boundary conditions

system.S = final(system.S);
system.S = addcl(system.S,[]);

%% Stiffness matrices and sollicitation vectors

% Source term f
f = 100;

% Stiffness matrix system.A and sollicitation vector system.b associated to mesh system.S
if israndom(system.S)
    system.A = [];
else
    system.A = calc_rigi(system.S);
end

if israndom(f)
    system.b = [];
else
    system.b = bodyload(system.S,[],'QN',f);
end

%% Sampling-based approach/method: L2 Projection, Least-squares minimization/Regression, Interpolation/Collocation

initPC = POLYCHAOS(RV,0,'typebase',1);

method = METHOD('type','leastsquares','display',true,'displayiter',true,...
    'basis','adaptive','initPC',initPC,'maxcoeff',Inf,...
    'algorithm','RMS','bulkparam',0.5,...
    'sampling','adaptive','initsample',2,'addsample',0.1,'maxsample',Inf,...
    'regul','','cv','leaveout','k',10,...
    'tol',1e-2,'tolstagn',1e-1,'toloverfit',1.1,'correction',false,...
    'decompKL',false,'tolKL',1e-12,'cvKL','leaveout','kKL',10);

%% Resolution of problem

fun = @(xi) solve_system(calc_system(randomeval_system(system,xi)));

[u,result] = solve_random(method,fun);

PC = getPC(u);

%% Save all variables

save(fullfile(PathName,'all.mat'));

%% Display domain, partition and mesh

% Display partition of mesh system.S
plot_partition(system.S);
mysaveas(PathName,'mesh_partition',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(PathName,'mesh_partition',renderer);

% Display mesh system.S
plot_model(system.S,'nolegend');
mysaveas(PathName,'mesh',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(PathName,'mesh',renderer);

%% Display multi-index set

for m=1:2:M
    plot_multi_index_set(PC,'dim',[m m+1],'nolegend')
    mysaveas(PathName,['multi_index_set_dim_' num2str(m) '_' num2str(m+1)],'fig');
    mymatlab2tikz(PathName,['multi_index_set_dim_' num2str(m) '_' num2str(m+1) '.tex']);
end

%% Display evolution of multi-index set

% if isfield(result,'PC_seq')
%     for m=1:2:M
%         video_indices(result.PC_seq,'dim',[m m+1],'filename','multi_index_set','pathname',PathName)
%     end
% end

%% Display evolution of cross-validation error indicator, dimension of stochastic space and number of samples w.r.t. number of iterations

if isfield(result,{'cv_error_indicator_seq','PC_seq','N_seq'})
    plot_adaptive_algorithm(result.cv_error_indicator_seq,result.PC_seq,result.N_seq);
    mysaveas(PathName,'adaptive_algorithm.fig','fig');
    mymatlab2tikz(PathName,'adaptive_algorithm.tex');
end

%% Display statistical outputs : mean, variance, standard deviation, Sobol and other sensitivity indices

% plot_stats(system.S,u);

plot_mean(system.S,u);
mysaveas(PathName,'mean_sol',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(PathName,'mean_sol',renderer);

plot_var(system.S,u);
mysaveas(PathName,'var_sol',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(PathName,'var_sol',renderer);

plot_std(system.S,u);
mysaveas(PathName,'std_sol',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(PathName,'std_sol',renderer);

M = getM(PC);
for m=1:M
    plot_sobol_indices(system.S,u,m);
    mysaveas(PathName,['sobol_indices_sol_var_' num2str(m)],{'fig','epsc2','pdf'},renderer);
    mysaveaspdf(PathName,['sobol_indices_sol_var_' num2str(m)],renderer);
end
for m=1:M
    plot_sensitivity_indices_max_var(system.S,u,m);
    mysaveas(PathName,['sensitivity_indices_sol_var_' num2str(m)],{'fig','epsc2','pdf'},renderer);
    mysaveaspdf(PathName,['sensitivity_indices_sol_var_' num2str(m)],renderer);
end

%% Display random evaluations of solution u

% nbsamples = 3;
% for s=1:nbsamples
%     xi = random(RANDVARS(PC));
%     u_xi = randomeval(u,xi);
%     S_xi = randomeval(system.S,xi);
%     plot_solution(S_xi,u_xi);
% end

% myparallel('stop');
