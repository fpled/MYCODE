%% Monoscale stochastic linear diffusion circular inclusions anisotropic %%
%%-----------------------------------------------------------------------%%
% [Beck, Nobile, Tamellini, Tempone 2011,2014]

% clc
clear all
close all

%% Input data
M = 4; % number of random variables
filename = ['monoscale_sto_lin_diff_' num2str(M) '_circ_inclusions_aniso'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'RESULTS',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
renderer = 'OpenGL';

% Parallel computing
myparallel('start');

%% Domain and mesh definition

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
system.S = gmshdomainwithinclusion(D,B,cl,cl,[pathname 'gmsh_circular_inclusions']);

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

save(fullfile(pathname,'all.mat'));

%% Display domain, partition and mesh

% Display partition of mesh system.S
plot_partition(system.S);
mysaveas(pathname,'mesh_partition',{'fig','epsc2'},renderer);

% Display mesh system.S
plot_model(system.S,'nolegend');
mysaveas(pathname,'mesh',{'fig','epsc2'},renderer);

%% Display multi-index set

for m=1:2:M
    plot_multi_index_set(PC,'dim',[m m+1],'nolegend')
    mysaveas(pathname,['multi_index_set_dim_' num2str(m) '_' num2str(m+1)],'fig');
    mymatlab2tikz(pathname,['multi_index_set_dim_' num2str(m) '_' num2str(m+1) '.tex']);
end

%% Display evolution of multi-index set

% if isfield(result,'PC_seq')
%     for m=1:2:M
%         video_indices(result.PC_seq,'dim',[m m+1],'filename','multi_index_set','pathname',pathname)
%     end
% end

%% Display evolution of cross-validation error indicator, dimension of stochastic space and number of samples w.r.t. number of iterations

% if isfield(result,{'cv_error_indicator_seq','PC_seq','N_seq'})
%     plot_adaptive_algorithm(result.cv_error_indicator_seq,result.PC_seq,result.N_seq);
%     mysaveas(pathname,'adaptive_algorithm.fig','fig');
%     mymatlab2tikz(pathname,'adaptive_algorithm.tex');
% end

%% Display statistical outputs : mean, variance, standard deviation, Sobol and other sensitivity indices

% plot_stats(system.S,u);

plot_mean(system.S,u);
mysaveas(pathname,'mean_sol',{'fig','epsc2'},renderer);

plot_var(system.S,u);
mysaveas(pathname,'var_sol',{'fig','epsc2'},renderer);

plot_std(system.S,u);
mysaveas(pathname,'std_sol',{'fig','epsc2'},renderer);

M = getM(PC);
for m=1:M
    plot_sobol_indices(system.S,u,m);
    mysaveas(pathname,['sobol_indices_sol_var_' num2str(m)],{'fig','epsc2'},renderer);
    
    plot_sensitivity_indices_max_var(system.S,u,m);
    mysaveas(pathname,['sensitivity_indices_sol_var_' num2str(m)],{'fig','epsc2'},renderer);
end

%% Display random evaluations of solution u

% nbsamples = 3;
% for s=1:nbsamples
%     xi = random(RANDVARS(PC));
%     u_xi = randomeval(u,xi);
%     S_xi = randomeval(system.S,xi);
%     plot_solution(S_xi,u_xi);
% end

myparallel('stop');
