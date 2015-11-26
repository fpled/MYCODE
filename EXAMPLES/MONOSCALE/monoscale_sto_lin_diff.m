%% Monoscale stochastic linear diffusion %%
%%---------------------------------------%%

% clc
clear all
close all

%% Input data
filename = 'monoscale_sto_lin_diff';
pathname = [getfemobjectoptions('path') 'MYCODE/RESULTS/' filename '/'];
if ~exist(pathname,'dir')
    mkdir(pathname);
end
set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
renderer = 'OpenGL';

% Parallel computing
myparallel('start');

%% Domain and mesh definition

D = DOMAIN(2,[0.0,0.0],[5.0,5.0]);
nbelem = [20,20];
system.S = build_model(D,'nbelem',nbelem);
% cl = 0.25;
% system.S = build_model(D,'cl',cl);

%% Random variables

rv = RVUNIFORM(1,2);
RV = RANDVARS(rv);
[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',1);
K = X{1};

%% Materials

% a(u,v) = int( K.grad(u).grad(v) )
mat = FOUR_ISOT('k',K); % uniform value
mat = setnumber(mat,1);
system.S = setmaterial(system.S,mat);

%% Finalization and application of Dirichlet boundary conditions

system.S = final(system.S);
system.S = addcl(system.S,[]);

%% Stiffness matrices and sollicitation vectors

% Source term f
f = 1;

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
    'tol',1e-12,'tolstagn',1e-1,'toloverfit',1.1,'correction',false,...
    'decompKL',false,'tolKL',1e-12,'cvKL','leaveout','kKL',10);

%% Resolution of problem

fun = @(xi) solve_system(calc_system(randomeval_system(system,xi)));

[u,result] = solve_random(method,fun);

PC = getPC(u);

%% Save all variables

save(fullfile(pathname,'all.mat'));

%% Display domain, partition and mesh

% Display domain
plot_domain(D);
mysaveas(pathname,'domain',{'fig','epsc2'},renderer);
mymatlab2tikz(pathname,'domain.tex');

% Display partition of mesh system.S
% plot_partition(system.S);
% mysaveas(pathname,'mesh_partition',{'fig','epsc2'},renderer);

% Display mesh system.S
plot_model(system.S,'nolegend');
mysaveas(pathname,'mesh',{'fig','epsc2'},renderer);

%% Display multi-index set

plot_multi_index_set(PC,'nolegend')
mysaveas(pathname,'multi_index_set','fig');
mymatlab2tikz(pathname,'multi_index_set.tex');

%% Display evolution of multi-index set

% if isfield(result,'PC_seq')
%     video_indices(result.PC_seq,'filename','multi_index_set','pathname',pathname)
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
