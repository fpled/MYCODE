%% Monoscale stochastic linear diffusion %%
%%---------------------------------------%%

% clc
clear all
close all

%% Input data

filename = 'monoscale_sto_lin_diff';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'RESULTS',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
% set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
renderer = 'OpenGL';

% Parallel computing
myparallel('start');

%% Domains and meshes

D = DOMAIN(2,[0.0,0.0],[1.0,1.0]);

% nbelem = [20,20];
% system.S = build_model(D,'nbelem',nbelem);
cl = 0.05;
system.S = build_model(D,'cl',cl,'filename',[pathname 'gmsh_domain']);

%% Random variables

rv = RVUNIFORM(0,1);
RV = RANDVARS(rv);
[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',1);

%% Materials

% Linear diffusion coefficient K
% K(xi) = 1 + xi
K = ones(1,1,PC) + X{1};
mat = FOUR_ISOT('k',K); % uniform value
system.S = setmaterial(system.S,mat);

%% Dirichlet boundary conditions

system.S = final(system.S);
system.S = addcl(system.S,[]);

%% Stiffness matrices and sollicitation vectors

if israndom(system.S)
    system.A = [];
else
    system.A = calc_rigi(system.S);
end

% Source term
f = 100;

if israndom(f)
    system.b = [];
else
    system.b = bodyload(system.S,[],'QN',f);
end

%% Resolution

initPC = POLYCHAOS(RV,0,'typebase',1);

method = METHOD('type','leastsquares','display',true,'displayiter',true,...
    'basis','adaptive','initPC',initPC,'maxcoeff',Inf,...
    'algorithm','RMS','bulkparam',0.5,...
    'sampling','adaptive','initsample',2,'addsample',0.1,'maxsample',Inf,...
    'regul','','cv','leaveout','k',10,...
    'tol',1e-12,'tolstagn',1e-1,'toloverfit',1.1,'correction',false,...
    'decompKL',false,'tolKL',1e-12,'cvKL','leaveout','kKL',10);

fun = @(xi) solve_system(calc_system(randomeval_system(system,xi)));
[u,result] = solve_random(method,fun);
PC = getPC(u);

%% Save variables

save(fullfile(pathname,'all.mat'));

%% Display domains and meshes

plot_domain(D);
mysaveas(pathname,'domain',{'fig','epsc2'},renderer);
mymatlab2tikz(pathname,'domain.tex');

% plot_partition(system.S,'nolegend');
% mysaveas(pathname,'mesh_partition',{'fig','epsc2'},renderer);

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

%% Display statistical outputs of solution

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

%% Display random evaluations of solution

% nbsamples = 3;
% for s=1:nbsamples
%     xi = random(RANDVARS(PC));
%     u_xi = randomeval(u,xi);
%     S_xi = randomeval(system.S,xi);
%     plot_solution(S_xi,u_xi);
% end

myparallel('stop');
