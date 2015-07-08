%% Monoscale stochastic nonlinear diffusion reaction %%
%%---------------------------------------------------%%

% clc
clear all
close all

%% Input data
FileName = 'monoscale_sto_nonlin_diff_reac';
PathName = [getfemobjectoptions('path') 'MYCODE/RESULTS/' FileName '/'];
if ~exist(PathName,'dir')
    dos(['mkdir ' PathName]);
end
set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
renderer = 'opengl';

% Parallel computing
% myparallel('start');

%% Domain and mesh definition

D = DOMAIN(2,[0.0,0.0],[5.0,5.0]);
nbelem = [20,20];
system.S = build_model(D,'nbelem',nbelem);
% cl = 0.25;
% system.S = build_model(D,'cl',cl);

%% Random variables

rv1 = RVUNIFORM(1,2);
rv2 = RVUNIFORM(0,1);
RV = RANDVARS(rv1,rv2);
[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',1);
K = X{1};
R = X{2};
% K2 = X{2};

%% Materials

% a(u,v) = int( (K+K2.u^2).grad(u).grad(v) + R.u^3.v )
mat = FOUR_ISOT('k',K,'r',R); % uniform value
% mat = FOUR_ISOT('k',K,'k2',K2); % uniform value
mat = setnumber(mat,1);
system.S = setmaterial(system.S,mat);

%% Finalization and application of Dirichlet boundary conditions

system.S = final(system.S);
system.S = addcl(system.S,[]);

%% Stiffness matrices and sollicitation vectors

% Source term f
f = 1;

% Stiffness operator system.A, tangent stiffness operator system.Atang and sollicitation vector system.b associated to mesh system.S
if israndom(system.S)
    system.A = [];
    system.Atang = [];
else
    system.A = @(u) calc_fint(system.S,u);
    system.Atang = @(u) calc_rigitang(system.S,u);
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
    'tol',1e-6,'tolstagn',1e-1,'toloverfit',1.1,'correction',false,...
    'decompKL',false,'tolKL',1e-12,'cvKL','leaveout','kKL',10);

%% Resolution of problem

system.solver = NEWTONSOLVER('type','full','increment',true,...
    'maxiter',100,'tol',1e-12,'tolreact',1e-1,'display',false,'stopini',true);

% u0 = solve_system(calc_system(randomeval(system,mean(RANDVARS(PC)))));
% fun = @(xi) solve_system(calc_system(randomeval(system,xi)),'inittype',u0);
fun = @(xi) solve_system(calc_system(randomeval_system(system,xi)));

[u,result] = solve_random(method,fun);

PC = getPC(u);

%% Save all variables

save(fullfile(PathName,'all.mat'));

%% Display domain, partition and mesh

% Display domain
plot_domain(D);
mysaveas(PathName,'domain',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(PathName,'domain',renderer);
mymatlab2tikz(PathName,'domain.tex');

% Display partition of mesh system.S
% plot_partition(system.S);
% mysaveas(PathName,'mesh_partition',{'fig','epsc2','pdf'},renderer);
% mysaveaspdf(PathName,'mesh_partition',renderer);

% Display mesh system.S
plot_model(system.S,'nolegend');
mysaveas(PathName,'mesh',{'fig','epsc2','pdf'},renderer);
mysaveaspdf(PathName,'mesh',renderer);

%% Display multi-index set

plot_multi_index_set(PC,'nolegend')
mysaveas(PathName,'multi_index_set','fig');
mymatlab2tikz(PathName,'multi_index_set.tex');

%% Display evolution of multi-index set

% if isfield(result,'PC_seq')
%     video_indices(result.PC_seq,'filename','multi_index_set','pathname',PathName)
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
