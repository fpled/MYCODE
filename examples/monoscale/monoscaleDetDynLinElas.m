%% Monoscale deterministic linear elasticity dynamics problem %%
%%------------------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

filename = 'dynLinElas';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','monoscaleDet',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    D = DOMAIN(2,[0.0,0.0],[2.0,0.3]);
    
    elemtype = 'TRI3';
    % elemtype = 'QUA4';
    % option = 'DEFO'; % plane strain
    option = 'CONT'; % plane stress
    % nbelem = [20,20];
    % pb.S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'option',option);
    cl = 0.03;
    pb.S = build_model(D,'cl',cl,'elemtype',elemtype,'option',option,'filename',[pathname 'gmsh_domain']);
    
    %% Materials
    % Poisson ratio
    NU = 0;
    % Thickness
    DIM3 = 1;
    % Density
    RHO = 1;
    % Young modulus
    E = 1;
    % Material
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
    mat = setnumber(mat,1);
    pb.S = setmaterial(pb.S,mat);
    
    %% Dirichlet boundary conditions
    L1 = LIGNE(getvertex(D,1),getvertex(D,4));
    L2 = LIGNE(getvertex(D,2),getvertex(D,3));
    
    pb.S = final(pb.S);
    pb.S = addcl(pb.S,L1);
    
    %% Initial conditions
    pb.u0 = zeros(getnbddlfree(pb.S),1);
    pb.v0 = zeros(getnbddlfree(pb.S),1);
    
    %% Time scheme
    t0 = 0;
    t1 = 2;
    nt = 50;
    T = TIMEMODEL(t0,t1,nt);
    
    % pb.N = NEWMARKSOLVER(T,'alpha',0.05,'display',false);
    pb.N = DGTIMESOLVER(T,1,'outputsplit',true,'display',false,'lu',true);
    
    %% Mass, stiffness and damping matrices and sollicitation vectors
    pb.M = calc_mass(pb.S);
    pb.K = calc_rigi(pb.S);
    f = surfload(pb.S,L2,'FX',-1);
    
    tc = get(T,'t1')/6;
    loadfun = @(N) rampe(N,t0,tc);
    % loadfun = @(N) dirac(N,t0,tc);
    % loadfun = @(N) one(N);
    pb.f = f*loadfun(pb.N);
    
    save(fullfile(pathname,'problem.mat'),'pb','D');
else
    load(fullfile(pathname,'problem.mat'),'pb','D');
end

%% Solution
if solveProblem
    t = tic;
    [ut,result,vt] = ddsolve(pb.N,pb.f,pb.M,pb.K,[],pb.u0,pb.v0);
    time = toc(t);
    
    et = calc_epsilon(pb.S,ut);
    st = calc_sigma(pb.S,ut);
    
    save(fullfile(pathname,'solution.mat'),'ut','result','vt','et','st','loadfun','time');
else
    load(fullfile(pathname,'solution.mat'),'ut','result','vt','et','st','loadfun','time');
end

%% Outputs
fprintf('\n');
fprintf(['load function : ' func2str(loadfun) '\n']);
fprintf(['spatial mesh  : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(pb.S));
fprintf('nb nodes    = %g\n',getnbnode(pb.S));
fprintf('nb dofs     = %g\n',getnbddl(pb.S));
fprintf('time solver : %s\n',class(pb.N));
fprintf('nb time steps = %g\n',getnt(pb.N));
fprintf('nb time dofs  = %g\n',getnbtimedof(pb.N));
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

%% Display
if displaySolution
    %% Display domains and meshes
    plotDomain(D,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    plotModel(pb.S,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    %% Display evolution of solution
    for i=1:2
        evolSolution(pb.S,ut,'displ',i,'filename',['evol_u_' num2str(i)],'pathname',pathname);
        evolSolution(pb.S,vt,'displ',i,'filename',['evol_v_' num2str(i)],'pathname',pathname);
    end
    
    for i=1:3
        evolSolution(pb.S,ut,'epsilon',i,'filename',['evol_eps_' num2str(i)],'pathname',pathname);
        evolSolution(pb.S,ut,'sigma',i,'filename',['evol_sig_' num2str(i)],'pathname',pathname);
    end
    
    evolSolution(pb.S,ut,'epsilon','mises','filename','evol_eps_von_mises','pathname',pathname);
    evolSolution(pb.S,ut,'sigma','mises','filename','evol_sig_von_mises','pathname',pathname);
end

% myparallel('stop');
