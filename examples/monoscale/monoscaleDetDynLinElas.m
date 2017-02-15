%% Monoscale deterministic linear elasticity dynamics problem %%
%%------------------------------------------------------------%%

% clc
% clear all
close all
% set(0,'DefaultFigureVisible','off');
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

filename = 'dynLinElas';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,...
    'results',filesep,'monoscaleDet',filesep,filename,filesep);
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
    
    % pb.N = NEWMARKSOLVER(T,'alpha',0.05);
    pb.N = DGTIMESOLVER(T,1);
    pb.N = setparam(pb.N,'display',true);
    
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

%% Newmark time scheme
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
fprintf(['mesh          : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(pb.S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(pb.S));
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
    
    %% Display evolution of stress field
    mov = VideoWriter(fullfile(pathname,['evol_sigma_' num2str(i)]));%,'Uncompressed AVI');
    mov.FrameRate = 30;
    mov.Quality = 100;
    i=1;
    stmin = min(min(st(i)));
    stmax = max(max(st(i)));
    ut = setevolparam(ut,'plotstep',10)
    pb.N = setevolparam(pb.N,'plotstep',1,'setaxis',false,'setcaxis',true,...
        'caxis',[stmin,stmax],'pausetime',1/get(gettimemodel(pb.N),'nt'),'colorbar',true);
    frame = evol(pb.N,st,pb.S,'compo','SXX');
    
    open(mov);
    writeVideo(mov,frame);
    close(mov);
    
    %% Display solution
    % ampl = 0;
%     ampl = getsize(S)/max(abs(u))/5;
%     
%     for i=1:2
%         plotSolution(S,ut,'displ',i,'ampl',ampl);
%         mysaveas(pathname,['ut_' num2str(i)],formats,renderer);
%     end
%     
%     for i=1:3
%         plotSolution(S,ut,'epsilon',i,'ampl',ampl);
%         mysaveas(pathname,['epst_' num2str(i)],formats,renderer);
%         
%         plotSolution(S,ut,'sigma',i,'ampl',ampl);
%         mysaveas(pathname,['sigt_' num2str(i)],formats,renderer);
%     end
end

% myparallel('stop');
