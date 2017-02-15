%% Monoscale deterministic transient linear advection-diffusion-reaction problem %%
%%-------------------------------------------------------------------------------%%

% clc
% clear all
close all
% set(0,'DefaultFigureVisible','off');
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

filename = 'transientLinAdvDiffReac';
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
    cl1 = 0.02;
    cl2 = 0.04;
    cl0 = 0.02;
    cltip = 0.01;
    pb.S = gmshcanister(cl1,cl2,cl0,cltip,[pathname 'gmsh_canister']);
    
    %% Materials
    % Linear diffusion coefficient
    K = 1;
    % Thermal capacity
    c = 1;
    % Advection velocity
    Sc = pb.S;
    mat = FOUR_ISOT('k',1);
    Sc = setmaterial(Sc,mat);
    P1 = POINT(getnode(getridge(Sc,1)));
    P5 = POINT(getnode(getridge(Sc,5)));
    P6 = POINT(getnode(getridge(Sc,6)));
    P15 = POINT(getnode(getridge(Sc,15)));
    P16 = POINT(getnode(getridge(Sc,16)));
    L1 = LIGNE(P5,P6);
    L2 = LIGNE(P15,P16);
    Sc = final(Sc);
    Sc = addcl(Sc,P1,'T',0);
    A = calc_rigi(Sc);
    b1 = surfload(Sc,L1,'QN',-1);
    b2 = surfload(Sc,L2,'QN',1);
    b = b1+b2;
    phi = A\b;
    v = FENODEFIELD(calc_sigma(Sc,phi,'node'));
    V = getvalue(v);
    V = {{FENODEFIELD(V(:,1)),FENODEFIELD(V(:,2))}};
    % Linear reaction parameter
    R1 = 0.1;
    R2 = 10;
    
    % Materials
    mat = MATERIALS();
    mat{1} = FOUR_ISOT('k',K,'c',c,'b',V,'r',R1);
    mat{2} = FOUR_ISOT('k',K,'c',c,'b',V,'r',R2);
    pb.S = setmaterial(pb.S,mat{1},1);
    pb.S = setmaterial(pb.S,mat{2},2:3);
    
    %% Dirichlet boundary conditions
    L1 = LIGNE(P5,P6);
    L2 = LIGNE(P15,P16);
    
    pb.S = final(pb.S);
    pb.S = addcl(pb.S,L1,'T',0);
    pb.S = addcl(pb.S,L2,'T',1);
    
    %% Initial conditions
    [~,numnode1] = intersect(pb.S,L1);
    pb.u0 = zeros(pb.S.nbnode,1);
    pb.u0(numnode1) = 1;
    
    %% Time scheme
    t0 = 0;
    t1 = 2;
    nt = 100;
    T = TIMEMODEL(t0,t1,nt);
    
    % pb.N = EULERTIMESOLVER(T,'eulertype','explicit');
    % pb.N = EULERTIMESOLVER(T,'eulertype','implicit');
    pb.N = DGTIMESOLVER(T,1);
    pb.N = setparam(pb.N,'display',true);
    
    %% Mass and stifness matrices and sollicitation vectors
    pb.M = calc_mass(pb.S);
    pb.K = calc_rigi(pb.S,'nofree');
    
    
    
    pb.f0 = -pb.K*pb.u0;
    pb.K = freematrix(pb.S,pb.K);
    pb.f0 = freevector(pb.S,pb.f0);
    
    pb.f = pb.f_init*one(pb.N);
    
    save(fullfile(pathname,'problem.mat'),'pb','D','v');
else
    load(fullfile(pathname,'problem.mat'),'pb','D','v');
end

%% Newmark time scheme
if solveProblem
    % Static solution at initial
    u = pb.K\pb.f;
    u = unfreevector(pb.S,u) + pb.u_init;
    
    % Dynamic solution
    t = tic;
    u0 = zeros(size(pb.M,1),1);
    [ut,result,vt] = dsolve(pb.N,pb.b,pb.M,pb.K,u0);
    time = toc(t);
    
    save(fullfile(pathname,'solution.mat'),'u','ut','result','vt','time');
else
    load(fullfile(pathname,'solution.mat'),'u','ut','result','vt','time');
end

%% Outputs
fprintf('\n');
fprintf('nb elements = %g\n',getnbelem(pb.S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(pb.S));
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

foutput = bodyload(keepgroupelem(pb.S,2),[],'QN',1,'free');

%% Display
if displaySolution
    %% Display domains and meshes
    figure('Name','Domain')
    clf;
    plot(create_boundary(pb.S));
    h1 = plot(pb.S,'selgroup',1,'FaceColor',getfacecolor(1),'EdgeColor','none');
    h2 = plot(pb.S,'selgroup',2,'FaceColor',getfacecolor(2),'EdgeColor','none');
    h3 = plot(get(pb.S,5),'selgroup',3,'FaceColor',getfacecolor(3),'EdgeColor','none');
    h4 = plot(pb.S,'selgroup',3,'FaceColor',getfacecolor(3),'EdgeColor','none');
    h5 = plot(pb.S,'selgroup',3,'FaceColor',getfacecolor(3),'EdgeColor','none');
    set(gca,'FontSize',16)
    l = legend([h1(1),h2(1),h3(1)],'$\Omega_1$','$\Omega_2$','$\Omega_0$');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    figure('Name','Advection velocity')
    clf
    plot(pb.phi,pb.S);
    colorbar
    set(gca,'FontSize',16)
    hold on
    ampl = 6;
    quiver(v,pb.S,ampl,'k');
    mysaveas(pathname,'advection',formats,renderer);

    plotModel(pb.S,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    %% Display solution
    plotSolution(pb.S,u,'surface');
    mysaveas(pathname,'solution',formats,renderer);
    
    %% Display evolution of solution
    utmin = min(min(ut(1)));
    utmax = max(max(ut(1)));
    pb.N = setevolparam(pb.N,'plotstep',1,'setaxis',false,'setcaxis',true,...
        'caxis',[stmin,stmax],'pausetime',1/get(T,'nt'),'colorbar',true);
    frame = evol(pb.N,ut,pb.S);
    mov = VideoWriter(fullfile(pathname,'evol_solution_time'));%,'Uncompressed AVI');
    mov.FrameRate = 30;
    mov.Quality = 100;
    open(mov);
    writeVideo(mov,frame);
    close(mov);
    
    % pb.N = setevolparam(pb.N,'step',3,'view',2,'setcaxis',true,'caxis',[smmin,smmax],...
    %     'pausetime',1/nt,'setaxis',false,'colormap',jet);
    % evol(pb.N,utx,pb.S)
    evol(ut,pb.S)
    
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
