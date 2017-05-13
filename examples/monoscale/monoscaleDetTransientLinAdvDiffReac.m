%% Monoscale deterministic transient linear advection-diffusion-reaction problem %%
%%-------------------------------------------------------------------------------%%
% [Pares, Diez, Huerta, 2008], [Nouy, 2010]

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

filename = 'transientLinAdvDiffReac';
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
    cl1 = 0.02;
    cl2 = 0.04;
    cl0 = 0.02;
    cltip = 0.01;
    pb.S = gmshcanister(cl1,cl2,cl0,cltip,fullfile(pathname,'gmsh_canister'));
    
    %% Materials
    % Linear diffusion coefficient
    K = 0.01;
    % Thermal capacity
    c = 1;
    % Advection velocity
    Sc = pb.S;
    mat = FOUR_ISOT('k',1);
    mat = setnumber(mat,1);
    Sc = setmaterial(Sc,mat);
    P = @(i) POINT(getnode(getridge(Sc,i)));
    L1 = LIGNE(P(5),P(6));
    L2 = LIGNE(P(15),P(16));
    Sc = final(Sc);
    Sc = addcl(Sc,P(1),'T',0);
    A = calc_rigi(Sc);
    b1 = surfload(Sc,L1,'QN',-1);
    b2 = surfload(Sc,L2,'QN',1);
    b = b1+b2;
    phi = A\b;
    v = FENODEFIELD(calc_sigma(Sc,phi,'node'));
    V = 2*getvalue(v);
    V = {{FENODEFIELD(V(:,1)),FENODEFIELD(V(:,2))}};
    % Linear reaction parameter
    R1 = 0.1;
    R2 = 10;
    
    % Materials
    mat = MATERIALS();
    mat{1} = FOUR_ISOT('k',K,'c',c,'b',V,'r',R1);
    mat{2} = FOUR_ISOT('k',K,'c',c,'b',V,'r',R2);
    mat{1} = setnumber(mat{1},1);
    mat{2} = setnumber(mat{2},2);
    pb.S = setmaterial(pb.S,mat{1},1);
    pb.S = setmaterial(pb.S,mat{2},2:3);
    
    %% Dirichlet boundary conditions
    numnode = getnumber(getnode(create_boundary(pb.S)));
    [~,numnode1] = intersect(pb.S,L1);
    [~,numnode2] = intersect(pb.S,L2);
    numnoderest = setdiff(setdiff(numnode,numnode1),numnode2);
    
    pb.S = final(pb.S);
    pb.S = addcl(pb.S,numnode1,'T',1);
    pb.S = addcl(pb.S,numnode2,'T',0);
    permeability = true; % if false, the walls of the canister are considered to be impermeable
    if permeability
        pb.S = addcl(pb.S,numnoderest,'T',0);
    end
    
    %% Initial conditions
    pb.u0 = calc_init_dirichlet(pb.S);
    
    %% Time scheme
    t0 = 0;
    t1 = 2;
    nt = 100;
    T = TIMEMODEL(t0,t1,nt);
    
    % pb.N = EULERTIMESOLVER(T,'eulertype','explicit','display',true);
    pb.N = EULERTIMESOLVER(T,'eulertype','implicit','display',true);
    % pb.N = DGTIMESOLVER(T,1,'outputsplit',true,'display',true,'lu',true);
    
    %% Mass and stifness matrices and sollicitation vectors
    pb.M = calc_mass(pb.S);
    [pb.K,pb.f0] = calc_rigi(pb.S);
    pb.f0 = -pb.f0;
    pb.f = pb.f0*one(pb.N);
    
    save(fullfile(pathname,'problem.mat'),'pb','Sc','v','phi');
else
    load(fullfile(pathname,'problem.mat'),'pb','Sc','v','phi');
end

%% Solution
if solveProblem
    % Static solution at initial time step
    u = pb.K\pb.f0;
    u = unfreevector(pb.S,u);
    
    % Dynamic solution
    t = tic;
    [ut,result,vt] = dsolve(pb.N,pb.f,pb.M,pb.K);
    ut = unfreevector(pb.S,ut);
    vt = unfreevector(pb.S,vt);
    time = toc(t);
    
    save(fullfile(pathname,'solution.mat'),'u','ut','result','vt','time');
else
    load(fullfile(pathname,'solution.mat'),'u','ut','result','vt','time');
end

%% Outputs
fprintf('\n');
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
    figure('Name','Domain')
    clf
    plot(create_boundary(pb.S));
    h1 = plot(pb.S,'selgroup',1,'FaceColor',getfacecolor(1),'EdgeColor','none');
    h2 = plot(pb.S,'selgroup',2,'FaceColor',getfacecolor(2),'EdgeColor','none');
    h3 = plot(pb.S,'selgroup',3,'FaceColor',getfacecolor(3),'EdgeColor','none');
    h4 = plotfacets(pb.S,5,'FaceColor',getfacecolor(4),'EdgeColor','none');
    h5 = plotfacets(pb.S,17,'FaceColor',getfacecolor(5),'EdgeColor','none');
    set(gca,'FontSize',16)
    l = legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'$\Omega_1$','$\Omega_2$','$\Omega_3$','$\Gamma_1$','$\Gamma_2$');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    figure('Name','Advection velocity')
    clf
    plot(phi,Sc);
    colorbar
    set(gca,'FontSize',16)
    hold on
    ampl = 6;
    quiver(v,Sc,ampl,'k');
    ylim([0,1.7])
    mysaveas(pathname,'advection_velocity',formats,renderer);

    plotModel(pb.S,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    %% Display static solution
    plotSolution(pb.S,u);
    mysaveas(pathname,'solution',formats,renderer);
    plotSolution(pb.S,u,'surface',true);
    mysaveas(pathname,'solution_surface',formats);
    
    %% Display evolution of solution
    evolSolution(pb.S,ut,'filename','evol_u','pathname',pathname);
    evolSolution(pb.S,vt,'rescale',false,'filename','evol_v','pathname',pathname);
    for i=1:2
        evolSolution(pb.S,ut,'epsilon',i,'filename',['evol_eps_' num2str(i)],'pathname',pathname);
        evolSolution(pb.S,ut,'sigma',i,'filename',['evol_sig_' num2str(i)],'pathname',pathname);
    end
    
    %% Display quantity of interest
    % boutput: concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) as a function of time
    % Ioutput: total concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) along the complete time evolution,
    %          corresponding to all the pollutant that the actual filter
    %          (group #1 in mesh) is not able to retain
    foutput = bodyload(keepgroupelem(pb.S,2),[],'QN',1,'nofree');
    boutput = foutput'*ut;
    
    figure('Name','Quantity of interest')
    clf
    plot(boutput,'-b','LineWidth',1);
    grid on
    box on
    set(gca,'FontSize',16)
    xlabel('Time (s)')
    ylabel('Quantity of interest')
    mysaveas(pathname,'quantity_of_interest',formats,renderer);
    
    Ioutput = integrate(boutput);
    fprintf('quantity of interest = %e\n',Ioutput);
end

% myparallel('stop');
