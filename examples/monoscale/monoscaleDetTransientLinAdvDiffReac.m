%% Monoscale deterministic transient linear advection-diffusion-reaction problem %%
%%-------------------------------------------------------------------------------%%
% [Pares, Diez, Huerta, 2008], [Nouy, 2010]

% clc
clearvars
close all

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
formats = {'fig','epsc'};
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
    Sadv = pb.S;
    mat = FOUR_ISOT('k',1);
    mat = setnumber(mat,1);
    Sadv = setmaterial(Sadv,mat);
    P = @(i) POINT(getnode(getridge(Sadv,i)));
    L1 = LIGNE(P(5),P(6));
    L2 = LIGNE(P(15),P(16));
    Sadv = final(Sadv);
    Sadv = addcl(Sadv,P(1),'T',0);
    A = calc_rigi(Sadv);
    b1 = surfload(Sadv,L1,'QN',-1);
    b2 = surfload(Sadv,L2,'QN',1);
    b = b1+b2;
    phi = A\b;
    v = 2*FENODEFIELD(calc_sigma(Sadv,phi,'node'));
    V = getvalue(v);
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
    % already taken into account by Dirichlet boundary conditions
    % pb.u0 = calc_init_dirichlet(pb.S);
    
    %% Time scheme
    t0 = 0;
    t1 = 2;
    nt = 100;
    T = TIMEMODEL(t0,t1,nt);
    
    % pb.N = EULERTIMESOLVER(T,'eulertype','explicit','display',true);
    pb.N = EULERTIMESOLVER(T,'eulertype','implicit','display',true);
    % pb.N = DGTIMESOLVER(T,1,'outputsplit',true,'display',true,'lu',true);
    
    pb.loadFunction = @(N) one(N);
    
    %% Mass and stifness matrices and sollicitation vectors
    pb.M = calc_mass(pb.S);
    [pb.A,pb.b0] = calc_rigi(pb.S);
    pb.b0 = -pb.b0;
    pb.b = pb.b0*pb.loadFunction(pb.N);
    
    save(fullfile(pathname,'problem.mat'),'pb','Sadv','v','phi');
else
    load(fullfile(pathname,'problem.mat'),'pb','Sadv','v','phi');
end

%% Solution
if solveProblem
    % Stationary solution
    t = tic;
    u = pb.A\pb.b0;
    % u = unfreevector(pb.S,u);
    time = toc(t);
    
    % Transient solution
    t = tic;
    [ut,result,vt] = dsolve(pb.N,pb.b,pb.M,pb.A);
    % ut = unfreevector(pb.S,ut);
    vt = unfreevector(pb.S,vt)-calc_init_dirichlet(pb.S);
    timet = toc(t);
    
    save(fullfile(pathname,'solution.mat'),'u','time');
    save(fullfile(pathname,'solution_transient.mat'),'ut','result','vt','timet');
else
    load(fullfile(pathname,'solution.mat'),'u','time');
    load(fullfile(pathname,'solution_transient.mat'),'ut','result','vt','timet');
end

%% Outputs
fprintf('\n');
fprintf('nb elements = %g\n',getnbelem(pb.S));
fprintf('nb nodes    = %g\n',getnbnode(pb.S));
fprintf('nb dofs     = %g\n',getnbddl(pb.S));
fprintf('time solver : %s\n',class(pb.N));
fprintf('nb time steps = %g\n',getnt(pb.N));
fprintf('nb time dofs  = %g\n',getnbtimedof(pb.N));
fprintf('elapsed time = %f s for stationary solution\n',time);
fprintf('elapsed time = %f s for transient solution\n',timet);

%% Display
if displaySolution
    %% Display domains and meshes
    figure('Name','Domain')
    clf
    plot(create_boundary(pb.S));
    hold on
    h1 = plot(pb.S,'selgroup',1,'FaceColor',getfacecolor(1),'EdgeColor','none');
    h2 = plot(pb.S,'selgroup',2,'FaceColor',getfacecolor(2),'EdgeColor','none');
    h3 = plot(pb.S,'selgroup',3,'FaceColor',getfacecolor(3),'EdgeColor','none');
    h4 = plotfacets(pb.S,5,'FaceColor',getfacecolor(4),'EdgeColor','none');
    h5 = plotfacets(pb.S,17,'FaceColor',getfacecolor(5),'EdgeColor','none');
    hold off
    set(gca,'FontSize',16)
    l = legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'$\Omega_1$','$\Omega_2$','$\Omega_3$','$\Gamma_D^1$','$\Gamma_D^2$','Location','NorthEast');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    figure('Name','Advection velocity')
    clf
    plot(phi,Sadv);
    colorbar
    set(gca,'FontSize',16)
    hold on
    ampl = 6;
    quiver(v,Sadv,ampl,'k');
    hold off
    ylim([0,1.7])
    mysaveas(pathname,'advection_velocity',formats,renderer);
    
    figure('Name','Mesh')
    clf
    h1 = plot(pb.S,'selgroup',1,'FaceColor',getfacecolor(1));
    hold on
    h2 = plot(pb.S,'selgroup',2,'FaceColor',getfacecolor(2));
    h3 = plot(pb.S,'selgroup',3,'FaceColor',getfacecolor(3));
    h4 = plotfacets(pb.S,5,'FaceColor',getfacecolor(4));
    h5 = plotfacets(pb.S,17,'FaceColor',getfacecolor(5));
    hold off
    set(gca,'FontSize',16)
    l = legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'$\Omega_1$','$\Omega_2$','$\Omega_3$','$\Gamma_D^1$','$\Gamma_D^2$','Location','NorthEast');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'mesh',formats,renderer);
    
    %% Display stationary solution
    plotSolution(pb.S,u);
    mysaveas(pathname,'solution',formats,renderer);
    plotSolution(pb.S,u,'surface',true);
    mysaveas(pathname,'solution_surface',formats);
    
    %% Display evolution of transient solution
    evolSolution(pb.S,ut,'filename','evol_solution','pathname',pathname);
    evolSolution(pb.S,ut,'surface',true,'filename','evol_solution_surface','pathname',pathname);
    
    evolSolution(pb.S,vt,'rescale',false,'filename','evol_velocity','pathname',pathname);
    evolSolution(pb.S,vt,'rescale',false,'surface',true,'filename','evol_velocity_surface','pathname',pathname);
    
%     for i=1:2
%         evolSolution(pb.S,ut,'epsilon',i,'filename',['evol_eps_' num2str(i)],'pathname',pathname);
%         evolSolution(pb.S,ut,'sigma',i,'filename',['evol_sig_' num2str(i)],'pathname',pathname);
%     end
    
    %% Display transient solution at differents instants
%     [t,rep] = gettevol(pb.N);
%     for k=1:floor(length(rep)/4):length(rep)
%         close all
%         uk = getmatrixatstep(ut,rep(k));
%         vk = getmatrixatstep(vt,rep(k));
%         
%         plotSolution(pb.S,uk);
%         mysaveas(pathname,['solution_t' num2str(k-1)],formats,renderer);
%         plotSolution(pb.S,uk,'surface',true);
%         mysaveas(pathname,['solution_t' num2str(k-1) '_surface'],formats);
%         
%         plotSolution(pb.S,vk);
%         mysaveas(pathname,['velocity_t' num2str(k-1)],formats,renderer);
%         plotSolution(pb.S,vk,'surface',true);
%         mysaveas(pathname,['velocity_t' num2str(k-1) '_surface'],formats);
%     end
    
    %% Display quantity of interest
    % boutput: concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) as a function of time
    % Ioutput: total concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) along the complete time evolution,
    %          corresponding to all the pollutant that the actual filter
    %          (group #1 in mesh) is not able to retain
    foutput = bodyload(keepgroupelem(pb.S,2),[],'QN',1,'nofree');
    ut = unfreevector(pb.S,ut);
    boutput = foutput'*ut;
    
    figure('Name','Quantity of interest')
    clf
    plot(boutput,'-b','LineWidth',1);
    grid on
    box on
    set(gca,'FontSize',16)
    xlabel('Time (s)')
    ylabel('Concentration of pollutant in trap domain')
    mysaveas(pathname,'quantity_of_interest',formats,renderer);
    mymatlab2tikz(pathname,'quantity_of_interest.tex');
    
    Ioutput = integrate(boutput);
    fprintf('quantity of interest = %e\n',Ioutput);
end
