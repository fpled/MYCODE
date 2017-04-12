%% Specimen - Determinsitic isotropic linear elasticity problem %%
%%--------------------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

filename = 'specimenDetLinElasIsot';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    L = 1;
    D = DOMAIN(2,[0.0,0.0],[L,L]);
    
    % elemtype = 'TRI3';
    elemtype = 'QUA4';
    option = 'DEFO'; % plane strain
    % option = 'CONT'; % plane stress
    nbelem = [50,50];
    S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'option',option);
    % cl = 0.05;
    % S = build_model(D,'cl',cl,'elemtype',elemtype,'option',option,'filename',[pathname 'gmsh_domain']);
    
    %% Materials
    % Poisson ratio
    NU = 0.3;
    % Thickness
    DIM3 = 1;
    % Density
    RHO = 1;
    % Young modulus
    E = 1;
    
    % Material
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
    mat = setnumber(mat,1);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    LU = LIGNE([0.0,L],[L,L]);
    LL = LIGNE([0.0,0.0],[L,0.0]);
    
    S = final(S);
    S = addcl(S,LL);
    
    % loading = 'Dirichlet'; % Imposed displacement
    loading = 'Neumann'; % Traction force density
    % degree = 'cst'; % constant loading
    % degree = 'lin'; % linear loading
    degree = 'qua'; % quadratic loading
    if strcmpi(loading,'dirichlet')
        udmax = -1e-2;
        switch lower(degree)
            case 'cst'
                ud = udmax;
            case 'lin'
                ud = @(x) udmax*(L-x(:,1));
            case 'qua'
                ud = @(x) 4*udmax/L^2*x(:,1).*(L-x(:,1));
        end
        S = addcl(S,LU,'UY',ud);
    end
    
    %% Stiffness matrices and sollicitation vectors
    switch lower(loading)
        case 'neumann'
            A = calc_rigi(S);
            fmax = -1;
            switch lower(degree)
                case 'cst'
                    f = fmax;
                case 'lin'
                    f = @(x) fmax*(L-x(:,1));
                case 'qua'
                    f = @(x) fmax*4/L^2*x(:,1).*(L-x(:,1));
            end
            b = surfload(S,LU,'FY',f);
        case 'dirichlet'
            [A,b] = calc_rigi(S);
            b = -b;
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'loading','elemtype','S','D','b');
else
    load(fullfile(pathname,'problem.mat'),'loading','elemtype','S','D','b');
end

%% Solution
if solveProblem
    t = tic;
    u = A\b;
    time = toc(t);
    
    e = calc_epsilon(S,u);
    s = calc_sigma(S,u);
    
    save(fullfile(pathname,'solution.mat'),'u','time','e','s');
else
    load(fullfile(pathname,'solution.mat'),'u','time','e','s');
end

%% Outputs
fprintf('\nSquare specimen\n');
fprintf(['load     : ' loading '\n']);
fprintf(['mesh     : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain(D,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 0.5;
    switch lower(loading)
        case 'neumann'
            [hN,legN] = vectorplot(S,'F',b,ampl,'r','LineWidth',1);
        case 'dirichlet'
            v = calc_init_dirichlet(S);
            [hN,legN] = vectorplot(S,'U',v,ampl,'r','LineWidth',1);
    end
    % legend([hD,hN],[legD,legN])
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    % plotModel(S,'legend',false);
    % mysaveas(pathname,'mesh',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    ampl = getsize(S)/max(abs(u))/5;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*unfreevector(S,u),'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    %% Display solution
    % ampl = 0;
    ampl = getsize(S)/max(abs(u))/5;
    
    for i=1:2
        plotSolution(S,u,'displ',i,'ampl',ampl);
        mysaveas(pathname,['u_' num2str(i)],formats,renderer);
    end
    
    for i=1:3
        plotSolution(S,u,'epsilon',i,'ampl',ampl);
        mysaveas(pathname,['eps_' num2str(i)],formats,renderer);
        
        plotSolution(S,u,'sigma',i,'ampl',ampl);
        mysaveas(pathname,['sig_' num2str(i)],formats,renderer);
    end
    
    plotSolution(S,u,'epsilon','mises','ampl',ampl);
    mysaveas(pathname,'eps_von_mises',formats,renderer);
    
    plotSolution(S,u,'sigma','mises','ampl',ampl);
    mysaveas(pathname,'sig_von_mises',formats,renderer);
    
    % u = unfreevector(S,u);
    %
    % figure('Name','Solution eps_xx')
    % clf
    % plot(e,S+ampl*u,'compo','EPXX')
    % colorbar
    % set(gca,'FontSize',16)
    % mysaveas(pathname,'eps_xx',formats,renderer);
    %
    % figure('Name','Solution eps_yy')
    % clf
    % plot(e,S+ampl*u,'compo','EPYY')
    % colorbar
    % set(gca,'FontSize',16)
    % mysaveas(pathname,'eps_yy',formats,renderer);
    %
    % figure('Name','Solution eps_xy')
    % clf
    % plot(e,S+ampl*u,'compo','EPXY')
    % colorbar
    % set(gca,'FontSize',16)
    % mysaveas(pathname,'eps_xy',formats,renderer);
    
    % figure('Name','Solution sig_xx')
    % clf
    % plot(s,S+ampl*u,'compo','SMXX')
    % colorbar
    % set(gca,'FontSize',16)
    % mysaveas(pathname,'sig_xx',formats,renderer);
    %
    % figure('Name','Solution sig_yy')
    % clf
    % plot(s,S+ampl*u,'compo','SMYY')
    % colorbar
    % set(gca,'FontSize',16)
    % mysaveas(pathname,'sig_yy',formats,renderer);
    %
    % figure('Name','Solution sig_xy')
    % clf
    % plot(s,S+ampl*u,'compo','SMXY')
    % colorbar
    % set(gca,'FontSize',16)
    % mysaveas(pathname,'sig_xy',formats,renderer);
    
end
