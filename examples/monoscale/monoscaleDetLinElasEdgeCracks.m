%% Monoscale deterministic linear elasticity problem with n edge cracks %%
%%----------------------------------------------------------------------%%

% clc
clearvars
close all
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

loading = 'Tension'; % 'Tension' or 'Shear'
filename = ['linElasEdgeCracks' loading];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','monoscaleDet',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
fontsize = 16;
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    L = 16;
    w = 7;
    a = w/2;
    D = DOMAIN(2,[0.0,-L/2],[w,L/2]);
    
    P = [a,0.0];
    B = LINE([0.0,0.0],P);
    
    option = 'DEFO'; % plane strain
    clD = 0.25;
    clB = 0.05;
    S = gmshDomainWithSingleEdgeCrack(D,B,clD,clB,fullfile(pathname,'gmsh_domain_single_edge_crack'));
    S = setoption(S,option);
    
    %% Materials
    % Young modulus
    E = 1;
    % Poisson ratio
    NU = 0.3;
    % Thickness
    DIM3 = 1;
    % Density
    RHO = 1;
    
    % Material
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
    mat = setnumber(mat,1);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    LU = LINE([0.0,L/2],[w,L/2]);
    LL = LINE([0.0,-L/2],[w,-L/2]);
    
    S = final(S,'duplicate');
    switch lower(loading)
        case 'tension'
            S = addcl(S,POINT([a,0.0]),'UY');
            % S = addcl(S,POINT([w,0.0]),'U');
            S = addclperiodic(S,LL,LU,'UX');
            % S = addclperiodic(S,POINT([0.0,-L/2]),POINT([0.0,L/2]),'UX');
        case 'shear'
            S = addcl(S,LL);
        otherwise
            error('Wrong loading case')
    end
    
    %% Stiffness matrices and sollicitation vectors
    % Traction force density
    f = 1;
    
    A = calc_rigi(S);
    switch lower(loading)
        case 'tension'
            b = surfload(S,LU,'FY',f);
            b = b + surfload(S,LL,'FY',-f);
        case 'shear'
            b = surfload(S,LU,'FX',f);
        otherwise
            error('Wrong loading case')
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','D','B','A','b');
else
    load(fullfile(pathname,'problem.mat'),'S','D','B','A','b');
end 

%% Solution
if solveProblem
    t = tic;
    u = A\b;
    time = toc(t);
    
    e = calc_epsilon(S,u);
    s = calc_sigma(S,u);
    
    save(fullfile(pathname,'solution.mat'),'u','e','s','time');
else
    load(fullfile(pathname,'solution.mat'),'u','e','s','time');
end

%% Outputs
fprintf('\n');
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('elapsed time = %f s\n',time);

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain({D,B},'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 0.5;
    [hN,legN] = vectorplot(S,'F',b,ampl,'r','LineWidth',1);
    % legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    % plotModel(S,'legend',false);
    % mysaveas(pathname,'mesh',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    ampl = getsize(S)/max(abs(u))/20;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*unfreevector(S,u),'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    %% Display solutions
    ampl = 0;
    % ampl = getsize(S)/max(abs(u))/20;
    
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
    % set(gca,'FontSize',fontsize)
    % mysaveas(pathname,'eps_xx',formats,renderer);
    %
    % figure('Name','Solution eps_yy')
    % clf
    % plot(e,S+ampl*u,'compo','EPYY')
    % colorbar
    % set(gca,'FontSize',fontsize)
    % mysaveas(pathname,'eps_yy',formats,renderer);
    %
    % figure('Name','Solution eps_xy')
    % clf
    % plot(e,S+ampl*u,'compo','EPXY')
    % colorbar
    % set(gca,'FontSize',fontsize)
    % mysaveas(pathname,'eps_xy',formats,renderer);
    
    % figure('Name','Solution sig_xx')
    % clf
    % plot(s,S+ampl*u,'compo','SMXX')
    % colorbar
    % set(gca,'FontSize',fontsize)
    % mysaveas(pathname,'sig_xx',formats,renderer);
    %
    % figure('Name','Solution sig_yy')
    % clf
    % plot(s,S+ampl*u,'compo','SMYY')
    % colorbar
    % set(gca,'FontSize',fontsize)
    % mysaveas(pathname,'sig_yy',formats,renderer);
    %
    % figure('Name','Solution sig_xy')
    % clf
    % plot(s,S+ampl*u,'compo','SMXY')
    % colorbar
    % set(gca,'FontSize',fontsize)
    % mysaveas(pathname,'sig_xy',formats,renderer);
end

% myparallel('stop');
