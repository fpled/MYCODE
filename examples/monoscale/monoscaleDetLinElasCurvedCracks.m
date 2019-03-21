%% Monoscale deterministic linear elasticity problem with curve crack %%
%%--------------------------------------------------------_-----------%%
% [Miehe, Hofacker, Welschinger, 2010, CMAME]
% [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]


% clc
clearvars
close all
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

loading = 'Shear'; % 'Pull' or 'Shear'
filename = ['linElasCurvedCracks' loading];
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
    L = 1e-3;
    a = L/2;
    D = DOMAIN(2,[0.0,0.0],[L,L]);
    
    P = [a,L/2];
    B = LIGNE([0.0,L/2],P);
    
    option = 'CONT'; % plane stress
    clD = 0.02e-3;
    clP = 6e-7;
    S = gmshdomainwithedgecrack(D,P,clD,clP,fullfile(pathname,'gmsh_domain_curved_crack'));
    S = setoption(S,option);
    
    %% Materials
    % Lame coefficients
    lambda = 121.15e9;
    mu = 80.77e9;
    % Young modulus
    % Poisson ratio
    switch option
        case 'DEFO'
            E = mu*(3*lambda+2*mu)/(lambda+mu);
            NU = lambda/(lambda+mu)/2;
        case 'CONT'
            E = 4*mu*(lambda+mu)/(lambda+2*mu);
            NU = lambda/(lambda+2*mu);
    end
    % Thickness
    DIM3 = 1;
    % Density
    RHO = 1;
    % Fracture toughness
    gc = 2700;
    
    % Material
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
    mat = setnumber(mat,1);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    LU = LIGNE([0.0,L],[L,L]);
    LL = LIGNE([0.0,0.0],[L,0.0]);
    LRight = LIGNE([L,0.0],[L,L]);
    LLeft = LIGNE([0.0,0.0],[0.0,L]);
    
    S = final(S,'duplicate');
    
    ud = 1e-5;
    switch lower(loading)
        case 'pull'
            S = addcl(S,LU,'UY',ud);
        case 'shear'
            S = addcl(S,LU,{'UX','UY'},[ud;0]);
            S = addcl(S,LRight,'UY');
            S = addcl(S,LLeft,'UY');
        otherwise
            error('Wrong loading case')
    end
    S = addcl(S,LL);
    
    %% Stiffness matrices and sollicitation vectors
    [A,b] = calc_rigi(S);
    b = -b;
    
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
    v = calc_init_dirichlet(S);
    [hN,legN] = vectorplot(S,'U',v,ampl,'r','LineWidth',1);
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
