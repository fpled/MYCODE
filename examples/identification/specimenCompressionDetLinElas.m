%% Specimen under compression - Determinsitic linear elasticity problem %%
%%----------------------------------------------------------------------%%

% clc
clearvars
close all

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

Dim = 2; % space dimension Dim = 2, 3
filename = ['specimenCompressionDetLinElas_' num2str(Dim) 'D'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

t = tic;

%% Problem
if setProblem
    tMesh = tic;
    %% Domains and meshes
    L = 1e-2; % [m]
    if Dim==2
        D = DOMAIN(2,[0.0,0.0],[L,L]);
        % elemtype = 'TRI3';
        elemtype = 'QUA4';
        % elemtype = 'TRI6';
    elseif Dim==3
        D = DOMAIN(3,[0.0,0.0,0.0],[L,L,L]);
        % elemtype = 'TET4';
        elemtype = 'CUB8';
        % elemtype = 'TET10';
    end
    % option = 'DEFO'; % plane strain
    option = 'CONT'; % plane stress
    nbelem = repmat(10,1,Dim);
    S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'option',option);
    % cl = L/20;
    % S = build_model(D,'cl',cl,'elemtype',elemtype,'option',option,'filename',fullfile(pathname,'gmsh_domain'));
    timeMesh = toc(tMesh);
    fprintf('Meshing: elapsed time = %f s\n',timeMesh);
    
    tAssemble = tic;
    %% Materials
    % Thickness
    DIM3 = 1; % [m]
    % Density
    RHO = 1; % [kg/m3]
    
    % Material
    % Material Symmetry
    materialSym = 'isotTrans';
    
    switch lower(materialSym)
        case 'isot'
            % Young modulus
            E = 1; % [Pa]
            % Poisson ratio
            NU = 0.3;
            % Material
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
        case 'isottrans'
            % Transverse Young modulus
            ET = 1; % [Pa]
            % Longitudinal shear modulus
            GL = 1.5; % [Pa]
            % Longitudinal Young modulus
            EL = 2; % [Pa]
            % Longitudinal Poisson ratio
            NUL = 0.2;
            % Transverse Poisson ratio
            NUT = 0.3;
            % Material
            mat = ELAS_ISOT_TRANS('AXISL',[0;0;1],'AXIST',[1;0;0],'EL',EL,'ET',ET,'NUL',NUL,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',DIM3);
        otherwise
            error('Wrong material symmetry !')
    end
    mat = setnumber(mat,1);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BU = LINE([0.0,L],[L,L]);
        BL = LINE([0.0,0.0],[L,0.0]);
    elseif Dim==3
        BU = PLANE([0.0,0.0,L],[L,0.0,L],[0.0,L,L]);
        BL = PLANE([0.0,0.0,0.0],[L,0.0,0.0],[0.0,L,0.0]);
    end
    
    S = final(S);
    S = addcl(S,BL);
    
    loading = 'Dirichlet'; % Imposed displacement
    % loading = 'Neumann'; % Traction force density
    degree = 'cst'; % constant loading
    % degree = 'lin'; % linear loading
    % degree = 'qua'; % quadratic loading
    if strcmpi(loading,'dirichlet')
        udmax = -1e-3; % [m]
        switch lower(degree)
            case 'cst'
                ud = udmax;
            case 'lin'
                ud = @(x) udmax*(L-x(:,1))/L;
            case 'qua'
                ud = @(x) 4*udmax*x(:,1).*(L-x(:,1))/L^2;
        end
        if Dim==2
            S = addcl(S,BU,'UY',ud);
        elseif Dim==3
            S = addcl(S,BU,'UZ',ud);
        end
    end
    
    %% Stiffness matrices and sollicitation vectors
    switch lower(loading)
        case 'neumann'
            A = calc_rigi(S);
            fmax = -5e-5/L^2; % surface load [N/m2]
            switch lower(degree)
                case 'cst'
                    f = fmax;
                case 'lin'
                    f = @(x) fmax*(L-x(:,1))/L;
                case 'qua'
                    f = @(x) 4*fmax*x(:,1).*(L-x(:,1))/L^2;
            end
            if Dim==2
                b = surfload(S,BU,'FY',f);
            elseif Dim==3
                b = surfload(S,BU,'FZ',f);
            end
        case 'dirichlet'
            [A,b] = calc_rigi(S);
            b = -b;
    end
    timeAssemble = toc(tAssemble);
    fprintf('Assembling: elapsed time = %f s\n',timeAssemble);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'loading','elemtype','S','D','A','b');
else
    load(fullfile(pathname,'problem.mat'),'loading','elemtype','S','D','A','b');
end

%% Solution
if solveProblem
    tSolve = tic;
    % u = A\b;
    R = chol(A);
    u = R\(R'\b);
    timeSolve = toc(tSolve);
    fprintf('Solving: elapsed time = %f s\n',timeSolve);
    
    tPostProcess = tic;
    e = calc_epsilon(S,u);
    s = calc_sigma(S,u);
    timePostProcess = toc(tPostProcess);
    fprintf('Post-processing: elapsed time = %f s\n',timePostProcess);
    
    save(fullfile(pathname,'solution.mat'),'u','e','s');
else
    load(fullfile(pathname,'solution.mat'),'u','e','s');
end

time = toc(t);

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
    % legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
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
    ampl = 0;
    % ampl = getsize(S)/max(abs(u))/5;
    
    for i=1:Dim
        plotSolution(S,u,'displ',i,'ampl',ampl);
        mysaveas(pathname,['u_' num2str(i)],formats,renderer);
    end
    
    for i=1:(Dim*(Dim+1)/2)
        plotSolution(S,u,'epsilon',i,'ampl',ampl);
        mysaveas(pathname,['eps_' num2str(i)],formats,renderer);
        
        plotSolution(S,u,'sigma',i,'ampl',ampl);
        mysaveas(pathname,['sig_' num2str(i)],formats,renderer);
    end
    
    plotSolution(S,u,'epsilon','mises','ampl',ampl);
    mysaveas(pathname,'eps_von_mises',formats,renderer);
    
    plotSolution(S,u,'sigma','mises','ampl',ampl);
    mysaveas(pathname,'sig_von_mises',formats,renderer);
    
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
