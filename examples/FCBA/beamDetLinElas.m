%% Beam deterministic linear elasticity %%
%%--------------------------------------%%

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

Dim = 3; % space dimension Dim = 2, 3
filename = ['beamDetLinElas_' num2str(Dim) 'D'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','FCBA',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Problem
if solveProblem
    %% Domains and meshes
    % Beams dimensions
    L = 1; % [m]
    b = 0.2;
    h = 0.5;
    
    % Points
    if Dim==2
        x = [0.0,0.0;
            L,0.0];
        x_load = [L/2,0.0];
    elseif Dim==3
        x = [0.0,0.0,0.0;
            L,0.0,0.0];
        x_load = [L/2,0.0,0.0];
    end
    P = {x(1,:),x_load,x(2,:)};
    
    % Beams meshes
    elemtype = 'BEAM';
    cl = L/100;
    
    S = gmshbeam(P,cl,fullfile(pathname,'gmsh'));
    S = concatgroupelem(S);
    if Dim==2
        S = convertelem(S,elemtype);
    elseif Dim==3
        S = convertelem(S,elemtype,'param',VECTEUR([0;1;0]));
    end
    
    %% Materials
    % Gravitational acceleration
    g = 9.81; % [m/s2]
    
    % Density
    RHO = 1; % [kg/m3]
    
    % Cross-section area
    Sec = b*h;
    % Planar second moment of area (or Planar area moment of inertia)
    IY = h*b^3/12;
    IZ = b*h^3/12;
    % Polar second moment of area (or Polar area moment of inertia)
    IX = IY+IZ;
    
    % Material symmetry
    materialSym = 'isot';
    
    switch lower(materialSym)
        case 'isot'
            % Young modulus
            E = 10e9; % [Pa]
            % Poisson ratio
            NU = 0.3;
            % Material
            mat = ELAS_BEAM('E',E,'NU',NU,'S',Sec,'IZ',IZ,'IY',IY,'IX',IX,'RHO',RHO);
            mat = setnumber(mat,1);
            S = setmaterial(S,mat);
        otherwise
            error('Wrong material symmetry !')
    end
    
    %% Neumann boundary conditions
    pl = RHO*g*Sec; % line load (body load for beams) [N/m]
    p = 1; % pointwise load, 200 [N]
    
    %% Dirichlet boundary conditions
    S = final(S);
    P = POINT(x);
    S = addcl(S,P,'U');
    if Dim==3
        S = addcl(S,P,{'RX','RY'});
    end
    
    %% Stiffness matrix and sollicitation vector
    A = calc_rigi(S);
    P_load = POINT(x_load);
    f = nodalload(S,P_load,'FY',-p);
    %f = f + bodyload(S,[],'FY',-pl);
    
    %% Solution
    t = tic;
    u = A\f;
    time = toc(t);
    
    u = unfreevector(S,u);
    
    % U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    Ux = u(findddl(S,'UX'),:);
    Uy = u(findddl(S,'UY'),:);
    if Dim==3
        Uz = u(findddl(S,'UZ'),:);
    end
    % R = u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    if Dim==3
        Rx = u(findddl(S,'RX'),:);
        Ry = u(findddl(S,'RY'),:);
    end
    Rz = u(findddl(S,'RZ'),:);
    
    x = getcoord(S.node);
    
    e = calc_epsilon(S,u,'smooth');
    s = calc_sigma(S,u,'smooth');
    
    Epsx = e(1);
    if Dim==2
        Gamz = e(2);
    elseif Dim==3
        Gamx = e(2);
        Gamy = e(3);
        Gamz = e(4);
    end
    N = s(1);
    if Dim==2
        Mz = s(2);
    elseif Dim==3
        Mx = s(2);
        My = s(3);
        Mz = s(4);
    end
    
    %% Test solution
    P = P_load;
    numnode = find(S.node==P);
    xP = x(numnode,:);
    
    ux = eval_sol(S,u,P,'UX');
    uy = eval_sol(S,u,P,'UY');
    if Dim==3
        uz = eval_sol(S,u,P,'UZ');
        rx = eval_sol(S,u,P,'RX');
        ry = eval_sol(S,u,P,'RY');
    end
    rz = eval_sol(S,u,P,'RZ');
    
    [~,~,numgroupelem] = findelemwithnode(S,numnode);
    n  = reshape(N{numgroupelem},[getnbnode(S),1]);
    if Dim==3
        mx = reshape(Mx{numgroupelem},[getnbnode(S),1]);
        my = reshape(My{numgroupelem},[getnbnode(S),1]);
    end
    mz = reshape(Mz{numgroupelem},[getnbnode(S),1]);
    epsx = reshape(Epsx{numgroupelem},[getnbnode(S),1]);
    if Dim==3
        gamx = reshape(Gamx{numgroupelem},[getnbnode(S),1]);
        gamy = reshape(Gamy{numgroupelem},[getnbnode(S),1]);
    end
    gamz = reshape(Gamz{numgroupelem},[getnbnode(S),1]);
    
    n = double(n(numnode));
    if Dim==3
        mx = double(mx(numnode));
        my = double(my(numnode));
    end
    mz = double(mz(numnode));
    epsx = double(epsx(numnode));
    if Dim==3
        gamx = double(gamx(numnode));
        gamy = double(gamy(numnode));
    end
    gamz = double(gamz(numnode));
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S',...
        'L','h','b','f');
    if Dim==2
        save(fullfile(pathname,'solution.mat'),'u','s','e','time',...
            'Ux','Uy','Rz',...
            'N','Mz','Epsx','Gamz');
        save(fullfile(pathname,'test_solution.mat'),'P',...
            'ux','uy','rz',...
            'n','mz','epsx','gamz');
    elseif Dim==3
        save(fullfile(pathname,'solution.mat'),'u','s','e','time',...
            'Ux','Uy','Uz','Rx','Ry','Rz',...
            'N','Mx','My','Mz','Epsx','Gamx','Gamy','Gamz');
        save(fullfile(pathname,'test_solution.mat'),'P',...
            'ux','uy','uz','rx','ry','rz',...
            'n','mx','my','mz','epsx','gamx','gamy','gamz');
    end
else
    load(fullfile(pathname,'problem.mat'),'S',...
        'L','h','b','f');
    if Dim==2
        load(fullfile(pathname,'solution.mat'),'u','s','e','time',...
            'Ux','Uy','Rz',...
            'N','Mz','Epsx','Gamz');
        load(fullfile(pathname,'test_solution.mat'),'P',...
            'ux','uy','rz',...
            'n','mz','epsx','gamz');
    elseif Dim==3
        load(fullfile(pathname,'solution.mat'),'u','s','e','time',...
            'Ux','Uy','Uz','Rx','Ry','Rz',...
            'N','Mx','My','Mz','Epsx','Gamx','Gamy','Gamz');
        load(fullfile(pathname,'test_solution.mat'),'P',...
            'ux','uy','uz','rx','ry','rz',...
            'n','mx','my','mz','epsx','gamx','gamy','gamz');
    end
end

%% Outputs
fprintf('\nBeam\n');
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

disp('Displacement u and rotation r at point'); disp(P);
fprintf('ux = %g m\n',ux);
fprintf('uy = %g m\n',uy);
if Dim==3
    fprintf('uz = %g m\n',uz);
    fprintf('rx = %g rad = %g deg\n',rx,rad2deg(rx));
    fprintf('ry = %g rad = %g deg\n',ry,rad2deg(ry));
end
fprintf('rz = %g rad = %g deg\n',rz,rad2deg(rz));
fprintf('\n');

if Dim==2
    disp('Force N and moment Mz at point'); disp(P);
elseif Dim==3
    disp('Force N and moments Mx, My, Mz at point'); disp(P);
end
fprintf('N  = %g N\n',n);
if Dim==3
    fprintf('Mx = %g N.m\n',mx);
    fprintf('My = %g N.m\n',my);
end
fprintf('Mz = %g N.m\n',mz);
fprintf('\n');

if Dim==2
    disp('Axial strain Epsx and bending strain (curvature) Gamz at point'); disp(P);
elseif Dim==3
    disp('Axial strain Epsx, torsion and bending strains (curvatures) Gamx, Gamy, Gamz at point'); disp(P);
end
fprintf('Epsx = %g\n',epsx);
if Dim==3
    fprintf('Gamx = %g\n',gamx);
    fprintf('Gamy = %g\n',gamy);
end
fprintf('Gamz = %g\n',gamz);
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
%     plotparamelem(S,'group')
%     plotparamelem(S,'material')
    
%     plotDomain(S,'legend',false);
%     mysaveas(pathname,'domain',formats,renderer);
%     mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'FaceColor','k','legend',false);
    ampl = getsize(S)/max(abs(f))/10;
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
    %legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','node',true,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    ampl = getsize(S)/max(abs(U))/20;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','node',true,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','node',true);
    plot(S+ampl*u,'Color','b','FaceColor','b','node',true);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    %% Display solution
    ampl = 0;
    % ampl = getsize(S)/max(abs(U))/20;
    options = {'solid',true};
    % options = {};
    
    % Displacements
    plotSolution(S,u,'displ',1,'ampl',ampl,options{:});
    mysaveas(pathname,'Ux',formats,renderer);
    
    plotSolution(S,u,'displ',2,'ampl',ampl,options{:});
    mysaveas(pathname,'Uy',formats,renderer);
    
    if Dim==3
        plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
        mysaveas(pathname,'Uz',formats,renderer);
    end
    
    % Rotations
    if Dim==2
        plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
        mysaveas(pathname,'Rz',formats,renderer);
    elseif Dim==3
        plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
        mysaveas(pathname,'Rx',formats,renderer);
        
        plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
        mysaveas(pathname,'Ry',formats,renderer);
        
        plotSolution(S,u,'rotation',3,'ampl',ampl,options{:});
        mysaveas(pathname,'Rz',formats,renderer);
    end
    
    % DO NOT WORK WITH BEAM ELEMENTS
    % Strains
    % plotSolution(S,u,'epsilon',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'Epsx',formats,renderer);
    % if Dim==2
    %     plotSolution(S,u,'epsilon',2,'ampl',ampl,options{:});
    %     mysaveas(pathname,'Gamz',formats,renderer);
    % elseif Dim==3
    %     plotSolution(S,u,'epsilon',2,'ampl',ampl,options{:});
    %     mysaveas(pathname,'Gamx',formats,renderer);
    %     
    %     plotSolution(S,u,'epsilon',3,'ampl',ampl,options{:});
    %     mysaveas(pathname,'Gamy',formats,renderer);
    %     
    %     plotSolution(S,u,'epsilon',4,'ampl',ampl,options{:});
    %     mysaveas(pathname,'Gamz',formats,renderer);
    % end
    %
    % Stresses
    % plotSolution(S,u,'sigma',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'N',formats,renderer);
    % if Dim==2
    %     plotSolution(S,u,'sigma',2,'ampl',ampl,options{:});
    %     mysaveas(pathname,'Mz',formats,renderer);
    % elseif Dim==3
    %     plotSolution(S,u,'sigma',2,'ampl',ampl,options{:});
    %     mysaveas(pathname,'Mx',formats,renderer);
    %     
    %     plotSolution(S,u,'sigma',3,'ampl',ampl,options{:});
    %     mysaveas(pathname,'My',formats,renderer);
    %     
    %     plotSolution(S,u,'sigma',4,'ampl',ampl,options{:});
    %     mysaveas(pathname,'Mz',formats,renderer);
    % end
    
    % Strains
    figure('Name','Solution eps_x')
    clf
    plot(e,S+ampl*u,'compo','EPSX')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Epsx',formats,renderer);
    
    if Dim==3
        figure('Name','Solution gam_x')
        clf
        plot(e,S+ampl*u,'compo','GAMX')
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,'Gamx',formats,renderer);
        
        figure('Name','Solution gam_y')
        clf
        plot(e,S+ampl*u,'compo','GAMY')
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,'Gamy',formats,renderer);
    end
    
    figure('Name','Solution gam_z')
    clf
    plot(e,S+ampl*u,'compo','GAMZ')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Gamz',formats,renderer);
    
    % Stresses
    figure('Name','Solution N')
    clf
    plot(s,S+ampl*u,'compo','EFFX')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'N',formats,renderer);
    
    if Dim==3
        figure('Name','Solution Mx')
        clf
        plot(s,S+ampl*u,'compo','MOMX')
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,'Mx',formats,renderer);
        
        figure('Name','Solution My')
        clf
        plot(s,S+ampl*u,'compo','MOMY')
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,'My',formats,renderer);
    end
    
    figure('Name','Solution Mz')
    clf
    plot(s,S+ampl*u,'compo','MOMZ')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Mz',formats,renderer);
end
