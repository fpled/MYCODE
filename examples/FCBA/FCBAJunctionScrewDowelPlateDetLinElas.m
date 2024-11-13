%% FCBA junction screw dowel plate deterministic linear elasticity %%
%%-----------------------------------------------------------------%%

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

filename = 'FCBAJunctionScrewDowelPlateDetLinElas';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','FCBA',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

junction = false; % junction modeling
selfweight = false; % self-weight

%% Problem
if solveProblem
    %% Domains and meshes
    % Plates dimensions
    L1 = 142.5e-3; % [m]
    L2 = 67.5e-3;
    b = 113e-3;
    h = 15e-3;
    
    Q1 = QUADRANGLE([0,0,0],[0,b,0],...
        [0,b,L1],[0,0,L1]);
    Q2 = QUADRANGLE([0,b,L1],[0,0,L1],...
        [L2,0,L1],[L2,b,L1]);
    
    % Plates meshes
    elemtype = 'DKT';
    cl = h/5;
    S1 = build_model(Q1,'cl',cl,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_plate_1_' elemtype]));
    S2 = build_model(Q2,'cl',cl,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_plate_2_' elemtype]));
    
    %% Materials
    % Gravitational acceleration
    g = 9.81; % [m/s2]
    
    % Density
    RHO = 707.1384; % [kg/m3]
    Vol_total = h*(L1+L2)*b;
    Mass_total = Vol_total*RHO; % [kg]
    
    % Data
    filenameAna = 'data_ET_GL.mat';
    filenameNum = 'data_EL_NUL.mat';
    pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','identification','materialParticleBoard');
    load(fullfile(pathnameIdentification,filenameAna));
    load(fullfile(pathnameIdentification,filenameNum));
    
    % Material symmetry
    materialSym = 'isotTrans';
    
    switch lower(materialSym)
        case 'isot'
            % Young modulus
            E = mean(mean_ET_data)*1e6; % [Pa]
            % Shear modulus
            %G = mean(mean_GL_data)*1e6*13; % [Pa]
            % Poisson ratio
            %NU = E./(2*G)-1;
            NU = 0.25;
            % Material
            mat = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',h,'k',5/6);
        case 'isottrans'
            % Transverse Young modulus
            ET = mean(mean_ET_data)*1e6; % [Pa]
            % Longitudinal shear modulus
            GL = mean(mean_GL_data)*1e6; % [Pa]
            % Longitudinal Young modulus
            % EL = mean(mean_EL_data)*1e6; % [Pa]
            % Longitudinal Poisson ratio
            % NUL = mean(mean_NUL_data);
            % Transverse Poisson ratio
            NUT = 0.25;
            % Material
            mat = ELAS_SHELL_ISOT_TRANS('ET',ET,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',h,'k',5/6);
        otherwise
            error('Wrong material symmetry !')
    end
    
    mat = setnumber(mat,1);
    if junction
        S = union(S1,S2,'duplicate');
    else
        S = union(S1,S2);
    end
    S = setmaterial(S,mat);
    
    %% Neumann boundary conditions
    if selfweight
        p_plate = RHO*g*h; % surface load (body load for plates) [N/m2]
    else
        p_plate = 0;
    end
    
    L_load = LIGNE([L2,0,L1],[L2,b,L1]);
    junction_type = 'S1';
    
    switch lower(junction_type)
        case 's1'
            p = [24 65 114 163 200 252 291];
        case 's2'
            p = [14 54 113 159 199 249 299];
        case 's3'
            p = [15 65 111 153 203 243 295 341];
        case 's4'
            p = [23 34 91 141 180 229 290];
        case 'd1'
            p = [40 46 101 146];
        case 'd2'
            p = [6 63 110 175];
    end
    p = 100; % pointwise load [N]
    % p = 130;
    % p = 160;
    p = p(1)/b; % line load (surface load for plates) [N/m]
    
    %% Dirichlet boundary conditions
    L1 = getedge(Q1,1);
    L3 = getedge(Q1,3);
    if junction
        S = final(S,'duplicate');
        [~,numnode3] = intersect(S,L3);
        numnode = cell(getnbgroupelem(S),1);
        xnode = cell(getnbgroupelem(S),1);
        for i=1:getnbgroupelem(S)
            numnode{i} = double(getnumnodeingroupelem(S,i));
            numnode{i} = intersect(numnode3,numnode{i});
            xnode{i} = double(getcoord(getnode(S),numnode{i}));
        end
        for i=1:size(xnode{1},2)
            [xnode{1},I1] = sortrows(xnode{1},i);
            numnode{1} = numnode{1}(I1);
            [xnode{2},I2] = sortrows(xnode{2},i);
            numnode{2} = numnode{2}(I2);
        end
        S = addclperiodic(S,numnode{1},numnode{2},{'U','RX','RZ'});
    else
        S = final(S);
    end
    [~,numnode1] = intersect(S,L1);
    S = addcl(S,numnode1);
    
    %% Stiffness matrix and sollicitation vector
    A = calc_rigi(S);
    f = surfload(S,L_load,'FZ',-p);
    f = f + bodyload(S,[],'FZ',-p_plate);
    
    if junction
        k = 1e2/length(numnode{1}); % additonal junction rotational stiffness
        numddl1 = findddl(S,'RY',numnode{1},'free');
        numddl2 = findddl(S,'RY',numnode{2},'free');
        numddl = [numddl1 numddl2];
        A_add = [k -k;-k k];
        for i=1:size(numddl,1)
            A(numddl(i,:),numddl(i,:)) = A(numddl(i,:),numddl(i,:)) + A_add;
        end
    end
    
    %% Solution
    t = tic;
    u = A\f;
    time = toc(t);
    
    u = unfreevector(S,u);
    
    % U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    Ux = u(findddl(S,'UX'),:);
    Uy = u(findddl(S,'UY'),:);
    Uz = u(findddl(S,'UZ'),:);
    
    % R = u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    Rx = u(findddl(S,'RX'),:);
    Ry = u(findddl(S,'RY'),:);
    Rz = u(findddl(S,'RZ'),:);
    
    x = getcoord(S.node);
    
    e = calc_epsilon(S,u,'smooth');
    s = calc_sigma(S,u,'smooth');
    
    Exx = e(1);
    Eyy = e(2);
    Exy = e(3);
    Gxx = e(4);
    Gyy = e(5);
    Gxy = e(6);
    if strcmp(elemtype,'DST') || strcmp(elemtype,'DSQ') || strcmp(elemtype,'COQ4')
        Gzx = e(7);
        Gzy = e(8);
    else
        Gzx = 0;
        Gzy = 0;
    end
    
    Nxx = s(1);
    Nyy = s(2);
    Nxy = s(3);
    Mxx = s(4);
    Myy = s(5);
    Mxy = s(6);
    if strcmp(elemtype,'DST') || strcmp(elemtype,'DSQ') || strcmp(elemtype,'COQ4')
        Qx = s(7);
        Qy = s(8);
    else
        Qx = 0;
        Qy = 0;
    end
    
    %% Test solution
    [~,numnode3] = intersect(S,L3);
    % xP3 = x(numnode3(2),:);
    % P3 = POINT(xP3);
    P3 = getcenter(L3);
    
    ux = eval_sol(S,u,P3,'UX');
    uy = eval_sol(S,u,P3,'UY');
    uz = eval_sol(S,u,P3,'UZ');
    rx = eval_sol(S,u,P3,'RX');
    rz = eval_sol(S,u,P3,'RZ');
    if junction
        % ux = eval_sol(S.groupelem{1},S.node,u,P3,'UX'); % ux = eval_sol(S.groupelem{2},S.node,u,P3,'UX');
        % uy = eval_sol(S.groupelem{1},S.node,u,P3,'UY'); % uy = eval_sol(S.groupelem{2},S.node,u,P3,'UY');
        % uz = eval_sol(S.groupelem{1},S.node,u,P3,'UZ'); % uz = eval_sol(S.groupelem{2},S.node,u,P3,'UZ');
        % rx = eval_sol(S.groupelem{1},S.node,u,P3,'RX'); % rx = eval_sol(S.groupelem{2},S.node,u,P3,'RX');
        % rz = eval_sol(S.groupelem{1},S.node,u,P3,'RZ'); % rz = eval_sol(S.groupelem{2},S.node,u,P3,'RZ');
        ry = [eval_sol(S.groupelem{1},S.node,u,P3,'RY') ...
            eval_sol(S.groupelem{2},S.node,u,P3,'RY')];
        ry = double(ry);
    else
        ry = eval_sol(S,u,P3,'RY');
    end
    
    nxx = max(abs(Nxx));
    nyy = max(abs(Nyy));
    nxy = max(abs(Nxy));
    mxx = max(abs(Mxx));
    myy = max(abs(Myy));
    mxy = max(abs(Mxy));
    qx = max(abs(Qx));
    qy = max(abs(Qy));
    
    exx = max(abs(Exx));
    eyy = max(abs(Eyy));
    exy = max(abs(Exy));
    gxx = max(abs(Gxx));
    gyy = max(abs(Gyy));
    gxy = max(abs(Gxy));
    gzx = max(abs(Gzx));
    gzy = max(abs(Gzy));
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','elemtype',...
        'L1','L2','b','h',...
        'f','p','junction','junction_type');
    save(fullfile(pathname,'solution.mat'),'u','s','e','time',...
        'Ux','Uy','Uz','Rx','Ry','Rz',...
        'Nxx','Nyy','Nxy','Mxx','Myy','Mxy','Qx','Qy',...
        'Exx','Eyy','Exy','Gxx','Gyy','Gxy','Gzx','Gzy');
    save(fullfile(pathname,'test_solution.mat'),'P3',...
        'ux','uy','uz','rx','ry','rz',...
        'nxx','nyy','nxy','mxx','myy','mxy','qx','qy',...
        'exx','eyy','exy','gxx','gyy','gxy','gzx','gzy');
else
    load(fullfile(pathname,'problem.mat'),'S','elemtype',...
        'L1','L2','b','h',...
        'f','p','junction','junction_type');
    load(fullfile(pathname,'solution.mat'),'u','s','e','time',...
        'Ux','Uy','Uz','Rx','Ry','Rz',...
        'Nxx','Nyy','Nxy','Mxx','Myy','Mxy','Qx','Qy',...
        'Exx','Eyy','Exy','Gxx','Gyy','Gxy','Gzx','Gzy');
    load(fullfile(pathname,'test_solution.mat'),'P3',...
        'ux','uy','uz','rx','ry','rz',...
        'nxx','nyy','nxy','mxx','myy','mxy','qx','qy',...
        'exx','eyy','exy','gxx','gyy','gxy','gzx','gzy');
end

%% Outputs
fprintf('\nJunction %s\n',junction_type);
fprintf(['mesh : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

fprintf('Displacement u and rotation r at point (%g,%g,%g) m\n',double(P3));
fprintf('ux    = %g m\n',ux);
fprintf('uy    = %g m\n',uy);
fprintf('uz    = %g m\n',uz);
fprintf('rx    = %g rad = %g deg\n',rx,rad2deg(rx))
if junction
    fprintf('ry(1) = %g rad = %g deg\n',ry(1),rad2deg(ry(1)));
    fprintf('ry(2) = %g rad = %g deg\n',ry(2),rad2deg(ry(2)));
    fprintf('|ry(1) - ry(2)| = %g rad = %g deg\n',abs(ry(1)-ry(2)),rad2deg(abs(ry(1)-ry(2))));
else
    fprintf('ry    = %g rad = %g deg\n',ry,rad2deg(ry));
end
fprintf('rz    = %g rad = %g deg\n',rz,rad2deg(rz));
fprintf('\n');

if strcmp(elemtype,'DST') || strcmp(elemtype,'DSQ') || strcmp(elemtype,'COQ4')
    disp('Maximum membrane forces Nxx, Nyy, Nxy, moments Mxx, Myy, Mxy and shear forces Qx, Qy');
else
    disp('Maximum forces Nxx, Nyy, Nxy and moments Mxx, Myy, Mxy');
end
fprintf('Nxx = %g N/m\n',nxx);
fprintf('Nyy = %g N/m\n',nyy);
fprintf('Nxy = %g N/m\n',nxy);
fprintf('Mxx = %g N\n',mxx);
fprintf('Myy = %g N\n',myy);
fprintf('Mxy = %g N\n',mxy);
if strcmp(elemtype,'DST') || strcmp(elemtype,'DSQ') || strcmp(elemtype,'COQ4')
    fprintf('Qx = %g N/m\n',qx);
    fprintf('Qy = %g N/m\n',qy);
end
fprintf('\n');

if strcmp(elemtype,'DST') || strcmp(elemtype,'DSQ') || strcmp(elemtype,'COQ4')
    disp('Maximum membrane strains Exx, Eyy, Exy, bending strains (curvatures) Gxx, Gyy, Gxy and shear strains Gzx, Gzy');
else
    disp('Maximum membrane strains Exx, Eyy, Exy and bending strains (curvatures) Gxx, Gyy, Gxy');
end
fprintf('Exx = %g\n',exx);
fprintf('Eyy = %g\n',eyy);
fprintf('Exy = %g\n',exy);
fprintf('Gxx = %g\n',gxx);
fprintf('Gyy = %g\n',gyy);
fprintf('Gxy = %g\n',gxy);
if strcmp(elemtype,'DST') || strcmp(elemtype,'DSQ') || strcmp(elemtype,'COQ4')
    fprintf('Gzx = %g\n',gzx);
    fprintf('Gzy = %g\n',gzy);
end
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain(S,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 3;
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
    hP = plot(P3,'g+');
    legend([hD,hN,hP],[legD,legN,'measure'],'Location','NorthEastOutside')
    %legend([hD,hN,hP],[legD,legN,'mesure'],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    ampl = getsize(S)/max(abs(U))/20;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*u,'Color','b','FaceColor','b','FaceAlpha',0.1);
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
    
    plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
    mysaveas(pathname,'Uz',formats,renderer);
    
    % Rotations
    plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
    mysaveas(pathname,'Rx',formats,renderer);
    
    plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
    mysaveas(pathname,'Ry',formats,renderer);
    
    plotSolution(S,u,'rotation',3,'ampl',ampl,options{:});
    mysaveas(pathname,'Rz',formats,renderer);
    
%     % Membrane strains
%     plotSolution(S,u,'epsilon',1,'ampl',ampl,options{:});
%     mysaveas(pathname,'Exx',formats,renderer);
%     
%     plotSolution(S,u,'epsilon',2,'ampl',ampl,options{:});
%     mysaveas(pathname,'Eyy',formats,renderer);
%     
%     plotSolution(S,u,'epsilon',3,'ampl',ampl,options{:});
%     mysaveas(pathname,'Exy',formats,renderer);
% 
%     % Bending strains (curvatures)
%     plotSolution(S,u,'epsilon',4,'ampl',ampl,options{:});
%     mysaveas(pathname,'Gxx',formats,renderer);
%     
%     plotSolution(S,u,'epsilon',5,'ampl',ampl,options{:});
%     mysaveas(pathname,'Gyy',formats,renderer);
%     
%     plotSolution(S,u,'epsilon',6,'ampl',ampl,options{:});
%     mysaveas(pathname,'Gxy',formats,renderer);
%     
%     if strcmp(elemtype,'DST') || strcmp(elemtype,'DSQ') || strcmp(elemtype,'COQ4')
%         % Shear strains
%         plotSolution(S,u,'epsilon',7,'ampl',ampl,options{:});
%         mysaveas(pathname,'Gzx',formats,renderer);
%         plotSolution(S,u,'epsilon',8,'ampl',ampl,options{:});
%         mysaveas(pathname,'Gzy',formats,renderer);
%     end
% 
%     % Membrane forces
%     plotSolution(S,u,'sigma',1,'ampl',ampl,options{:});
%     mysaveas(pathname,'Nxx',formats,renderer);
%     
%     plotSolution(S,u,'sigma',2,'ampl',ampl,options{:});
%     mysaveas(pathname,'Nyy',formats,renderer);
%     
%     plotSolution(S,u,'sigma',3,'ampl',ampl,options{:});
%     mysaveas(pathname,'Nxy',formats,renderer);
%     
%     % Moments
%     plotSolution(S,u,'sigma',4,'ampl',ampl,options{:});
%     mysaveas(pathname,'Mxx',formats,renderer);
%     
%     plotSolution(S,u,'sigma',5,'ampl',ampl,options{:});
%     mysaveas(pathname,'Myy',formats,renderer);
%     
%     plotSolution(S,u,'sigma',6,'ampl',ampl,options{:});
%     mysaveas(pathname,'Mxy',formats,renderer);
% 
%     if strcmp(elemtype,'DST') || strcmp(elemtype,'DSQ') || strcmp(elemtype,'COQ4')
%         % Shear forces
%         plotSolution(S,u,'sigma',7,'ampl',ampl,options{:});
%         mysaveas(pathname,'Qx',formats,renderer);
%         plotSolution(S,u,'sigma',8,'ampl',ampl,options{:});
%         mysaveas(pathname,'Qy',formats,renderer);
%     end
    
    % Membrane strains
    figure('Name','Solution Exx')
    clf
    plot(e,S+ampl*u,'compo','EXX',options{:})
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Exx',formats,renderer);
    
    figure('Name','Solution Eyy')
    clf
    plot(e,S+ampl*u,'compo','EYY',options{:})
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Eyy',formats,renderer);
    
    figure('Name','Solution Exy')
    clf
    plot(e,S+ampl*u,'compo','EXY',options{:})
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Exy',formats,renderer);
    
    % Bending strains (curvatures)
    figure('Name','Solution Gxx')
    clf
    plot(e,S+ampl*u,'compo','GXX',options{:})
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Gxx',formats,renderer);
    
    figure('Name','Solution Gyy')
    clf
    plot(e,S+ampl*u,'compo','GYY',options{:})
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Gyy',formats,renderer);
    
    figure('Name','Solution Gxy')
    clf
    plot(e,S+ampl*u,'compo','GXY',options{:})
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Gxy',formats,renderer);
    
    if strcmp(elemtype,'DST') || strcmp(elemtype,'DSQ') || strcmp(elemtype,'COQ4')
        % Shear strains
        figure('Name','Solution Gzx')
        clf
        plot(e,S+ampl*u,'compo','GZX',options{:})
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,'Gzx',formats,renderer);
        
        figure('Name','Solution Gzy')
        clf
        plot(e,S+ampl*u,'compo','GZY',options{:})
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,'Gzy',formats,renderer);
    end
    
    % Membrane forces
    figure('Name','Solution Nxx')
    clf
    plot(s,S+ampl*u,'compo','NXX',options{:})
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Nxx',formats,renderer);
    
    figure('Name','Solution Nyy')
    clf
    plot(s,S+ampl*u,'compo','NYY',options{:})
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Nyy',formats,renderer);
    
    figure('Name','Solution Nxy')
    clf
    plot(s,S+ampl*u,'compo','NXY',options{:})
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Nxy',formats,renderer);
    
    % Moments
    figure('Name','Solution Mxx')
    clf
    plot(s,S+ampl*u,'compo','MXX',options{:})
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Mxx',formats,renderer);
    
    figure('Name','Solution Myy')
    clf
    plot(s,S+ampl*u,'compo','MYY',options{:})
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Myy',formats,renderer);
    
    figure('Name','Solution Mxy')
    clf
    plot(s,S+ampl*u,'compo','MXY',options{:})
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Mxy',formats,renderer);
    
    if strcmp(elemtype,'DST') || strcmp(elemtype,'DSQ') || strcmp(elemtype,'COQ4')
        % Shear forces
        figure('Name','Solution Qx')
        clf
        plot(s,S+ampl*u,'compo','QX',options{:})
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,'Qx',formats,renderer);
        
        figure('Name','Solution Qy')
        clf
        plot(s,S+ampl*u,'compo','QY',options{:})
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,'Qy',formats,renderer);
    end
end
