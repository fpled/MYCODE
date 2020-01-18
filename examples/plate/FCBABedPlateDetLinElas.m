%% FCBA bed plate deterministic linear elasticity %%
%%-----------------------------------------------%%

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

% tests = {'StaticVert'}; % test under static vertical load
% tests = {'StaticHoriIn'}; % test under static horizontal inward load
tests = {'StaticHoriOut'}; % test under static horizontal outward load
% tests = {'StaticVert','StaticHoriIn','StaticHoriOut'};

for it=1:length(tests)
    test = tests{it};
    
filename = ['FCBABedPlateDetLinElas' test];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','plate',filename);
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
    % Beams and plates dimensions
    H = 1940e-3; % m
    L = 1990e-3;
    L1 = 445e-3;
    l = 990e-3;
    h = 410e-3;
    H1 = 620e-3;
    H2 = 1530e-3;
    b1 = 20e-3;
    b2 = 27e-3;
    h1 = 90e-3;
    h2 = 65e-3;
    c = 50e-3;
    d = 70e-3;
    e = 10e-3;
    
    % Points
    x_bot = [0.0,0.0,0.0;
        L-c,0.0,0.0;
        L-c,l-c,0.0;
        0.0,l-c,0.0];
    x_top = [0.0,0.0,H;
        L-c,0.0,H;
        L-c,l-c,H;
        0.0,l-c,H];
    x_mid = [0.0,0.0,H1;
        0.0,0.0,H1+h1;
        L-c,0.0,H1;
        L-c,0.0,H1+h1;
        L-c,l-c,H1;
        L-c,l-c,H1+h1;
        0.0,l-c,H1;
        0.0,l-c,H1+h1];
    x_bed = [0.0,0.0,H2;
        0.0,0.0,H2+h1;
        L-L1-c/2-h2,0.0,H2;
        L-L1-c/2,0.0,H2;
        L-L1-c/2-h2/2,0.0,H2+h1;
        L-L1-c/2-h2/2,0.0,H2+h1;
        L-c,0.0,H2;
        L-c,0.0,H2+h1;
        L-c,l-c,H2;
        L-c,l-c,H2+h1;
        0.0,l-c,H2;
        0.0,l-c,H2+h1];
    x_barrier1 = [0.0,0.0,H2+h1+d;
        0.0,0.0,H2+h1+d+h1;
        L-L1-c/2-h2/2,0.0,H2+h1+d;
        L-L1-c/2-h2/2,0.0,H2+h1+d+h1;
        L-c,0.0,H2+h1+d;
        L-c,0.0,H2+h1+d+h1;
        L-c,l-c,H2+h1+d;
        L-c,l-c,H2+h1+d+h1;
        0.0,l-c,H2+h1+d+h1/2];
    x_barrier2 = [0.0,0.0,H2+2*(h1+d);
        0.0,0.0,H2+2*(h1+d)+h1;
        L-L1-c/2-h2/2,0.0,H2+2*(h1+d);
        L-L1-c/2-h2/2,0.0,H2+2*(h1+d)+h1;
        L-c,0.0,H2+2*(h1+d);
        L-c,0.0,H2+2*(h1+d)+h1;
        L-c,l-c,H2+2*(h1+d);
        L-c,l-c,H2+2*(h1+d)+h1;
        0.0,l-c,H2+2*(h1+d);
        0.0,l-c,H2+2*(h1+d)+h1];
    x_lath = zeros(14*4,3);
    for i=1:14
        x_lath(4*i-3,:) = [b2/2+e+d/2+(i-1)*2*d,0.0,H2];
        x_lath(4*i-2,:) = [b2/2+e+d/2+(i-1)*2*d,0.0,H2+h1];
        x_lath(4*i-1,:) = [b2/2+e+d/2+(i-1)*2*d,l-c,H2];
        x_lath(4*i,:) = [b2/2+e+d/2+(i-1)*2*d,l-c,H2+h1];
    end
    
    x_load = [(L-L1-c/2-h2/2)/2,0.0,H2+2*(h1+d)+h1/2];
    
    P_post{1} = {x_bot(1,:),x_mid(1,:),x_mid(2,:),x_bed(1,:),x_bed(2,:),x_barrier1(1,:),x_barrier1(2,:),x_barrier2(1,:),x_barrier2(2,:),x_top(1,:)};
    P_post{2} = {x_bot(2,:),x_mid(3,:),x_mid(4,:),x_bed(7,:),x_bed(8,:),x_barrier1(7,:),x_barrier1(8,:),x_barrier2(7,:),x_barrier2(8,:),x_top(2,:)};
    P_post{3} = {x_bot(3,:),x_mid(5,:),x_mid(6,:),x_bed(9,:),x_bed(10,:),x_barrier1(9,:),x_barrier1(10,:),x_barrier2(9,:),x_barrier2(10,:),x_top(3,:)};
    P_post{4} = {x_bot(4,:),x_mid(7,:),x_mid(8,:),x_bed(11,:),x_bed(12,:),x_barrier1(11,:),x_barrier1(12,:),x_barrier2(11,:),x_barrier2(12,:),x_top(4,:)};
    
    Q_plate{1} = QUADRANGLE(x_mid(3,:),x_mid(4,:),x_mid(6,:),x_mid(5,:));
    Q_plate{2} = QUADRANGLE(x_mid(5,:),x_mid(6,:),x_mid(8,:),x_mid(7,:));
    Q_plate{3} = QUADRANGLE(x_mid(7,:),x_mid(8,:),x_mid(2,:),x_mid(1,:));
    
    Q_bed{1} = QUADRANGLE(x_bed(1,:),x_bed(2,:),x_bed(4,:),x_bed(3,:));
    Q_bed{2} = QUADRANGLE(x_bed(5,:),x_bed(6,:),x_bed(8,:),x_bed(7,:));
    Q_bed{3} = QUADRANGLE(x_bed(9,:),x_bed(10,:),x_bed(12,:),x_bed(11,:));
    Q_bed{4} = QUADRANGLE(x_bed(7,:),x_bed(8,:),x_bed(10,:),x_bed(9,:));
    Q_bed{5} = QUADRANGLE(x_bed(11,:),x_bed(12,:),x_bed(2,:),x_bed(1,:));
    
    Q_lath = cell(1,14);
    for i=1:14
        Q_lath{i} = QUADRANGLE(x_lath(4*i-3,:),x_bed(4*i-2,:),x_bed(4*i,:),x_bed(4*i-1,:));
    end
    
    Q_barrier1{1} = QUADRANGLE(x_barrier1(1,:),x_barrier1(2,:),x_barrier1(4,:),x_barrier1(3,:));
    Q_barrier1{2} = QUADRANGLE(x_barrier1(7,:),x_barrier1(8,:),x_barrier1(10,:),x_barrier1(9,:));
    Q_barrier1{3} = QUADRANGLE(x_barrier1(9,:),x_barrier1(10,:),x_barrier1(12,:),x_barrier1(11,:));
    Q_barrier1{4} = QUADRANGLE(x_barrier1(11,:),x_barrier1(12,:),x_barrier1(2,:),x_barrier1(1,:));
    
    Q_barrier2{1} = QUADRANGLE(x_barrier2(1,:),x_barrier2(2,:),x_barrier2(4,:),x_barrier2(3,:));
    Q_barrier2{2} = QUADRANGLE(x_barrier2(7,:),x_barrier2(8,:),x_barrier2(10,:),x_barrier2(9,:));
    Q_barrier2{3} = QUADRANGLE(x_barrier2(9,:),x_barrier2(10,:),x_barrier2(12,:),x_barrier2(11,:));
    Q_barrier2{4} = QUADRANGLE(x_barrier2(11,:),x_barrier2(12,:),x_barrier2(2,:),x_barrier2(1,:));
    
    Q_frame{1} = QUADRANGLE(x_bed(3,:),x_bed(4,:),x_bed(6,:),x_bed(5,:));
    Q_frame{2} = QUADRANGLE(x_bed(4,:),x_barrier1(3,:),x_barrier1(5,:),x_bed(6,:));
    Q_frame{3} = QUADRANGLE(x_barrier1(3,:),x_barrier1(4,:),x_barrier1(6,:),x_barrier1(5,:));
    Q_frame{4} = QUADRANGLE(x_barrier1(4,:),x_barrier2(3,:),x_barrier2(5,:),x_barrier1(6,:));
    Q_frame{5} = QUADRANGLE(x_barrier2(3,:),x_barrier2(4,:),x_barrier2(6,:),x_barrier2(5,:));
    
    % Beams meshes
    cl = b1;
    S_post = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_post_' num2str(n)])),P_post,num2cell(1:length(P_post)),'UniformOutput',false);
    S_post = cellfun(@(S) concatgroupelem(S),S_post,'UniformOutput',false);
    S_post = union(S_post{:});
    S_post = concatgroupelem(S_post);
    S_post = convertelem(S_post,elemtype,'param',VECTEUR([1;0;0]));
    
    % Plate meshes
    elemtype = 'DKT';
    S_plate = cellfun(@(Q,n) build_model(Q,'cl',cl,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_plate_' num2str(n)]),Q_plate,num2cell(1:length(Q_plate)),'UniformOutput',false);
    S_plate = union(S_plate{:});
    S_plate = concatgroupelem(S_plate);
    
    S_bed = cell(1,5);
    for n=1:3
        S_bed{n} = gmshbeam(P,cl,fullfile(pathname,['gmsh_bed_' num2str(n)]))
    end
    for n=4:5
        S_bed{n} = build_model(Q_bed{n},'cl',cl,'elemtype',elemtype,...
            'filename',fullfile(pathname,['gmsh_bed_' num2str(n) '_elemtype_' elemtype]));
    end
    
    S_bed = union(S_bed{:});
    
    S_beam1 = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_plate1_' num2str(n)])),P_beam1,num2cell(1:length(P_beam1)),'UniformOutput',false);
    S_beam2 = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_plate2_' num2str(n)])),P_beam2,num2cell(1:length(P_beam2)),'UniformOutput',false);
    S_beam3 = gmshbeam(P_beam3,cl,fullfile(pathname,['gmsh_beam3_cl_' num2str(cl)]));
    S_lath = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_lath_' num2str(n) '_cl_' num2str(cl)])),P_lath,num2cell(1:length(P_lath)),'UniformOutput',false);
    
    
    S_beam1 = cellfun(@(S) concatgroupelem(S),S_beam1,'UniformOutput',false);
    S_beam2 = cellfun(@(S) concatgroupelem(S),S_beam2,'UniformOutput',false);
    S_lath = cellfun(@(S) concatgroupelem(S),S_lath,'UniformOutput',false);
    
    S_beam1 = union(S_beam1{:});
    S_beam2 = union(S_beam2{:});
    S_lath = union(S_lath{:});
    
    S_beam1 = concatgroupelem(S_beam1);
    S_beam2 = concatgroupelem(S_beam2);
    S_beam3 = concatgroupelem(S_beam3);
    S_lath = concatgroupelem(S_lath);
    
    
    S_beam1 = convertelem(S_beam1,elemtype);
    S_beam2 = convertelem(S_beam2,elemtype);
    S_beam3 = convertelem(S_beam3,elemtype);
    S_lath = convertelem(S_lath,elemtype,'param',VECTEUR([1;0;0]));
    
    %% Materials
    % Gravitational acceleration
    g = 9.81; % m/s2
    
    % Density
    RHO = 800; % kg/m3
    
    % Cross-section area
    Sec0 = c^2;
    Sec1 = b1*h1;
    Sec2 = b2*h1;
    Sec3 = b2*h2;
    Sec4 = b1*d;
    % Planar second moment of area (or Planar area moment of inertia)
    IY0 = c^4/12;
    IY1 = h1*b1^3/12;
    IY2 = h1*b2^3/12;
    IY3 = h2*b2^3/12;
    IY4 = d*b1^3/12;
    IZ0 = IY0;
    IZ1 = b1*h1^3/12;
    IZ2 = b2*h1^3/12;
    IZ3 = b2*h2^3/12;
    IZ4 = b1*d^3/12;
    % Polar second moment of area (or Polar area moment of inertia)
    IX0 = IY0+IZ0;
    IX1 = IY1+IZ1;
    IX2 = IY2+IZ2;
    IX3 = IY3+IZ3;
    IX4 = IY4+IZ4;
    
    % Material symmetry
    materialSym = 'isot';
    
    switch lower(materialSym)
        case 'isot'
            % Young modulus
            E = 2e9; % Pa
            % Poisson ratio
            NU = 0.3;
            % Material
            mat_0 = ELAS_BEAM('E',E,'NU',NU,'S',Sec0,'IZ',IZ0,'IY',IY0,'IX',IX0,'RHO',RHO);
            mat_0 = setnumber(mat_0,1);
            mat_1 = ELAS_BEAM('E',E,'NU',NU,'S',Sec1,'IZ',IZ1,'IY',IY1,'IX',IX1,'RHO',RHO);
            mat_1 = setnumber(mat_1,2);
            mat_2 = ELAS_BEAM('E',E,'NU',NU,'S',Sec2,'IZ',IZ2,'IY',IY2,'IX',IX2,'RHO',RHO);
            mat_2 = setnumber(mat_2,3);
            mat_3 = ELAS_BEAM('E',E,'NU',NU,'S',Sec3,'IZ',IZ3,'IY',IY3,'IX',IX3,'RHO',RHO);
            mat_3 = setnumber(mat_3,4);
            mat_4 = ELAS_BEAM('E',E,'NU',NU,'S',Sec4,'IZ',IZ4,'IY',IY4,'IX',IX4,'RHO',RHO);
            mat_4 = setnumber(mat_4,5);
            S_post = setmaterial(S_post,mat_0);
            S_beam1 = setmaterial(S_beam1,mat_1);
            S_beam2 = setmaterial(S_beam2,mat_2);
            S_beam3 = setmaterial(S_beam3,mat_3);
            S_lath = setmaterial(S_lath,mat_4);
        otherwise
            error('Wrong material symmetry !')
    end
    
    S = union(S_post,S_beam1,S_beam2,S_beam3,S_lath);
    
    %% Neumann boundary conditions
    p0 = RHO*g*Sec0; % line load (body load for posts)
    p1 = RHO*g*Sec1; % line load (body load for beams 1)
    p2 = RHO*g*Sec2; % line load (body load for beams 2)
    p3 = RHO*g*Sec3; % line load (body load for beams 3)
    p4 = RHO*g*Sec4; % line load (body load for laths)
    p = 500; % pointwise load, 500N
    
    %% Dirichlet boundary conditions
    S = final(S);
    S = addcl(S,P_bot);
    
    %% Stiffness matrix and sollicitation vector
    A = calc_rigi(S);
    f = nodalload(S,P_load,'FY',-p);
    f = f + bodyload(keepgroupelem(S,1),[],'FZ',-p0);
    f = f + bodyload(keepgroupelem(S,2),[],'FZ',-p1);
    f = f + bodyload(keepgroupelem(S,3),[],'FZ',-p2);
    f = f + bodyload(keepgroupelem(S,4),[],'FZ',-p3);
    f = f + bodyload(keepgroupelem(S,5),[],'FZ',-p4);
    
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
    
    Epsx = e(1);
    Gamx = e(2);
    Gamy = e(3);
    Gamz = e(4);
    N = s(1);
    Mx = s(2);
    My = s(3);
    Mz = s(4);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S',...
        'H','L','L1','l','h','H1','H2','h1','h2','b1','b2','c','d','e',...
        'f');
    save(fullfile(pathname,'solution.mat'),'u','s','e','time',...
        'Ux','Uy','Uz','Rx','Ry','Rz','N','Mx','My','Mz',...
        'Epsx','Gamx','Gamy','Gamz');
else
    load(fullfile(pathname,'problem.mat'),'S',...
        'H','L','L1','l','h','H1','H2','h1','h2','b1','b2','c','d','e',...
        'f');
    load(fullfile(pathname,'solution.mat'),'u','s','e','time',...
        'Ux','Uy','Uz','Rx','Ry','Rz','N','Mx','My','Mz',...
        'Epsx','Gamx','Gamy','Gamz');
end

%% Outputs
fprintf('\nDesk\n');
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    figure('Name','Materials')
    clf
    h1 = plot(S,'selgroup',1,'color','r');
    hold on
    h2 = plot(S,'selgroup',2,'color','m');
    h3 = plot(S,'selgroup',3,'color','b');
    h4 = plot(S,'selgroup',4,'color','k');
    h5 = plot(S,'selgroup',5,'color','g');
    hold off
    mysaveas(pathname,'material',formats,renderer);
    
%     plotDomain(S,'legend',false);
%     mysaveas(pathname,'domain',formats,renderer);
%     mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'FaceColor','k','legend',false);
    ampl = 5;
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
    legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
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
    
    plotSolution(S,u,'displ',1,'ampl',ampl,options{:});
    mysaveas(pathname,'Ux',formats,renderer);
    
    plotSolution(S,u,'displ',2,'ampl',ampl,options{:});
    mysaveas(pathname,'Uy',formats,renderer);
    
    plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
    mysaveas(pathname,'Uz',formats,renderer);
    
    plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
    mysaveas(pathname,'Rx',formats,renderer);
    
    plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
    mysaveas(pathname,'Ry',formats,renderer);
    
    plotSolution(S,u,'rotation',3,'ampl',ampl,options{:});
    mysaveas(pathname,'Rz',formats,renderer);
    
    % DO NOT WORK WITH BEAM ELEMENTS
    % plotSolution(S,u,'epsilon',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'eps_x',formats,renderer);
    %
    % plotSolution(S,u,'epsilon',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'gam_z',formats,renderer);
    %
    % plotSolution(S,u,'sigma',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'eff_x',formats,renderer);
    %
    % plotSolution(S,u,'sigma',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'mom_z',formats,renderer);
    
    figure('Name','Solution eps_x')
    clf
    plot(e,S+ampl*u,'compo','EPSX')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Epsx',formats,renderer);
    
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
    
    figure('Name','Solution gam_z')
    clf
    plot(e,S+ampl*u,'compo','GAMZ')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Gamz',formats,renderer);
    
    figure('Name','Solution N')
    clf
    plot(s,S+ampl*u,'compo','EFFX')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'N',formats,renderer);
    
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
    
    figure('Name','Solution Mz')
    clf
    plot(s,S+ampl*u,'compo','MOMZ')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Mz',formats,renderer);
end

end
