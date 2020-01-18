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
        L-L1-c/2+h2/2,0.0,H2+h1+d;
        L-L1-c/2+h2/2,0.0,H2+h1+d+h1;
        L-c,0.0,H2+h1+d;
        L-c,0.0,H2+h1+d+h1;
        L-c,l-c,H2+h1+d;
        L-c,l-c,H2+h1+d+h1;
        0.0,l-c,H2+h1+d;
        0.0,l-c,H2+h1+d+h1];
    x_barrier2 = [0.0,0.0,H2+2*(h1+d);
        0.0,0.0,H2+2*(h1+d)+h1;
        L-L1-c/2-h2/2,0.0,H2+2*(h1+d);
        L-L1-c/2-h2/2,0.0,H2+2*(h1+d)+h1;
        L-L1-c/2+h2/2,0.0,H2+2*(h1+d);
        L-L1-c/2+h2/2,0.0,H2+2*(h1+d)+h1;
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
        Q_lath{i} = QUADRANGLE(x_lath(4*i-3,:),x_lath(4*i-2,:),x_lath(4*i,:),x_lath(4*i-1,:));
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
    S_post = convertelem(S_post,'BEAM','param',VECTEUR([1;0;0]));
    
    % Plate meshes
    elemtype = 'DKT';
    S_plate = cellfun(@(Q,n) build_model(Q,'cl',cl,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_plate_' num2str(n)])),Q_plate,num2cell(1:length(Q_plate)),'UniformOutput',false);
    S_plate = union(S_plate{:});
    S_plate = concatgroupelem(S_plate);
    
    S_bed = cell(1,5);
    for n=1:3
        S_bed{n} = gmsh(Q_bed{n},P,cl,fullfile(pathname,['gmsh_bed_' num2str(n)]));
        S_bed{n} = convertelem(S_bed{n},elemtype);
    end
    for n=4:5
        S_bed{n} = build_model(Q_bed{n},'cl',cl,'elemtype',elemtype,...
            'filename',fullfile(pathname,['gmsh_bed_' num2str(n) '_elemtype_' elemtype]));
    end
    S_bed = union(S_bed{:});
    
    S_barrier1 = cellfun(@(Q,n) build_model(Q,'cl',cl,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_barrier1_' num2str(n)])),Q_barrier1,num2cell(1:length(Q_barrier1)),'UniformOutput',false);
    S_barrier2 = cell(1,4);
    for n=1:4
        if n==1
            S_barrier2{n} = build_model(Q_barrier2{n},'cl',cl,'elemtype',elemtype,...
                'filename',fullfile(pathname,['gmsh_barrier2_' num2str(n)]),'points',x_load);
        else
            S_barrier2{n} = build_model(Q_barrier2{n},'cl',cl,'elemtype',elemtype,...
                'filename',fullfile(pathname,['gmsh_barrier2_' num2str(n)]));
        end
    end
    S_barrier = union(S_barrier1{:},S_barrier2{:});
    S_barrier = concatgroupelem(S_barrier);
    S_barrier = convertelem(S_barrier,elemtype);
    
    S_frame = gmshbeam(P_frame,cl,fullfile(pathname,'gmsh_frame'));
    S_frame = concatgroupelem(S_frame);
    S_frame = convertelem(S_frame,elemtype);
    
    S_lath = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_lath_' num2str(n)])),P_lath,num2cell(1:length(P_lath)),'UniformOutput',false);
    S_lath = union(S_lath{:});
    S_lath = concatgroupelem(S_lath);
    S_lath = convertelem(S_lath,elemtype);
    
    %% Materials
    % Gravitational acceleration
    g = 9.81; % m/s2
    
    % Density
    RHO = 800; % kg/m3
    
    % Cross-section area
    Sec0 = c^2;
    % Planar second moment of area (or Planar area moment of inertia)
    IY0 = c^4/12;
    IZ0 = IY0;
    % Polar second moment of area (or Polar area moment of inertia)
    IX0 = IY0+IZ0;
    
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
            mat_1 = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',b1,'k',5/6);
            mat_1 = setnumber(mat_1,2);
            mat_2 = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',b2,'k',5/6);
            mat_2 = setnumber(mat_2,3);
            mat_3 = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',b1+b2,'k',5/6);
            mat_3 = setnumber(mat_3,4);
            S_post = setmaterial(S_post,mat_0);
            S_beam = setmaterial(S_beam,mat_1);
            S_bed = setmaterial(S_bed,mat_1,1:3);
            S_bed = setmaterial(S_bed,mat_2,4:5);
            S_barrier = setmaterial(S_barrier,mat_1);
            S_frame = setmaterial(S_frame,mat_3,[1,3,5]);
            S_frame = setmaterial(S_frame,mat_2,[2,4]);
            S_lath = setmaterial(S_lath,mat_1);
        otherwise
            error('Wrong material symmetry !')
    end
    
    S_bed = concatgroupelem(S_bed);
    S_bed = union(S_bed,S_lath);
    S = union(S_post,S_beam,S_bed,S_barrier,S_frame);
    
    %% Neumann boundary conditions
    p0 = RHO*g*Sec0; % line load (body load for posts)
    p1 = RHO*g*Sec1; % line load (body load for beams 1)
    p2 = RHO*g*Sec2; % line load (body load for beams 2)
    p3 = RHO*g*Sec3; % line load (body load for beams 3)
    p4 = RHO*g*Sec4; % line load (body load for laths)
    switch lower(test)
        case 'staticvert'
            p = 200; % pointwise load, 200N
        case {'statichoriin','statichoriout'}
            p = 500; % pointwise load, 500N
    end
    
    %% Dirichlet boundary conditions
    S = final(S);
    P_bot = POINT(x_bot);
    S = addcl(S,P_bot);
    
    %% Stiffness matrix and sollicitation vector
    A = calc_rigi(S);
    P_load = POINT(x_load);
    switch lower(test)
        case 'staticvert'
            f = nodalload(S,P_load,'FZ',p);
        case 'statichoriin'
            f = nodalload(S,P_load,'FY',p);
        case 'statichoriout'
            f = nodalload(S,P_load,'FY',-p);
    end
    f = f + bodyload(keepgroupelem(S,getnumgroupelemwithfield(S,'material',mat_0)),[],'FZ',-p0);
    f = f + bodyload(keepgroupelem(S,getnumgroupelemwithfield(S,'material',mat_1)),[],'FZ',-p1);
    f = f + bodyload(keepgroupelem(S,getnumgroupelemwithfield(S,'material',mat_2)),[],'FZ',-p2);
    f = f + bodyload(keepgroupelem(S,getnumgroupelemwithfield(S,'material',mat_3)),[],'FZ',-p3);
    f = f + bodyload(keepgroupelem(S,getnumgroupelemwithfield(S,'material',mat_4)),[],'FZ',-p4);
    
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
    
    %% Test solution
    P = P_load;
    numnode = find(S.node==P);
    xP = x(numnode,:);
    
    ux = eval_sol(S,u,P,'UX');
    uy = eval_sol(S,u,P,'UY');
    uz = eval_sol(S,u,P,'UZ');
    rx = eval_sol(S,u,P,'RX');
    ry = eval_sol(S,u,P,'RY');
    rz = eval_sol(S,u,P,'RZ');
    
    n  = reshape(N{6},[getnbnode(S),1]);
    mx = reshape(Mx{6},[getnbnode(S),1]);
    my = reshape(My{6},[getnbnode(S),1]);
    mz = reshape(Mz{6},[getnbnode(S),1]);
    epsx = reshape(Epsx{6},[getnbnode(S),1]);
    gamx = reshape(Gamx{6},[getnbnode(S),1]);
    gamy = reshape(Gamy{6},[getnbnode(S),1]);
    gamz = reshape(Gamz{6},[getnbnode(S),1]);
    n = double(n(numnode));
    mx = double(mx(numnode));
    my = double(my(numnode));
    mz = double(mz(numnode));
    epsx = double(epsx(numnode));
    gamx = double(gamx(numnode));
    gamy = double(gamy(numnode));
    gamz = double(gamz(numnode));
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','test',...
        'H','L','L1','l','h','H1','H2','h1','h2','b1','b2','c','d','e',...
        'f');
    save(fullfile(pathname,'solution.mat'),'u','s','e','time',...
        'Ux','Uy','Uz','Rx','Ry','Rz','N','Mx','My','Mz',...
        'Epsx','Gamx','Gamy','Gamz');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz','rx','ry','rz','n','mx','my','mz',...
        'epsx','gamx','gamy','gamz');
else
    load(fullfile(pathname,'problem.mat'),'S','test',...
        'H','L','L1','l','h','H1','H2','h1','h2','b1','b2','c','d','e',...
        'f');
    load(fullfile(pathname,'solution.mat'),'u','s','e','time',...
        'Ux','Uy','Uz','Rx','Ry','Rz',...
        'Epsx','Gamx','Gamy','Gamz',...
        'N','Mx','My','Mz');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz','rx','ry','rz',...
        'epsx','gamx','gamy','gamz',...
        'n','mx','my','mz');
end

%% Outputs
fprintf('\nDesk\n');
fprintf(['test : ' test '\n']);
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('elapsed time = %f s\n',time);

disp('Displacement u and rotation r at point'); disp(P);
fprintf('ux = %g m\n',ux);
fprintf('uy = %g m\n',uy);
fprintf('uz = %g m\n',uz);
fprintf('rx = %g rad = %g deg\n',rx,rad2deg(rx));
fprintf('ry = %g rad = %g deg\n',ry,rad2deg(ry));
fprintf('rz = %g rad = %g deg\n',rz,rad2deg(rz));
fprintf('\n');

disp('Axial strain Epsx, torsion and bending strains (curvatures) Gamx, Gamy, Gamz at point'); disp(P);
fprintf('Epsx = %g\n',epsx);
fprintf('Gamx = %g\n',gamx);
fprintf('Gamy = %g\n',gamy);
fprintf('Gamz = %g\n',gamz);
fprintf('\n');

disp('Force N and moments Mx, My, Mz at point'); disp(P);
fprintf('N  = %g N\n',n);
fprintf('Mx = %g N.m\n',mx);
fprintf('My = %g N.m\n',my);
fprintf('Mz = %g N.m\n',mz);
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    figure('Name','Domain')
    clf
    h1 = plot(S,'selgroup',1,'color','r');
    hold on
    h2 = plot(S,'selgroup',2,'color','g');
    h3 = plot(S,'selgroup',3:5,'color','b');
    h4 = plot(S,'selgroup',6,'color','k');
    h5 = plot(S,'selgroup',7,'color','m');
    hold off
    set(gca,'FontSize',16)
    l = legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'post','beam','bed','barrier','frame','Location','NorthEastOutside');
    %set(l,'Interpreter','latex')
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
%     plot(S,'group')
%     plot(S,'mat')
    
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
