%% FCBA bed beam deterministic linear elasticity %%
%%-----------------------------------------------%%

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

% tests = {'StaticVertUp'}; % test under static vertical upward load
% tests = {'StaticVertDown'}; % test under static vertical downward load
% tests = {'StaticHoriIn'}; % test under static horizontal inward load
tests = {'StaticHoriOut'}; % test under static horizontal outward load
% tests = {'StaticVertUp','StaticVertDown','StaticHoriIn','StaticHoriOut'};
% tests = {'StaticVertDown','StaticHoriOut'};

for it=1:length(tests)
    test = tests{it};
    
filename = ['FCBABedBeamDetLinElas' test];
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
    % Beams dimensions
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
    x_botrail = [0.0,0.0,H1+h1/2;
        L-c,0.0,H1+h1/2;
        L-c,l-c,H1+h1/2;
        0.0,l-c,H1+h1/2];
    x_toprail = [0.0,0.0,H2+h1/2;
        L-L1-c/2-h2/2,0.0,H2+h1/2;
        L-c,0.0,H2+h1/2;
        L-c,l-c,H2+h1/2;
        0.0,l-c,H2+h1/2];
    x_botguardrail = [0.0,0.0,H2+h1+d+h1/2;
        L-L1-c/2-h2/2,0.0,H2+h1+d+h1/2;
        L-c,0.0,H2+h1+d+h1/2;
        L-c,l-c,H2+h1+d+h1/2;
        0.0,l-c,H2+h1+d+h1/2];
    x_topguardrail = [0.0,0.0,H2+2*(h1+d)+h1/2;
        L-L1-c/2-h2/2,0.0,H2+2*(h1+d)+h1/2;
        L-c,0.0,H2+2*(h1+d)+h1/2;
        L-c,l-c,H2+2*(h1+d)+h1/2;
        0.0,l-c,H2+2*(h1+d)+h1/2];
    x_slat = zeros(14*2,3);
    for i=1:14
        x_slat(2*i-1,:) = [b2/2+e+d/2+(i-1)*2*d,0.0,H2+h1/2];
        x_slat(2*i,:) = [b2/2+e+d/2+(i-1)*2*d,l-c,H2+h1/2];
    end
    
    x_load = [(L-L1-c/2-h2/2)/2,0.0,H2+2*(h1+d)+h1/2];
    
    P_leg{1} = {x_bot(1,:),x_botrail(1,:),x_toprail(1,:),x_botguardrail(1,:),x_topguardrail(1,:),x_top(1,:)};
    P_leg{2} = {x_bot(2,:),x_botrail(2,:),x_toprail(3,:),x_botguardrail(3,:),x_topguardrail(3,:),x_top(2,:)};
    P_leg{3} = {x_bot(3,:),x_botrail(3,:),x_toprail(4,:),x_botguardrail(4,:),x_topguardrail(4,:),x_top(3,:)};
    P_leg{4} = {x_bot(4,:),x_botrail(4,:),x_toprail(5,:),x_botguardrail(5,:),x_topguardrail(5,:),x_top(4,:)};
    
    P_botrail{1} = {x_botrail(2,:),x_botrail(3,:)};
    P_botrail{2} = {x_botrail(3,:),x_botrail(4,:)};
    P_botrail{3} = {x_botrail(4,:),x_botrail(1,:)};
    
    P_siderail{1} = {x_toprail(1,:),x_slat(1,:),x_slat(3,:),x_slat(5,:),x_slat(7,:),x_slat(9,:),...
        x_slat(11,:),x_slat(13,:),x_slat(15,:),x_slat(17,:),x_slat(19,:),x_slat(21,:),...
        x_toprail(2,:),x_slat(23,:),x_slat(25,:),x_slat(27,:),x_toprail(3,:)};
    P_siderail{2} = {x_toprail(5,:),x_slat(2,:),x_slat(4,:),x_slat(6,:),x_slat(8,:),x_slat(10,:),...
        x_slat(12,:),x_slat(14,:),x_slat(16,:),x_slat(18,:),x_slat(20,:),x_slat(22,:),...
        x_slat(24,:),x_slat(26,:),x_slat(28,:),x_toprail(4,:)};
    
    P_endrail{1} = {x_toprail(3,:),x_toprail(4,:)};
    P_endrail{2} = {x_toprail(5,:),x_toprail(1,:)};
    
    P_slat = cell(1,14);
    for i=1:14
        P_slat{i} = {x_slat(2*i-1,:),x_slat(2*i,:)};
    end
    
    P_botguardrail{1} = {x_botguardrail(1,:),x_botguardrail(2,:)};
    P_botguardrail{2} = {x_botguardrail(3,:),x_botguardrail(4,:)};
    P_botguardrail{3} = {x_botguardrail(4,:),x_botguardrail(5,:)};
    P_botguardrail{4} = {x_botguardrail(5,:),x_botguardrail(1,:)};
    
    P_topguardrail{1} = {x_topguardrail(1,:),x_load,x_topguardrail(2,:)};
    P_topguardrail{2} = {x_topguardrail(3,:),x_topguardrail(4,:)};
    P_topguardrail{3} = {x_topguardrail(4,:),x_topguardrail(5,:)};
    P_topguardrail{4} = {x_topguardrail(5,:),x_topguardrail(1,:)};
    
    P_guardrailsupport = {x_toprail(2,:),x_botguardrail(2,:),x_topguardrail(2,:)};
    
    % Beams meshes
    elemtype = 'BEAM';
    cl = b1/2;
    S_leg = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_leg_' num2str(n)])),P_leg,num2cell(1:length(P_leg)),'UniformOutput',false);
    S_leg = cellfun(@(S) concatgroupelem(S),S_leg,'UniformOutput',false);
    S_leg = union(S_leg{:});
    S_leg = concatgroupelem(S_leg);
    S_leg = convertelem(S_leg,elemtype,'param',VECTEUR([1;0;0]));
    
    S_botrail = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_botrail_' num2str(n)])),P_botrail,num2cell(1:length(P_botrail)),'UniformOutput',false);
    S_botrail = union(S_botrail{:});
    S_botrail = concatgroupelem(S_botrail);
    S_botrail = convertelem(S_botrail,elemtype,'param',VECTEUR([0;0;1]));
    
    S_siderail = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_siderail_' num2str(n)])),P_siderail,num2cell(1:length(P_siderail)),'UniformOutput',false);
    S_siderail = cellfun(@(S) concatgroupelem(S),S_siderail,'UniformOutput',false);
    S_siderail = union(S_siderail{:});
    S_siderail = concatgroupelem(S_siderail);
    S_siderail = convertelem(S_siderail,elemtype,'param',VECTEUR([0;0;1]));
    
    S_endrail = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_endrail_' num2str(n)])),P_endrail,num2cell(1:length(P_endrail)),'UniformOutput',false);
    S_endrail = union(S_endrail{:});
    S_endrail = concatgroupelem(S_endrail);
    S_endrail = convertelem(S_endrail,elemtype,'param',VECTEUR([0;0;1]));
    
    S_botguardrail = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_botguardrail_' num2str(n)])),P_botguardrail,num2cell(1:length(P_botguardrail)),'UniformOutput',false);
    S_topguardrail = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_topguardrail_' num2str(n)])),P_topguardrail,num2cell(1:length(P_topguardrail)),'UniformOutput',false);
    S_guardrail = union(S_botguardrail{:},S_topguardrail{:});
    S_guardrail = concatgroupelem(S_guardrail);
    S_guardrail = convertelem(S_guardrail,elemtype,'param',VECTEUR([0;0;1]));
    
    S_guardrailsupport = gmshbeam(P_guardrailsupport,cl,fullfile(pathname,'gmsh_guardrailsupport'));
    S_guardrailsupport = concatgroupelem(S_guardrailsupport);
    S_guardrailsupport = convertelem(S_guardrailsupport,elemtype,'param',VECTEUR([1;0;0]));
    
    S_slat = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_slat_' num2str(n)])),P_slat,num2cell(1:length(P_slat)),'UniformOutput',false);
    S_slat = union(S_slat{:});
    S_slat = concatgroupelem(S_slat);
    S_slat = convertelem(S_slat,elemtype,'param',VECTEUR([1;0;0]));
    
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
            E = 10e9; % Pa
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
            S_leg = setmaterial(S_leg,mat_0);
            S_botrail = setmaterial(S_botrail,mat_1);
            S_siderail = setmaterial(S_siderail,mat_1);
            S_endrail = setmaterial(S_endrail,mat_2);
            S_guardrail = setmaterial(S_guardrail,mat_1);
            S_guardrailsupport = setmaterial(S_guardrailsupport,mat_3);
            S_slat = setmaterial(S_slat,mat_4);
        otherwise
            error('Wrong material symmetry !')
    end
    S = union(S_leg,S_botrail,S_siderail,S_endrail,S_guardrail,S_guardrailsupport,S_slat);
    
    %% Neumann boundary conditions
    p0 = RHO*g*Sec0; % line load (body load for legs)
    p1 = RHO*g*Sec1; % line load (body load for bottom rails, side rails and guard rails)
    p2 = RHO*g*Sec2; % line load (body load for end rails)
    p3 = RHO*g*Sec3; % line load (body load for guardrail support)
    p4 = RHO*g*Sec4; % line load (body load for slats)
    switch lower(test)
        case {'staticvertup','staticvertdown'}
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
        case 'staticvertup'
            f = nodalload(S,P_load,'FZ',p);
        case 'staticvertdown'
            f = nodalload(S,P_load,'FZ',-p);
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
    
    [~,~,numgroupelem] = findelemwithnode(S,numnode);
    n  = reshape(N{numgroupelem},[getnbnode(S),1]);
    mx = reshape(Mx{numgroupelem},[getnbnode(S),1]);
    my = reshape(My{numgroupelem},[getnbnode(S),1]);
    mz = reshape(Mz{numgroupelem},[getnbnode(S),1]);
    epsx = reshape(Epsx{numgroupelem},[getnbnode(S),1]);
    gamx = reshape(Gamx{numgroupelem},[getnbnode(S),1]);
    gamy = reshape(Gamy{numgroupelem},[getnbnode(S),1]);
    gamz = reshape(Gamz{numgroupelem},[getnbnode(S),1]);
    
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
        'Ux','Uy','Uz','Rx','Ry','Rz',...
        'N','Mx','My','Mz','Epsx','Gamx','Gamy','Gamz');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz','rx','ry','rz',...
        'n','mx','my','mz','epsx','gamx','gamy','gamz');
else
    load(fullfile(pathname,'problem.mat'),'S','test',...
        'H','L','L1','l','h','H1','H2','h1','h2','b1','b2','c','d','e',...
        'f');
    load(fullfile(pathname,'solution.mat'),'u','s','e','time',...
        'Ux','Uy','Uz','Rx','Ry','Rz',...
        'N','Mx','My','Mz','Epsx','Gamx','Gamy','Gamz');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz','rx','ry','rz',...
        'n','mx','my','mz','epsx','gamx','gamy','gamz');
end

%% Outputs
fprintf('\nBed\n');
fprintf(['test : ' test '\n']);
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

disp('Displacement u and rotation r at point'); disp(P);
fprintf('ux = %g m\n',ux);
fprintf('uy = %g m\n',uy);
fprintf('uz = %g m\n',uz);
fprintf('rx = %g rad = %g deg\n',rx,rad2deg(rx));
fprintf('ry = %g rad = %g deg\n',ry,rad2deg(ry));
fprintf('rz = %g rad = %g deg\n',rz,rad2deg(rz));
fprintf('\n');

disp('Force N and moments Mx, My, Mz at point'); disp(P);
fprintf('N  = %g N\n',n);
fprintf('Mx = %g N.m\n',mx);
fprintf('My = %g N.m\n',my);
fprintf('Mz = %g N.m\n',mz);
fprintf('\n');

disp('Axial strain Epsx, torsion and bending strains (curvatures) Gamx, Gamy, Gamz at point'); disp(P);
fprintf('Epsx = %g\n',epsx);
fprintf('Gamx = %g\n',gamx);
fprintf('Gamy = %g\n',gamy);
fprintf('Gamz = %g\n',gamz);
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    figure('Name','Domain')
    clf
    h1 = plot(S,'selgroup',1,'EdgeColor','k');
    hold on
    h2 = plot(S,'selgroup',2,'EdgeColor','c');
    h3 = plot(S,'selgroup',3,'EdgeColor','r');
    h4 = plot(S,'selgroup',4,'EdgeColor',[1 0.5 0]);
    h5 = plot(S,'selgroup',5,'EdgeColor','b');
    h6 = plot(S,'selgroup',6,'EdgeColor','m');
    h7 = plot(S,'selgroup',7,'EdgeColor','g');
    hold off
    set(gca,'FontSize',16)
    l = legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1),h7(1)],'leg','bottom rail','side rail','end rail','guard rail','guard rail support','slat','Location','NorthEastOutside');
    %set(l,'Interpreter','latex')
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
%     plotparamelem(S,'group')
%     plotparamelem(S,'material')
    
%     plotDomain(S,'legend',false);
%     mysaveas(pathname,'domain',formats,renderer);
%     mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'FaceColor','k','legend',false);
    ampl = 5;
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
    
    plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
    mysaveas(pathname,'Uz',formats,renderer);
    
    % Rotations
    plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
    mysaveas(pathname,'Rx',formats,renderer);
    
    plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
    mysaveas(pathname,'Ry',formats,renderer);
    
    plotSolution(S,u,'rotation',3,'ampl',ampl,options{:});
    mysaveas(pathname,'Rz',formats,renderer);
    
    % DO NOT WORK WITH BEAM ELEMENTS
    % Strains
    % plotSolution(S,u,'epsilon',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'Epsx',formats,renderer);
    %
    % plotSolution(S,u,'epsilon',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'Gamx',formats,renderer);
    %
    % plotSolution(S,u,'epsilon',3,'ampl',ampl,options{:});
    % mysaveas(pathname,'Gamy',formats,renderer);
    %
    % plotSolution(S,u,'epsilon',4,'ampl',ampl,options{:});
    % mysaveas(pathname,'Gamz',formats,renderer);
    %
    % Stresses
    % plotSolution(S,u,'sigma',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'N',formats,renderer);
    %
    % plotSolution(S,u,'sigma',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'Mx',formats,renderer);
    %
    % plotSolution(S,u,'sigma',3,'ampl',ampl,options{:});
    % mysaveas(pathname,'My',formats,renderer);
    %
    % plotSolution(S,u,'sigma',4,'ampl',ampl,options{:});
    % mysaveas(pathname,'Mz',formats,renderer);
    
    % Strains
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
    
    % Stresses
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
