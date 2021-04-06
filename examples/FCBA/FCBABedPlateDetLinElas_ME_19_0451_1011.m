%% FCBA bed plate deterministic linear elasticity ME 19-0451-1011 %%
%%----------------------------------------------------------------%%

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

junction = false; % junction modeling
materialSym = 'isot'; % isotropic material symmetry class
% materialSym = 'isotTrans'; % transversely isotropic material symmetry class
slat = true; % slat modeling

for it=1:length(tests)
    test = tests{it};
    
filename = ['FCBABedPlateDetLinElas_ME_19_0451_1011_' test '_' materialSym '_slat_' num2str(slat)];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','FCBA',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc','png'};
renderer = 'OpenGL';

%% Problem
if solveProblem
    %% Domains and meshes
    % Beams and plates dimensions
    H = 1940e-3; % [m]
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
    d_load = 100e-3;
    
    % Points
    x_bot = [0.0,0.0,0.0;
        L-c,0.0,0.0;
        L-c,l-c,0.0;
        0.0,l-c,0.0];
    x_top = [0.0,0.0,H;
        L-c,0.0,H;
        L-c,l-c,H;
        0.0,l-c,H];
    x_botrail = [0.0,0.0,H1;
        0.0,0.0,H1+h1;
        L-c,0.0,H1;
        L-c,0.0,H1+h1;
        L-c,l-c,H1;
        L-c,l-c,H1+h1;
        0.0,l-c,H1;
        0.0,l-c,H1+h1];
    x_toprail = [0.0,0.0,H2;
        0.0,0.0,H2+h1;
        L-c,0.0,H2;
        L-c,0.0,H2+h1;
        L-c,l-c,H2;
        L-c,l-c,H2+h1;
        0.0,l-c,H2;
        0.0,l-c,H2+h1];
    x_botguardrail = [0.0,0.0,H2+h1+d;
        0.0,0.0,H2+h1+d+h1;
        L-L1-c/2,0.0,H2+h1+d;
        L-L1-c/2,0.0,H2+h1+d+h1;
        L-c,0.0,H2+h1+d;
        L-c,0.0,H2+h1+d+h1;
        L-c,l-c,H2+h1+d;
        L-c,l-c,H2+h1+d+h1;
        0.0,l-c,H2+h1+d;
        0.0,l-c,H2+h1+d+h1];
    x_topguardrail = [0.0,0.0,H2+2*(h1+d);
        0.0,0.0,H2+2*(h1+d)+h1;
        L-L1-c/2,0.0,H2+2*(h1+d);
        L-L1-c/2,0.0,H2+2*(h1+d)+h1;
        L-c,0.0,H2+2*(h1+d);
        L-c,0.0,H2+2*(h1+d)+h1;
        L-c,l-c,H2+2*(h1+d);
        L-c,l-c,H2+2*(h1+d)+h1;
        0.0,l-c,H2+2*(h1+d);
        0.0,l-c,H2+2*(h1+d)+h1];
    x_guardrailsupport = [L-L1-c/2-h2,0.0,H2;
        L-L1-c/2-h2,0.0,H2+h1;
        L-L1-c/2-h2,0.0,H2+h1+d;
        L-L1-c/2-h2,0.0,H2+h1+d+h1;
        L-L1-c/2-h2,0.0,H2+2*(h1+d);
        L-L1-c/2-h2,0.0,H2+2*(h1+d)+h1;
        L-L1-c/2,0.0,H2;
        L-L1-c/2,0.0,H2+h1];
    x_slat = zeros(14*4,3);
    for i=1:14
        x_slat(4*i-3,:) = [b2/2+e+(i-1)*2*d,0.0,H2+h1/2];
        x_slat(4*i-2,:) = [b2/2+e+d+(i-1)*2*d,0.0,H2+h1/2];
        x_slat(4*i-1,:) = [b2/2+e+(i-1)*2*d,l-c,H2+h1/2];
        x_slat(4*i,:) = [b2/2+e+d+(i-1)*2*d,l-c,H2+h1/2];
    end
    x_cross = [L-L1-c/2-h2,0.0,H2+h1/2];
    x_load = [(L-c)/2,l-c,H2+2*(h1+d)+h1/2;
        d_load-c/2,l-c,H2+2*(h1+d)+h1/2];
    x_measure = [0.0,l-c,H2+2*(h1+d)+h1/2];
    
    P_leg{1} = {x_bot(1,:),x_botrail(1,:),x_botrail(2,:),x_toprail(1,:),x_toprail(2,:),x_botguardrail(1,:),x_botguardrail(2,:),x_topguardrail(1,:),x_top(1,:)};
    P_leg{2} = {x_bot(2,:),x_botrail(3,:),x_botrail(4,:),x_toprail(3,:),x_toprail(4,:),x_botguardrail(5,:),x_botguardrail(6,:),x_topguardrail(5,:),x_top(2,:)};
    P_leg{3} = {x_bot(3,:),x_botrail(5,:),x_botrail(6,:),x_toprail(5,:),x_toprail(6,:),x_botguardrail(7,:),x_botguardrail(8,:),x_topguardrail(7,:),x_top(3,:)};
    P_leg{4} = {x_bot(4,:),x_botrail(7,:),x_botrail(8,:),x_toprail(7,:),x_toprail(8,:),x_botguardrail(9,:),x_botguardrail(10,:),x_measure,x_topguardrail(9,:),x_top(4,:)};
    
    Q_botrail{1} = QUADRANGLE(x_botrail(3,:),x_botrail(4,:),x_botrail(6,:),x_botrail(5,:));
    Q_botrail{2} = QUADRANGLE(x_botrail(5,:),x_botrail(6,:),x_botrail(8,:),x_botrail(7,:));
    Q_botrail{3} = QUADRANGLE(x_botrail(7,:),x_botrail(8,:),x_botrail(2,:),x_botrail(1,:));
    
    Q_siderail{1} = QUADRANGLE(x_toprail(1,:),x_toprail(2,:),x_guardrailsupport(2,:),x_guardrailsupport(1,:));
    Q_siderail{2} = QUADRANGLE(x_guardrailsupport(7,:),x_guardrailsupport(8,:),x_toprail(4,:),x_toprail(3,:));
    Q_siderail{3} = QUADRANGLE(x_toprail(5,:),x_toprail(6,:),x_toprail(8,:),x_toprail(7,:));
    
    Q_endrail{1} = QUADRANGLE(x_toprail(3,:),x_toprail(4,:),x_toprail(6,:),x_toprail(5,:));
    Q_endrail{2} = QUADRANGLE(x_toprail(7,:),x_toprail(8,:),x_toprail(2,:),x_toprail(1,:));
    
    Q_slat = cell(1,14);
    for i=[1:10,12:14]
        Q_slat{i} = QUADRANGLE(x_slat(4*i-3,:),x_slat(4*i-2,:),x_slat(4*i,:),x_slat(4*i-1,:));
    end
    Q_slat{11} = POLYGON(x_slat(4*11-3,:),x_cross,x_slat(4*11-2,:),x_slat(4*11,:),x_slat(4*11-1,:));
    
    Q_botguardrail{1} = QUADRANGLE(x_botguardrail(1,:),x_botguardrail(2,:),x_guardrailsupport(4,:),x_guardrailsupport(3,:));
    Q_botguardrail{2} = QUADRANGLE(x_botguardrail(5,:),x_botguardrail(6,:),x_botguardrail(8,:),x_botguardrail(7,:));
    Q_botguardrail{3} = QUADRANGLE(x_botguardrail(7,:),x_botguardrail(8,:),x_botguardrail(10,:),x_botguardrail(9,:));
    Q_botguardrail{4} = QUADRANGLE(x_botguardrail(9,:),x_botguardrail(10,:),x_botguardrail(2,:),x_botguardrail(1,:));
    
    Q_topguardrail{1} = QUADRANGLE(x_topguardrail(1,:),x_topguardrail(2,:),x_guardrailsupport(6,:),x_guardrailsupport(5,:));
    Q_topguardrail{2} = QUADRANGLE(x_topguardrail(5,:),x_topguardrail(6,:),x_topguardrail(8,:),x_topguardrail(7,:));
    Q_topguardrail{3} = POLYGON(x_topguardrail(7,:),x_topguardrail(8,:),x_topguardrail(10,:),x_measure,x_topguardrail(9,:));
    Q_topguardrail{4} = POLYGON(x_topguardrail(9,:),x_measure,x_topguardrail(10,:),x_topguardrail(2,:),x_topguardrail(1,:));
    
    Q_guardrailsupport{1} = QUADRANGLE(x_guardrailsupport(1,:),x_guardrailsupport(2,:),x_guardrailsupport(8,:),x_guardrailsupport(7,:));
    Q_guardrailsupport{2} = QUADRANGLE(x_guardrailsupport(2,:),x_guardrailsupport(3,:),x_botguardrail(3,:),x_guardrailsupport(8,:));
    Q_guardrailsupport{3} = QUADRANGLE(x_guardrailsupport(3,:),x_guardrailsupport(4,:),x_botguardrail(4,:),x_botguardrail(3,:));
    Q_guardrailsupport{4} = QUADRANGLE(x_guardrailsupport(4,:),x_guardrailsupport(5,:),x_topguardrail(3,:),x_botguardrail(4,:));
    Q_guardrailsupport{5} = QUADRANGLE(x_guardrailsupport(5,:),x_guardrailsupport(6,:),x_topguardrail(4,:),x_topguardrail(3,:));
    
    % Beams meshes
    cl = b1/2;
    S_leg = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_leg_' num2str(n)])),P_leg,num2cell(1:length(P_leg)),'UniformOutput',false);
    S_leg = cellfun(@(S) concatgroupelem(S),S_leg,'UniformOutput',false);
    S_leg = union(S_leg{:});
    S_leg = convertelem(S_leg,'BEAM','param',VECTEUR([1;0;0]));
    
    % Plate meshes
    elemtype = 'DKT';
    S_botrail = cellfun(@(Q,n) build_model(Q,'cl',cl,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_botrail_' num2str(n) '_elemtype_' elemtype])),Q_botrail,num2cell(1:length(Q_botrail)),'UniformOutput',false);
    S_botrail = union(S_botrail{:});
    
    L_slat1 = cellfun(@(Q) getedge(Q,1),Q_slat(1:11),'UniformOutput',false);
    L_slat2 = cellfun(@(Q) getedge(Q,1),Q_slat(12:14),'UniformOutput',false);
    L_slat3 = cellfun(@(Q) getedge(Q,3),Q_slat,'UniformOutput',false);
    L_slat3{11} = getedge(Q_slat{11},4);
    S_siderail = cell(1,3);
    if ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4')
        S_siderail{1} = gmshFCBAbedsiderail(Q_siderail{1},L_slat1,cl,cl,fullfile(pathname,['gmsh_siderail_1_elemtype_' elemtype]),3);
        S_siderail{2} = gmshdomainwithinclusion(Q_siderail{2},L_slat2,cl,cl,fullfile(pathname,['gmsh_siderail_2_elemtype_' elemtype]),3);
        S_siderail{3} = gmshdomainwithinclusion(Q_siderail{3},L_slat3,cl,cl,fullfile(pathname,['gmsh_siderail_3_elemtype_' elemtype]),3);
    else
        S_siderail{1} = gmshFCBAbedsiderail(Q_siderail{1},L_slat1,cl,cl,fullfile(pathname,['gmsh_siderail_1_elemtype_' elemtype]),3,'recombine');
        S_siderail{2} = gmshdomainwithinclusion(Q_siderail{2},L_slat2,cl,cl,fullfile(pathname,['gmsh_siderail_2_elemtype_' elemtype]),3,'recombine');
        S_siderail{3} = gmshdomainwithinclusion(Q_siderail{3},L_slat3,cl,cl,fullfile(pathname,['gmsh_siderail_3_elemtype_' elemtype]),3,'recombine');
    end
    S_siderail = union(S_siderail{:});
    S_siderail = convertelem(S_siderail,elemtype);
    
    S_endrail = cellfun(@(Q,n) build_model(Q,'cl',cl,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_endrail_' num2str(n) '_elemtype_' elemtype])),Q_endrail,num2cell(1:length(Q_endrail)),'UniformOutput',false);
    S_endrail = union(S_endrail{:});
    
    S_botguardrail = cellfun(@(Q,n) build_model(Q,'cl',cl,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_botguardrail_' num2str(n) '_elemtype_' elemtype])),Q_botguardrail,num2cell(1:length(Q_botguardrail)),'UniformOutput',false);
    S_topguardrail = cell(1,4);
    S_topguardrail([1,2,4]) = cellfun(@(Q,n) build_model(Q,'cl',cl,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_topguardrail_' num2str(n) '_elemtype_' elemtype])),Q_topguardrail([1,2,4]),num2cell([1,2,4]),'UniformOutput',false);
    S_topguardrail{3} = build_model(Q_topguardrail{3},'cl',cl,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_topguardrail_3_elemtype_' elemtype]),'points',{x_load(1,:),x_load(2,:)});
    S_guardrail = union(S_botguardrail{:},S_topguardrail{:});
    
    L_slat = getedge(Q_slat{11},2);
    S_guardrailsupport = cell(1,5);
    if ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4')
        S_guardrailsupport{1} = gmshFCBAbedguardrailsupport(Q_guardrailsupport{1},L_slat,cl,cl,fullfile(pathname,['gmsh_guardrailsupport_1_elemtype_' elemtype]),3);
    else
        S_guardrailsupport{1} = gmshFCBAbedguardrailsupport(Q_guardrailsupport{1},L_slat,cl,cl,fullfile(pathname,['gmsh_guardrailsupport_1_elemtype_' elemtype]),3,'recombine');
    end
    S_guardrailsupport{1} = convertelem(S_guardrailsupport{1},elemtype);
    S_guardrailsupport(2:5) = cellfun(@(Q,n) build_model(Q,'cl',cl,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_guardrailsupport_' num2str(n) '_elemtype_' elemtype])),Q_guardrailsupport(2:5),num2cell(2:5),'UniformOutput',false);
    S_guardrailsupport = union(S_guardrailsupport{:});
    
    if slat
        S_slat = cellfun(@(Q,n) build_model(Q,'cl',cl,'elemtype',elemtype,...
            'filename',fullfile(pathname,['gmsh_slat_' num2str(n) '_elemtype_' elemtype])),Q_slat,num2cell(1:14),'UniformOutput',false);
        S_slat = union(S_slat{:});
    end
    
    %% Materials
    % Gravitational acceleration
    g = 9.81; % [m/s2]
    
    % Density
    RHO = 500; % [kg/m3]
    
    % Cross-section area
    Sec0 = c^2;
    % Planar second moment of area (or Planar area moment of inertia)
    IY0 = c^4/12;
    IZ0 = IY0;
    % Polar second moment of area (or Polar area moment of inertia)
    IX0 = IY0+IZ0;
    
    % Material
    switch lower(materialSym)
        case 'isot'
            % Young modulus
            E = 11e9; % [Pa]
            % Poisson ratio
            NU = 0.2;
            % Material
            mat_0 = ELAS_BEAM('E',E,'NU',NU,'S',Sec0,'IZ',IZ0,'IY',IY0,'IX',IX0,'RHO',RHO);
            mat_0 = setnumber(mat_0,1);
            mat_1 = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',b1,'k',5/6);
            mat_1 = setnumber(mat_1,2);
            mat_2 = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',b2,'k',5/6);
            mat_2 = setnumber(mat_2,3);
            mat_3 = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',b1+b2,'k',5/6);
            mat_3 = setnumber(mat_3,4);
            S_leg = setmaterial(S_leg,mat_0);
            S_botrail = setmaterial(S_botrail,mat_1);
            S_siderail = setmaterial(S_siderail,mat_1);
            S_endrail = setmaterial(S_endrail,mat_2);
            S_guardrail = setmaterial(S_guardrail,mat_1);
            S_guardrailsupport = setmaterial(S_guardrailsupport,mat_2,[2,4]);
            S_guardrailsupport = setmaterial(S_guardrailsupport,mat_3,[1,3,5]);
            if slat
                S_slat = setmaterial(S_slat,mat_1);
            end
        case 'isottrans'
            % Transverse Young modulus
            ET = 11e9; % [Pa]
            % Longitudinal shear modulus
            GL = 500e6; % [Pa]
            % Longitudinal Young modulus
            % EL = ET; % [Pa]
            % Longitudinal Poisson ratio
            % NUL = 0.2;
            % Transverse Poisson ratio
            NUT = 0.2;
            % Material
            mat_0 = ELAS_BEAM_ISOT_TRANS('EL',ET,'NUT',NUT,'GL',GL,'S',Sec0,'IZ',IZ0,'IY',IY0,'IX',IX0,'RHO',RHO);
            mat_0 = setnumber(mat_0,1);
            mat_1 = ELAS_SHELL_ISOT_TRANS('ET',ET,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',b1,'k',5/6);
            mat_1 = setnumber(mat_1,2);
            mat_2 = ELAS_SHELL_ISOT_TRANS('ET',ET,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',b2,'k',5/6);
            mat_2 = setnumber(mat_2,3);
            mat_3 = ELAS_SHELL_ISOT_TRANS('ET',ET,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',b1+b2,'k',5/6);
            mat_3 = setnumber(mat_3,4);
            S_leg = setmaterial(S_leg,mat_0);
            S_botrail = setmaterial(S_botrail,mat_1);
            S_siderail = setmaterial(S_siderail,mat_1);
            S_endrail = setmaterial(S_endrail,mat_2);
            S_guardrail = setmaterial(S_guardrail,mat_1);
            S_guardrailsupport = setmaterial(S_guardrailsupport,mat_2,[2,4]);
            S_guardrailsupport = setmaterial(S_guardrailsupport,mat_3,[1,3,5]);
            if slat
                S_slat = setmaterial(S_slat,mat_1);
            end
        otherwise
            error('Wrong material symmetry !')
    end
    if slat
        S_plate = union(S_botrail,S_siderail,S_endrail,S_guardrail,S_guardrailsupport,S_slat);
    else
        S_plate = union(S_botrail,S_siderail,S_endrail,S_guardrail,S_guardrailsupport);
    end
    S = union(S_leg,S_plate);
    
    %% Neumann boundary conditions
    p0 = RHO*g*Sec0; % line load (body load for legs) [N/m]
    p1 = RHO*g*b1; % surface load (body load for bottom rails, side rails, guard rails and slats) [N/m2]
    p2 = RHO*g*b2; % surface load (body load for end rails and parts 2 and 4 of guard rail support) [N/m2]
    p3 = RHO*g*(b1+b2); % surface load (body load for parts 1, 3 and 5 of guard rail support) [N/m2]
    p = 500; % pointwise load [N]
    
    %% Dirichlet boundary conditions
    S = final(S);
    P_bot = POINT(x_bot);
    S = addcl(S,P_bot);
    
    %% Stiffness matrix and sollicitation vector
    A = calc_rigi(S);
    switch lower(test)
        case 'staticvertup'
            P_load = POINT(x_load(1,:));
            f = nodalload(S,P_load,'FZ',p);
        case 'staticvertdown'
            P_load = POINT(x_load(1,:));
            f = nodalload(S,P_load,'FZ',-p);
        case 'statichoriin'
            P_load = POINT(x_load(2,:));
            f = nodalload(S,P_load,'FY',-p);
        case 'statichoriout'
            P_load = POINT(x_load(2,:));
            f = nodalload(S,P_load,'FY',p);
    end
    f = f + bodyload(keepgroupelem(S,getnumgroupelemwithfield(S,'material',mat_0)),[],'FZ',-p0);
    f = f + bodyload(keepgroupelem(S,getnumgroupelemwithfield(S,'material',mat_1)),[],'FZ',-p1);
    f = f + bodyload(keepgroupelem(S,getnumgroupelemwithfield(S,'material',mat_2)),[],'FZ',-p2);
    f = f + bodyload(keepgroupelem(S,getnumgroupelemwithfield(S,'material',mat_3)),[],'FZ',-p3);
    
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
    
    S_beam = final(S_leg);
    P_beam = calc_P(S,S_beam);
    u_beam = P_beam*u; % u_beam = transfer(S,S_beam,u);
    
    S_plate = final(S_plate);
    P_plate = calc_P(S,S_plate);
    u_plate = P_plate*u; % u_plate = transfer(S,S_plate,u);
    
    e_beam = calc_epsilon(S_beam,u_beam,'smooth');
    s_beam = calc_sigma(S_beam,u_beam,'smooth');
    
    e_plate = calc_epsilon(S_plate,u_plate,'smooth');
    s_plate = calc_sigma(S_plate,u_plate,'smooth');
    
    Epsx = e_beam(1);
    Gamx = e_beam(2);
    Gamy = e_beam(3);
    Gamz = e_beam(4);
    N = s_beam(1);
    Mx = s_beam(2);
    My = s_beam(3);
    Mz = s_beam(4);
    
    Exx = e_plate(1);
    Eyy = e_plate(2);
    Exy = e_plate(3);
    Gxx = e_plate(4);
    Gyy = e_plate(5);
    Gxy = e_plate(6);
    Nxx = s_plate(1);
    Nyy = s_plate(2);
    Nxy = s_plate(3);
    Mxx = s_plate(4);
    Myy = s_plate(5);
    Mxy = s_plate(6);
    
    %% Test solution
    tol = getfemobjectoptions('tolerancepoint');
    P = POINT(x_measure);
    % numnode = find(S.node==P);
    numnode = find(distance(S.node,P)<tol);
    xP = x(numnode,:);
    
    ux = eval_sol(S,u,P,'UX');
    uy = eval_sol(S,u,P,'UY');
    uz = eval_sol(S,u,P,'UZ');
    rx = eval_sol(S,u,P,'RX');
    ry = eval_sol(S,u,P,'RY');
    rz = eval_sol(S,u,P,'RZ');
    
    % numnode_beam = find(S_beam.node==P);
    numnode_beam = find(distance(S_beam.node,P)<tol);
    [~,~,numgroupelem_beam] = findelemwithnode(S_beam,numnode_beam);
    n = 0;
    mx = 0;
    my = 0;
    mz = 0;
    epsx = 0;
    gamx = 0;
    gamy = 0;
    gamz = 0;
    for i=1:length(numgroupelem_beam)
        Ni  = reshape(N{numgroupelem_beam(i)},[getnbnode(S_beam),1]);
        Mxi = reshape(Mx{numgroupelem_beam(i)},[getnbnode(S_beam),1]);
        Myi = reshape(My{numgroupelem_beam(i)},[getnbnode(S_beam),1]);
        Mzi = reshape(Mz{numgroupelem_beam(i)},[getnbnode(S_beam),1]);
        Epsxi = reshape(Epsx{numgroupelem_beam(i)},[getnbnode(S_beam),1]);
        Gamxi = reshape(Gamx{numgroupelem_beam(i)},[getnbnode(S_beam),1]);
        Gamyi = reshape(Gamy{numgroupelem_beam(i)},[getnbnode(S_beam),1]);
        Gamzi = reshape(Gamz{numgroupelem_beam(i)},[getnbnode(S_beam),1]);
        ni = abs(double(Ni(numnode_beam)));
        mxi = abs(double(Mxi(numnode_beam)));
        myi = abs(double(Myi(numnode_beam)));
        mzi = abs(double(Mzi(numnode_beam)));
        epsxi = abs(double(Epsxi(numnode_beam)));
        gamxi = abs(double(Gamxi(numnode_beam)));
        gamyi = abs(double(Gamyi(numnode_beam)));
        gamzi = abs(double(Gamzi(numnode_beam)));
        n = max(n,ni);
        mx = max(mx,mxi);
        my = max(my,myi);
        mz = max(mz,mzi);
        epsx = max(epsx,epsxi);
        gamx = max(gamx,gamxi);
        gamy = max(gamy,gamyi);
        gamz = max(gamz,gamzi);
    end
    
    % numnode_plate = find(S_plate.node==P);
    numnode_plate = find(distance(S_plate.node,P)<tol);
    [~,~,numgroupelem_plate] = findelemwithnode(S_plate,numnode_plate);
    nxx = 0;
    nyy = 0;
    nxy = 0;
    mxx = 0;
    myy = 0;
    mxy = 0;
    exx = 0;
    eyy = 0;
    exy = 0;
    gxx = 0;
    gyy = 0;
    gxy = 0;
    for i=1:length(numgroupelem_plate)
        Nxxi = reshape(Nxx{numgroupelem_plate(i)},[getnbnode(S_plate),1]);
        Nyyi = reshape(Nyy{numgroupelem_plate(i)},[getnbnode(S_plate),1]);
        Nxyi = reshape(Nxy{numgroupelem_plate(i)},[getnbnode(S_plate),1]);
        Mxxi = reshape(Mxx{numgroupelem_plate(i)},[getnbnode(S_plate),1]);
        Myyi = reshape(Myy{numgroupelem_plate(i)},[getnbnode(S_plate),1]);
        Mxyi = reshape(Mxy{numgroupelem_plate(i)},[getnbnode(S_plate),1]);
        Exxi = reshape(Exx{numgroupelem_plate(i)},[getnbnode(S_plate),1]);
        Eyyi = reshape(Eyy{numgroupelem_plate(i)},[getnbnode(S_plate),1]);
        Exyi = reshape(Exy{numgroupelem_plate(i)},[getnbnode(S_plate),1]);
        Gxxi = reshape(Gxx{numgroupelem_plate(i)},[getnbnode(S_plate),1]);
        Gyyi = reshape(Gyy{numgroupelem_plate(i)},[getnbnode(S_plate),1]);
        Gxyi = reshape(Gxy{numgroupelem_plate(i)},[getnbnode(S_plate),1]);
        nxxi = abs(double(Nxxi(numnode_plate)));
        nyyi = abs(double(Nyyi(numnode_plate)));
        nxyi = abs(double(Nxyi(numnode_plate)));
        mxxi = abs(double(Mxxi(numnode_plate)));
        myyi = abs(double(Myyi(numnode_plate)));
        mxyi = abs(double(Mxyi(numnode_plate)));
        exxi = abs(double(Exxi(numnode_plate)));
        eyyi = abs(double(Eyyi(numnode_plate)));
        exyi = abs(double(Exyi(numnode_plate)));
        gxxi = abs(double(Gxxi(numnode_plate)));
        gyyi = abs(double(Gyyi(numnode_plate)));
        gxyi = abs(double(Gxyi(numnode_plate)));
        nxx = max(nxx,nxxi);
        nyy = max(nyy,nyyi);
        nxy = max(nxy,nxyi);
        mxx = max(mxx,nxxi);
        myy = max(myy,nyyi);
        mxy = max(mxy,nxyi);
        exx = max(exx,exxi);
        eyy = max(eyy,eyyi);
        exy = max(exy,exyi);
        gxx = max(gxx,gxxi);
        gyy = max(gyy,gyyi);
        gxy = max(gxy,gxyi);
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','S_beam','S_plate','test','elemtype',...
        'H','L','L1','l','h','H1','H2','h1','h2','b1','b2','c','d','e',...
        'f');
    save(fullfile(pathname,'solution.mat'),'u','s_beam','e_beam','s_plate','e_plate','time',...
        'Ux','Uy','Uz','Rx','Ry','Rz',...
        'N','Mx','My','Mz',...
        'Epsx','Gamx','Gamy','Gamz',...
        'Nxx','Nyy','Nxy','Mxx','Myy','Mxy',...
        'Exx','Eyy','Exy','Gxx','Gyy','Gxy');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz','rx','ry','rz',...
        'n','mx','my','mz','epsx','gamx','gamy','gamz',...
        'nxx','nyy','nxy','mxx','myy','mxy',...
        'exx','eyy','exy','gxx','gyy','gxy');
else
    load(fullfile(pathname,'problem.mat'),'S','S_beam','S_plate','test','elemtype',...
        'H','L','L1','l','h','H1','H2','h1','h2','b1','b2','c','d','e',...
        'f');
    load(fullfile(pathname,'solution.mat'),'u','s_beam','e_beam','s_plate','e_plate','time',...
        'Ux','Uy','Uz','Rx','Ry','Rz',...
        'N','Mx','My','Mz',...
        'Epsx','Gamx','Gamy','Gamz',...
        'Nxx','Nyy','Nxy','Mxx','Myy','Mxy',...
        'Exx','Eyy','Exy','Gxx','Gyy','Gxy');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz','rx','ry','rz',...
        'n','mx','my','mz','epsx','gamx','gamy','gamz',...
        'nxx','nyy','nxy','mxx','myy','mxy',...
        'exx','eyy','exy','gxx','gyy','gxy');
end

%% Outputs
fprintf('\nBed\n');
fprintf(['test : ' test '\n']);
fprintf(['mesh : ' elemtype ' elements\n']);
fprintf('nb elements       = %g\n',getnbelem(S));
fprintf('nb beam elements  = %g\n',getnbelem(S_beam));
fprintf('nb plate elements = %g\n',getnbelem(S_plate));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

fprintf('Displacement u and rotation r at point (%g,%g,%g) m\n',double(P));
fprintf('ux = %g m\n',ux);
fprintf('uy = %g m\n',uy);
fprintf('uz = %g m\n',uz);
fprintf('rx = %g rad = %g deg\n',rx,rad2deg(rx));
fprintf('ry = %g rad = %g deg\n',ry,rad2deg(ry));
fprintf('rz = %g rad = %g deg\n',rz,rad2deg(rz));
fprintf('\n');

fprintf('Maximum force N and moments Mx, My, Mz at point (%g,%g,%g) m\n',double(P));
fprintf('N  = %g N\n',n);
fprintf('Mx = %g N.m\n',mx);
fprintf('My = %g N.m\n',my);
fprintf('Mz = %g N.m\n',mz);
fprintf('\n');

fprintf('Maximum axial strain Epsx, torsion and bending strains (curvatures) Gamx, Gamy, Gamz at point (%g,%g,%g) m\n',double(P));
fprintf('Epsx = %g\n',epsx);
fprintf('Gamx = %g\n',gamx);
fprintf('Gamy = %g\n',gamy);
fprintf('Gamz = %g\n',gamz);
fprintf('\n');

fprintf('Maximum forces Nxx, Nyy, Nxy and moments Mxx, Myy, Mxy at point (%g,%g,%g) m\n',double(P));
fprintf('Nxx = %g N/m\n',nxx);
fprintf('Nyy = %g N/m\n',nyy);
fprintf('Nxy = %g N/m\n',nxy);
fprintf('Mxx = %g N\n',mxx);
fprintf('Myy = %g N\n',myy);
fprintf('Mxy = %g N\n',mxy);
fprintf('\n');

fprintf('Maximum membrane strains Exx, Eyy, Exy and bending strains (curvatures) Gxx, Gyy, Gxy at point (%g,%g,%g) m\n',double(P));
fprintf('Exx = %g\n',exx);
fprintf('Eyy = %g\n',eyy);
fprintf('Exy = %g\n',exy);
fprintf('Gxx = %g\n',gxx);
fprintf('Gyy = %g\n',gyy);
fprintf('Gxy = %g\n',gxy);
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    figure('Name','Domain')
    clf
    h1 = plot(S,'selgroup',getnumgroupelemwithfield(S,'material',mat_0),'EdgeColor','k');
    hold on
    numelem1 = getnumgroupelemwithfield(S,'material',mat_1);
    numelem2 = getnumgroupelemwithfield(S,'material',mat_2);
    numelem3 = getnumgroupelemwithfield(S,'material',mat_3);
    h2 = plot(S,'selgroup',numelem1(1:3),'EdgeColor','c','FaceColor','c','FaceAlpha',0.1);
    h3 = plot(S,'selgroup',numelem1(4:6),'EdgeColor','r','FaceColor','r','FaceAlpha',0.1);
    h4 = plot(S,'selgroup',numelem2(1:2),'EdgeColor',[1 0.5 0],'FaceColor',[1 0.5 0],'FaceAlpha',0.1);
    h5 = plot(S,'selgroup',numelem1(7:14),'EdgeColor','b','FaceColor','b','FaceAlpha',0.1);
    h6 = plot(S,'selgroup',numelem2(3:4),'EdgeColor','m','FaceColor','m','FaceAlpha',0.1);
    if slat
        h7 = plot(S,'selgroup',numelem1(15:28),'EdgeColor','g','FaceColor','g','FaceAlpha',0.1);
    end
    h8 = plot(S,'selgroup',numelem3(1),'EdgeColor',([1 0 1]+[1 0 0])/2,'FaceColor',([1 0 1]+[1 0 0])/2,'FaceAlpha',0.1);
    h9 = plot(S,'selgroup',numelem3(2:3),'EdgeColor',([1 0 1]+[0 0 1])/2,'FaceColor',([1 0 1]+[0 0 1])/2,'FaceAlpha',0.1);
    hold off
    set(gca,'FontSize',16)
    if slat
        l = legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1),h7(1),h8(1),h9(1)],'leg','bottom rail','side rail','end rail','guard rail','guard rail support','slat',...
            'junction guard rail support - side rail','junction guard rail support - guard rail','Location','NorthEastOutside');
    else
        l = legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1),h8(1),h9(1)],'leg','bottom rail','side rail','end rail','guard rail','guard rail support',...
            'junction guard rail support - side rail','junction guard rail support - guard rail','Location','NorthEastOutside');
    end
    %set(l,'Interpreter','latex')
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
%     figure('Name','Group of elements')
%     plotparamelem(S,'group')
%     mysaveas(pathname,'groupelem',formats,renderer);
    
%     figure('Name','Materials')
%     plotparamelem(S,'material')
%     mysaveas(pathname,'material',formats,renderer);
    
%     plotDomain(S,'legend',false);
%     mysaveas(pathname,'domain',formats,renderer);
%     mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'FaceColor','k','legend',false);
    ampl = 5;
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
    hP = plot(P,'g+');
    legend([hD,hN,hP],[legD,legN,'measure'],'Location','NorthEastOutside')
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
%     plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
%     mysaveas(pathname,'Rx',formats,renderer);
%     
%     plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
%     mysaveas(pathname,'Ry',formats,renderer);
%     
%     plotSolution(S,u,'rotation',3,'ampl',ampl,options{:});
%     mysaveas(pathname,'Rz',formats,renderer);
    
    % DO NOT WORK WITH BEAM ELEMENTS
    % Beams strains
    % plotSolution(S_beam,u_beam,'epsilon',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'Epsx',formats,renderer);
    %
    % plotSolution(S_beam,u_beam,'epsilon',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'Gamx',formats,renderer);
    %
    % plotSolution(S_beam,u_beam,'epsilon',3,'ampl',ampl,options{:});
    % mysaveas(pathname,'Gamy',formats,renderer);
    %
    % plotSolution(S_beam,u_beam,'epsilon',4,'ampl',ampl,options{:});
    % mysaveas(pathname,'Gamz',formats,renderer);
    %
    % Beams stresses
    % plotSolution(S_beam,u_beam,'sigma',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'N',formats,renderer);
    %
    % plotSolution(S_beam,u_beam,'sigma',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'Mx',formats,renderer);
    %
    % plotSolution(S_beam,u_beam,'sigma',3,'ampl',ampl,options{:});
    % mysaveas(pathname,'My',formats,renderer);
    %
    % plotSolution(S_beam,u_beam,'sigma',4,'ampl',ampl,options{:});
    % mysaveas(pathname,'Mz',formats,renderer);
    
    % Beams strains
%     figure('Name','Solution Epsx')
%     clf
%     plot(e_beam,S_beam+ampl*u_beam,'compo','EPSX')
%     colorbar
%     set(gca,'FontSize',fontsize)
%     mysaveas(pathname,'Epsx',formats,renderer);
%     
%     figure('Name','Solution Gamx')
%     clf
%     plot(e_beam,S_beam+ampl*u_beam,'compo','GAMX')
%     colorbar
%     set(gca,'FontSize',fontsize)
%     mysaveas(pathname,'Gamx',formats,renderer);
%     
%     figure('Name','Solution Gamy')
%     clf
%     plot(e_beam,S_beam+ampl*u_beam,'compo','GAMY')
%     colorbar
%     set(gca,'FontSize',fontsize)
%     mysaveas(pathname,'Gamy',formats,renderer);
%     
%     figure('Name','Solution Gamz')
%     clf
%     plot(e_beam,S_beam+ampl*u_beam,'compo','GAMZ')
%     colorbar
%     set(gca,'FontSize',fontsize)
%     mysaveas(pathname,'Gamz',formats,renderer);
    
    % Beams stresses
    figure('Name','Solution N')
    clf
    plot(s_beam,S_beam+ampl*u_beam,'compo','EFFX')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'N',formats,renderer);
    
    figure('Name','Solution Mx')
    clf
    plot(s_beam,S_beam+ampl*u_beam,'compo','MOMX')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Mx',formats,renderer);
    
    figure('Name','Solution My')
    clf
    plot(s_beam,S_beam+ampl*u_beam,'compo','MOMY')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'My',formats,renderer);
    
    figure('Name','Solution Mz')
    clf
    plot(s_beam,S_beam+ampl*u_beam,'compo','MOMZ')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Mz',formats,renderer);
    
    % Plates strains
%     plotSolution(S_plate,u_plate,'epsilon',1,'ampl',ampl,options{:});
%     mysaveas(pathname,'Exx',formats,renderer);
%     
%     plotSolution(S_plate,u_plate,'epsilon',2,'ampl',ampl,options{:});
%     mysaveas(pathname,'Eyy',formats,renderer);
%     
%     plotSolution(S_plate,u_plate,'epsilon',3,'ampl',ampl,options{:});
%     mysaveas(pathname,'Exy',formats,renderer);
%     
%     plotSolution(S_plate,u_plate,'epsilon',4,'ampl',ampl,options{:});
%     mysaveas(pathname,'Gxx',formats,renderer);
%     
%     plotSolution(S_plate,u_plate,'epsilon',5,'ampl',ampl,options{:});
%     mysaveas(pathname,'Gyy',formats,renderer);
%     
%     plotSolution(S_plate,u_plate,'epsilon',6,'ampl',ampl,options{:});
%     mysaveas(pathname,'Gxy',formats,renderer);
    
    % Plates stresses
    plotSolution(S_plate,u_plate,'sigma',1,'ampl',ampl,options{:});
    mysaveas(pathname,'Nxx',formats,renderer);
    
    plotSolution(S_plate,u_plate,'sigma',2,'ampl',ampl,options{:});
    mysaveas(pathname,'Nyy',formats,renderer);
    
    plotSolution(S_plate,u_plate,'sigma',3,'ampl',ampl,options{:});
    mysaveas(pathname,'Nxy',formats,renderer);
    
    plotSolution(S_plate,u_plate,'sigma',4,'ampl',ampl,options{:});
    mysaveas(pathname,'Mxx',formats,renderer);
    
    plotSolution(S_plate,u_plate,'sigma',5,'ampl',ampl,options{:});
    mysaveas(pathname,'Myy',formats,renderer);
    
    plotSolution(S_plate,u_plate,'sigma',6,'ampl',ampl,options{:});
    mysaveas(pathname,'Mxy',formats,renderer);
end

end
