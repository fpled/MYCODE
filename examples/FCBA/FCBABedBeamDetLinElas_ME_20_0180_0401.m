%% FCBA bed beam deterministic linear elasticity ME 20-0180-0401 %%
%%---------------------------------------------------------------%%

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
    
filename = ['FCBABedBeamDetLinElas_ME_20_0180_0401_' test '_' materialSym '_slat_' num2str(slat)];
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
    % Beams dimensions
    L = 1990e-3; % [m]
    l = 1495e-3;
    H = 1950e-3;
    L1 = 445e-3;
    L2 = 1890e-3;
    l1 = 425e-3;
    l2 = 490e-3;
    h = 410e-3;
    H1 = 620e-3;
    H2 = 1530e-3;
    b1 = 21e-3;
    b2 = 28e-3;
    b4 = 15e-3;
    b5 = 40e-3;
    h1 = 90e-3;
    h2 = 65e-3;
    h5 = 50e-3;
    c = 50e-3;
    d = 68e-3;
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
    x_botrail = [0.0,0.0,H1+h1/2;
        L-c,0.0,H1+h1/2;
        L-c,l-c,H1+h1/2;
        0.0,l-c,H1+h1/2];
    x_toprail = [0.0,0.0,H2+h1/2;
        L-L1-c/2-h2/2,0.0,H2+h1/2;
        L-c,0.0,H2+h1/2;
        L-c,l-c,H2+h1/2;
        0.0,l-c,H2+h1/2];
    x_botguardrail = [0.0,0.0,(H+H2)/2;
        L-L1-c/2-h2/2,0.0,(H+H2)/2;
        L-c,0.0,(H+H2)/2;
        L-c,l-c,(H+H2)/2;
        0.0,l-c,(H+H2)/2];
    x_topguardrail = [0.0,0.0,H-h1/2;
        L-L1-c/2-h2/2,0.0,H-h1/2;
        L-c,0.0,H-h1/2;
        L-c,l-c,H-h1/2;
        0.0,l-c,H-h1/2];
    x_slat = zeros(14*2,3);
    for i=1:14
        x_slat(2*i-1,:) = [(L-c)/2-(14-1)*d+(i-1)*2*d,0.0,H2+h1/2];
        x_slat(2*i,:) = [(L-c)/2-(14-1)*d+(i-1)*2*d,l-c,H2+h1/2];
    end
    x_longslat = [(L-c-L2)/2,(l-c)/2-l2/2-h5,H2+h1/2;
        (L-c+L2)/2,(l-c)/2-l2/2-h5,H2+h1/2;
        (L-c+L2)/2,(l-c)/2+l2/2+h5,H2+h1/2;
        (L-c-L2)/2,(l-c)/2+l2/2+h5,H2+h1/2];
    x_slatcross = zeros(14*2,3);
    for i=1:14
        x_slatcross(2*i-1,:) = [(L-c)/2-(14-1)*d+(i-1)*2*d,(l-c)/2-l2/2-h5,H2+h1/2];
        x_slatcross(2*i,:) = [(L-c)/2-(14-1)*d+(i-1)*2*d,(l-c)/2+l2/2+h5,H2+h1/2];
    end
    x_load = [(L-c)/2,l-c,H-h1/2;
        d_load-c/2,l-c,H-h1/2];
    
    P_leg{1} = {x_bot(1,:),x_botrail(1,:),x_toprail(1,:),x_botguardrail(1,:),x_topguardrail(1,:),x_top(1,:)};
    P_leg{2} = {x_bot(2,:),x_botrail(2,:),x_toprail(3,:),x_botguardrail(3,:),x_topguardrail(3,:),x_top(2,:)};
    P_leg{3} = {x_bot(3,:),x_botrail(3,:),x_toprail(4,:),x_botguardrail(4,:),x_topguardrail(4,:),x_top(3,:)};
    P_leg{4} = {x_bot(4,:),x_botrail(4,:),x_toprail(5,:),x_botguardrail(5,:),x_topguardrail(5,:),x_top(4,:)};
    
    P_botrail{1} = {x_botrail(2,:),x_botrail(3,:)};
    P_botrail{2} = {x_botrail(3,:),x_botrail(4,:)};
    P_botrail{3} = {x_botrail(4,:),x_botrail(1,:)};
    
    P_siderail{1} = {x_toprail(1,:),x_slat(2*1-1,:),x_slat(2*2-1,:),x_slat(2*3-1,:),x_slat(2*4-1,:),x_slat(2*5-1,:),...
        x_slat(2*6-1,:),x_slat(2*7-1,:),x_slat(2*8-1,:),x_slat(2*9-1,:),x_slat(2*10-1,:),x_slat(2*11-1,:),...
        x_toprail(2,:),x_slat(2*12-1,:),x_slat(2*13-1,:),x_slat(2*14-1,:),x_toprail(3,:)};
    P_siderail{2} = {x_toprail(5,:),x_slat(2*1,:),x_slat(2*2,:),x_slat(2*3,:),x_slat(2*4,:),x_slat(2*5,:),...
        x_slat(2*6,:),x_slat(2*7,:),x_slat(2*8,:),x_slat(2*9,:),x_slat(2*10,:),x_slat(2*11,:),...
        x_slat(2*12,:),x_slat(2*13,:),x_slat(2*14,:),x_toprail(4,:)};
    
    P_endrail{1} = {x_toprail(3,:),x_toprail(4,:)};
    P_endrail{2} = {x_toprail(5,:),x_toprail(1,:)};
    
    if slat
        P_slat = cell(1,14);
        for i=1:14
            P_slat{i} = {x_slat(2*i-1,:),x_slatcross(2*i-1,:),x_slatcross(2*i,:),x_slat(2*i,:)};
        end
        P_longslat{1} = {x_longslat(1,:),x_slatcross(2*1-1,:),x_slatcross(2*2-1,:),x_slatcross(2*3-1,:),x_slatcross(2*4-1,:),x_slatcross(2*5-1,:),...
            x_slatcross(2*6-1,:),x_slatcross(2*7-1,:),x_slatcross(2*8-1,:),x_slatcross(2*9-1,:),x_slatcross(2*10-1,:),x_slatcross(2*11-1,:),...
            x_slatcross(2*12-1,:),x_slatcross(2*13-1,:),x_slatcross(2*14-1,:),x_longslat(2,:)};
        P_longslat{2} = {x_longslat(4,:),x_slatcross(2*1,:),x_slatcross(2*2,:),x_slatcross(2*3,:),x_slatcross(2*4,:),x_slatcross(2*5,:),...
            x_slatcross(2*6,:),x_slatcross(2*7,:),x_slatcross(2*8,:),x_slatcross(2*9,:),x_slatcross(2*10,:),x_slatcross(2*11,:),...
            x_slatcross(2*12,:),x_slatcross(2*13,:),x_slatcross(2*14,:),x_longslat(3,:)};
    end
    
    P_botguardrail{1} = {x_botguardrail(1,:),x_botguardrail(2,:)};
    P_botguardrail{2} = {x_botguardrail(3,:),x_botguardrail(4,:)};
    P_botguardrail{3} = {x_botguardrail(4,:),x_botguardrail(5,:)};
    P_botguardrail{4} = {x_botguardrail(5,:),x_botguardrail(1,:)};
    
    P_topguardrail{1} = {x_topguardrail(1,:),x_topguardrail(2,:)};
    P_topguardrail{2} = {x_topguardrail(3,:),x_topguardrail(4,:)};
    P_topguardrail{3} = {x_topguardrail(4,:),x_load(1,:),x_load(2,:),x_topguardrail(5,:)};
    P_topguardrail{4} = {x_topguardrail(5,:),x_topguardrail(1,:)};
    
    P_guardrailsupport = {x_toprail(2,:),x_botguardrail(2,:),x_topguardrail(2,:)};
    
    % Beams meshes
    elemtype = 'BEAM';
    cl = b1/2;
    S_leg = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_leg_' num2str(n)])),P_leg,num2cell(1:length(P_leg)),'UniformOutput',false);
    % S_leg = cellfun(@(S) concatgroupelem(S),S_leg,'UniformOutput',false);
    S_leg = union(S_leg{:});
    S_leg = convertelem(S_leg,elemtype,'param',VECTEUR([1;0;0]));
    
    S_botrail = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_botrail_' num2str(n)])),P_botrail,num2cell(1:length(P_botrail)),'UniformOutput',false);
    S_botrail = union(S_botrail{:});
    S_botrail = convertelem(S_botrail,elemtype,'param',VECTEUR([0;0;1]));
    
    S_siderail = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_siderail_' num2str(n)])),P_siderail,num2cell(1:length(P_siderail)),'UniformOutput',false);
    % S_siderail = cellfun(@(S) concatgroupelem(S),S_siderail,'UniformOutput',false);
    S_siderail = union(S_siderail{:});
    S_siderail = convertelem(S_siderail,elemtype,'param',VECTEUR([0;0;1]));
    
    S_endrail = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_endrail_' num2str(n)])),P_endrail,num2cell(1:length(P_endrail)),'UniformOutput',false);
    S_endrail = union(S_endrail{:});
    S_endrail = convertelem(S_endrail,elemtype,'param',VECTEUR([0;0;1]));
    
    S_botguardrail = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_botguardrail_' num2str(n)])),P_botguardrail,num2cell(1:length(P_botguardrail)),'UniformOutput',false);
    S_topguardrail = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_topguardrail_' num2str(n)])),P_topguardrail,num2cell(1:length(P_topguardrail)),'UniformOutput',false);
    % S_topguardrail = cellfun(@(S) concatgroupelem(S),S_topguardrail,'UniformOutput',false);
    S_guardrail = union(S_botguardrail{:},S_topguardrail{:});
    S_guardrail = convertelem(S_guardrail,elemtype,'param',VECTEUR([0;0;1]));
    
    S_guardrailsupport = gmshbeam(P_guardrailsupport,cl,fullfile(pathname,'gmsh_guardrailsupport'));
    % S_guardrailsupport = concatgroupelem(S_guardrailsupport);
    S_guardrailsupport = convertelem(S_guardrailsupport,elemtype,'param',VECTEUR([1;0;0]));
    
    if slat
        S_slat = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_slat_' num2str(n)])),P_slat,num2cell(1:length(P_slat)),'UniformOutput',false);
        S_slat = union(S_slat{:});
        S_slat = convertelem(S_slat,elemtype,'param',VECTEUR([1;0;0]));
        S_longslat = cellfun(@(P,n) gmshbeam(P,cl,fullfile(pathname,['gmsh_longslat_' num2str(n)])),P_longslat,num2cell(1:length(P_longslat)),'UniformOutput',false);
        S_longslat = union(S_longslat{:});
        S_longslat = convertelem(S_longslat,elemtype,'param',VECTEUR([0;1;0]));
    end
    
    %% Materials
    % Gravitational acceleration
    g = 9.81; % [m/s2]
    
    % Density
    RHO = 500; % [kg/m3]
    
    % Cross-section area
    Sec0 = c^2;
    Sec1 = b1*h1;
    Sec2 = b2*h1;
    Sec3 = b2*h2;
    if slat
        Sec4 = b4*d;
        Sec5 = b5*h5;
    end
    % Planar second moment of area (or Planar area moment of inertia)
    IY0 = c^4/12;
    IY1 = h1*b1^3/12;
    IY2 = h1*b2^3/12;
    IY3 = h2*b2^3/12;
    if slat
        IY4 = d*b4^3/12;
        IY5 = h5*b5^3/12;
    end
    IZ0 = IY0;
    IZ1 = b1*h1^3/12;
    IZ2 = b2*h1^3/12;
    IZ3 = b2*h2^3/12;
    if slat
        IZ4 = b4*d^3/12;
        IZ5 = b5*h5^3/12;
    end
    % Polar second moment of area (or Polar area moment of inertia)
    IX0 = IY0+IZ0;
    IX1 = IY1+IZ1;
    IX2 = IY2+IZ2;
    IX3 = IY3+IZ3;
    if slat
        IX4 = IY4+IZ4;
        IX5 = IY5+IZ5;
    end
    
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
            mat_1 = ELAS_BEAM('E',E,'NU',NU,'S',Sec1,'IZ',IZ1,'IY',IY1,'IX',IX1,'RHO',RHO);
            mat_1 = setnumber(mat_1,2);
            mat_2 = ELAS_BEAM('E',E,'NU',NU,'S',Sec2,'IZ',IZ2,'IY',IY2,'IX',IX2,'RHO',RHO);
            mat_2 = setnumber(mat_2,3);
            mat_3 = ELAS_BEAM('E',E,'NU',NU,'S',Sec3,'IZ',IZ3,'IY',IY3,'IX',IX3,'RHO',RHO);
            mat_3 = setnumber(mat_3,4);
            if slat
                mat_4 = ELAS_BEAM('E',E,'NU',NU,'S',Sec4,'IZ',IZ4,'IY',IY4,'IX',IX4,'RHO',RHO);
                mat_4 = setnumber(mat_4,5);
                mat_5 = ELAS_BEAM('E',E,'NU',NU,'S',Sec5,'IZ',IZ5,'IY',IY5,'IX',IX5,'RHO',RHO);
                mat_5 = setnumber(mat_5,6);
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
            mat_1 = ELAS_BEAM_ISOT_TRANS('EL',ET,'NUT',NUT,'GL',GL,'S',Sec1,'IZ',IZ1,'IY',IY1,'IX',IX1,'RHO',RHO);
            mat_1 = setnumber(mat_1,2);
            mat_2 = ELAS_BEAM_ISOT_TRANS('EL',ET,'NUT',NUT,'GL',GL,'S',Sec2,'IZ',IZ2,'IY',IY2,'IX',IX2,'RHO',RHO);
            mat_2 = setnumber(mat_2,3);
            mat_3 = ELAS_BEAM_ISOT_TRANS('EL',ET,'NUT',NUT,'GL',GL,'S',Sec3,'IZ',IZ3,'IY',IY3,'IX',IX3,'RHO',RHO);
            mat_3 = setnumber(mat_3,4);
            if slat
                mat_4 = ELAS_BEAM_ISOT_TRANS('EL',ET,'NUT',NUT,'GL',GL,'S',Sec4,'IZ',IZ4,'IY',IY4,'IX',IX4,'RHO',RHO);
                mat_4 = setnumber(mat_4,5);
                mat_5 = ELAS_BEAM_ISOT_TRANS('EL',ET,'NUT',NUT,'GL',GL,'S',Sec5,'IZ',IZ5,'IY',IY5,'IX',IX5,'RHO',RHO);
                mat_5 = setnumber(mat_5,6);
            end
        otherwise
            error('Wrong material symmetry !')
    end
    S_leg = setmaterial(S_leg,mat_0);
    S_botrail = setmaterial(S_botrail,mat_1);
    S_siderail = setmaterial(S_siderail,mat_1);
    S_endrail = setmaterial(S_endrail,mat_2);
    S_guardrail = setmaterial(S_guardrail,mat_1);
    S_guardrailsupport = setmaterial(S_guardrailsupport,mat_3);
    if slat
        S_slat = setmaterial(S_slat,mat_4);
        S_longslat = setmaterial(S_longslat,mat_5);
        S = union(S_leg,S_botrail,S_siderail,S_endrail,S_guardrail,S_guardrailsupport,S_slat,S_longslat);
    else
        S = union(S_leg,S_botrail,S_siderail,S_endrail,S_guardrail,S_guardrailsupport);
    end
    
    %% Neumann boundary conditions
    p0 = RHO*g*Sec0; % line load (body load for legs) [N/m]
    p1 = RHO*g*Sec1; % line load (body load for bottom rails, side rails and guard rails) [N/m]
    p2 = RHO*g*Sec2; % line load (body load for end rails) [N/m]
    p3 = RHO*g*Sec3; % line load (body load for guardrail support) [N/m]
    if slat
        p4 = RHO*g*Sec4; % line load (body load for slats) [N/m]
        p5 = RHO*g*Sec5; % line load (body load for long slats) [N/m]
    end
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
    if slat
        f = f + bodyload(keepgroupelem(S,getnumgroupelemwithfield(S,'material',mat_4)),[],'FZ',-p4);
        f = f + bodyload(keepgroupelem(S,getnumgroupelemwithfield(S,'material',mat_5)),[],'FZ',-p5);
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
    
    Epsx = e(1);
    Gamx = e(2);
    Gamy = e(3);
    Gamz = e(4);
    N = s(1);
    Mx = s(2);
    My = s(3);
    Mz = s(4);
    
    %% Test solution
    tol = getfemobjectoptions('tolerancepoint');
    x_measure = x_topguardrail(5,:);
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
    
    [~,~,numgroupelem] = findelemwithnode(S,numnode);
    n = 0;
    mx = 0;
    my = 0;
    mz = 0;
    epsx = 0;
    gamx = 0;
    gamy = 0;
    gamz = 0;
    for i=1:length(numgroupelem)
        Ni  = reshape(N{numgroupelem(i)},[getnbnode(S),1]);
        Mxi = reshape(Mx{numgroupelem(i)},[getnbnode(S),1]);
        Myi = reshape(My{numgroupelem(i)},[getnbnode(S),1]);
        Mzi = reshape(Mz{numgroupelem(i)},[getnbnode(S),1]);
        Epsxi = reshape(Epsx{numgroupelem(i)},[getnbnode(S),1]);
        Gamxi = reshape(Gamx{numgroupelem(i)},[getnbnode(S),1]);
        Gamyi = reshape(Gamy{numgroupelem(i)},[getnbnode(S),1]);
        Gamzi = reshape(Gamz{numgroupelem(i)},[getnbnode(S),1]);
        ni = abs(double(Ni(numnode)));
        mxi = abs(double(Mxi(numnode)));
        myi = abs(double(Myi(numnode)));
        mzi = abs(double(Mzi(numnode)));
        epsxi = abs(double(Epsxi(numnode)));
        gamxi = abs(double(Gamxi(numnode)));
        gamyi = abs(double(Gamyi(numnode)));
        gamzi = abs(double(Gamzi(numnode)));
        n = max(n,ni);
        mx = max(mx,mxi);
        my = max(my,myi);
        mz = max(mz,mzi);
        epsx = max(epsx,epsxi);
        gamx = max(gamx,gamxi);
        gamy = max(gamy,gamyi);
        gamz = max(gamz,gamzi);
    end
    
    [uxmax,numnodeUxmax] = max(abs(full(Ux)));
    [uymax,numnodeUymax] = max(abs(full(Uy)));
    [uzmax,numnodeUzmax] = max(abs(full(Uz)));
    [rxmax,numnodeRxmax] = max(abs(full(Rx)));
    [rymax,numnodeRymax] = max(abs(full(Ry)));
    [rzmax,numnodeRzmax] = max(abs(full(Rz)));
    
    [nmax,numgroupelemnodeNmax] = max(abs(N));
    [mxmax,numgroupelemnodeMxmax] = max(abs(Mx));
    [mymax,numgroupelemnodeMymax] = max(abs(My));
    [mzmax,numgroupelemnodeMzmax] = max(abs(Mz));
    [epsxmax,numgroupelemnodeEpsxmax] = max(abs(Epsx));
    [gamxmax,numgroupelemnodeGamxmax] = max(abs(Gamx));
    [gamymax,numgroupelemnodeGamymax] = max(abs(Gamy));
    [gamzmax,numgroupelemnodeGamzmax] = max(abs(Gamz));
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','test',...
        'L','l','H','L1','L2','l1','l2','h','H1','H2','b1','b2','b4','b5','h1','h2','h5','c','d','d_load',...
        'f');
    save(fullfile(pathname,'solution.mat'),'u','s','e','time',...
        'Ux','Uy','Uz','Rx','Ry','Rz',...
        'N','Mx','My','Mz','Epsx','Gamx','Gamy','Gamz');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz','rx','ry','rz',...
        'n','mx','my','mz','epsx','gamx','gamy','gamz',...
        'uxmax','uymax','uzmax','rxmax','rymax','rzmax',...
        'numnodeUxmax','numnodeUymax','numnodeUzmax','numnodeRxmax','numnodeRymax','numnodeRzmax',... 
        'nmax','mxmax','mymax','mzmax','epsxmax','gamxmax','gamymax','gamzmax',...
        'numgroupelemnodeNmax','numgroupelemnodeMxmax','numgroupelemnodeMymax','numgroupelemnodeMzmax','numgroupelemnodeEpsxmax','numgroupelemnodeGamxmax','numgroupelemnodeGamymax','numgroupelemnodeGamzmax');
else
    load(fullfile(pathname,'problem.mat'),'S','test',...
        'L','l','H','L1','L2','l1','l2','h','H1','H2','b1','b2','b4','b5','h1','h2','h5','c','d','d_load',...
        'f');
    load(fullfile(pathname,'solution.mat'),'u','s','e','time',...
        'Ux','Uy','Uz','Rx','Ry','Rz',...
        'N','Mx','My','Mz','Epsx','Gamx','Gamy','Gamz');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz','rx','ry','rz',...
        'n','mx','my','mz','epsx','gamx','gamy','gamz',...
        'uxmax','uymax','uzmax','rxmax','rymax','rzmax',...
        'numnodeUxmax','numnodeUymax','numnodeUzmax','numnodeRxmax','numnodeRymax','numnodeRzmax',... 
        'nmax','mxmax','mymax','mzmax','epsxmax','gamxmax','gamymax','gamzmax',...
        'numgroupelemnodeNmax','numgroupelemnodeMxmax','numgroupelemnodeMymax','numgroupelemnodeMzmax','numgroupelemnodeEpsxmax','numgroupelemnodeGamxmax','numgroupelemnodeGamymax','numgroupelemnodeGamzmax');
end

%% Outputs
fprintf('\nBed\n');
fprintf(['test : ' test '\n']);
fprintf('nb elements = %g\n',getnbelem(S));
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

fprintf('Maximum displacement u and rotation r\n');
Puxmax = POINT(getcoord(S.node(numnodeUxmax)));
Puymax = POINT(getcoord(S.node(numnodeUymax)));
Puzmax = POINT(getcoord(S.node(numnodeUzmax)));
Prxmax = POINT(getcoord(S.node(numnodeRxmax)));
Prymax = POINT(getcoord(S.node(numnodeRymax)));
Przmax = POINT(getcoord(S.node(numnodeRzmax)));
fprintf('ux = %g m at point (%g,%g,%g) m\n',uxmax,double(Puxmax));
fprintf('uy = %g m at point (%g,%g,%g) m\n',uymax,double(Puymax));
fprintf('uz = %g m at point (%g,%g,%g) m\n',uzmax,double(Puzmax));
fprintf('rx = %g rad = %g deg at point (%g,%g,%g) m\n',rxmax,rad2deg(rxmax),double(Prxmax));
fprintf('ry = %g rad = %g deg at point (%g,%g,%g) m\n',rymax,rad2deg(rymax),double(Prymax));
fprintf('rz = %g rad = %g deg at point (%g,%g,%g) m\n',rzmax,rad2deg(rzmax),double(Przmax));
fprintf('\n');

fprintf('Maximum force N and moments Mx, My, Mz\n');
PNmax = POINT(getcoord(S.node(numgroupelemnodeNmax(2))));
PMxmax = POINT(getcoord(S.node(numgroupelemnodeMxmax(2))));
PMymax = POINT(getcoord(S.node(numgroupelemnodeMymax(2))));
PMzmax = POINT(getcoord(S.node(numgroupelemnodeMzmax(2))));
fprintf('N  = %g N in groupelem #%d at point (%g,%g,%g) m\n',nmax,numgroupelemnodeNmax(1),double(PNmax));
fprintf('Mx = %g N.m in groupelem #%d at point (%g,%g,%g) m\n',mxmax,numgroupelemnodeMxmax(1),double(PMxmax));
fprintf('My = %g N.m in groupelem #%d at point (%g,%g,%g) m\n',mymax,numgroupelemnodeMymax(1),double(PMymax));
fprintf('Mz = %g N.m in groupelem #%d at point (%g,%g,%g) m\n',mzmax,numgroupelemnodeMzmax(1),double(PMzmax));
fprintf('\n');

fprintf('Maximum axial strain Epsx, torsion and bending strains (curvatures) Gamx, Gamy, Gamz\n');
PEpsxmax = POINT(getcoord(S.node(numgroupelemnodeEpsxmax(2))));
PGamxmax = POINT(getcoord(S.node(numgroupelemnodeGamxmax(2))));
PGamymax = POINT(getcoord(S.node(numgroupelemnodeGamymax(2))));
PGamzmax = POINT(getcoord(S.node(numgroupelemnodeGamzmax(2))));
fprintf('Epsx = %g in groupelem #%d at point (%g,%g,%g) m\n',epsxmax,numgroupelemnodeEpsxmax(1),double(PEpsxmax));
fprintf('Gamx = %g in groupelem #%d at point (%g,%g,%g) m\n',gamxmax,numgroupelemnodeGamxmax(1),double(PGamxmax));
fprintf('Gamy = %g in groupelem #%d at point (%g,%g,%g) m\n',gamymax,numgroupelemnodeGamymax(1),double(PGamymax));
fprintf('Gamz = %g in groupelem #%d at point (%g,%g,%g) m\n',gamzmax,numgroupelemnodeGamzmax(1),double(PGamzmax));
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    figure('Name','Domain')
    clf
    h1 = plot(S,'selgroup',getnumgroupelemwithfield(S,'material',mat_0),'EdgeColor','k');
    hold on
    numelem = getnumgroupelemwithfield(S,'material',mat_1);
    h2 = plot(S,'selgroup',numelem(1:3),'EdgeColor','c');
    h3 = plot(S,'selgroup',numelem(4:34),'EdgeColor','r');
    h4 = plot(S,'selgroup',getnumgroupelemwithfield(S,'material',mat_2),'EdgeColor',[1 0.5 0]);
    h5 = plot(S,'selgroup',numelem(end-10+1:end),'EdgeColor','b');
    h6 = plot(S,'selgroup',getnumgroupelemwithfield(S,'material',mat_3),'EdgeColor','m');
    if slat
        h7 = plot(S,'selgroup',getnumgroupelemwithfield(S,'material',mat_4),'EdgeColor','g');
        h8 = plot(S,'selgroup',getnumgroupelemwithfield(S,'material',mat_5),'EdgeColor',[0 0.5 0]);
    end
    hold off
    set(gca,'FontSize',16)
    if slat
        l = legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1),h7(1),h8(1)],'leg','bottom rail','side rail','end rail','guard rail','guard rail support','slat','long slat','Location','NorthEastOutside');
    else
        l = legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1)],'leg','bottom rail','side rail','end rail','guard rail','guard rail support','Location','NorthEastOutside');
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
    hPNmax = plot(PNmax,'r+');
    hPMxmax = plot(PMxmax,'g+');
    hPMymax = plot(PMymax,'b+');
    hPMzmax = plot(PMzmax,'m+');
    legend([hD,hN,hPNmax,hPMxmax,hPMymax,hPMzmax],[legD,legN,'N max','Mx max','My max','Mz max'],'Location','NorthEastOutside')
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
%     figure('Name','Solution Epsx')
%     clf
%     plot(e,S+ampl*u,'compo','EPSX')
%     colorbar
%     set(gca,'FontSize',fontsize)
%     mysaveas(pathname,'Epsx',formats,renderer);
%     
%     figure('Name','Solution Gamx')
%     clf
%     plot(e,S+ampl*u,'compo','GAMX')
%     colorbar
%     set(gca,'FontSize',fontsize)
%     mysaveas(pathname,'Gamx',formats,renderer);
%     
%     figure('Name','Solution Gamy')
%     clf
%     plot(e,S+ampl*u,'compo','GAMY')
%     colorbar
%     set(gca,'FontSize',fontsize)
%     mysaveas(pathname,'Gamy',formats,renderer);
%     
%     figure('Name','Solution Gamz')
%     clf
%     plot(e,S+ampl*u,'compo','GAMZ')
%     colorbar
%     set(gca,'FontSize',fontsize)
%     mysaveas(pathname,'Gamz',formats,renderer);
    
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
