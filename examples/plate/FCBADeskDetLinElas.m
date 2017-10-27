%% FCBA desk deterministic linear elasticity %%
%%-------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');

%% Input data
solveProblem = true;
displaySolution = true;

% test = 'Stability'; % stability test under vertical load
% test = 'StaticHori1'; % test under static horizontal load 1
% test = 'StaticHori2'; % test under static horizontal load 2
% test = 'StaticHori3'; % test under static horizontal load 3 (soulèvement)
% test = 'StaticHori4'; % test under static horizontal load 4 (soulèvement)
test = 'StaticVert'; % test under static vertical load
% test = 'Fatigue1'; % fatigue test under horizontal load 1
% test = 'Fatigue2'; % fatigue test under horizontal load 2
% test = 'Fatigue3'; % fatigue test under horizontal load 3 (soulèvement)
% test = 'Fatigue4'; % fatigue test under horizontal load 4 (soulèvement)
% test = 'Impact'; % vertical impact test
% test = 'Drop'; % drop test

formats = {'fig','epsc2'};
renderer = 'OpenGL';

filename = ['FCBADeskDetLinElas' test];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','plate',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Problem
if solveProblem
    %% Domains and meshes
    % Plates Dimensions
    a12 = 750e-3; % m
    b12 = 396e-3;
    a3 = 1006e-3;
    b3 = 501e-3;
    a5 = 940e-3;
    b5 = 113e-3;
    % Plates 1 and 2
    c = 74e-3;
    d = 147e-3;
    e = 63e-3;
    f = 30e-3;
    % Plate 5
    a = 20e-3;
    b = 30e-3;
    % Thickness
    % Same thickness for all the plates
    h = 15e-3;
    %
    x_S1 = (a5+h)/2;
    y14_S1 = -b12+c;
    y23_S1 = c;
    z12_S1 = 0;
    z34_S1 = a12+h/2;
    S1 = QUADRANGLE([x_S1,y14_S1,z12_S1],[x_S1,y23_S1,z12_S1],[x_S1,y23_S1,z34_S1],[x_S1,y14_S1,z34_S1]);
    x_S2 = -(a5+h)/2;
    y14_S2 = -b12+c;
    y23_S2 = c;
    z12_S2 = 0;
    z34_S2 = a12+h/2;
    S2 = QUADRANGLE([x_S2,y14_S2,z12_S2],[x_S2,y23_S2,z12_S2],[x_S2,y23_S2,z34_S2],[x_S2,y14_S2,z34_S2]);
    x14_S3 = -a3/2;
    x23_S3 = a3/2;
    y12_S3 = (-b12-b3+2*c)/2;
    y34_S3 = c+(b3-b12)/2;
    z_S3 = a12+h/2;
    S3 = QUADRANGLE([x14_S3,y12_S3,z_S3],[x23_S3,y12_S3,z_S3],[x23_S3,y34_S3,z_S3],[x14_S3,y34_S3,z_S3]);
    x14_S5a = -a5/2-h/2;
    x23_S5a = a5/2+h/2;
    y_S5a = 0;
    z12_S5a = a12-d-e-b;
    z34_S5a = a12-d+a;
    S5a = QUADRANGLE([x14_S5a,y_S5a,z12_S5a],[x23_S5a,y_S5a,z12_S5a],[x23_S5a,y_S5a,z34_S5a],[x14_S5a,y_S5a,z34_S5a]);
    x14_S5b = -a5/2-h/2;
    x23_S5b = a5/2+h/2;
    y_S5b = 0;
    z12_S5b = f-b;
    z34_S5b = f-b+b5;
    S5b = QUADRANGLE([x14_S5b,y_S5b,z12_S5b],[x23_S5b,y_S5b,z12_S5b],[x23_S5b,y_S5b,z34_S5b],[x14_S5b,y_S5b,z34_S5b]);
    
    % Points
    L_S3 = getedges(S3);
    x_load_hori = {double(getcenter(L_S3{2})),double(getcenter(L_S3{4})),...
        double(getcenter(L_S3{3})),double(getcenter(L_S3{1}))};
    x_load_vert = double(getcenter(S3));
    x_load_fati = {[x23_S3,y12_S3+50e-3,z_S3],[x14_S3,y12_S3+50e-3,z_S3],[x23_S3-50e-3,y12_S3,z_S3],[x23_S3-50e-3,y34_S3,z_S3]};
    x_load_stab = double(POINT([x_load_hori{4}(1),x_load_hori{4}(2)+50e-3,x_load_hori{4}(3)]));
    x_load = [x_load_hori,x_load_vert,x_load_fati,x_load_stab];
    P_load_hori = cellfun(@(x) POINT(x),x_load_hori,'UniformOutput',false);
    P_load_vert = POINT(x_load_vert);
    P_load_fati = cellfun(@(x) POINT(x),x_load_fati,'UniformOutput',false);
    P_load_stab = POINT(x_load_stab);
    P_load = cellfun(@(x) POINT(x),x_load,'UniformOutput',false);
    
    % Plates meshes
    elemtype = 'DST';
    cl_plate12 = b12/10;
    cl_plate3 = b3/10;
    cl_plate5 = b5/4;
    r_masse = 100e-3;
    C_masse = CIRCLE(0.0,y12_S3+1/2*b3,z_S3,r_masse);
    %
    vertices_S1 = getvertices(S1);
    PbQ_S1 = {vertices_S1{1},[x23_S5b,y_S5b,z12_S5b],vertices_S1{2},vertices_S1{3},vertices_S1{4}};
    PL_S1 = {[x23_S5a,y_S5a,z12_S5a],[x23_S5a,y_S5a,z34_S5a],[x23_S5b,y_S5b,z34_S5b]};
    S1_plate = gmshFCBAdesk12(S1,PL_S1,PbQ_S1,cl_plate12,cl_plate5,cl_plate12,fullfile(pathname,['gmsh_desk_1_' elemtype '_cl_' num2str(cl_plate12)]),3);
    S1_plate = convertelem(S1_plate,elemtype);
    %
    vertices_S2 = getvertices(S2);
    PbQ_S2 = {vertices_S2{1},[x14_S5b,y_S5b,z12_S5b],vertices_S2{2},vertices_S2{3},vertices_S2{4}};
    PL_S2 = {[x14_S5a,y_S5a,z12_S5a],[x14_S5a,y_S5a,z34_S5a],[x14_S5b,y_S5b,z34_S5b]};
    S2_plate = gmshFCBAdesk12(S2,PL_S2,PbQ_S2,cl_plate12,cl_plate5,cl_plate12,fullfile(pathname,['gmsh_desk_2_' elemtype '_cl_' num2str(cl_plate12)]),3);
    S2_plate = convertelem(S2_plate,elemtype);
    %
    vertices_S3 = getvertices(S3);
    PbQ_S3 = {vertices_S3{1},x_load_hori{4},x_load_fati{3},vertices_S3{2},...
        x_load_fati{1},x_load_hori{1},vertices_S3{3},...
        x_load_fati{4},x_load_hori{3},vertices_S3{4},...
        x_load_hori{2},x_load_fati{2}};
    PL_S3 = {[x_S1,y23_S1,z34_S1],[x_S1,y14_S1,z34_S1],[x_S2,y23_S2,z34_S2],[x_S2,y14_S2,z34_S2]};
    PiQeI = x_load_stab;
    PiI = double(getcoord(getcenter(C_masse)));
    S3_plate = gmshFCBAdesk3(S3,C_masse,PL_S3,PbQ_S3,PiQeI,PiI,cl_plate3,cl_plate3,cl_plate12,cl_plate3,cl_plate3,cl_plate3,fullfile(pathname,['gmsh_desk_3_' elemtype '_cl_' num2str(cl_plate3)]),3);
    S3_plate = convertelem(S3_plate,elemtype);
    %
    S5a_plate = build_model(S5a,'cl',cl_plate5,'elemtype',elemtype,'filename',fullfile(pathname,['gmsh_desk_5a_' elemtype '_cl_' num2str(cl_plate5)]));
    S5a_plate = convertelem(S5a_plate,elemtype);
    %
    S5b_plate = build_model(S5b,'cl',cl_plate5,'elemtype',elemtype,'filename',fullfile(pathname,['gmsh_desk_5b_' elemtype '_cl_' num2str(cl_plate5)]));
    S5b_plate = convertelem(S5b_plate,elemtype);
    
    %% Materials
    % Gravitational acceleration
    g = 9.81;
    
    % Plate
    MaterialModel = 'isot';
    % Young modulus
    ET = 1.9e9;
    % Poisson ratio
    NUT = 0.25;
    % Shear modulus
    switch MaterialModel
        case 'isot'
            GL = ET/2/(1+NUT);
        case 'isotTrans'
            GL = 0.13e9;
        otherwise
            error('Wrong material model !')
    end
    
    % Density
    Mass_plates_total = 13.9; % kg
    Volum_plates_total = h*(a12*b12*2+a3*b3+a5*b5*2);
    RHO = Mass_plates_total/(Volum_plates_total);
    % Material
    switch MaterialModel
        case 'isot'
            mat_plate = ELAS_SHELL('E',ET,'NU',NUT,'RHO',RHO,'DIM3',h,'k',5/6);
        case 'isotTrans'
            mat_plate = ELAS_SHELL_ISOT_TRANS('ET',ET,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',h,'k',5/6);
    end
    mat_plate = setnumber(mat_plate,1);
    S1_plate = setmaterial(S1_plate,mat_plate);
    S2_plate = setmaterial(S2_plate,mat_plate);
    S3_plate = setmaterial(S3_plate,mat_plate);
    S5a_plate = setmaterial(S5a_plate,mat_plate);
    S5b_plate = setmaterial(S5b_plate,mat_plate);
    
    S = union(S1_plate,S2_plate,S3_plate,S5a_plate,S5b_plate);
    
    %% Neumann boundary conditions
    p_plate = RHO*g*h;
    switch lower(test)
        case 'stability'
            p = 400;
        case {'statichori1','statichori2','statichori3','statichori4'}
            masse = 50.5;
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse;
            p = 100; % F1 F2: 100N 200N; F3 F4: 100N
            slope = 0;
        case 'staticvert'
            p = 300; % 300N, 400N, 500N
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            masse = 50.5;
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse;
            p = 100;
        case 'impact'
            H = 180e-3;
        case 'drop'
            H = 100e-3;
    end
    
    %% Dirichlet boundary conditions
    P1_S1 = POINT([x_S1,y14_S1,z12_S1]);
    P1_S2 = POINT([x_S2,y14_S2,z12_S2]);
    L1_S1 = getedge(S1,1);
    L1_S2 = getedge(S2,1);
    L1_S5b = getedge(S5b,1);
    [~,numnode1] = intersect(S,L1_S1);
    [~,numnode2] = intersect(S,L1_S2);
    [~,numnode5b] = intersect(S,L1_S5b);
    
    S = final(S);
    switch lower(test)
        case 'stability'
            % S = addcl(S,P1_S1);
            % S = addcl(S,P1_S2,'U');
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
        case {'statichori1','statichori2'}
            S = addcl(S,numnode2);
            S = addcl(S,union(numnode1,numnode5b),{'UY','UZ'});
        case {'statichori3','statichori4'}
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
        case 'staticvert'
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
        case {'impact','drop'}
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
    end
    
    %% Stiffness matrices and sollicitation vectors
    A = calc_rigi(S);
    
    switch lower(test)
        case 'stability'
            f = nodalload(S,P_load_stab,'FZ',-p);
            if isempty(ispointin(P_load_stab,POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        case {'statichori1','statichori2','statichori3','statichori4'}
            if strcmpi(test,'statichori1')
                f = nodalload(S,P_load_hori{1},{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
                if isempty(ispointin(P_load_hori{1},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            elseif strcmpi(test,'statichori2')
                f = nodalload(S,P_load_hori{2},{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
                if isempty(ispointin(P_load_hori{2},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            elseif strcmpi(test,'statichori3')
                f = nodalload(S,P_load_hori{3},{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
                if isempty(ispointin(P_load_hori{3},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            elseif strcmpi(test,'statichori4')
                f = nodalload(S,P_load_hori{4},{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
                if isempty(ispointin(P_load_hori{4},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            end
            f = f + bodyload(keepgroupelem(S,4),[],'FZ',-p_masse);
        case 'staticvert'
            f = nodalload(S,P_load_vert,'FZ',-p);
            if isempty(ispointin(P_load_vert,POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            if strcmpi(test,'fatigue1')
                f = nodalload(S,P_load_fati{1},'FX',-p);
                if isempty(ispointin(P_load_fati{1},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            elseif strcmpi(test,'fatigue2')
                f = nodalload(S,P_load_fati{2},'FX',p);
                if isempty(ispointin(P_load_fati{2},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            elseif strcmpi(test,'fatigue3')
                f = nodalload(S,P_load_fati{3},'FY',p);
                if isempty(ispointin(P_load_fati{3},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            elseif strcmpi(test,'fatigue4')
                f = nodalload(S,P_load_fati{4},'FY',-p);
                if isempty(ispointin(P_load_fati{4},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            end
            f = f + bodyload(keepgroupelem(S,4),[],'FZ',-p_masse);
        case {'impact','drop'}
            error('Not implemented')
    end
    f = f + bodyload(keepgroupelem(S,[1,2,3,4,5,6]),[],'FZ',-p_plate);
    
    %% Solution
    t = tic;
    u = A\f;
    time = toc(t);
    
    u = unfreevector(S,u);
    
    U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    Ux = u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,u,S.node,'UX')));
    Uy = u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,u,S.node,'UY')));
    Uz = u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,u,S.node,'UZ')));
    
    R = u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    Rx = u(findddl(S,'RX'),:); % Rx = double(squeeze(eval_sol(S,u,S.node,'RX')));
    Ry = u(findddl(S,'RY'),:); % Ry = double(squeeze(eval_sol(S,u,S.node,'RY')));
    Rz = u(findddl(S,'RZ'),:); % Rz = double(squeeze(eval_sol(S,u,S.node,'RZ')));
    
    %% Test solution
    P_ver = getcenter(S3);
    P_stability = POINT([x_load_hori{4}(1),x_load_hori{4}(2)+50e-3,x_load_hori{4}(3)]);
    P_hori1 = getcenter(L_S3{2});
    P_hori2 = getcenter(L_S3{4});
    P_hori3 = getcenter(L_S3{3});
    P_hori4 = getcenter(L_S3{1});
    P_fati1 = POINT([x_load_fati{1}]);
    P_fati2 = POINT([x_load_fati{2}]);
    P_fati3 = POINT([x_load_fati{3}]);
    P_fati4 = POINT([x_load_fati{4}]);
    
    ux_P_ver = eval_sol(S,u,P_ver,'UX');
    uy_P_ver = eval_sol(S,u,P_ver,'UY');
    uz_P_ver = eval_sol(S,u,P_ver,'UZ');
    ux_P_stability = eval_sol(S,u,P_stability,'UX');
    uy_P_stability = eval_sol(S,u,P_stability,'UY');
    uz_P_stability = eval_sol(S,u,P_stability,'UZ');
    ux_P_hori1 = eval_sol(S,u,P_hori1,'UX');
    uy_P_hori1 = eval_sol(S,u,P_hori1,'UY');
    uz_P_hori1 = eval_sol(S,u,P_hori1,'UZ');
    ux_P_hori2 = eval_sol(S,u,P_hori2,'UX');
    uy_P_hori2 = eval_sol(S,u,P_hori2,'UY');
    uz_P_hori2 = eval_sol(S,u,P_hori2,'UZ');
    ux_P_hori3 = eval_sol(S,u,P_hori3,'UX');
    uy_P_hori3 = eval_sol(S,u,P_hori3,'UY');
    uz_P_hori3 = eval_sol(S,u,P_hori3,'UZ');
    ux_P_hori4 = eval_sol(S,u,P_hori4,'UX');
    uy_P_hori4 = eval_sol(S,u,P_hori4,'UY');
    uz_P_hori4 = eval_sol(S,u,P_hori4,'UZ');
    ux_P_fati1 = eval_sol(S,u,P_fati1,'UX');
    uy_P_fati1 = eval_sol(S,u,P_fati1,'UY');
    uz_P_fati1 = eval_sol(S,u,P_fati1,'UZ');
    ux_P_fati2 = eval_sol(S,u,P_fati2,'UX');
    uy_P_fati2 = eval_sol(S,u,P_fati2,'UY');
    uz_P_fati2 = eval_sol(S,u,P_fati2,'UZ');
    ux_P_fati3 = eval_sol(S,u,P_fati3,'UX');
    uy_P_fati3 = eval_sol(S,u,P_fati3,'UY');
    uz_P_fati3 = eval_sol(S,u,P_fati3,'UZ');
    ux_P_fati4 = eval_sol(S,u,P_fati4,'UX');
    uy_P_fati4 = eval_sol(S,u,P_fati4,'UY');
    uz_P_fati4 = eval_sol(S,u,P_fati4,'UZ');
    
    rx_P_ver = eval_sol(S,u,P_ver,'RX');
    ry_P_ver = eval_sol(S,u,P_ver,'RY');
    rz_P_ver = eval_sol(S,u,P_ver,'RZ');
    rx_P_stability = eval_sol(S,u,P_stability,'RX');
    ry_P_stability = eval_sol(S,u,P_stability,'RY');
    rz_P_stability = eval_sol(S,u,P_stability,'RZ');
    rx_P_hori1 = eval_sol(S,u,P_hori1,'RX');
    ry_P_hori1 = eval_sol(S,u,P_hori1,'RY');
    rz_P_hori1 = eval_sol(S,u,P_hori1,'RZ');
    rx_P_hori2 = eval_sol(S,u,P_hori2,'RX');
    ry_P_hori2 = eval_sol(S,u,P_hori2,'RY');
    rz_P_hori2 = eval_sol(S,u,P_hori2,'RZ');
    rx_P_hori3 = eval_sol(S,u,P_hori3,'RX');
    ry_P_hori3 = eval_sol(S,u,P_hori3,'RY');
    rz_P_hori3 = eval_sol(S,u,P_hori3,'RZ');
    rx_P_hori4 = eval_sol(S,u,P_hori4,'RX');
    ry_P_hori4 = eval_sol(S,u,P_hori4,'RY');
    rz_P_hori4 = eval_sol(S,u,P_hori4,'RZ');
    rx_P_fati1 = eval_sol(S,u,P_fati1,'RX');
    ry_P_fati1 = eval_sol(S,u,P_fati1,'RY');
    rz_P_fati1 = eval_sol(S,u,P_fati1,'RZ');
    rx_P_fati2 = eval_sol(S,u,P_fati2,'RX');
    ry_P_fati2 = eval_sol(S,u,P_fati2,'RY');
    rz_P_fati2 = eval_sol(S,u,P_fati2,'RZ');
    rx_P_fati3 = eval_sol(S,u,P_fati3,'RX');
    ry_P_fati3 = eval_sol(S,u,P_fati3,'RY');
    rz_P_fati3 = eval_sol(S,u,P_fati3,'RZ');
    rx_P_fati4 = eval_sol(S,u,P_fati4,'RX');
    ry_P_fati4 = eval_sol(S,u,P_fati4,'RY');
    rz_P_fati4 = eval_sol(S,u,P_fati4,'RZ');
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','S1_plate','S2_plate',...
        'S3_plate','S5a_plate','S5b_plate','elemtype','a12','b12',...
        'a3','b3','a5','b5','h','f','p');
    save(fullfile(pathname,'solution.mat'),'u','time',...
        'U','Ux','Uy','Uz',...
        'R','Rx','Ry','Rz');
    save(fullfile(pathname,'test_solution.mat'),'P_ver','P_stability',...
        'P_hori1','P_hori2','P_hori3','P_hori4',...
        'P_fati1','P_fati2','P_fati3','P_fati4',...
        'ux_P_ver','uy_P_ver','uz_P_ver',...
        'ux_P_stability','uy_P_stability','uz_P_stability',...
        'ux_P_hori1','uy_P_hori1','uz_P_hori1',...
        'ux_P_hori2','uy_P_hori2','uz_P_hori2',...
        'ux_P_hori3','uy_P_hori3','uz_P_hori3',...
        'ux_P_hori4','uy_P_hori4','uz_P_hori4',...        
        'ux_P_fati1','uy_P_fati1','uz_P_fati1',...
        'ux_P_fati2','uy_P_fati2','uz_P_fati2',...
        'ux_P_fati3','uy_P_fati3','uz_P_fati3',...
        'ux_P_fati4','uy_P_fati4','uz_P_fati4',...       
        'rx_P_ver','ry_P_ver','rz_P_ver',...
        'rx_P_stability','ry_P_stability','rz_P_stability',...
        'rx_P_hori1','ry_P_hori1','rz_P_hori1',...
        'rx_P_hori2','ry_P_hori2','rz_P_hori2',...
        'rx_P_hori3','ry_P_hori3','rz_P_hori3',...
        'rx_P_hori4','ry_P_hori4','rz_P_hori4',...
        'rx_P_fati1','ry_P_fati1','rz_P_fati1',...
        'rx_P_fati2','ry_P_fati2','rz_P_fati2',...
        'rx_P_fati3','ry_P_fati3','rz_P_fati3',...
        'rx_P_fati4','ry_P_fati4','rz_P_fati4');
else
    load(fullfile(pathname,'problem.mat'),'S','S1_plate','S2_plate',...
        'S3_plate','S5a_plate','S5b_plate','elemtype','a12','b12',...
        'a3','b3','a5','b5','h','f','p');
    load(fullfile(pathname,'solution.mat'),'u','time',...
        'U','Ux','Uy','Uz',...
        'R','Rx','Ry','Rz');
    load(fullfile(pathname,'test_solution.mat'),'P_ver','P_stability',...
        'P_hori1','P_hori2','P_hori3','P_hori4',...
        'P_fati1','P_fati2','P_fati3','P_fati4',...
        'ux_P_ver','uy_P_ver','uz_P_ver',...
        'ux_P_stability','uy_P_stability','uz_P_stability',...
        'ux_P_hori1','uy_P_hori1','uz_P_hori1',...
        'ux_P_hori2','uy_P_hori2','uz_P_hori2',...
        'ux_P_hori3','uy_P_hori3','uz_P_hori3',...
        'ux_P_hori4','uy_P_hori4','uz_P_hori4',...        
        'ux_P_fati1','uy_P_fati1','uz_P_fati1',...
        'ux_P_fati2','uy_P_fati2','uz_P_fati2',...
        'ux_P_fati3','uy_P_fati3','uz_P_fati3',...
        'ux_P_fati4','uy_P_fati4','uz_P_fati4',...       
        'rx_P_ver','ry_P_ver','rz_P_ver',...
        'rx_P_stability','ry_P_stability','rz_P_stability',...
        'rx_P_hori1','ry_P_hori1','rz_P_hori1',...
        'rx_P_hori2','ry_P_hori2','rz_P_hori2',...
        'rx_P_hori3','ry_P_hori3','rz_P_hori3',...
        'rx_P_hori4','ry_P_hori4','rz_P_hori4',...
        'rx_P_fati1','ry_P_fati1','rz_P_fati1',...
        'rx_P_fati2','ry_P_fati2','rz_P_fati2',...
        'rx_P_fati3','ry_P_fati3','rz_P_fati3',...
        'rx_P_fati4','ry_P_fati4','rz_P_fati4');
end

%% Outputs
fprintf('\nDesk\n');
fprintf(['test : ' test '\n']);
fprintf(['mesh : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('span-to-thickness ratio of plates 1 and 2 = %g\n',min(a12,b12)/h);
fprintf('span-to-thickness ratio of plate 3 = %g\n',min(a3,b3)/h);
fprintf('span-to-thickness ratio of plates 5a and 5b = %g\n',min(a5,b5)/h);
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

switch lower(test)
    case 'staticvert'
        disp('Displacement u at point'); disp(P_ver);
        fprintf('ux = %g\n',ux_P_ver);
        fprintf('uy = %g\n',uy_P_ver);
        fprintf('uz = %g\n',uz_P_ver);
        if p == 300
            uz_exp_start = -0.69*1e-3;
            uz_exp_end = -[10.10 9.88 9.64 9.88 9.94 9.79 9.92 9.93 9.82 9.95]*1e-3;
        elseif p == 400
            uz_exp_start = -0.75*1e-3;
            uz_exp_end = -[13.45 13.52 13.56 13.64 13.65 13.74 13.75 13.44 13.74 13.53]*1e-3;
        elseif p == 500
            uz_exp_start = -0.78*1e-3;
            uz_exp_end = -[16.66 16.57 16.59 16.78 16.55 16.69 16.75 16.59 16.73 16.76]*1e-3;
        end
        uz_exp = mean(uz_exp_end - uz_exp_start);
        err_uz = norm(uz_P_ver-uz_exp)/norm(uz_exp);
        fprintf('uz_exp   = %g, error    = %.3e\n',uz_exp,err_uz);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_ver);
        fprintf('rx = %g\n',rx_P_ver);
        fprintf('ry = %g\n',ry_P_ver);
        fprintf('rz = %g\n',rz_P_ver);
        fprintf('\n');
    case 'stability'
        disp('Displacement u at point'); disp(P_stability);
        fprintf('ux = %g\n',ux_P_stability);
        fprintf('uy = %g\n',uy_P_stability);
        fprintf('uz = %g\n',uz_P_stability);
        uz_exp_start = -1.93*1e-3;
        uz_exp_end = -[18.46 18.44 18.53 18.58 18.59 18.7 18.77 18.73 18.85 18.76]*1e-3;
        uz_exp = mean(uz_exp_end - uz_exp_start);
        err_uz = norm(uz_P_stability-uz_exp)/norm(uz_exp);
        fprintf('uz_exp   = %g, error    = %.3e\n',uz_exp,err_uz);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_stability);
        fprintf('rx = %g\n',rx_P_stability);
        fprintf('ry = %g\n',ry_P_stability);
        fprintf('rz = %g\n',rz_P_stability);
        fprintf('\n');
    case 'statichori1'
        disp('Displacement u at point'); disp(P_hori2);
        fprintf('ux = %g\n',ux_P_hori2);
        fprintf('uy = %g\n',uy_P_hori2);
        fprintf('uz = %g\n',uz_P_hori2);
        if p==100
        ux_exp_start = -6.88*1e-3;
        ux_exp_end = -[10.5 10.51 10.44 10.8 10.72 10.62 10.67 10.65 10.66 10.87 10.86]*1e-3;
        elseif p==200
        ux_exp_start = -6.16*1e-3;
        ux_exp_end = -[16.78 16.74 16.72 17.13 17 16.8 16.87 16.78 17.04 16.82 16.71 17.17]*1e-3;
        end
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux_P_hori2-ux_exp)/norm(ux_exp);
        fprintf('ux_exp   = %g, error    = %.3e\n',ux_exp,err_ux);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_hori2);
        fprintf('rx = %g\n',rx_P_hori2);
        fprintf('ry = %g\n',ry_P_hori2);
        fprintf('rz = %g\n',rz_P_hori2);
        fprintf('\n');
    case 'statichori2'
        disp('Displacement u at point'); disp(P_hori1);
        fprintf('ux = %g\n',ux_P_hori1);
        fprintf('uy = %g\n',uy_P_hori1);
        fprintf('uz = %g\n',uz_P_hori1);
        if p==100
        ux_exp_start = 2.12*1e-3;
        ux_exp_end = [6.22 6.17 6.26 6.31 6.33 6.24 6.26 6.4 6.26 6.49 6.48 6.42 6.36 6.56 6.37 6.39]*1e-3;
        elseif p==200
        ux_exp_start = 1.91*1e-3;
        ux_exp_end = [12.45 12.68 12.66 12.65 12.71 12.64 12.82 12.73 12.89 12.86 12.79 12.86]*1e-3;    
        end    
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux_P_hori1-ux_exp)/norm(ux_exp);
        fprintf('ux_exp   = %g, error    = %.3e\n',ux_exp,err_ux);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_hori1);
        fprintf('rx = %g\n',rx_P_hori1);
        fprintf('ry = %g\n',ry_P_hori1);
        fprintf('rz = %g\n',rz_P_hori1);
        fprintf('\n'); 
    case 'statichori3'
        disp('Displacement u at point'); disp(P_hori4);
        fprintf('ux = %g\n',ux_P_hori4);
        fprintf('uy = %g\n',uy_P_hori4);
        fprintf('uz = %g\n',uz_P_hori4);
        uy_exp_start = -3.77*1e-3;
        uy_exp_end = -[4.71 4.73 4.69 4.56 4.47 4.73]*1e-3;   
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy_P_hori4-uy_exp)/norm(uy_exp);
        fprintf('uy_exp   = %g, error    = %.3e\n',uy_exp,err_uy);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_hori4);
        fprintf('rx = %g\n',rx_P_hori4);
        fprintf('ry = %g\n',ry_P_hori4);
        fprintf('rz = %g\n',rz_P_hori4);
        fprintf('\n'); 
    case 'statichori4'
        disp('Displacement u at point'); disp(P_hori3);
        fprintf('ux = %g\n',ux_P_hori3);
        fprintf('uy = %g\n',uy_P_hori3);
        fprintf('uz = %g\n',uz_P_hori3);
        uy_exp_start = 9.71*1e-3;
        uy_exp_end = [12.21 12.2 12.2 12.23 12.2 12.19 12.21]*1e-3;   
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy_P_hori3-uy_exp)/norm(uy_exp);
        fprintf('uy_exp   = %g, error    = %.3e\n',uy_exp,err_uy);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_hori3);
        fprintf('rx = %g\n',rx_P_hori3);
        fprintf('ry = %g\n',ry_P_hori3);
        fprintf('rz = %g\n',rz_P_hori3);
        fprintf('\n');
    case 'fatigue1'
        disp('Displacement u at point'); disp(P_fati2);
        fprintf('ux = %g\n',ux_P_fati2);
        fprintf('uy = %g\n',uy_P_fati2);
        fprintf('uz = %g\n',uz_P_fati2);
        ux_exp_start = -4.42*1e-3;
        ux_exp_end = -[8.4 8.3 8.37 8.41 8.54 8.39 8.56 8.48 8.46 8.49 8.49 8.43 8.55 8.52]*1e-3;   
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux_P_fati2-ux_exp)/norm(ux_exp);
        fprintf('ux_exp   = %g, error    = %.3e\n',ux_exp,err_ux);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_fati2);
        fprintf('rx = %g\n',rx_P_fati2);
        fprintf('ry = %g\n',ry_P_fati2);
        fprintf('rz = %g\n',rz_P_fati2);
        fprintf('\n');
    case 'fatigue2'
        disp('Displacement u at point'); disp(P_fati1);
        fprintf('ux = %g\n',ux_P_fati1);
        fprintf('uy = %g\n',uy_P_fati1);
        fprintf('uz = %g\n',uz_P_fati1);
        ux_exp_start = 3.48*1e-3;
        ux_exp_end = [7.89 7.85 8.1 8.4 8.36 8.55 8.27 8.27 8.47 8.49 8.64 8.35 8.5 8.63 8.73]*1e-3;   
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux_P_fati1-ux_exp)/norm(ux_exp);
        fprintf('ux_exp   = %g, error    = %.3e\n',ux_exp,err_ux);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_fati1);
        fprintf('rx = %g\n',rx_P_fati1);
        fprintf('ry = %g\n',ry_P_fati1);
        fprintf('rz = %g\n',rz_P_fati1);
        fprintf('\n'); 
    case 'fatigue3'
        disp('Displacement u at point'); disp(P_fati4);
        fprintf('ux = %g\n',ux_P_fati4);
        fprintf('uy = %g\n',uy_P_fati4);
        fprintf('uz = %g\n',uz_P_fati4);
        uy_exp_start = 3.35*1e-3;
        uy_exp_end = [6.16 5.76 5.97 5.81 5.84 5.61 5.86 5.64 5.62 5.68]*1e-3;   
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy_P_fati1-uy_exp)/norm(uy_exp);
        fprintf('uy_exp   = %g, error    = %.3e\n',uy_exp,err_uy);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_fati4);
        fprintf('rx = %g\n',rx_P_fati4);
        fprintf('ry = %g\n',ry_P_fati4);
        fprintf('rz = %g\n',rz_P_fati4);
        fprintf('\n');
     case 'fatigue4'
        disp('Displacement u at point'); disp(P_fati3);
        fprintf('ux = %g\n',ux_P_fati3);
        fprintf('uy = %g\n',uy_P_fati3);
        fprintf('uz = %g\n',uz_P_fati3);
        uy_exp_start = -3.75*1e-3;
        uy_exp_end = -[3.89 3.88 3.89 3.88 3.89]*1e-3;   
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy_P_fati1-uy_exp)/norm(uy_exp);
        fprintf('uy_exp   = %g, error    = %.3e\n',uy_exp,err_uy);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_fati3);
        fprintf('rx = %g\n',rx_P_fati3);
        fprintf('ry = %g\n',ry_P_fati3);
        fprintf('rz = %g\n',rz_P_fati3);
        fprintf('\n');
    
end

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain(S,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 8;
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
    % legend([hD,hN],[legD,legN])
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'node',true,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    ampl = getsize(S)/max(abs(u))/10;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'node',true,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'node',true);
    plot(S+ampl*u,'Color','b','FaceColor','b','FaceAlpha',0.1,'node',true);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    %% Display solution
    ampl = 0;
%     ampl = getsize(S)/max(abs(u))/10;
    options = {'solid',true};
    % options = {};
    
    switch lower(test)
        case 'stability'
            plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'Uz',formats,renderer);
        case {'statichori1','statichori2'}
            plotSolution(S,u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'Ux',formats,renderer);
        case {'statichori3','statichori4'}
            plotSolution(S,u,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'Uy',formats,renderer);
        case 'staticvert'
            plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'Uz',formats,renderer);
        case {'fatigue1','fatigue2'}
            plotSolution(S,u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'Ux',formats,renderer);
        case {'fatigue3','fatigue4'}
            plotSolution(S,u,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'Uy',formats,renderer);
        case 'impact'
            plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'Uz',formats,renderer);
        case 'drop'
            plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'Uz',formats,renderer);
    end
    
    % plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'Rx',formats,renderer);
    %
    % plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'Ry',formats,renderer);
end
