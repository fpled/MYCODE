%% FCBA desk plate junction deterministic linear elasticity %%
%%----------------------------------------------------------%%

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

% tests = {'StaticHori1'}; % strength test under static horizontal load 1
% tests = {'StaticHori2'}; % strength test under static horizontal load 2
% tests = {'StaticHori3'}; % strength test under static horizontal load 3 (lifting)
% tests = {'StaticHori4'}; % strength test under static horizontal load 4 (lifting)
tests = {'StaticVert'}; % strength test under static vertical load
% tests = {'DurabilityHori1'}; % durability test under horizontal load 1
% tests = {'DurabilityHori2'}; % durability test under horizontal load 2
% tests = {'DurabilityHori3'}; % durability test under horizontal load 3 (lifting)
% tests = {'DurabilityHori4'}; % durability test under horizontal load 4 (lifting)
% tests = {'StabilityVert'}; % stability test under vertical load
% tests = {'Impact'}; % vertical impact test
% tests = {'Drop'}; % drop test
% tests = {'StaticHori1','StaticHori2',...
%     'StaticVert',...
%     'DurabilityHori1','DurabilityHori2',...
%     'StabilityVert'};

pointLoad = false; % point load

junction = true; % junction modeling

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

for it=1:length(tests)
    test = tests{it};
if pointLoad
    filename = ['FCBADeskPlateJunctionDetLinElas' test 'PointLoad'];
else
    filename = ['FCBADeskPlateJunctionDetLinElas' test];
end
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','FCBA',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Problem
if solveProblem
    %% Domains and meshes
    % Plates dimensions
    a12 = 750e-3; % [m]
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
    x1 = (a5+h)/2;
    y1_14 = -b12+c;
    y1_23 = c;
    z1_12 = 0;
    z1_34 = a12+h/2;
    Q1 = QUADRANGLE([x1,y1_14,z1_12],[x1,y1_23,z1_12],...
                    [x1,y1_23,z1_34],[x1,y1_14,z1_34]);
    x2 = -(a5+h)/2;
    y2_14 = -b12+c;
    y2_23 = c;
    z2_12 = 0;
    z2_34 = a12+h/2;
    Q2 = QUADRANGLE([x2,y2_14,z2_12],[x2,y2_23,z2_12],...
                    [x2,y2_23,z2_34],[x2,y2_14,z2_34]);
    x3_14 = -a3/2;
    x3_23 = a3/2;
    y3_12 = c-(b3+b12)/2;
    y3_34 = c+(b3-b12)/2;
    z3 = a12+h/2;
    Q3 = QUADRANGLE([x3_14,y3_12,z3],[x3_23,y3_12,z3],...
                    [x3_23,y3_34,z3],[x3_14,y3_34,z3]);
    x5a_14 = -(a5+h)/2;
    x5a_23 = (a5+h)/2;
    y5a = 0;
    z5a_12 = a12-d-e-b;
    z5a_34 = a12-d+a;
    Q5a = QUADRANGLE([x5a_14,y5a,z5a_12],[x5a_23,y5a,z5a_12],...
                     [x5a_23,y5a,z5a_34],[x5a_14,y5a,z5a_34]);
    x5b_14 = -(a5+h)/2;
    x5b_23 = (a5+h)/2;
    y5b = 0;
    z5b_12 = f-b;
    z5b_34 = f-b+b5;
    Q5b = QUADRANGLE([x5b_14,y5b,z5b_12],[x5b_23,y5b,z5b_12],...
                     [x5b_23,y5b,z5b_34],[x5b_14,y5b,z5b_34]);
    
    % Points
    L3 = getedges(Q3);
    x_hori = {double(getcenter(L3{2})),double(getcenter(L3{4})),...
              double(getcenter(L3{3})),double(getcenter(L3{1}))};
    x_vert = double(getcenter(Q3));
    x_dura = {[x3_23,y3_12+50e-3,z3],[x3_14,y3_12+50e-3,z3],...
              [x3_23-50e-3,y3_12,z3],[x3_23-50e-3,y3_34,z3]};
    x_stab = double(getcenter(L3{1}))+[0.0,50e-3,0.0];
    x_meas = [x_hori,x_vert,x_dura,x_stab];
    P_hori = cellfun(@(x) POINT(x),x_hori,'UniformOutput',false);
    P_vert = POINT(x_vert);
    P_dura = cellfun(@(x) POINT(x),x_dura,'UniformOutput',false);
    P_stab = POINT(x_stab);
    P_meas = cellfun(@(x) POINT(x),x_meas,'UniformOutput',false);
    
    % Plates meshes
    elemtype = 'DKT';
    cl = h;
    cl_12 = cl;
    cl_3 = cl;
    cl_5 = cl;
    r_load = 40e-3;
    r_masse = 100e-3;
    C_masse = CIRCLE(0.0,y3_12+b3/2,z3,r_masse);
    x_masse = double(getcenter(C_masse));
    %
    L1_a = LINE([x5a_23,y5a,z5a_12],[x5a_23,y5a,z5a_34]);
    L1_b = LINE([x5b_23,y5b,z5b_12],[x5b_23,y5b,z5b_34]);
    if ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4')
        S1 = gmshFCBAdesk12(Q1,L1_a,L1_b,cl_12,cl_5,cl_5,...
            fullfile(pathname,['gmsh_desk_1_' elemtype]),3);
    else
        S1 = gmshFCBAdesk12(Q1,L1_a,L1_b,cl_12,cl_5,cl_5,...
            fullfile(pathname,['gmsh_desk_1_' elemtype]),3,'recombine');
    end
    S1 = convertelem(S1,elemtype);
    %
    L2_a = LINE([x5a_14,y5a,z5a_12],[x5a_14,y5a,z5a_34]);
    L2_b = LINE([x5b_14,y5b,z5b_12],[x5b_14,y5b,z5b_34]);
    if ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4')
        S2 = gmshFCBAdesk12(Q2,L2_a,L2_b,cl_12,cl_5,cl_5,...
            fullfile(pathname,['gmsh_desk_2_' elemtype]),3);
    else
        S2 = gmshFCBAdesk12(Q2,L2_a,L2_b,cl_12,cl_5,cl_5,...
            fullfile(pathname,['gmsh_desk_2_' elemtype]),3,'recombine');
    end
    S2 = convertelem(S2,elemtype);
    %
    L3_1 = LINE([x1,y1_23,z1_34],[x1,y1_14,z1_34]);
    L3_2 = LINE([x2,y2_23,z2_34],[x2,y2_14,z2_34]);
    if pointLoad
        PbQ3 = {x_hori{4},x_dura{3},x_dura{1},x_hori{1},...
                x_dura{4},x_hori{3},x_hori{2},x_dura{2}};
        if ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4')
            S3 = gmshFCBAdesk3pointload(Q3,C_masse,L3_1,L3_2,PbQ3,x_stab,x_masse,...
                cl_3,cl_3,cl_12,cl_12,cl_3,cl_3,cl_3,...
                fullfile(pathname,['gmsh_desk_3_' elemtype]),3);
        else
            S3 = gmshFCBAdesk3pointload(Q3,C_masse,L3_1,L3_2,PbQ3,x_stab,x_masse,...
                cl_3,cl_3,cl_12,cl_12,cl_3,cl_3,cl_3,...
                fullfile(pathname,['gmsh_desk_3_' elemtype]),3,'recombine');
        end
    else
        L_hori{1} = LINE(x_hori{1}+[0,-r_load,0],x_hori{1}+[0,r_load,0]);
        L_hori{2} = LINE(x_hori{2}+[0,r_load,0],x_hori{2}+[0,-r_load,0]);
        L_hori{3} = LINE(x_hori{3}+[r_load,0,0],x_hori{3}+[-r_load,0,0]);
        L_hori{4} = LINE(x_hori{4}+[-r_load,0,0],x_hori{4}+[r_load,0,0]);
        L_dura{1} = LINE(x_dura{1}+[0,-r_load,0],x_dura{1}+[0,r_load,0]);
        L_dura{2} = LINE(x_dura{2}+[0,r_load,0],x_dura{2}+[0,-r_load,0]);
        L_dura{3} = LINE(x_dura{3}+[-r_load,0,0],x_dura{3}+[r_load,0,0]);
        L_dura{4} = LINE(x_dura{4}+[r_load,0,0],x_dura{4}+[-r_load,0,0]);
        LbQ3 = {L_hori{4},L_dura{3},L_dura{1},L_hori{1},...
                L_dura{4},L_hori{3},L_hori{2},L_dura{2}};
        C_vert = CIRCLE(x_vert(1),x_vert(2),x_vert(3),r_load);
        C_stab = CIRCLE(x_stab(1),x_stab(2),x_stab(3),r_load);
        if ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4')
            S3 = gmshFCBAdesk3(Q3,C_masse,L3_1,L3_2,LbQ3,C_stab,C_vert,...
                cl_3,cl_3,cl_12,cl_12,cl_3,cl_3,cl_3,...
                fullfile(pathname,['gmsh_desk_3_' elemtype]),3);
        else
            S3 = gmshFCBAdesk3(Q3,C_masse,L3_1,L3_2,LbQ3,C_stab,C_vert,...
                cl_3,cl_3,cl_12,cl_12,cl_3,cl_3,cl_3,...
                fullfile(pathname,['gmsh_desk_3_' elemtype]),3,'recombine');
        end
    end
    S3 = convertelem(S3,elemtype);
    %
    S5a = build_model(Q5a,'cl',cl_5,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_desk_5a_' elemtype]));
    %
    S5b = build_model(Q5b,'cl',cl_5,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_desk_5b_' elemtype]));
    
    %% Materials
    % Gravitational acceleration
    g = 9.81; % [m/s2]
    
    % Density
    Mass_total = 13.9; % [kg]
    Vol_total = h*(a12*b12*2+a3*b3+a5*b5*2);
    RHO = Mass_total/(Vol_total); % [kg/m3]
    
    % Data
    filenameAna = 'data_ET_GL.mat';
    filenameNum = 'data_EL_NUL.mat';
    filenameS = 'data_KS.mat';
    filenameD = 'data_KD.mat';
    pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','identification','materialParticleBoard');
    load(fullfile(pathnameIdentification,filenameAna));
    load(fullfile(pathnameIdentification,filenameNum));
    load(fullfile(pathnameIdentification,filenameS));
    load(fullfile(pathnameIdentification,filenameD));
    
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
            NU = 0.2;
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
            %NUL = mean(mean_NUL_data);
            NUL = 0.045;
            % Transverse Poisson ratio
            NUT = 0.2;
            % Material
            mat = ELAS_SHELL_ISOT_TRANS('ET',ET,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',h,'k',5/6);
        otherwise
            error('Wrong material symmetry !')
    end
    mat = setnumber(mat,1);
    if junction
        S = union(S1,S2,S3,S5a,S5b,'duplicate');
    else
        S = union(S1,S2,S3,S5a,S5b);
    end
    S = setmaterial(S,mat);
    
    %% Neumann boundary conditions
    p_plate = RHO*g*h; % surface load (body load for plates) [N/m2]
    L_hori_dura = 2*r_load;
    Sec_vert_stab = pi*r_load^2;
    switch lower(test)
        case {'statichori1','statichori2','statichori3','statichori4'}
            masse = 50.5; % [kg]
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse; % surface load (body load for plates) [N/m2]
            p = 100; % point load, F1=F2=100, 200 [N], F3=F4=100 [N]
            if ~pointLoad
                p = p/L_hori_dura; % line load (surface load for plates) [N/m]
            end
            slope = 0;
        case 'staticvert'
            p = 300; % point load, F=300, 400, 500 [N]
            if ~pointLoad
                p = p/Sec_vert_stab; % surface load (body load for plates) [N/m2]
            end
        case {'durabilityhori1','durabilityhori2','durabilityhori3','durabilityhori4'}
            masse = 50.5; % [kg]
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse; % surface load (body load for plates) [N/m2]
            p = 100; % point load, F=100 [N]
            if ~pointLoad
                p = p/L_hori_dura; % line load (surface load for plates) [N/m]
            end
        case 'stabilityvert'
            p = 400; % point load, V=400 [N]
            if ~pointLoad
                p = p/Sec_vert_stab; % surface load (body load for plates) [N/m2]
            end
        case 'impact'
            H = 180e-3; % [m]
        case 'drop'
            H = 100e-3; % [m]
    end
    
    %% Dirichlet boundary conditions
    if junction
        S = final(S,'duplicate');
    else
        S = final(S);
    end
    
    L1 = getedge(Q1,1);
    L2 = getedge(Q2,1);
    L5b = getedge(Q5b,1);
    [~,numnodeL1] = intersect(S,L1);
    [~,numnodeL2] = intersect(S,L2);
    [~,numnodeL5b] = intersect(S,L5b);
    switch lower(test)
        case {'statichori1','statichori2'}
            S = addcl(S,numnodeL2);
            S = addcl(S,union(numnodeL1,numnodeL5b),{'UY','UZ'});
        case {'statichori3','statichori4'}
            S = addcl(S,union(numnodeL1,numnodeL2));
            S = addcl(S,numnodeL5b,'UZ');
        case 'staticvert'
            S = addcl(S,union(numnodeL1,numnodeL2));
            S = addcl(S,numnodeL5b,'UZ');
        case {'durabilityhori1','durabilityhori2','durabilityhori3','durabilityhori4'}
            S = addcl(S,union(numnodeL1,numnodeL2));
            S = addcl(S,numnodeL5b,'UZ');
        case 'stabilityvert'
            S = addcl(S,union(numnodeL1,numnodeL2));
            S = addcl(S,numnodeL5b,'UZ');
        case {'impact','drop'}
            S = addcl(S,union(numnodeL1,numnodeL2));
            S = addcl(S,numnodeL5b,'UZ');
    end
    
    if junction
        L13 = getedge(Q1,3);
        L23 = getedge(Q2,3);
        L15a = getedge(Q5a,2);
        L15b = getedge(Q5b,2);
        L25a = getedge(Q5a,4);
        L25b = getedge(Q5b,4);
        % intersection points between plates 1 and 3
        [~,numnodeL13] = intersect(S,L13);
        numnode13 = cell(2,1);
        xnode13 = cell(2,1);
        for i=1:2
            numnode13{i} = double(getnumnodeingroupelem(S,2*i-1));
            numnode13{i} = intersect(numnodeL13,numnode13{i});
            xnode13{i} = double(getcoord(getnode(S),numnode13{i}));
        end
        for i=1:size(xnode13{1},2)
            [xnode13{1},I1] = sortrows(xnode13{1},i);
            numnode13{1} = numnode13{1}(I1);
            [xnode13{2},I2] = sortrows(xnode13{2},i);
            numnode13{2} = numnode13{2}(I2);
        end
        S = addclperiodic(S,numnode13{1},numnode13{2},{'U','RX','RZ'});
        % intersection points between plates 2 and 3
        [~,numnodeL23] = intersect(S,L23);
        numnode23 = cell(2,1);
        xnode23 = cell(2,1);
        for i=1:2
            numnode23{i} = double(getnumnodeingroupelem(S,i+1));
            numnode23{i} = intersect(numnodeL23,numnode23{i});
            xnode23{i} = double(getcoord(getnode(S),numnode23{i}));
        end
        for i=1:size(xnode23{1},2)
            [xnode23{1},I1] = sortrows(xnode23{1},i);
            numnode23{1} = numnode23{1}(I1);
            [xnode23{2},I2] = sortrows(xnode23{2},i);
            numnode23{2} = numnode23{2}(I2);
        end
        S = addclperiodic(S,numnode23{1},numnode23{2},{'U','RX','RZ'});
        % intersection points between plates 1 and 5a
        [~,numnodeL15a] = intersect(S,L15a);
        numnode15a = cell(2,1);
        xnode15a = cell(2,1);
        for i=1:2
            if pointLoad
                numnode15a{i} = double(getnumnodeingroupelem(S,4*i-3));
            else
                numnode15a{i} = double(getnumnodeingroupelem(S,6*i-5));
            end
            numnode15a{i} = intersect(numnodeL15a,numnode15a{i});
            xnode15a{i} = double(getcoord(getnode(S),numnode15a{i}));
        end
        for i=1:size(xnode15a{1},2)
            [xnode15a{1},I1] = sortrows(xnode15a{1},i);
            numnode15a{1} = numnode15a{1}(I1);
            [xnode15a{2},I2] = sortrows(xnode15a{2},i);
            numnode15a{2} = numnode15a{2}(I2);
        end
        S = addclperiodic(S,numnode15a{1},numnode15a{2},{'U','RX','RY'});
        % intersection points between plates 1 and 5b
        [~,numnodeL15b] = intersect(S,L15b);
        numnode15b = cell(2,1);
        xnode15b = cell(2,1);
        for i=1:2
            if pointLoad
                numnode15b{i} = double(getnumnodeingroupelem(S,5*i-4));
            else
                numnode15b{i} = double(getnumnodeingroupelem(S,7*i-6));
            end
            numnode15b{i} = intersect(numnodeL15b,numnode15b{i});
            xnode15b{i} = double(getcoord(getnode(S),numnode15b{i}));
        end
        for i=1:size(xnode15b{1},2)
            [xnode15b{1},I1] = sortrows(xnode15b{1},i);
            numnode15b{1} = numnode15b{1}(I1);
            [xnode15b{2},I2] = sortrows(xnode15b{2},i);
            numnode15b{2} = numnode15b{2}(I2);
        end
        S = addclperiodic(S,numnode15b{1},numnode15b{2},{'U','RX','RY'});
        % intersection points between plates 2 and 5a
        [~,numnodeL25a] = intersect(S,L25a);
        numnode25a = cell(2,1);
        xnode25a = cell(2,1);
        for i=1:2
            if pointLoad
                numnode25a{i} = double(getnumnodeingroupelem(S,3*i-1));
            else
                numnode25a{i} = double(getnumnodeingroupelem(S,5*i-3));
            end
            numnode25a{i} = intersect(numnodeL25a,numnode25a{i});
            xnode25a{i} = double(getcoord(getnode(S),numnode25a{i}));
        end
        for i=1:size(xnode25a{1},2)
            [xnode25a{1},I1] = sortrows(xnode25a{1},i);
            numnode25a{1} = numnode25a{1}(I1);
            [xnode25a{2},I2] = sortrows(xnode25a{2},i);
            numnode25a{2} = numnode25a{2}(I2);
        end
        S = addclperiodic(S,numnode25a{1},numnode25a{2},{'U','RX','RY'});
        % intersection points between plates 2 and 5b
        [~,numnodeL25b] = intersect(S,L25b);
        numnode25b = cell(2,1);
        xnode25b = cell(2,1);
        for i=1:2
            if pointLoad
                numnode25b{i} = double(getnumnodeingroupelem(S,4*i-2));
            else
                numnode25b{i} = double(getnumnodeingroupelem(S,6*i-4));
            end
            numnode25b{i} = intersect(numnodeL25b,numnode25b{i});
            xnode25b{i} = double(getcoord(getnode(S),numnode25b{i}));
        end
        for i=1:size(xnode25b{1},2)
            [xnode25b{1},I1] = sortrows(xnode25b{1},i);
            numnode25b{1} = numnode25b{1}(I1);
            [xnode25b{2},I2] = sortrows(xnode25b{2},i);
            numnode25b{2} = numnode25b{2}(I2);
        end
        S = addclperiodic(S,numnode25b{1},numnode25b{2},{'U','RX','RY'});
    end
    
    %% Stiffness matrix and sollicitation vector
    A = calc_rigi(S);
    
    switch lower(test)
        case {'statichori1','statichori2','statichori3','statichori4'}
            if strcmpi(test,'statichori1')
                if pointLoad
                    f = nodalload(S,P_hori{1},{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
                    if isempty(ispointin(P_hori{1},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{1},{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
                end
            elseif strcmpi(test,'statichori2')
                if pointLoad
                    f = nodalload(S,P_hori{2},{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
                    if isempty(ispointin(P_hori{2},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{2},{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
                end
            elseif strcmpi(test,'statichori3')
                if pointLoad
                    f = nodalload(S,P_hori{3},{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
                    if isempty(ispointin(P_hori{3},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{3},{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
                end
            elseif strcmpi(test,'statichori4')
                if pointLoad
                    f = nodalload(S,P_hori{4},{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
                    if isempty(ispointin(P_hori{4},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{4},{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
                end
            end
            if pointLoad
                f = f + bodyload(keepgroupelem(S,4),[],'FZ',-p_masse);
            else
                f = f + bodyload(keepgroupelem(S,[4,5]),[],'FZ',-p_masse);
            end
        case 'staticvert'
            if pointLoad
                f = nodalload(S,P_vert,'FZ',-p);
                if isempty(ispointin(P_vert,POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            else
                f = bodyload(keepgroupelem(S,4),[],'FZ',-p);
            end
        case {'durabilityhori1','durabilityhori2','durabilityhori3','durabilityhori4'}
            if strcmpi(test,'durabilityhori1')
                if pointLoad
                    f = nodalload(S,P_dura{1},'FX',-p);
                    if isempty(ispointin(P_dura{1},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_dura{1},'FX',-p);
                end
            elseif strcmpi(test,'durabilityhori2')
                if pointLoad
                    f = nodalload(S,P_dura{2},'FX',p);
                    if isempty(ispointin(P_dura{2},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_dura{2},'FX',p);
                end
            elseif strcmpi(test,'durabilityhori3')
                if pointLoad
                    f = nodalload(S,P_dura{3},'FY',p);
                    if isempty(ispointin(P_dura{3},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_dura{3},'FY',p);
                end
            elseif strcmpi(test,'durabilityhori4')
                if pointLoad
                    f = nodalload(S,P_dura{4},'FY',-p);
                    if isempty(ispointin(P_dura{4},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_dura{4},'FY',-p);
                end
            end
            if pointLoad
                f = f + bodyload(keepgroupelem(S,4),[],'FZ',-p_masse);
            else
                f = f + bodyload(keepgroupelem(S,[4,5]),[],'FZ',-p_masse);
            end
        case 'stabilityvert'
            if pointLoad
                f = nodalload(S,P_stab,'FZ',-p);
                if isempty(ispointin(P_stab,POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            else
                f = bodyload(keepgroupelem(S,6),[],'FZ',-p);
            end
        case {'impact','drop'}
            error('Not implemented')
    end
    f = f + bodyload(S,[],'FZ',-p_plate);
    
    if junction
        % junction dofs
        numddl131 = findddl(S,'RY',numnode13{1},'free');
        numddl132 = findddl(S,'RY',numnode13{2},'free');
        numddl13 = [numddl131 numddl132];
        numddl231 = findddl(S,'RY',numnode23{1},'free');
        numddl232 = findddl(S,'RY',numnode23{2},'free');
        numddl23 = [numddl231 numddl232];
        numddl15a1 = findddl(S,'RZ',numnode15a{1},'free');
        numddl15a2 = findddl(S,'RZ',numnode15a{2},'free');
        numddl15a = [numddl15a1 numddl15a2];
        numddl15b1 = findddl(S,'RZ',numnode15b{1},'free');
        numddl15b2 = findddl(S,'RZ',numnode15b{2},'free');
        numddl15b = [numddl15b1 numddl15b2];
        numddl25a1 = findddl(S,'RZ',numnode25a{1},'free');
        numddl25a2 = findddl(S,'RZ',numnode25a{2},'free');
        numddl25a = [numddl25a1 numddl25a2];
        numddl25b1 = findddl(S,'RZ',numnode25b{1},'free');
        numddl25b2 = findddl(S,'RZ',numnode25b{2},'free');
        numddl25b = [numddl25b1 numddl25b2];
        kS = mean(mean_KS_data); % additonal bending stiffness for screw junction
        kD = mean(mean_KD_data); % additonal bending stiffness for dowel junction
        AD_add = [kD -kD;-kD kD];
        AS_add = [kS -kS;-kS kS];
        for i=1:size(numddl13,1)
            A(numddl13(i,:),numddl13(i,:)) = A(numddl13(i,:),numddl13(i,:)) + AD_add;
        end
        for i=1:size(numddl23,1)
            A(numddl23(i,:),numddl23(i,:)) = A(numddl23(i,:),numddl23(i,:)) + AD_add;
        end
        for i=1:size(numddl15a,1)
            A(numddl15a(i,:),numddl15a(i,:)) = A(numddl15a(i,:),numddl15a(i,:)) + AS_add;
        end
        for i=1:size(numddl15b,1)
            A(numddl15b(i,:),numddl15b(i,:)) = A(numddl15b(i,:),numddl15b(i,:)) + AS_add;
        end
        for i=1:size(numddl25a,1)
            A(numddl25a(i,:),numddl25a(i,:)) = A(numddl25a(i,:),numddl25a(i,:)) + AS_add;
        end
        for i=1:size(numddl25b,1)
            A(numddl25b(i,:),numddl25b(i,:)) = A(numddl25b(i,:),numddl25b(i,:)) + AS_add;
        end
    end
    
    %% Solution
    t = tic;
    u = A\f;
    time = toc(t);
    
    u = unfreevector(S,u);
    
    % U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    % Ux = u(findddl(S,'UX'),:);
    % Uy = u(findddl(S,'UY'),:);
    % Uz = u(findddl(S,'UZ'),:);
    
    % R = u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    % Rx = u(findddl(S,'RX'),:);
    % Ry = u(findddl(S,'RY'),:);
    % Rz = u(findddl(S,'RZ'),:);
    
    %% Test solution
    switch lower(test)
        case 'statichori1'
            P = P_hori{2};
        case 'statichori2'
            P = P_hori{1};
        case 'statichori3'
            P = P_hori{4};
        case 'statichori4'
            P = P_hori{3};
        case 'staticvert'
            P = P_vert;
        case 'durabilityhori1'
            P = P_dura{2};
        case 'durabilityhori2'
            P = P_dura{1};
        case 'durabilityhori3'
            P = P_dura{4};
        case 'durabilityhori4'
            P = P_dura{3};
        case 'stabilityvert'
            P = P_stab;
    end
    ux = eval_sol(S,u,P,'UX');
    uy = eval_sol(S,u,P,'UY');
    uz = eval_sol(S,u,P,'UZ');
    rx = eval_sol(S,u,P,'RX');
    ry = eval_sol(S,u,P,'RY');
    rz = eval_sol(S,u,P,'RZ');
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','elemtype',...
        'a12','b12','a3','b3','a5','b5','h','L_hori_dura','Sec_vert_stab',...
        'f','p','pointLoad');
    save(fullfile(pathname,'solution.mat'),'u','time');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz',...
        'rx','ry','rz');
else
    load(fullfile(pathname,'problem.mat'),'S','elemtype',...
        'a12','b12','a3','b3','a5','b5','h','L_hori_dura','Sec_vert_stab',...
        'f','p','pointLoad');
    load(fullfile(pathname,'solution.mat'),'u','time');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz',...
        'rx','ry','rz');
end

%% Outputs
fprintf('2D Plate Desk\n');
fprintf('\n');
fprintf(['test : ' test '\n']);
fprintf(['mesh : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('span-to-thickness ratio of plates 1 and 2 = %g\n',min(a12,b12)/h);
fprintf('span-to-thickness ratio of plate 3 = %g\n',min(a3,b3)/h);
fprintf('span-to-thickness ratio of plates 5a and 5b = %g\n',min(a5,b5)/h);
fprintf('elapsed time = %f s\n',time);

switch lower(test)
    case 'statichori1'
        if (pointLoad && p==100) || (~pointLoad && p==100/L_hori_dura)
            ux_exp_start = -6.88*1e-3;
            ux_exp_end = -[10.5 10.51 10.44 10.8 10.72 10.62 10.67 10.65 10.66 10.87 10.86]*1e-3;
        elseif (pointLoad && p==200) || (~pointLoad && p==200/L_hori_dura)
            ux_exp_start = -6.16*1e-3;
            ux_exp_end = -[16.78 16.74 16.72 17.13 17 16.8 16.87 16.78 17.04 16.82 16.71 17.17]*1e-3;
        end
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux-ux_exp)/norm(ux_exp);
        
        fprintf('\n');
        fprintf('Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf('ux     = %g m\n',ux);
        fprintf('ux_exp = %g m, error = %.3e\n',ux_exp,err_ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uz     = %g m\n',uz);
    case 'statichori2'
        if (pointLoad && p==100) || (~pointLoad && p==100/L_hori_dura)
            ux_exp_start = 2.12*1e-3;
            ux_exp_end = [6.22 6.17 6.26 6.31 6.33 6.24 6.26 6.4 6.26 6.49 6.48 6.42 6.36 6.56 6.37 6.39]*1e-3;
        elseif (pointLoad && p==200) || (~pointLoad && p==200/L_hori_dura)
            ux_exp_start = 1.91*1e-3;
            ux_exp_end = [12.45 12.68 12.66 12.65 12.71 12.64 12.82 12.73 12.89 12.86 12.79 12.86]*1e-3;
        end
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux-ux_exp)/norm(ux_exp);
        
        fprintf('\n');
        fprintf('Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf('ux     = %g m\n',ux);
        fprintf('ux_exp = %g m, error = %.3e\n',ux_exp,err_ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uz     = %g m\n',uz);
    case 'statichori3'
        uy_exp_start = -3.77*1e-3;
        uy_exp_end = -[4.71 4.73 4.69 4.56 4.47 4.73]*1e-3;
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy-uy_exp)/norm(uy_exp);
        
        fprintf('\n');
        fprintf('Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf('ux     = %g m\n',ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uy_exp = %g m, error = %.3e\n',uy_exp,err_uy);
        fprintf('uz     = %g m\n',uz);
    case 'statichori4'
        uy_exp_start = 9.71*1e-3;
        uy_exp_end = [12.21 12.2 12.2 12.23 12.2 12.19 12.21]*1e-3;
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy-uy_exp)/norm(uy_exp);
        
        fprintf('\n');
        fprintf('Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf('ux     = %g m\n',ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uy_exp = %g m, error = %.3e\n',uy_exp,err_uy);
        fprintf('uz     = %g m\n',uz);
    case 'staticvert'
        if (pointLoad && p==300) || (~pointLoad && p==300/Sec_vert_stab)
            uz_exp_start = -0.69*1e-3;
            uz_exp_end = -[10.10 9.88 9.64 9.88 9.94 9.79 9.92 9.93 9.82 9.95]*1e-3;
        elseif (pointLoad && p==400) || (~pointLoad && p==400/Sec_vert_stab)
            uz_exp_start = -0.75*1e-3;
            uz_exp_end = -[13.45 13.52 13.56 13.64 13.65 13.74 13.75 13.44 13.74 13.53]*1e-3;
        elseif (pointLoad && p==500) || (~pointLoad && p==500/Sec_vert_stab)
            uz_exp_start = -0.78*1e-3;
            uz_exp_end = -[16.66 16.57 16.59 16.78 16.55 16.69 16.75 16.59 16.73 16.76]*1e-3;
        end
        uz_exp = mean(uz_exp_end - uz_exp_start);
        err_uz = norm(uz-uz_exp)/norm(uz_exp);
        
        fprintf('\n');
        fprintf('Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf('ux     = %g m\n',ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uz     = %g m\n',uz);
        fprintf('uz_exp = %g m, error = %.3e\n',uz_exp,err_uz);
    case 'durabilityhori1'
        ux_exp_start = -4.42*1e-3;
        ux_exp_end = -[8.4 8.3 8.37 8.41 8.54 8.39 8.56 8.48 8.46 8.49 8.49 8.43 8.55 8.52]*1e-3;
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux-ux_exp)/norm(ux_exp);
        
        fprintf('\n');
        fprintf('Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf('ux     = %g m\n',ux);
        fprintf('ux_exp = %g m, error = %.3e\n',ux_exp,err_ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uz     = %g m\n',uz);
    case 'durabilityhori2'
        ux_exp_start = 3.48*1e-3;
        ux_exp_end = [7.89 7.85 8.1 8.4 8.36 8.55 8.27 8.27 8.47 8.49 8.64 8.35 8.5 8.63 8.73]*1e-3;
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux-ux_exp)/norm(ux_exp);
        
        fprintf('\n');
        fprintf('Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf('ux     = %g m\n',ux);
        fprintf('ux_exp = %g m, error = %.3e\n',ux_exp,err_ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uz     = %g m\n',uz);
    case 'durabilityhori3'
        uy_exp_start = 3.35*1e-3;
        uy_exp_end = [6.16 5.76 5.97 5.81 5.84 5.61 5.86 5.64 5.62 5.68]*1e-3;
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy-uy_exp)/norm(uy_exp);
        
        fprintf('\n');
        fprintf('Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf('ux     = %g m\n',ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uy_exp = %g m, error = %.3e\n',uy_exp,err_uy);
        fprintf('uz     = %g m\n',uz);
    case 'durabilityhori4'
        uy_exp_start = -3.75*1e-3;
        uy_exp_end = -[3.89 3.88 3.89 3.88 3.89]*1e-3;
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy-uy_exp)/norm(uy_exp);
        
        fprintf('\n');
        fprintf('Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf('ux     = %g m\n',ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uy_exp = %g m, error = %.3e\n',uy_exp,err_uy);
        fprintf('uz     = %g m\n',uz);
    case 'stabilityvert'
        uz_exp_start = -1.93*1e-3;
        uz_exp_end = -[18.46 18.44 18.53 18.58 18.59 18.7 18.77 18.73 18.85 18.76]*1e-3;
        uz_exp = mean(uz_exp_end - uz_exp_start);
        err_uz = norm(uz-uz_exp)/norm(uz_exp);
        
        fprintf('\n');
        fprintf('Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf('ux     = %g m\n',ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uz     = %g m\n',uz);
        fprintf('uz_exp = %g m, error = %.3e\n',uz_exp,err_uz);
end

fprintf('\n');
fprintf('Rotation r at point (%g,%g,%g) m\n',double(P));
fprintf('rx     = %g rad = %g deg\n',rx,rad2deg(rx));
fprintf('ry     = %g rad = %g deg\n',ry,rad2deg(ry));
fprintf('rz     = %g rad = %g deg\n',rz,rad2deg(rz));

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain(S,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 20;
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',linewidth);
    hP = plot(P,'g+');
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
    % ampl = 0;
    ampl = getsize(S)/max(abs(U))/20;
    options = {'solid',true};
    % options = {};
    
    switch lower(test)
        case {'statichori1','statichori2','durabilityhori1','durabilityhori2'}
            plotSolution(S,u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'Ux',formats,renderer);
        case {'statichori3','statichori4','durabilityhori3','durabilityhori4'}
            plotSolution(S,u,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'Uy',formats,renderer);
        case {'staticvert','stabilityvert','impact','drop'}
            plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'Uz',formats,renderer);
    end
    
    % plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'Rx',formats,renderer);
    %
    % plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'Ry',formats,renderer);
    
    % plotSolution(S,u,'epsilon','mises','ampl',ampl,options{:});
    % mysaveas(pathname,'EpsVM',formats,renderer);
    %
    % plotSolution(S,u,'sigma','mises','ampl',ampl,options{:});
    % mysaveas(pathname,'SigVM',formats,renderer);
end

end
