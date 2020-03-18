%% FCBA desk plate deterministic linear elasticity %%
%%-------------------------------------------------%%

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

% tests = {'Stability'}; % stability test under vertical load
% tests = {'StaticHori1'}; % test under static horizontal load 1
% tests = {'StaticHori2'}; % test under static horizontal load 2
% tests = {'StaticHori3'}; % test under static horizontal load 3 (lifting)
% tests = {'StaticHori4'}; % test under static horizontal load 4 (lifting)
tests = {'StaticVert'}; % test under static vertical load
% tests = {'Fatigue1'}; % fatigue test under horizontal load 1
% tests = {'Fatigue2'}; % fatigue test under horizontal load 2
% tests = {'Fatigue3'}; % fatigue test under horizontal load 3 (lifting)
% tests = {'Fatigue4'}; % fatigue test under horizontal load 4 (lifting)
% tests = {'Impact'}; % vertical impact test
% tests = {'Drop'}; % drop test
% tests = {'Stability','StaticVert',...
%     'StaticHori1','StaticHori2',...
%     'Fatigue1','Fatigue2'};

pointwiseLoading = false; % pointwise loading

formats = {'fig','epsc'};
renderer = 'OpenGL';

for it=1:length(tests)
    test = tests{it};
if pointwiseLoading
    filename = ['FCBADeskPlateDetLinElas' test 'PointwiseLoading'];
else
    filename = ['FCBADeskPlateDetLinElas' test];
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
    x_fati = {[x3_23,y3_12+50e-3,z3],[x3_14,y3_12+50e-3,z3],...
              [x3_23-50e-3,y3_12,z3],[x3_23-50e-3,y3_34,z3]};
    x_stab = double(getcenter(L3{1}))+[0.0,50e-3,0.0];
    x_meas = [x_hori,x_vert,x_fati,x_stab];
    P_hori = cellfun(@(x) POINT(x),x_hori,'UniformOutput',false);
    P_vert = POINT(x_vert);
    P_fati = cellfun(@(x) POINT(x),x_fati,'UniformOutput',false);
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
    x_masse = double(getcoord(getcenter(C_masse)));
    if pointwiseLoading
        PbQ3 = {x_hori{4},x_fati{3},x_fati{1},x_hori{1},...
                x_fati{4},x_hori{3},x_hori{2},x_fati{2}};
        if ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4')
            S = gmshFCBAdesksimplified(Q1,Q2,Q3,Q5a,Q5b,C_masse,PbQ3,x_stab,x_masse,...
                cl_12,cl_12,cl_3,cl_5,cl_5,cl_3,cl_3,cl_3,cl_3,...
                fullfile(pathname,['gmsh_desk_' elemtype]),3);
        else
            S = gmshFCBAdesksimplified(Q1,Q2,Q3,Q5a,Q5b,C_masse,PbQ3,x_stab,x_masse,...
                cl_12,cl_12,cl_3,cl_5,cl_5,cl_3,cl_3,cl_3,cl_3,...
                fullfile(pathname,['gmsh_desk_' elemtype]),3,'recombine');
        end
    else
        L_hori{1} = LIGNE(x_hori{1}+[0,-r_load,0],x_hori{1}+[0,r_load,0]);
        L_hori{2} = LIGNE(x_hori{2}+[0,r_load,0],x_hori{2}+[0,-r_load,0]);
        L_hori{3} = LIGNE(x_hori{3}+[r_load,0,0],x_hori{3}+[-r_load,0,0]);
        L_hori{4} = LIGNE(x_hori{4}+[-r_load,0,0],x_hori{4}+[r_load,0,0]);
        L_fati{1} = LIGNE(x_fati{1}+[0,-r_load,0],x_fati{1}+[0,r_load,0]);
        L_fati{2} = LIGNE(x_fati{2}+[0,r_load,0],x_fati{2}+[0,-r_load,0]);
        L_fati{3} = LIGNE(x_fati{3}+[-r_load,0,0],x_fati{3}+[r_load,0,0]);
        L_fati{4} = LIGNE(x_fati{4}+[r_load,0,0],x_fati{4}+[-r_load,0,0]);
        LbQ3 = {L_hori{4},L_fati{3},L_fati{1},L_hori{1},...
                L_fati{4},L_hori{3},L_hori{2},L_fati{2}};
        C_vert = CIRCLE(x_vert(1),x_vert(2),x_vert(3),r_load);
        C_stab = CIRCLE(x_stab(1),x_stab(2),x_stab(3),r_load);
        if ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4')
            S = gmshFCBAdesk(Q1,Q2,Q3,Q5a,Q5b,C_masse,LbQ3,C_stab,C_vert,...
                cl_12,cl_12,cl_3,cl_5,cl_5,cl_3,cl_3,cl_3,cl_3,...
                fullfile(pathname,['gmsh_desk_' elemtype]),3);
        else
            S = gmshFCBAdesk(Q1,Q2,Q3,Q5a,Q5b,C_masse,LbQ3,C_stab,C_vert,...
                cl_12,cl_12,cl_3,cl_5,cl_5,cl_3,cl_3,cl_3,cl_3,...
                fullfile(pathname,['gmsh_desk_' elemtype]),3,'recombine');
        end
    end
    S = convertelem(S,elemtype);
    
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
            %E = 2e9; % [Pa]
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
    S = setmaterial(S,mat);
    
    %% Neumann boundary conditions
    p_plate = RHO*g*h; % surface load (body load for plates) [N/m2]
    Sec_stab_vert = pi*r_load^2;
    L_hori_fati = 2*r_load;
    switch lower(test)
        case 'stability'
            p = 400; % pointwise load [N]
            if ~pointwiseLoading
                p = p/Sec_stab_vert; % surface load (body load for plates) [N/m2]
            end
        case {'statichori1','statichori2','statichori3','statichori4'}
            masse = 50.5; % [kg]
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse; % surface load (body load for plates) [N/m2]
            p = 100; % pointwise load, F1=F2=100, 200 [N], F3=F4=100 [N]
            if ~pointwiseLoading
                p = p/L_hori_fati; % line load (surface load for plates) [N/m]
            end
            slope = 0;
        case 'staticvert'
            p = 300; % pointwise load, 300, 400, 500 [N]
            if ~pointwiseLoading
                p = p/Sec_stab_vert; % surface load (body load for plates) [N/m2]
            end
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            masse = 50.5; % [kg]
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse; % surface load (body load for plates) [N/m2]
            p = 100; % pointwise load [N]
            if ~pointwiseLoading
                p = p/L_hori_fati; % line load (surface load for plates) [N/m]
            end
        case 'impact'
            H = 180e-3; % [m]
        case 'drop'
            H = 100e-3; % [m]
    end
    
    %% Dirichlet boundary conditions
    S = final(S);
    L1 = getedge(Q1,1);
    L2 = getedge(Q2,1);
    L5b = getedge(Q5b,1);
    [~,numnode1] = intersect(S,L1);
    [~,numnode2] = intersect(S,L2);
    [~,numnode5b] = intersect(S,L5b);
    switch lower(test)
        case 'stability'
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
    
    %% Stiffness matrix and sollicitation vector
    A = calc_rigi(S);
    
    switch lower(test)
        case 'stability'
            if pointwiseLoading
                f = nodalload(S,P_stab,'FZ',-p);
                if isempty(ispointin(P_stab,POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            else
                f = bodyload(keepgroupelem(S,6),[],'FZ',-p);
            end
        case {'statichori1','statichori2','statichori3','statichori4'}
            if strcmpi(test,'statichori1')
                if pointwiseLoading
                    f = nodalload(S,P_hori{1},{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
                    if isempty(ispointin(P_hori{1},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{1},{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
                end
            elseif strcmpi(test,'statichori2')
                if pointwiseLoading
                    f = nodalload(S,P_hori{2},{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
                    if isempty(ispointin(P_hori{2},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{2},{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
                end
            elseif strcmpi(test,'statichori3')
                if pointwiseLoading
                    f = nodalload(S,P_hori{3},{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
                    if isempty(ispointin(P_hori{3},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{3},{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
                end
            elseif strcmpi(test,'statichori4')
                if pointwiseLoading
                    f = nodalload(S,P_hori{4},{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
                    if isempty(ispointin(P_hori{4},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{4},{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
                end
            end
            if pointwiseLoading
                f = f + bodyload(keepgroupelem(S,4),[],'FZ',-p_masse);
            else
                f = f + bodyload(keepgroupelem(S,[4,5]),[],'FZ',-p_masse);
            end
        case 'staticvert'
            if pointwiseLoading
                f = nodalload(S,P_vert,'FZ',-p);
                if isempty(ispointin(P_vert,POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            else
                f = bodyload(keepgroupelem(S,4),[],'FZ',-p);
            end
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            if strcmpi(test,'fatigue1')
                if pointwiseLoading
                    f = nodalload(S,P_fati{1},'FX',-p);
                    if isempty(ispointin(P_fati{1},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_fati{1},'FX',-p);
                end
            elseif strcmpi(test,'fatigue2')
                if pointwiseLoading
                    f = nodalload(S,P_fati{2},'FX',p);
                    if isempty(ispointin(P_fati{2},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_fati{2},'FX',p);
                end
            elseif strcmpi(test,'fatigue3')
                if pointwiseLoading
                    f = nodalload(S,P_fati{3},'FY',p);
                    if isempty(ispointin(P_fati{3},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_fati{3},'FY',p);
                end
            elseif strcmpi(test,'fatigue4')
                if pointwiseLoading
                    f = nodalload(S,P_fati{4},'FY',-p);
                    if isempty(ispointin(P_fati{4},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_fati{4},'FY',-p);
                end
            end
            if pointwiseLoading
                f = f + bodyload(keepgroupelem(S,4),[],'FZ',-p_masse);
            else
                f = f + bodyload(keepgroupelem(S,[4,5]),[],'FZ',-p_masse);
            end
        case {'impact','drop'}
            error('Not implemented')
    end
    f = f + bodyload(S,[],'FZ',-p_plate);
    
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
        case 'staticvert'
            P = P_vert;
        case 'stability'
            P = P_stab;
        case 'statichori1'
            P = P_hori{2};
        case 'statichori2'
            P = P_hori{1};
        case 'statichori3'
            P = P_hori{4};
        case 'statichori4'
            P = P_hori{3};
        case 'fatigue1'
            P = P_fati{2};
        case 'fatigue2'
            P = P_fati{1};
        case 'fatigue3'
            P = P_fati{4};
        case 'fatigue4'
            P = P_fati{3};
    end
    ux = eval_sol(S,u,P,'UX');
    uy = eval_sol(S,u,P,'UY');
    uz = eval_sol(S,u,P,'UZ');
    rx = eval_sol(S,u,P,'RX');
    ry = eval_sol(S,u,P,'RY');
    rz = eval_sol(S,u,P,'RZ');
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','elemtype',...
        'a12','b12','a3','b3','a5','b5','h','Sec_stab_vert','L_hori_fati',...
        'f','p','pointwiseLoading');
    save(fullfile(pathname,'solution.mat'),'u','time');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz',...
        'rx','ry','rz');
else
    load(fullfile(pathname,'problem.mat'),'S','elemtype',...
        'a12','b12','a3','b3','a5','b5','h','Sec_stab_vert','L_hori_fati',...
        'f','p','pointwiseLoading');
    load(fullfile(pathname,'solution.mat'),'u','time');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz',...
        'rx','ry','rz');
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
        if (pointwiseLoading && p==300) || (~pointwiseLoading && p==300/Sec_stab_vert)
            uz_exp_start = -0.69*1e-3;
            uz_exp_end = -[10.10 9.88 9.64 9.88 9.94 9.79 9.92 9.93 9.82 9.95]*1e-3;
        elseif (pointwiseLoading && p==400) || (~pointwiseLoading && p==400/Sec_stab_vert)
            uz_exp_start = -0.75*1e-3;
            uz_exp_end = -[13.45 13.52 13.56 13.64 13.65 13.74 13.75 13.44 13.74 13.53]*1e-3;
        elseif (pointwiseLoading && p==500) || (~pointwiseLoading && p==500/Sec_stab_vert)
            uz_exp_start = -0.78*1e-3;
            uz_exp_end = -[16.66 16.57 16.59 16.78 16.55 16.69 16.75 16.59 16.73 16.76]*1e-3;
        end
        uz_exp = mean(uz_exp_end - uz_exp_start);
        err_uz = norm(uz-uz_exp)/norm(uz_exp);
        
        disp('Displacement u at point'); disp(P);
        fprintf('ux     = %g m\n',ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uz     = %g m\n',uz);
        fprintf('uz_exp = %g m, error = %.3e\n',uz_exp,err_uz);
        fprintf('\n');
    case 'stability'
        uz_exp_start = -1.93*1e-3;
        uz_exp_end = -[18.46 18.44 18.53 18.58 18.59 18.7 18.77 18.73 18.85 18.76]*1e-3;
        uz_exp = mean(uz_exp_end - uz_exp_start);
        err_uz = norm(uz-uz_exp)/norm(uz_exp);
        
        disp('Displacement u at point'); disp(P);
        fprintf('ux     = %g m\n',ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uz     = %g m\n',uz);
        fprintf('uz_exp = %g m, error = %.3e\n',uz_exp,err_uz);
        fprintf('\n');
    case 'statichori1'
        if (pointwiseLoading && p==100) || (~pointwiseLoading && p==100/L_hori_fati)
            ux_exp_start = -6.88*1e-3;
            ux_exp_end = -[10.5 10.51 10.44 10.8 10.72 10.62 10.67 10.65 10.66 10.87 10.86]*1e-3;
        elseif (pointwiseLoading && p==200) || (~pointwiseLoading && p==200/L_hori_fati)
            ux_exp_start = -6.16*1e-3;
            ux_exp_end = -[16.78 16.74 16.72 17.13 17 16.8 16.87 16.78 17.04 16.82 16.71 17.17]*1e-3;
        end
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux-ux_exp)/norm(ux_exp);
        
        disp('Displacement u at point'); disp(P);
        fprintf('ux     = %g m\n',ux);
        fprintf('ux_exp = %g m, error = %.3e\n',ux_exp,err_ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uz     = %g m\n',uz);
        fprintf('\n');
    case 'statichori2'
        if (pointwiseLoading && p==100) || (~pointwiseLoading && p==100/L_hori_fati)
            ux_exp_start = 2.12*1e-3;
            ux_exp_end = [6.22 6.17 6.26 6.31 6.33 6.24 6.26 6.4 6.26 6.49 6.48 6.42 6.36 6.56 6.37 6.39]*1e-3;
        elseif (pointwiseLoading && p==200) || (~pointwiseLoading && p==200/L_hori_fati)
            ux_exp_start = 1.91*1e-3;
            ux_exp_end = [12.45 12.68 12.66 12.65 12.71 12.64 12.82 12.73 12.89 12.86 12.79 12.86]*1e-3;
        end
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux-ux_exp)/norm(ux_exp);
        
        disp('Displacement u at point'); disp(P);
        fprintf('ux     = %g m\n',ux);
        fprintf('ux_exp = %g m, error = %.3e\n',ux_exp,err_ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uz     = %g m\n',uz);
        fprintf('\n');
    case 'statichori3'
        uy_exp_start = -3.77*1e-3;
        uy_exp_end = -[4.71 4.73 4.69 4.56 4.47 4.73]*1e-3;   
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy-uy_exp)/norm(uy_exp);
        
        disp('Displacement u at point'); disp(P);
        fprintf('ux     = %g m\n',ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uy_exp = %g m, error = %.3e\n',uy_exp,err_uy);
        fprintf('uz     = %g m\n',uz);
        fprintf('\n'); 
    case 'statichori4'
        uy_exp_start = 9.71*1e-3;
        uy_exp_end = [12.21 12.2 12.2 12.23 12.2 12.19 12.21]*1e-3;   
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy-uy_exp)/norm(uy_exp);
        
        disp('Displacement u at point'); disp(P);
        fprintf('ux     = %g m\n',ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uy_exp = %g m, error = %.3e\n',uy_exp,err_uy);
        fprintf('uz     = %g m\n',uz);
        fprintf('\n');
    case 'fatigue1'
        ux_exp_start = -4.42*1e-3;
        ux_exp_end = -[8.4 8.3 8.37 8.41 8.54 8.39 8.56 8.48 8.46 8.49 8.49 8.43 8.55 8.52]*1e-3;   
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux-ux_exp)/norm(ux_exp);
        
        disp('Displacement u at point'); disp(P);
        fprintf('ux     = %g m\n',ux);
        fprintf('ux_exp = %g m, error = %.3e\n',ux_exp,err_ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uz     = %g m\n',uz);
        fprintf('\n');
    case 'fatigue2'
        ux_exp_start = 3.48*1e-3;
        ux_exp_end = [7.89 7.85 8.1 8.4 8.36 8.55 8.27 8.27 8.47 8.49 8.64 8.35 8.5 8.63 8.73]*1e-3;   
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux-ux_exp)/norm(ux_exp);
        
        disp('Displacement u at point'); disp(P);
        fprintf('ux     = %g m\n',ux);
        fprintf('ux_exp = %g m, error = %.3e\n',ux_exp,err_ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uz     = %g m\n',uz);
        fprintf('\n');
    case 'fatigue3'
        uy_exp_start = 3.35*1e-3;
        uy_exp_end = [6.16 5.76 5.97 5.81 5.84 5.61 5.86 5.64 5.62 5.68]*1e-3;   
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy-uy_exp)/norm(uy_exp);
        
        disp('Displacement u at point'); disp(P);
        fprintf('ux     = %g m\n',ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uy_exp = %g m, error = %.3e\n',uy_exp,err_uy);
        fprintf('uz     = %g m\n',uz);
        fprintf('\n');
    case 'fatigue4'
        uy_exp_start = -3.75*1e-3;
        uy_exp_end = -[3.89 3.88 3.89 3.88 3.89]*1e-3;
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy-uy_exp)/norm(uy_exp);
        
        disp('Displacement u at point'); disp(P);
        fprintf('ux     = %g m\n',ux);
        fprintf('uy     = %g m\n',uy);
        fprintf('uy_exp = %g m, error = %.3e\n',uy_exp,err_uy);
        fprintf('uz     = %g m\n',uz);
        fprintf('\n');
end

disp('Rotation r at point'); disp(P);
fprintf('rx     = %g rad = %g deg\n',rx,rad2deg(rx));
fprintf('ry     = %g rad = %g deg\n',ry,rad2deg(ry));
fprintf('rz     = %g rad = %g deg\n',rz,rad2deg(rz));
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain(S,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 20;
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
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
        case {'stability','staticvert','impact','drop'}
            plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'Uz',formats,renderer);
        case {'statichori1','statichori2','fatigue1','fatigue2'}
            plotSolution(S,u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'Ux',formats,renderer);
        case {'statichori3','statichori4','fatigue3','fatigue4'}
            plotSolution(S,u,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'Uy',formats,renderer);
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
