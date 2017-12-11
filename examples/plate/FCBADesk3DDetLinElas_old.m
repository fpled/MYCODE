%% FCBA desk 3D deterministic linear elasticity %%
%%----------------------------------------------%%

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

% test = 'Stability'; % stability test under vertical load
% test = 'StaticHori1'; % test under static horizontal load 1
% test = 'StaticHori2'; % test under static horizontal load 2
% test = 'StaticHori3'; % test under static horizontal load 3 (lifting)
% test = 'StaticHori4'; % test under static horizontal load 4 (lifting)
test = 'StaticVert'; % test under static vertical load
% test = 'Fatigue1'; % fatigue test under horizontal load 1
% test = 'Fatigue2'; % fatigue test under horizontal load 2
% test = 'Fatigue3'; % fatigue test under horizontal load 3 (lifting)
% test = 'Fatigue4'; % fatigue test under horizontal load 4 (lifting)
% test = 'Impact'; % vertical impact test
% test = 'Drop'; % drop test

formats = {'fig','epsc2'};
renderer = 'OpenGL';

filename = ['FCBADesk3DDetLinElas' test];
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
    x1 = (a5+h)/2;
    y1_14 = -b12+c;
    y1_23 = c;
    z1_12 = 0;
    z1_34 = a12+h/2;
    D1 = DOMAIN(3,[x1-h/2,y1_14,z1_12],[x1+h/2,y1_23,z1_34-h/2]);
    %
    x2 = -(a5+h)/2;
    y2_14 = -b12+c;
    y2_23 = c;
    z2_12 = 0;
    z2_34 = a12+h/2;
    D2 = DOMAIN(3,[x2-h/2,y2_14,z2_12],[x2+h/2,y2_23,z2_34-h/2]);
    %
    x3_14 = -a3/2;
    x3_23 = a3/2;
    y3_12 = c-(b3+b12)/2;
    y3_34 = c+(b3-b12)/2;
    z3 = a12+h/2;
    D3 = DOMAIN(3,[x3_14,y3_12,z3-h/2],[x3_23,y3_34,z3+h/2]);
    %
    x5a_14 = -(a5+h)/2;
    x5a_23 = (a5+h)/2;
    y5a = 0;
    z5a_12 = a12-d-e-b;
    z5a_34 = a12-d+a;
    D5a = DOMAIN(3,[x5a_14+h/2,y5a-h/2,z5a_12],[x5a_23-h/2,y5a+h/2,z5a_34]);
    %
    x5b_14 = -(a5+h)/2;
    x5b_23 = (a5+h)/2;
    y5b = 0;
    z5b_12 = f-b;
    z5b_34 = f-b+b5;
    D5b = DOMAIN(3,[x5b_14+h/2,y5b-h/2,z5b_12],[x5b_23-h/2,y5b+h/2,z5b_34]);
    
    % Points
    x_hori = {[x3_23,(y3_12+y3_34)/2,z3],[x3_14,(y3_12+y3_34)/2,z3],...
        [(x3_23+x3_14)/2,y3_34,z3],[(x3_23+x3_14)/2,y3_12,z3]};
    x_vert = double(getcenter(D3))-[0.0,0.0,h/2];
    x_fati = {[x3_23,y3_12+50e-3,z3-h/2],[x3_14,y3_12+50e-3,z3-h/2],[x3_23-50e-3,y3_12,z3-h/2],[x3_23-50e-3,y3_34,z3-h/2]};
    x_stab = double(POINT([x_hori{4}(1),x_hori{4}(2)+50e-3,x_hori{4}(3)-h/2]));
    x_meas = [x_hori,x_vert,x_fati,x_stab];
    P_hori = cellfun(@(x) POINT(x),x_hori,'UniformOutput',false);
    P_vert = POINT(x_vert);
    P_fati = cellfun(@(x) POINT(x),x_fati,'UniformOutput',false);
    P_stab = POINT(x_stab);
    P_meas = cellfun(@(x) POINT(x),x_meas,'UniformOutput',false);
    
    % Plates meshes
    elemtype = 'TET4';
    cl = h;
    cl_12 = cl;
    cl_3 = cl;
    cl_5 = cl;
    r_masse = 100e-3;
    C_masse = CIRCLE(0.0,y3_12+b3/2,z3+h/2,r_masse);
    r_load = 40e-3;
    C_vert = CIRCLE(x_vert(1),x_vert(2),x_vert(3)+h,r_load);
    C_stab = CIRCLE(x_stab(1),x_stab(2),x_stab(3)+h,r_load);
    %
    Q1_a = QUADRANGLE([x5a_23-h/2,y5a+h/2,z5a_12],[x5a_23-h/2,y5a-h/2,z5a_12],[x5a_23-h/2,y5a-h/2,z5a_34],[x5a_23-h/2,y5a+h/2,z5a_34]);
    Q1_b = QUADRANGLE([x5b_23-h/2,y5b+h/2,z5b_12],[x5b_23-h/2,y5b-h/2,z5b_12],[x5b_23-h/2,y5b-h/2,z5b_34],[x5b_23-h/2,y5b+h/2,z5b_34]);
    if ~strcmp(elemtype,'CUB8')
        S1 = gmshFCBAdesk3D1(D1,Q1_a,Q1_b,cl_12,cl_5,cl_5,...
            fullfile(pathname,['gmsh_desk_1_' elemtype '_cl_' num2str(cl_12)]),3);
    else
        S1 = gmshFCBAdesk3D1(D1,Q1_a,Q1_b,cl_12,cl_5,cl_5,...
            fullfile(pathname,['gmsh_desk_1_' elemtype '_cl_' num2str(cl_12)]),3,'recombine');
    end
    S1 = convertelem(S1,elemtype);
    %
    Q2_a = QUADRANGLE([x5a_14+h/2,y5a-h/2,z5a_12],[x5a_14+h/2,y5a+h/2,z5a_12],[x5a_14+h/2,y5a+h/2,z5a_34],[x5a_14+h/2,y5a-h/2,z5a_34]);
    Q2_b = QUADRANGLE([x5b_14+h/2,y5b-h/2,z5b_12],[x5b_14+h/2,y5b+h/2,z5b_12],[x5b_14+h/2,y5b+h/2,z5b_34],[x5b_14+h/2,y5b-h/2,z5b_34]);
    if ~strcmp(elemtype,'CUB8')
        S2 = gmshFCBAdesk3D2(D2,Q2_a,Q2_b,cl_12,cl_5,cl_5,...
            fullfile(pathname,['gmsh_desk_2_' elemtype '_cl_' num2str(cl_12)]),3);
    else
        S2 = gmshFCBAdesk3D2(D2,Q2_a,Q2_b,cl_12,cl_5,cl_5,...
            fullfile(pathname,['gmsh_desk_2_' elemtype '_cl_' num2str(cl_12)]),3,'recombine');
    end    
    S2 = convertelem(S2,elemtype);
    %
    LbD3 = {LIGNE(x_hori{4}+[-r_load,0,-h/2],x_hori{4}+[-r_load,0,h/2]),...
        LIGNE(x_hori{4}+[r_load,0,-h/2],x_hori{4}+[r_load,0,h/2]),...
        LIGNE(x_fati{3}+[-r_load,0,-h/2],x_fati{3}+[-r_load,0,h/2]),...
        LIGNE(x_fati{3}+[r_load,0,-h/2],x_fati{3}+[r_load,0,h/2]),...
        LIGNE(x_fati{1}+[0,-r_load,-h/2],x_fati{1}+[0,-r_load,h/2]),...
        LIGNE(x_fati{1}+[0,r_load,-h/2],x_fati{1}+[0,r_load,h/2]),...
        LIGNE(x_hori{1}+[0,-r_load,-h/2],x_hori{1}+[0,-r_load,h/2]),...
        LIGNE(x_hori{1}+[0,r_load,-h/2],x_hori{1}+[0,r_load,h/2]),...
        LIGNE(x_fati{4}+[r_load,0,-h/2],x_fati{4}+[r_load,0,h/2]),...
        LIGNE(x_fati{4}+[-r_load,0,-h/2],x_fati{4}+[-r_load,0,h/2]),...
        LIGNE(x_hori{3}+[r_load,0,-h/2],x_hori{3}+[r_load,0,h/2]),...
        LIGNE(x_hori{3}+[-r_load,0,-h/2],x_hori{3}+[-r_load,0,h/2]),...
        LIGNE(x_hori{2}+[0,r_load,-h/2],x_hori{2}+[0,r_load,h/2]),...
        LIGNE(x_hori{2}+[0,-r_load,-h/2],x_hori{2}+[0,-r_load,h/2]),...
        LIGNE(x_fati{2}+[0,r_load,-h/2],x_fati{2}+[0,r_load,h/2]),...
        LIGNE(x_fati{2}+[0,-r_load,-h/2],x_fati{2}+[0,-r_load,h/2])};
    Q3_1 = QUADRANGLE([x1-h/2,y1_14,z1_34-h/2],[x1+h/2,y1_14,z1_34-h/2],[x1+h/2,y1_23,z1_34-h/2],[x1-h/2,y1_23,z1_34-h/2]);
    Q3_2 = QUADRANGLE([x2-h/2,y2_14,z2_34-h/2],[x2+h/2,y2_14,z2_34-h/2],[x2+h/2,y2_23,z2_34-h/2],[x2-h/2,y2_23,z2_34-h/2]);
    CiD3eI = C_stab;
    CiI = C_vert;
    if ~strcmp(elemtype,'CUB8')
        S3 = gmshFCBAdesk3D3(D3,C_masse,Q3_1,Q3_2,LbD3,CiD3eI,CiI,cl_3,cl_3,cl_12,cl_12,cl_3,cl_3,cl_3,...
            fullfile(pathname,['gmsh_desk_3_' elemtype '_cl_' num2str(cl_3)]),3);
    else
        S3 = gmshFCBAdesk3D3(D3,C_masse,Q3_1,Q3_2,LbD3,CiD3eI,CiI,cl_3,cl_3,cl_12,cl_12,cl_3,cl_3,cl_3,...
            fullfile(pathname,['gmsh_desk_3_' elemtype '_cl_' num2str(cl_3)]),3,'recombine');
    end
    S3 = convertelem(S3,elemtype);
    %
    S5a = build_model(D5a,'cl',cl_5,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_desk_5a_' elemtype '_cl_' num2str(cl_5)]));
    %
    S5b = build_model(D5b,'cl',cl_5,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_desk_5b_' elemtype '_cl_' num2str(cl_5)]));
    
    %% Materials
    % Gravitational acceleration
    g = 9.81;
    
    % Density
    Mass_total = 13.9; % kg
    Vol_total = h*(a12*b12*2+a3*b3+a5*b5*2);
    RHO = Mass_total/(Vol_total);
    
    % Data
    filenameAna = 'data_ET_GL.mat';
    filenameNum = 'data_EL_NUL.mat';
    pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','identification','materialParticleBoard');
    load(fullfile(pathnameIdentification,filenameAna));
    load(fullfile(pathnameIdentification,filenameNum));
    
    % Sample number
    sampleNum = 'B13';
    
    % Material symmetry
    materialSym = 'isot';
    
    switch lower(materialSym)
        case 'isot'
           % Young modulus
            E = eval(['mean_ET_' sampleNum '_data;'])*1e9; % Pa
            %E = 1.7e9; % Pa
            % Shear modulus
            G = eval(['mean_GL_' sampleNum '_data;'])*13*1e6; % Pa
            % Poisson ratio
            NU = E./(2*G)-1;
            %NU = 0.25
            % Material
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO);
        case 'isottrans'
           % Transverse Young modulus
            ET = eval(['mean_ET_' sampleNum '_data;'])*1e9; % Pa
            % Longitudinal shear modulus
            GL = eval(['mean_GL_' sampleNum '_data;'])*1e6; % Pa
            % Longitudinal Young modulus
            EL = eval(['mean_EL_' sampleNum '_data;'])*1e6; % Pa
            % Longitudinal Poisson ratio
            NUL = eval(['mean_NUL_' sampleNum '_data;']);
            % Transverse Poisson ratio
            NUT = 0.25;
            % Material
            mat = ELAS_ISOT_TRANS('EL',EL,'ET',ET,'NUL',NUL,'NUT',NUT,'GL',GL,'RHO',RHO);
        otherwise
            error('Wrong material symmetry !')
    end
    mat = setnumber(mat,1);
    S = union(S1,S2,S3,S5a,S5b);
    S = setmaterial(S,mat);
    
    %% Neumann boundary conditions
    p_plate = RHO*g; % body load
    Sec_stab_vert = pi*r_load^2;
    Sec_hori_fati = 2*r_load*h;
    switch lower(test)
        case 'stability'
            p = 400/Sec_stab_vert; % surface load
        case {'statichori1','statichori2','statichori3','statichori4'}
            masse = 50.5;
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse; % surface load
            p = 100/Sec_hori_fati; % surface load, F1=F2=100N or 200N, F3=F4=100N
            slope = 0;
        case 'staticvert'
            p = 300/Sec_stab_vert; % surface load, 300N, 400N or 500N
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            masse = 50.5;
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse; % surface load 
            p = 100/Sec_hori_fati; % surface load
        case 'impact'
            H = 180e-3;
        case 'drop'
            H = 100e-3;
    end
    
    %% Dirichlet boundary conditions
    S = final(S);
    S1_5 = getfacet(S1,5);
    S2_5 = getfacet(S2,5);
    S5b_1 = getfacet(S5b,1);
    [~,numnode1] = intersect(S,S1_5);
    [~,numnode2] = intersect(S,S2_5);
    [~,numnode5b] = intersect(S,S5b_1);
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
    
    %% Stiffness matrices and sollicitation vectors
    A = calc_rigi(S);
    
    switch lower(test)
        case 'stability'
            S3_26 = getfacet(S3,26);
            f = surfload(S,S3_26,'FZ',-p);
        case {'statichori1','statichori2','statichori3','statichori4'}
            if strcmpi(test,'statichori1')
                S3_9 = getfacet(S3,9);
                f = surfload(S,S3_9,{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
            elseif strcmpi(test,'statichori2')
                S3_17 = getfacet(S3,17);
                f = surfload(S,S3_17,{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
            elseif strcmpi(test,'statichori3')
                S3_14 = getfacet(S3,14);
                f = surfload(S,S3_14,{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
            elseif strcmpi(test,'statichori4')
                S3_2 = getfacet(S3,2);
                f = surfload(S,S3_2,{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
            end
            S3_masse = union(getfacet(S3,24),getfacet(S3,25));
            f = f + surfload(S,S3_masse,'FZ',-p_masse);
        case 'staticvert'
            S3_24 = getfacet(S3,24);
            f = surfload(S,S3_24,'FZ',-p);
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            if strcmpi(test,'fatigue1')
                S3_7 = getfacet(S3,7);
                f = surfload(S,S3_7,'FX',-p);
            elseif strcmpi(test,'fatigue2')
                S3_19 = getfacet(S3,19);
                f = surfload(S,S3_19,'FX',p);
            elseif strcmpi(test,'fatigue3')
                S3_4 = getfacet(S3,4);
                f = surfload(S,S3_4,'FY',p);
            elseif strcmpi(test,'fatigue4')
                S3_12 = getfacet(S3,12);
                f = surfload(S,S3_12,'FY',-p);
            end
            S3_masse = union(getfacet(S3,24),getfacet(S3,25));
            f = f + surfload(S,S3_masse,'FZ',-p_masse);
        case {'impact','drop'}
            error('Not implemented')
    end
    f = f + bodyload(S,[],'FZ',-p_plate);
    
    %% Solution
    t = tic;
    u = A\f;
    time = toc(t);
    
    u = unfreevector(S,u);
    
    U = u(findddl(S,DDL(DDLVECT('U',S.syscoord))),:);
    Ux = u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,u,S.node,'UX')));
    Uy = u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,u,S.node,'UY')));
    Uz = u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,u,S.node,'UZ')));
    
    %% Test solution
    ux_P_vert = eval_sol(S,u,P_vert,'UX');
    uy_P_vert = eval_sol(S,u,P_vert,'UY');
    uz_P_vert = eval_sol(S,u,P_vert,'UZ');
    ux_P_stab = eval_sol(S,u,P_stab,'UX');
    uy_P_stab = eval_sol(S,u,P_stab,'UY');
    uz_P_stab = eval_sol(S,u,P_stab,'UZ');
    ux_P_hori(1) = eval_sol(S,u,P_hori{1},'UX');
    uy_P_hori(1) = eval_sol(S,u,P_hori{1},'UY');
    uz_P_hori(1) = eval_sol(S,u,P_hori{1},'UZ');
    ux_P_hori(2) = eval_sol(S,u,P_hori{2},'UX');
    uy_P_hori(2) = eval_sol(S,u,P_hori{2},'UY');
    uz_P_hori(2) = eval_sol(S,u,P_hori{2},'UZ');
    ux_P_hori(3) = eval_sol(S,u,P_hori{3},'UX');
    uy_P_hori(3) = eval_sol(S,u,P_hori{3},'UY');
    uz_P_hori(3) = eval_sol(S,u,P_hori{3},'UZ');
    ux_P_hori(4) = eval_sol(S,u,P_hori{4},'UX');
    uy_P_hori(4) = eval_sol(S,u,P_hori{4},'UY');
    uz_P_hori(4) = eval_sol(S,u,P_hori{4},'UZ');
    ux_P_fati(1) = eval_sol(S,u,P_fati{1},'UX');
    uy_P_fati(1) = eval_sol(S,u,P_fati{1},'UY');
    uz_P_fati(1) = eval_sol(S,u,P_fati{1},'UZ');
    ux_P_fati(2) = eval_sol(S,u,P_fati{2},'UX');
    uy_P_fati(2) = eval_sol(S,u,P_fati{2},'UY');
    uz_P_fati(2) = eval_sol(S,u,P_fati{2},'UZ');
    ux_P_fati(3) = eval_sol(S,u,P_fati{3},'UX');
    uy_P_fati(3) = eval_sol(S,u,P_fati{3},'UY');
    uz_P_fati(3) = eval_sol(S,u,P_fati{3},'UZ');
    ux_P_fati(4) = eval_sol(S,u,P_fati{4},'UX');
    uy_P_fati(4) = eval_sol(S,u,P_fati{4},'UY');
    uz_P_fati(4) = eval_sol(S,u,P_fati{4},'UZ');
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','elemtype',...
        'a12','b12','a3','b3','a5','b5','h','Sec_stab_vert','Sec_hori_fati',...
        'f','p');
    save(fullfile(pathname,'solution.mat'),'u','time',...
        'U','Ux','Uy','Uz');
    save(fullfile(pathname,'test_solution.mat'),'P_vert','P_stab',...
        'P_hori','P_fati',...
        'ux_P_vert','uy_P_vert','uz_P_vert',...
        'ux_P_stab','uy_P_stab','uz_P_stab',...
        'ux_P_hori','uy_P_hori','uz_P_hori',...        
        'ux_P_fati','uy_P_fati','uz_P_fati');
else
    load(fullfile(pathname,'problem.mat'),'S','elemtype',...
        'a12','b12','a3','b3','a5','b5','h','f','p','Sec_stab_vert','Sec_hori_fati');
    load(fullfile(pathname,'solution.mat'),'u','time',...
        'U','Ux','Uy','Uz');
    load(fullfile(pathname,'test_solution.mat'),'P_vert','P_stab',...
        'P_hori','P_fati',...
        'ux_P_vert','uy_P_vert','uz_P_vert',...
        'ux_P_stab','uy_P_stab','uz_P_stab',...
        'ux_P_hori','uy_P_hori','uz_P_hori',...        
        'ux_P_fati','uy_P_fati','uz_P_fati');
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
        disp('Displacement u at point'); disp(P_vert);
        fprintf('ux = %g\n',ux_P_vert);
        fprintf('uy = %g\n',uy_P_vert);
        fprintf('uz = %g\n',uz_P_vert);
        if p == 300/Sec_stab_vert
            uz_exp_start = -0.69*1e-3;
            uz_exp_end = -[10.10 9.88 9.64 9.88 9.94 9.79 9.92 9.93 9.82 9.95]*1e-3;
        elseif p == 400/Sec_stab_vert
            uz_exp_start = -0.75*1e-3;
            uz_exp_end = -[13.45 13.52 13.56 13.64 13.65 13.74 13.75 13.44 13.74 13.53]*1e-3;
        elseif p == 500/Sec_stab_vert
            uz_exp_start = -0.78*1e-3;
            uz_exp_end = -[16.66 16.57 16.59 16.78 16.55 16.69 16.75 16.59 16.73 16.76]*1e-3;
        end
        uz_exp = mean(uz_exp_end - uz_exp_start);
        err_uz = norm(uz_P_vert-uz_exp)/norm(uz_exp);
        fprintf('uz_exp   = %g, error    = %.3e\n',uz_exp,err_uz);
        fprintf('\n');
    case 'stability'
        disp('Displacement u at point'); disp(P_stab);
        fprintf('ux = %g\n',ux_P_stab);
        fprintf('uy = %g\n',uy_P_stab);
        fprintf('uz = %g\n',uz_P_stab);
        uz_exp_start = -1.93*1e-3;
        uz_exp_end = -[18.46 18.44 18.53 18.58 18.59 18.7 18.77 18.73 18.85 18.76]*1e-3;
        uz_exp = mean(uz_exp_end - uz_exp_start);
        err_uz = norm(uz_P_stab-uz_exp)/norm(uz_exp);
        fprintf('uz_exp   = %g, error    = %.3e\n',uz_exp,err_uz);
        fprintf('\n');
    case 'statichori1'
        disp('Displacement u at point'); disp(P_hori{2});
        fprintf('ux = %g\n',ux_P_hori(2));
        fprintf('uy = %g\n',uy_P_hori(2));
        fprintf('uz = %g\n',uz_P_hori(2));
        if p==100/Sec_hori_fati
            ux_exp_start = -6.88*1e-3;
            ux_exp_end = -[10.5 10.51 10.44 10.8 10.72 10.62 10.67 10.65 10.66 10.87 10.86]*1e-3;
        elseif p==200/Sec_hori_fati
            ux_exp_start = -6.16*1e-3;
            ux_exp_end = -[16.78 16.74 16.72 17.13 17 16.8 16.87 16.78 17.04 16.82 16.71 17.17]*1e-3;
        end
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux_P_hori(2)-ux_exp)/norm(ux_exp);
        fprintf('ux_exp   = %g, error    = %.3e\n',ux_exp,err_ux);
        fprintf('\n');
    case 'statichori2'
        disp('Displacement u at point'); disp(P_hori{1});
        fprintf('ux = %g\n',ux_P_hori(1));
        fprintf('uy = %g\n',uy_P_hori(1));
        fprintf('uz = %g\n',uz_P_hori(1));
        if p==100/Sec_hori_fati
            ux_exp_start = 2.12*1e-3;
            ux_exp_end = [6.22 6.17 6.26 6.31 6.33 6.24 6.26 6.4 6.26 6.49 6.48 6.42 6.36 6.56 6.37 6.39]*1e-3;
        elseif p==200/Sec_hori_fati
            ux_exp_start = 1.91*1e-3;
            ux_exp_end = [12.45 12.68 12.66 12.65 12.71 12.64 12.82 12.73 12.89 12.86 12.79 12.86]*1e-3;
        end
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux_P_hori(1)-ux_exp)/norm(ux_exp);
        fprintf('ux_exp   = %g, error    = %.3e\n',ux_exp,err_ux);
        fprintf('\n'); 
    case 'statichori3'
        disp('Displacement u at point'); disp(P_hori{4});
        fprintf('ux = %g\n',ux_P_hori(4));
        fprintf('uy = %g\n',uy_P_hori(4));
        fprintf('uz = %g\n',uz_P_hori(4));
        uy_exp_start = -3.77*1e-3;
        uy_exp_end = -[4.71 4.73 4.69 4.56 4.47 4.73]*1e-3;   
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy_P_hori(4)-uy_exp)/norm(uy_exp);
        fprintf('uy_exp   = %g, error    = %.3e\n',uy_exp,err_uy);
        fprintf('\n'); 
    case 'statichori4'
        disp('Displacement u at point'); disp(P_hori{3});
        fprintf('ux = %g\n',ux_P_hori(3));
        fprintf('uy = %g\n',uy_P_hori(3));
        fprintf('uz = %g\n',uz_P_hori(3));
        uy_exp_start = 9.71*1e-3;
        uy_exp_end = [12.21 12.2 12.2 12.23 12.2 12.19 12.21]*1e-3;   
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy_P_hori(3)-uy_exp)/norm(uy_exp);
        fprintf('uy_exp   = %g, error    = %.3e\n',uy_exp,err_uy);
        fprintf('\n');
    case 'fatigue1'
        disp('Displacement u at point'); disp(P_fati{2});
        fprintf('ux = %g\n',ux_P_fati(2));
        fprintf('uy = %g\n',uy_P_fati(2));
        fprintf('uz = %g\n',uz_P_fati(2));
        ux_exp_start = -4.42*1e-3;
        ux_exp_end = -[8.4 8.3 8.37 8.41 8.54 8.39 8.56 8.48 8.46 8.49 8.49 8.43 8.55 8.52]*1e-3;   
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux_P_fati(2)-ux_exp)/norm(ux_exp);
        fprintf('ux_exp   = %g, error    = %.3e\n',ux_exp,err_ux);
        fprintf('\n');
    case 'fatigue2'
        disp('Displacement u at point'); disp(P_fati{1});
        fprintf('ux = %g\n',ux_P_fati(1));
        fprintf('uy = %g\n',uy_P_fati(1));
        fprintf('uz = %g\n',uz_P_fati(1));
        ux_exp_start = 3.48*1e-3;
        ux_exp_end = [7.89 7.85 8.1 8.4 8.36 8.55 8.27 8.27 8.47 8.49 8.64 8.35 8.5 8.63 8.73]*1e-3;   
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(ux_P_fati(1)-ux_exp)/norm(ux_exp);
        fprintf('ux_exp   = %g, error    = %.3e\n',ux_exp,err_ux);
        fprintf('\n'); 
    case 'fatigue3'
        disp('Displacement u at point'); disp(P_fati{4});
        fprintf('ux = %g\n',ux_P_fati(4));
        fprintf('uy = %g\n',uy_P_fati(4));
        fprintf('uz = %g\n',uz_P_fati(4));
        uy_exp_start = 3.35*1e-3;
        uy_exp_end = [6.16 5.76 5.97 5.81 5.84 5.61 5.86 5.64 5.62 5.68]*1e-3;   
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy_P_fati(1)-uy_exp)/norm(uy_exp);
        fprintf('uy_exp   = %g, error    = %.3e\n',uy_exp,err_uy);
        fprintf('\n');
     case 'fatigue4'
        disp('Displacement u at point'); disp(P_fati{3});
        fprintf('ux = %g\n',ux_P_fati(3));
        fprintf('uy = %g\n',uy_P_fati(3));
        fprintf('uz = %g\n',uz_P_fati(3));
        uy_exp_start = -3.75*1e-3;
        uy_exp_end = -[3.89 3.88 3.89 3.88 3.89]*1e-3;   
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(uy_P_fati(1)-uy_exp)/norm(uy_exp);
        fprintf('uy_exp   = %g, error    = %.3e\n',uy_exp,err_uy);
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
    % legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    plotModel(S,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    ampl = getsize(S)/max(abs(u))/10;
    plotModelDeflection(S,u,'ampl',ampl,'FaceColor','b','legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S);
    plot(S+ampl*u,'FaceColor','b');
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
end
