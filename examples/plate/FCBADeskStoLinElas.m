%% FCBA desk stochastic linear elasticity %%
%%----------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');
myparallel('start');

%% Input data
solveProblem = true;
displaySolution = true;

% test = 'Stability'; % stability test under vertical load
% test = 'StaticHori1'; % test under static horizontal load 1
% test = 'StaticHori2'; % test under static horizontal load 2
% test = 'StaticHori3'; % test under static horizontal load 3
% test = 'StaticHori4'; % test under static horizontal load 4
% test = 'StaticVert'; % test under static vertical load
test = 'Fatigue1'; % fatigue test under horizontal load 1
% test = 'Fatigue2'; % fatigue test under horizontal load 2
% test = 'Fatigue3'; % fatigue test under horizontal load 3
% test = 'Fatigue4'; % fatigue test under horizontal load 4
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
    elemtype = 'DKT';
    cl_plate12 = b12/10;
    cl_plate3 = b3/10;
    cl_plate5 = b5/4;
    r_masse = 150e-3;
    C_masse = CIRCLE(0.0,y12_S3+1/2*b3,z_S3,r_masse);
    %
    vertices_S1 = getvertices(S1);
    PbQ_S1 = {vertices_S1{1},[x23_S5b,y_S5b,z12_S5b],vertices_S1{2},vertices_S1{3},vertices_S1{4}};
    PL_S1 = {[x23_S5a,y_S5a,z12_S5a],[x23_S5a,y_S5a,z34_S5a],[x23_S5b,y_S5b,z34_S5b]};
    S1_plate = gmshFCBAdesk12(S1,PL_S1,PbQ_S1,cl_plate12,cl_plate5,cl_plate12,[pathname 'gmsh_desk_1_' elemtype '_cl_' num2str(cl_plate12)],3);
    S1_plate = convertelem(S1_plate,elemtype);
    %
    vertices_S2 = getvertices(S2);
    PbQ_S2 = {vertices_S2{1},[x14_S5b,y_S5b,z12_S5b],vertices_S2{2},vertices_S2{3},vertices_S2{4}};
    PL_S2 = {[x14_S5a,y_S5a,z12_S5a],[x14_S5a,y_S5a,z34_S5a],[x14_S5b,y_S5b,z34_S5b]};
    S2_plate = gmshFCBAdesk12(S2,PL_S2,PbQ_S2,cl_plate12,cl_plate5,cl_plate12,[pathname 'gmsh_desk_2_' elemtype '_cl_' num2str(cl_plate12)],3);
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
    S3_plate = gmshFCBAdesk3(S3,C_masse,PL_S3,PbQ_S3,PiQeI,PiI,cl_plate3,cl_plate3,cl_plate12,cl_plate3,cl_plate3,cl_plate3,[pathname 'gmsh_desk_3_' elemtype '_cl_' num2str(cl_plate3)],3);
    S3_plate = convertelem(S3_plate,elemtype);
    %
    S5a_plate = build_model(S5a,'cl',cl_plate5,'elemtype',elemtype,'filename',[pathname 'gmsh_desk_5a_' elemtype '_cl_' num2str(cl_plate5)]);
    S5a_plate = convertelem(S5a_plate,elemtype);
    %
    S5b_plate = build_model(S5b,'cl',cl_plate5,'elemtype',elemtype,'filename',[pathname 'gmsh_desk_5b_' elemtype '_cl_' num2str(cl_plate5)]);
    S5b_plate = convertelem(S5b_plate,elemtype);
    
    %% Random variables
    % Experimental data
    samples_E = [4.211 4.057 3.685 3.921 3.839 3.845 3.795...
        3.406 3.389 3.299 3.485 3.319 3.267 3.349 3.307...
        4.684 4.245 4.076 4.407 4.283 4.054 4.226 4.041...
        4.104 4.075 3.556 3.319 3.848 3.707 3.664 3.493 3.550]*1e9;
    % samples_E = [2.7116 2.8945 2.9689 3.022 3.014 3.0069 2.9952 2.9711...
    %     2.7807 2.9541 2.9618 3.0014 2.9562 2.8909 2.8874 2.8923...
    %     2.6636 3.2635 3.1073 3.1657 3.1424 3.1731 3.1416 3.155...
    %     2.7356 2.8797 2.7230 2.7851 2.8312 2.8018 2.7592 2.743 2.7241 2.7055 2.7073...
    %     3.8267 3.2860 3.3753 3.1742 3.3089 3.2461 3.1331 3.0817	3.0093 3.0528]*1e9;
    % Parameters for Gamma distribution
    phat = gamfit(samples_E);
    % Number of samples
    N = 1e3;
    % Sample set
    e = gamrnd(phat(1),phat(2),N,1);
    
    %% Materials
    % Gravitational acceleration
    g = 9.81;
    
    % Plate
    % Young modulus
    E = mean(e);
    % Poisson ratio
    NU = 0.3;
    % Density
    Mass_plates_total = 13.9; % kg
    Volum_plates_total = h*(a12*b12*2+a3*b3+a5*b5*2);
    RHO = Mass_plates_total/(Volum_plates_total);
    % Material
    mat_plate = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',h,'k',5/6);
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
            masse = 50;
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse;
            p = 400;
            slope = 0;
        case 'staticvert'
            p = 1200;
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            masse = 50;
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse;
            p = 300;
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
            S = addcl(S,P1_S1);
            S = addcl(S,P1_S2,'U');
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
    
    %% Stiffness matrices
    A = cell(N,1);
    for i=1:N
        % Young modulus
        Ei = e(i);
        % Material
        mat_platei = setparam(mat_plate,'E',Ei);
        Si = setmaterial(S,mat_platei);
        % Stiffness matrix
        A{i} = calc_rigi(Si);
    end
    
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
    parfor i=1:N
        u(:,i) = A{i}\f;
    end
    time = toc(t);
    
    mean_u = mean(u,2);
    mean_u = unfreevector(S,mean_u);
    
    std_u = std(u,0,2);
    std_u = unfreevector(S,std_u);
    
    mean_U = mean_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    mean_Ux = mean_u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,mean_u,S.node,'UX')),:);
    mean_Uy = mean_u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,mean_u,S.node,'UY')),:);
    mean_Uz = mean_u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,mean_u,S.node,'UZ')),:);
    
    mean_R = mean_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    mean_Rx = mean_u(findddl(S,'RX'),:); % Rx = double(squeeze(eval_sol(S,mean_u,S.node,'RX')),:);
    mean_Ry = mean_u(findddl(S,'RY'),:); % Ry = double(squeeze(eval_sol(S,mean_u,S.node,'RY')),:);
    mean_Rz = mean_u(findddl(S,'RZ'),:); % Rz = double(squeeze(eval_sol(S,mean_u,S.node,'RZ')),:);
    
    std_U = std_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    std_Ux = std_u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,std_u,S.node,'UX')),:);
    std_Uy = std_u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,std_u,S.node,'UY')),:);
    std_Uz = std_u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,std_u,S.node,'UZ')),:);
    
    std_R = std_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    std_Rx = std_u(findddl(S,'RX'),:); % Rx = double(squeeze(eval_sol(S,std_u,S.node,'RX')),:);
    std_Ry = std_u(findddl(S,'RY'),:); % Ry = double(squeeze(eval_sol(S,std_u,S.node,'RY')),:);
    std_Rz = std_u(findddl(S,'RZ'),:); % Rz = double(squeeze(eval_sol(S,std_u,S.node,'RZ')),:);
    
    %% Test solution
    P = getcenter(S3);
    
    mean_ux = eval_sol(S,mean_u,P,'UX');
    mean_uy = eval_sol(S,mean_u,P,'UY');
    mean_uz = eval_sol(S,mean_u,P,'UZ');
    
    mean_rx = eval_sol(S,mean_u,P,'RX');
    mean_ry = eval_sol(S,mean_u,P,'RY');
    mean_rz = eval_sol(S,mean_u,P,'RZ');
    
    std_ux = eval_sol(S,std_u,P,'UX');
    std_uy = eval_sol(S,std_u,P,'UY');
    std_uz = eval_sol(S,std_u,P,'UZ');
    
    std_rx = eval_sol(S,std_u,P,'RX');
    std_ry = eval_sol(S,std_u,P,'RY');
    std_rz = eval_sol(S,std_u,P,'RZ');
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','S1_plate','S2_plate','S3_plate','S5a_plate','S5b_plate','elemtype','a12','b12','a3','b3','a5','b5','h','f');
    save(fullfile(pathname,'solution.mat'),'u','mean_u','std_u','time',...
        'mean_U','mean_Ux','mean_Uy','mean_Uz',...
        'mean_R','mean_Rx','mean_Ry','mean_Rz',...
        'std_U','std_Ux','std_Uy','std_Uz',...
        'std_R','std_Rx','std_Ry','std_Rz');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'mean_ux','mean_uy','mean_uz',...
        'mean_rx','mean_ry','mean_rz',...
        'std_ux','std_uy','std_uz',...
        'std_rx','std_ry','std_rz');
else
    load(fullfile(pathname,'problem.mat'),'S','S1_plate','S2_plate','S3_plate','S5a_plate','S5b_plate','elemtype','a12','b12','a3','b3','a5','b5','h','f');
    load(fullfile(pathname,'solution.mat'),'u','mean_u','std_u','time',...
        'mean_U','mean_Ux','mean_Uy','mean_Uz',...
        'mean_R','mean_Rx','mean_Ry','mean_Rz',...
        'std_U','std_Ux','std_Uy','std_Uz',...
        'std_R','std_Rx','std_Ry','std_Rz');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'mean_ux','mean_uy','mean_uz',...
        'mean_rx','mean_ry','mean_rz',...
        'std_ux','std_uy','std_uz',...
        'std_rx','std_ry','std_rz');
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

disp('Displacement u at point'); disp(P);
fprintf('mean(ux) = %g, std(ux) = %g\n',mean_ux,std_ux);
fprintf('mean(uy) = %g, std(uy) = %g\n',mean_uy,std_uy);
fprintf('mean(uz) = %g, std(uz) = %g\n',mean_uz,std_uz);
if strcmpi(test,'staticvert')
    uz_exp = -2.35e-3;
    err_uz = norm(mean_uz-uz_exp)/norm(uz_exp);
    fprintf('uz_exp   = %g, error    = %.3e\n',uz_exp,err_uz);
end
fprintf('\n');

disp('Rotation r at point'); disp(P);
fprintf('mean(rx) = %g, std(rx) = %g\n',mean_rx,std_rx);
fprintf('mean(ry) = %g, std(ry) = %g\n',mean_ry,std_ry);
fprintf('mean(rz) = %g, std(rz) = %g\n',mean_rz,std_rz);
fprintf('\n');

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
    
    ampl = getsize(S)/max(abs(mean_u))/10;
    plotModelDeflection(S,mean_u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'node',true,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'node',true);
    plot(S+ampl*mean_u,'Color','b','FaceColor','b','FaceAlpha',0.1,'node',true);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    %% Display solution
    % ampl = 0;
    ampl = getsize(S)/max(abs(mean_u))/10;
    options = {'solid',true};
    % options = {};
    
    switch lower(test)
        case 'stability'
            plotSolution(S,mean_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uz',formats,renderer);
            
            plotSolution(S,std_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uz',formats,renderer);
        case {'statichori1','statichori2','statichori3','statichori4'}
            plotSolution(S,mean_u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Ux',formats,renderer);
            
            plotSolution(S,std_u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Ux',formats,renderer);
        case 'staticvert'
            plotSolution(S,mean_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uz',formats,renderer);
            
            plotSolution(S,std_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uz',formats,renderer);
        case {'fatigue1','fatigue2'}
            plotSolution(S,mean_u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Ux',formats,renderer);
            
            plotSolution(S,std_u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Ux',formats,renderer);
        case {'fatigue3','fatigue4'}
            plotSolution(S,mean_u,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uy',formats,renderer);
            
            plotSolution(S,std_u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uy',formats,renderer);
        case 'impact'
            plotSolution(S,mean_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uz',formats,renderer);
            
            plotSolution(S,std_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uz',formats,renderer);
        case 'drop'
            plotSolution(S,mean_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uz',formats,renderer);
            
            plotSolution(S,std_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uz',formats,renderer);
    end
    
    % plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'Rx',formats,renderer);
    %
    % plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'Ry',formats,renderer);
end
%% Display convergence Monte-Carlo
if displayCv
    means_u = arrayfun(@(x) norm(mean(u(:,1:x),2)),1:N);
    stds_u = arrayfun(@(x) norm(std(u(:,1:x),0,2)),1:N);
    
    figure('Name','Convergence empirical mean')
    clf
    plot(1:N,means_u,'-b','LineWidth',1)
    grid on
    box on
    set(gca,'FontSize',16)
    % xlabel('Nombre de r\''ealisations','Interpreter','latex')
    % ylabel('Moyenne empirique','Interpreter','latex')
    xlabel('Number of samples','Interpreter','latex')
    ylabel('Empirical mean','Interpreter','latex')
    mysaveas(pathname,'convergence_empirical_mean','fig');
    mymatlab2tikz(pathname,'convergence_empirical_mean.tex');
    
    figure('Name','Convergence empirical standard deviation')
    clf
    plot(1:N,stds_u,'-r','LineWidth',1)
    grid on
    box on
    set(gca,'FontSize',16)
    % xlabel('Nombre de r\''ealisations','Interpreter','latex')
    % ylabel('Ecart-type empirique','Interpreter','latex')
    xlabel('Number of samples','Interpreter','latex')
    ylabel('Empirical standard deviation','Interpreter','latex')
    mysaveas(pathname,'convergence_empirical_std','fig');
    mymatlab2tikz(pathname,'convergence_empirical_std.tex');
end

% myparallel('stop');
