%% FCBA table circular stochastic linear elasticity %%
%%--------------------------------------------------%%

% clc
clearvars
close all
% rng('default');
myparallel('start');

%% Input data
solveProblem = true;
displaySolution = true;
displayCv = true;

% tests = {'Stability1'}; % stability test under vertical load 1
% tests = {'Stability2'}; % stability test under vertical load 2
% tests = {'Stability3'}; % stability test under vertical load 3
% tests = {'Stability4'}; % stability test under vertical load 4
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
% tests = {'Stability1','Stability2','Stability3','Stability4','StaticVert',...
%     'StaticHori1','StaticHori2','StaticHori3','StaticHori4',...
%     'Fatigue1','Fatigue2','Fatigue3','Fatigue4'};

belt = true; % belt modeling

formats = {'fig','epsc'};
renderer = 'OpenGL';

for it=1:length(tests)
    test = tests{it};
filename = ['FCBATableCircStoLinElas' test];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','plate',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Problem
if solveProblem
    %% Domains and meshes
    % Plate
    % Radius
    r = 600e-3;
    % Thickness
    h = 25e-3;
    C = CIRCLE(0.0,0.0,0.0,r);
    
    % Beams and Belt
    % Cross-section base
    b_beam_top = 48e-3;
    b_beam_bot = 38e-3;
    b_beam = (b_beam_top+b_beam_bot)/2;
    b_belt = 30e-3;
    % Cross-section height
    h_beam_top = 48e-3;
    h_beam_bot = 38e-3;
    h_beam = (h_beam_top+h_beam_bot)/2;
    h_belt = 80e-3;
    % Length
    l = 710e-3+h/2;
    a = 800e-3-h_beam_top;
    b = 800e-3-b_beam_top;
    L_beam{1} = LIGNE([-a/2,-b/2,0.0],[-a/2,-b/2,-l]);
    L_beam{2} = LIGNE([a/2,-b/2,0.0],[a/2,-b/2,-l]);
    L_beam{3} = LIGNE([a/2,b/2,0.0],[a/2,b/2,-l]);
    L_beam{4} = LIGNE([-a/2,b/2,0.0],[-a/2,b/2,-l]);
    Q_belt = QUADRANGLE([-a/2,-b/2,0.0],[a/2,-b/2,0.0],[a/2,b/2,0.0],[-a/2,b/2,0.0]);
    % L_belt{1} = LIGNE([-a/2,-b/2,0.0],[a/2,-b/2,0.0]);
    % L_belt{2} = LIGNE([a/2,-b/2,0.0],[a/2,b/2,0.0]);
    % L_belt{3} = LIGNE([a/2,b/2,0.0],[-a/2,b/2,0.0]);
    % L_belt{4} = LIGNE([-a/2,b/2,0.0],[-a/2,-b/2,0.0]);
    L_belt = getedges(Q_belt);
    
    % Points
    x_beam = cellfun(@(L) getvertex(L,1),L_beam,'UniformOutput',false);
    x_load_stab = {[-r+50e-3,0.0,0.0],[0.0,-r+50e-3,0.0],[r-50e-3,0.0,0.0],[0.0,r-50e-3,0.0]}; % stability test under vertical load
    x_load_hori = {getvertex(C,4),getvertex(C,2),getvertex(C,3),getvertex(C,1)}; % test under static horizontal load
    x_load_vert = double(getcenter(C)); % test under static vertical load
    x_load_fati = {getvertex(C,4),getvertex(C,2),[-str2double(num2str(sqrt(r^2-(r-50e-3)^2))),r-50e-3,0.0],[str2double(num2str(sqrt(r^2-(r-50e-3)^2))),r-50e-3,0.0]}; % fatigue test under horizontal load
    x_load = [x_load_stab,x_load_hori,x_load_vert,x_load_fati];
    P_beam = cellfun(@(x) POINT(x),x_beam,'UniformOutput',false);
    P_load_stab = cellfun(@(x) POINT(x),x_load_stab,'UniformOutput',false);
    P_load_hori = cellfun(@(x) POINT(x),x_load_hori,'UniformOutput',false);
    P_load_vert = POINT(x_load_vert);
    P_load_fati = cellfun(@(x) POINT(x),x_load_fati,'UniformOutput',false);
    P_load = cellfun(@(x) POINT(x),x_load,'UniformOutput',false);
    
    % Plate mesh
    cl_plate = r/10;
    cl_belt = cl_plate;
    elemtype = 'DKT';
    r_masse = 150e-3;
    C_masse = CIRCLE(0.0,0.0,0.0,r_masse);
    Pb = {getvertex(C,1),getvertex(C,2),getvertex(C,3),x_load_fati{4},getvertex(C,4),x_load_fati{3}};
    Pe = x_load_stab;
    Pi = double(getcoord(getcenter(C)));
    if ~strcmp(elemtype,'QUA4') && ~strcmp(elemtype,'CUB8') && ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4') && ~strcmp(elemtype,'STOKES')
        S_plate = gmshFCBAtablecirc(C,Q_belt,C_masse,Pb,Pe,[],Pi,cl_plate,cl_belt,cl_plate,cl_plate,cl_plate,cl_plate,cl_plate,fullfile(pathname,['gmsh_plate_circ_' elemtype '_cl_' num2str(cl_plate)]),3);
    else
        S_plate = gmshFCBAtablecirc(C,Q_belt,C_masse,Pb,Pe,[],Pi,cl_plate,cl_belt,cl_plate,cl_plate,cl_plate,cl_plate,cl_plate,fullfile(pathname,['gmsh_plate_circ_' elemtype '_cl_' num2str(cl_plate)]),3,'recombine');
    end
    S_plate = convertelem(S_plate,elemtype);
    
    % Beams meshes
    nbelem_beam = 80;
    S_beam = cellfun(@(L) build_model(L,'nbelem',nbelem_beam,'elemtype','BEAM'),L_beam,'UniformOutput',false);
    % cl_beam = l/80;
    % S_beam = cellfun(@(L,n) build_model(L,'cl',cl_beam,'elemtype','BEAM','filename',fullfile(pathname,['gmsh_beam_' num2str(n) '_cl_' num2str(cl_beam)])),L_beam,num2cell(1:length(L_beam)),'UniformOutput',false);
    
    % Belt mesh
    S_belt = cellfun(@(L,n) build_model(L,'cl',cl_belt,'elemtype','BEAM','filename',fullfile(pathname,['gmsh_belt_' num2str(n) '_cl_' num2str(cl_belt)])),L_belt,num2cell(1:length(L_belt)),'UniformOutput',false);
    
    %% Random variables
    % Data
    E_data = [4.211 4.057 3.685 3.921 3.839 3.845 3.795...
        3.406 3.389 3.299 3.485 3.319 3.267 3.349 3.307...
        4.684 4.245 4.076 4.407 4.283 4.054 4.226 4.041...
        4.104 4.075 3.556 3.319 3.848 3.707 3.664 3.493 3.550]*1e9; % GPa
    % E_data = [2.7116 2.8945 2.9689 3.022 3.014 3.0069 2.9952 2.9711...
    %     2.7807 2.9541 2.9618 3.0014 2.9562 2.8909 2.8874 2.8923...
    %     2.6636 3.2635 3.1073 3.1657 3.1424 3.1731 3.1416 3.155...
    %     2.7356 2.8797 2.7230 2.7851 2.8312 2.8018 2.7592 2.743 2.7241 2.7055 2.7073...
    %     3.8267 3.2860 3.3753 3.1742 3.3089 3.2461 3.1331 3.0817 3.0093 3.0528]*1e9; % GPa
    
    % Maximum likelihood estimation
    % Parameters for Gamma distribution
    phat = gamfit(E_data);
    a = phat(1);
    b = phat(2);
    
    % Number of samples
    N = 1e3;
    % Sample set
    E_sample = gamrnd(a,b,N,1);
    
    %% Materials
    % Gravitational acceleration
    g = 9.81;
    
    % Plate
    % Young modulus
    E = mean(E_sample);
    % Poisson ratio
    NU = 0.3;
    % Density
    mass_plate = 18.54;
    Sec_plate = pi*r^2;
    RHO = mass_plate/(Sec_plate*h);
    % Material
    mat_plate = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',h,'k',5/6);
    mat_plate = setnumber(mat_plate,1);
    S_plate = setmaterial(S_plate,mat_plate);
    
    % Beams and Belt
    % Young modulus
    E_beam = 15e9;
    % Poisson ratio
    NU_beam = 0.3;
    % Cross-section area
    Sec_beam_top = b_beam_top*h_beam_top;
    Sec_beam_bot = b_beam_bot*h_beam_bot;
    Sec_beam = b_beam*h_beam;
    Sec_belt = b_belt*h_belt;
    % Density
    Vol_beam = (l-h/2)*(Sec_beam_top+Sec_beam_bot+sqrt(Sec_beam_top*Sec_beam_bot))/3;
    Vol_belt = 2*(a-h_beam_top)*b_belt*h_belt + 2*(b-b_beam_top)*b_belt*h_belt;
    mass_beams = 8.48;
    RHO_beam = mass_beams/(length(L_beam)*Vol_beam + Vol_belt);
    % Planar second moment of area (or Planar area moment of inertia)
    IY = h_beam*b_beam^3/12;
    IZ = b_beam*h_beam^3/12;
    IY_belt = h_belt*b_belt^3/12;
    IZ_belt = b_belt*h_belt^3/12;
    % Polar second moment of area (or Polar area moment of inertia)
    IX = IY+IZ;
    IX_belt = IY_belt+IZ_belt;
    % Material
    mat_beam = ELAS_BEAM('E',E_beam,'NU',NU_beam,'S',Sec_beam,'IZ',IZ,'IY',IY,'IX',IX,'RHO',RHO_beam);
    mat_belt{1} = ELAS_BEAM('E',E_beam,'NU',NU_beam,'S',Sec_belt,'IZ',IZ_belt,'IY',IY_belt,'IX',IX_belt,'RHO',RHO_beam);
    mat_belt{2} = ELAS_BEAM('E',E_beam,'NU',NU_beam,'S',Sec_belt,'IZ',IY_belt,'IY',IZ_belt,'IX',IX_belt,'RHO',RHO_beam);
    mat_beam = setnumber(mat_beam,2);
    mat_belt{1} = setnumber(mat_belt{1},3);
    mat_belt{2} = setnumber(mat_belt{2},4);
    S_beam = cellfun(@(S) setmaterial(S,mat_beam),S_beam,'UniformOutput',false);
    S_belt([1,3]) = cellfun(@(S) setmaterial(S,mat_belt{1}),S_belt([1,3]),'UniformOutput',false);
    S_belt([2,4]) = cellfun(@(S) setmaterial(S,mat_belt{2}),S_belt([2,4]),'UniformOutput',false);
    
    S = union(S_plate,S_beam{:});
    if belt
        S = union(S,S_belt{:});
    end
    
    %% Neumann boundary conditions
    p_plate = RHO*g*h; % surface load (body load for plates)
    p_beam = RHO_beam*g*Sec_beam; % line load (body load for beams)
    p_belt = RHO_beam*g*Sec_belt; % line load (body load for beams)
    switch lower(test)
        case {'stability1','stability2','stability3','stability4'}
            p = 400; % pointwise load, pmin = 668N
        case {'statichori1','statichori2','statichori3','statichori4'}
            masse = 50.5;
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse; % surface load (body load for plates)
            p = 400; % pointwise load
            slope = 0;
        case 'staticvert'
            p = 1200; % pointwise load
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            masse = 50.5;
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse; % surface load (body load for plates)
            p = 300; % pointwise load
        case 'impact'
            H = 180e-3;
        case 'drop'
            H = 100e-3;
    end
    
    %% Dirichlet boundary conditions
    x_support = cellfun(@(L) getvertex(L,2)',L_beam,'UniformOutput',false);
    x_support = [x_support{:}]';
    P_support = POINT(x_support);
    
    S = final(S);
    switch lower(test)
        case 'stability1'
            S = addcl(S,P_support([1;4]));
            pmin = (p_plate*(-integral(@(x) 2*sqrt(r^2-x.^2).*(x-a/2),a/2,r,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(a/2-x),0,a/2,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(a/2+x),0,r,'RelTol',eps,'AbsTol',eps))...
                +p_belt*a*(a+b)+2*p_beam*a*l)/(-x_load_stab{1}(1)-a/2);
            if p < pmin % no lifting
                S = addcl(S,P_support([2;3]),'UZ');
            end
        case 'stability2'
            S = addcl(S,P_support([1;2]));
            pmin = (p_plate*(-integral(@(x) 2*sqrt(r^2-x.^2).*(x-b/2),b/2,r,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(b/2-x),0,b/2,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(b/2+x),0,r,'RelTol',eps,'AbsTol',eps))...
                +p_belt*b*(b+a)+2*p_beam*b*l)/(-x_load_stab{2}(2)-b/2);
            if p < pmin % no lifting
                S = addcl(S,P_support([3;4]),'UZ');
            end
        case 'stability3'
            S = addcl(S,P_support([2;3]));
            pmin = (p_plate*(-integral(@(x) 2*sqrt(r^2-x.^2).*(x-a/2),a/2,r,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(a/2-x),0,a/2,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(a/2+x),0,r,'RelTol',eps,'AbsTol',eps))...
                +p_belt*a*(a+b)+2*p_beam*a*l)/(x_load_stab{3}(1)-a/2);
            if p < pmin % no lifting
                S = addcl(S,P_support([1;4]),'UZ');
            end
        case 'stability4'
            S = addcl(S,P_support([3;4]));
            pmin = (p_plate*(-integral(@(x) 2*sqrt(r^2-x.^2).*(x-b/2),b/2,r,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(b/2-x),0,b/2,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(b/2+x),0,r,'RelTol',eps,'AbsTol',eps))...
                +p_belt*b*(b+a)+2*p_beam*b*l)/(x_load_stab{4}(2)-b/2);
            if p < pmin % no lifting
                S = addcl(S,P_support([1;2]),'UZ');
            end
        case {'statichori1','statichori2'}
            S = addcl(S,P_support([3;4]));
            S = addcl(S,P_support([1;2]),'UZ');
        case {'statichori3','statichori4'}
            S = addcl(S,P_support([1;4]));
            S = addcl(S,P_support([2;3]),'UZ');
        case 'staticvert'
            S = addcl(S,P_support,'U');
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            S = addcl(S,P_support);
        case {'impact','drop'}
            S = addcl(S,P_support);
    end
    
    %% Sollicitation vector
    switch lower(test)
        case 'stability1'
            f = nodalload(S,P_load_stab{1},'FZ',-p);
            if isempty(ispointin(P_load_stab{1},POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        case 'stability2'
            f = nodalload(S,P_load_stab{2},'FZ',-p);
            if isempty(ispointin(P_load_stab{2},POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        case 'stability3'
            f = nodalload(S,P_load_stab{3},'FZ',-p);
            if isempty(ispointin(P_load_stab{3},POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        case 'stability4'
            f = nodalload(S,P_load_stab{4},'FZ',-p);
            if isempty(ispointin(P_load_stab{4},POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        case {'statichori1','statichori2','statichori3','statichori4'}
            if strcmpi(test,'statichori1')
                f = nodalload(S,P_load_hori{1},{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
                if isempty(ispointin(P_load_hori{1},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            elseif strcmpi(test,'statichori2')
                f = nodalload(S,P_load_hori{2},{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
                if isempty(ispointin(P_load_hori{2},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            elseif strcmpi(test,'statichori3')
                f = nodalload(S,P_load_hori{3},{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
                if isempty(ispointin(P_load_hori{3},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            elseif strcmpi(test,'statichori4')
                f = nodalload(S,P_load_hori{4},{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
                if isempty(ispointin(P_load_hori{4},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            end
            f = f + bodyload(keepgroupelem(S,3),[],'FZ',-p_masse);
        case 'staticvert'
            f = nodalload(S,P_load_vert,'FZ',-p);
            if isempty(ispointin(P_load_vert,POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            if strcmpi(test,'fatigue1')
                f = nodalload(S,P_load_fati{1},'FY',-p);
                if isempty(ispointin(P_load_fati{1},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            elseif strcmpi(test,'fatigue2')
                f = nodalload(S,P_load_fati{2},'FY',p);
                if isempty(ispointin(P_load_fati{2},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            elseif strcmpi(test,'fatigue3')
                f = nodalload(S,P_load_fati{3},'FX',p);
                if isempty(ispointin(P_load_fati{3},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            elseif strcmpi(test,'fatigue4')
                f = nodalload(S,P_load_fati{4},'FX',-p);
                if isempty(ispointin(P_load_fati{4},POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            end
            f = f + bodyload(keepgroupelem(S,3),[],'FZ',-p_masse);
        case {'impact','drop'}
            error('Not implemented')
    end
    f = f + bodyload(keepgroupelem(S,[1,2,3]),[],'FZ',-p_plate);
    for k=1:length(L_beam)
        f = f + bodyload(keepgroupelem(S,3+k),[],'FZ',-p_beam);
    end
    f = f + bodyload(keepgroupelem(S,3+length(L_beam)+1),[],'FZ',-p_belt);
    f = f + bodyload(keepgroupelem(S,3+length(L_beam)+2),[],'FZ',-p_belt);
    f = f + bodyload(keepgroupelem(S,3+length(L_beam)+3),[],'FZ',-p_belt);
    f = f + bodyload(keepgroupelem(S,3+length(L_beam)+4),[],'FZ',-p_belt);
    
    %% Stiffness matrix and solution
    t = tic;
    u = sparse(getnbddlfree(S),N);
    parfor i=1:N
        % Young modulus
        Ei = E_sample(i);
        % Material
        mati_plate = setparam(mat_plate,'E',Ei);
        Si = setmaterial(S,mati_plate,[1,2,3]);
        % Stiffness matrix
        Ai = calc_rigi(Si);
        % Solution
        u(:,i) = Ai\f;
    end
    time = toc(t);
    
    %% Statistical outputs of solution
    x = getcoord(S.node);
    t = cart2pol(x(:,1),x(:,2),x(:,3));
    funr = @(x,y,theta) dot([cos(theta),sin(theta)],[x,y],2);
    funt = @(x,y,theta) dot([-sin(theta),cos(theta)],[x,y],2);
    funx = @(r,t,theta) dot([cos(theta),-sin(theta)],[r,t],2);
    funy = @(r,t,theta) dot([sin(theta),cos(theta)],[r,t],2);
    
    mean_u = mean(u,2);
    mean_u = unfreevector(S,mean_u);
    
    std_u = std(u,0,2);
    std_u = unfreevector(S,std_u);
    
    probs = [0.025 0.975];
    ci_u = quantile(u,probs,2);
    ci_u = unfreevector(S,ci_u);
    
    % mean_U = mean_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    % mean_Ux = mean_u(findddl(S,'UX'),:);
    % mean_Uy = mean_u(findddl(S,'UY'),:);
    % mean_Uz = mean_u(findddl(S,'UZ'),:);
    % mean_Ur = funr(mean_Ux,mean_Uy,t);
    % mean_Ut = funt(mean_Ux,mean_Uy,t);
    
    % mean_R = mean_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    % mean_Rx = mean_u(findddl(S,'RX'),:);
    % mean_Ry = mean_u(findddl(S,'RY'),:);
    % mean_Rz = mean_u(findddl(S,'RZ'),:);
    % mean_Rr = funr(mean_Rx,mean_Ry,t);
    % mean_Rt = funt(mean_Rx,mean_Ry,t);
    
    % std_U = std_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    % std_Ux = std_u(findddl(S,'UX'),:);
    % std_Uy = std_u(findddl(S,'UY'),:);
    % std_Uz = std_u(findddl(S,'UZ'),:);
    % std_Ur = funr(std_Ux,std_Uy,t);
    % std_Ut = funt(std_Ux,std_Uy,t);
    
    % std_R = std_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    % std_Rx = std_u(findddl(S,'RX'),:);
    % std_Ry = std_u(findddl(S,'RY'),:);
    % std_Rz = std_u(findddl(S,'RZ'),:);
    % std_Rr = funr(std_Rx,std_Ry,t);
    % std_Rt = funt(std_Rx,std_Ry,t);
    
    % ci_U = ci_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    % ci_Ux = ci_u(findddl(S,'UX'),:);
    % ci_Uy = ci_u(findddl(S,'UY'),:);
    % ci_Uz = ci_u(findddl(S,'UZ'),:);
    % ci_Ur = zeros(size(ci_Ux));
    % ci_Ut = zeros(size(ci_Ux));
    % for i=1:length(probs)
    %     ci_Ur(:,i) = funr(ci_Ux(:,i),ci_Uy(:,i),t);
    %     ci_Ut(:,i) = funt(ci_Ux(:,i),ci_Uy(:,i),t);
    % end
    
    % ci_R = ci_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    % ci_Rx = ci_u(findddl(S,'RX'),:);
    % ci_Ry = ci_u(findddl(S,'RY'),:);
    % ci_Rz = ci_u(findddl(S,'RZ'),:);
    % ci_Rr = zeros(size(ci_Rx));
    % ci_Rt = zeros(size(ci_Rx));
    % for i=1:length(probs)
    %     ci_Rr(:,i) = funr(ci_Rx(:,i),ci_Ry(:,i),t);
    %     ci_Rt(:,i) = funt(ci_Rx(:,i),ci_Ry(:,i),t);
    % end
    
    %% Test solution
    P = getcenter(C);
    xP = double(getcoord(P));
    tP = cart2pol(xP(:,1),xP(:,2),xP(:,3));
    
    mean_ux = eval_sol(S,mean_u,P,'UX');
    mean_uy = eval_sol(S,mean_u,P,'UY');
    mean_uz = eval_sol(S,mean_u,P,'UZ');
    mean_ur = funr(mean_ux,mean_uy,tP);
    mean_ut = funt(mean_ux,mean_uy,tP);
    
    mean_rx = eval_sol(S,mean_u,P,'RX');
    mean_ry = eval_sol(S,mean_u,P,'RY');
    mean_rz = eval_sol(S,mean_u,P,'RZ');
    mean_rr = funr(mean_rx,mean_ry,tP);
    mean_rt = funt(mean_rx,mean_ry,tP);
    
    std_ux = eval_sol(S,std_u,P,'UX');
    std_uy = eval_sol(S,std_u,P,'UY');
    std_uz = eval_sol(S,std_u,P,'UZ');
    std_ur = funr(std_ux,std_uy,tP);
    std_ut = funt(std_ux,std_uy,tP);
    
    std_rx = eval_sol(S,std_u,P,'RX');
    std_ry = eval_sol(S,std_u,P,'RY');
    std_rz = eval_sol(S,std_u,P,'RZ');
    std_rr = funr(std_rx,std_ry,tP);
    std_rt = funt(std_rx,std_ry,tP);
    
    ci_ux = eval_sol(S,ci_u,P,'UX');
    ci_uy = eval_sol(S,ci_u,P,'UY');
    ci_uz = eval_sol(S,ci_u,P,'UZ');
    ci_ur = zeros(size(ci_ux));
    ci_ut = zeros(size(ci_ux));
    for j=1:length(probs)
        ci_ur(:,j) = funr(ci_ux(:,j),ci_uy(:,j),tP);
        ci_ut(:,j) = funt(ci_ux(:,j),ci_uy(:,j),tP);
    end
    
    ci_rx = eval_sol(S,ci_u,P,'RX');
    ci_ry = eval_sol(S,ci_u,P,'RY');
    ci_rz = eval_sol(S,ci_u,P,'RZ');
    ci_rr = zeros(size(ci_rx));
    ci_rt = zeros(size(ci_rx));
    for j=1:length(probs)
        ci_rr(:,j) = funr(ci_rx(:,j),ci_ry(:,j),tP);
        ci_rt(:,j) = funt(ci_rx(:,j),ci_ry(:,j),tP);
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','S_plate','elemtype','C','L_beam','L_belt','r','h','f','N');
    save(fullfile(pathname,'solution.mat'),'u','mean_u','std_u','ci_u','probs','time');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'mean_ux','mean_uy','mean_uz','mean_ur','mean_ut',...
        'mean_rx','mean_ry','mean_rz','mean_rr','mean_rt',...
        'std_ux','std_uy','std_uz','std_ur','std_ut',...
        'std_rx','std_ry','std_rz','std_rr','std_rt',...
        'ci_ux','ci_uy','ci_uz','ci_ur','ci_ut',...
        'ci_rx','ci_ry','ci_rz','ci_rr','ci_rt');
else
    load(fullfile(pathname,'problem.mat'),'S','S_plate','elemtype','C','L_beam','L_belt','r','h','f','N');
    load(fullfile(pathname,'solution.mat'),'u','mean_u','std_u','ci_u','probs','time');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'mean_ux','mean_uy','mean_uz','mean_ur','mean_ut',...
        'mean_rx','mean_ry','mean_rz','mean_rr','mean_rt',...
        'std_ux','std_uy','std_uz','std_ur','std_ut',...
        'std_rx','std_ry','std_rz','std_rr','std_rt',...
        'ci','ci_uy','ci_uz','ci_ur','ci_ut',...
        'ci_rx','ci_ry','ci_rz','ci_rr','ci_rt');
end

%% Outputs
fprintf('\nCircular table\n');
fprintf(['test : ' test '\n']);
fprintf(['mesh : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(S_plate));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S_plate));
fprintf('span-to-thickness ratio = %g\n',r/h);
fprintf('nb samples = %g\n',N);
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

disp('Displacement u at point'); disp(P);
fprintf('mean(ux) = %g m, std(ux) = %g m, ci(ux) = [%g %g] m\n',mean_ux,std_ux,ci_ux(1),ci_ux(2));
fprintf('mean(uy) = %g m, std(uy) = %g m, ci(uy) = [%g %g] m\n',mean_uy,std_uy,ci_uy(1),ci_uy(2));
fprintf('mean(uz) = %g m, std(uz) = %g m, ci(uz) = [%g %g] m\n',mean_uz,std_uz,ci_uz(1),ci_uz(2));
if strcmpi(test,'staticvert')
    uz_exp = -2.35e-3;
    err_uz = norm(mean_uz-uz_exp)/norm(uz_exp);
    fprintf('uz_exp   = %g m, error    = %.3e\n',uz_exp,err_uz);
end
fprintf('mean(ur) = %g m, std(ur) = %g m, ci(ur) = [%g %g] m\n',mean_ur,std_ur,ci_ur(1),ci_ur(2));
fprintf('mean(ut) = %g m, std(ut) = %g m, ci(ut) = [%g %g] m\n',mean_ut,std_ut,ci_ut(1),ci_ut(2));
fprintf('\n');

disp('Rotation r at point'); disp(P);
fprintf('mean(rx) = %g rad = %g deg, std(rx) = %g rad = %g deg, ci(rx) = [%g %g] rad = [%g %g] deg\n',mean_rx,rad2deg(mean_rx),std_rx,rad2deg(std_rx),ci_rx(1),ci_rx(2),rad2deg(ci_rx(1)),rad2deg(ci_rx(2)));
fprintf('mean(ry) = %g rad = %g deg, std(ry) = %g rad = %g deg, ci(ry) = [%g %g] rad = [%g %g] deg\n',mean_ry,rad2deg(mean_ry),std_ry,rad2deg(std_ry),ci_ry(1),ci_ry(2),rad2deg(ci_ry(1)),rad2deg(ci_ry(2)));
fprintf('mean(rz) = %g rad = %g deg, std(rz) = %g rad = %g deg, ci(rz) = [%g %g] rad = [%g %g] deg\n',mean_rz,rad2deg(mean_rz),std_rz,rad2deg(std_rz),ci_rz(1),ci_rz(2),rad2deg(ci_rz(1)),rad2deg(ci_rz(2)));
fprintf('mean(rr) = %g rad = %g deg, std(rr) = %g rad = %g deg, ci(rr) = [%g %g] rad = [%g %g] deg\n',mean_rr,rad2deg(mean_rr),std_rr,rad2deg(std_rr),ci_rr(1),ci_rr(2),rad2deg(ci_rr(1)),rad2deg(ci_rr(2)));
fprintf('mean(rt) = %g rad = %g deg, std(rt) = %g rad = %g deg, ci(rt) = [%g %g] rad = [%g %g] deg\n',mean_rt,rad2deg(mean_rt),std_rt,rad2deg(std_rt),ci_rt(1),ci_rt(2),rad2deg(ci_rt(1)),rad2deg(ci_rt(2)));
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    if ~belt
        plotDomain(C,L_beam,'legend',false);
    else
        plotDomain(C,[L_beam,L_belt],'legend',false);
    end
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 5;
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
    hP = plot(P,'g+');
    % legend([hD,hN,hP],[legD,legN,'measure'],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    mean_U = mean_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    ampl = getsize(S)/max(abs(mean_U))/10;
    plotModelDeflection(S,mean_u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*mean_u,'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    %% Display solution
    % ampl = 0;
    ampl = getsize(S)/max(abs(mean_U))/10;
    options = {'solid',true};
    % options = {};
    
    switch lower(test)
        case {'stability1','stability2','stability3','stability4',...
                'staticvert','impact','drop'}
            plotSolution(S,mean_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uz',formats,renderer);
            
            plotSolution(S,std_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uz',formats,renderer);
        case {'statichori1','statichori2','fatigue1','fatigue2'}
            plotSolution(S,mean_u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Ux',formats,renderer);
            
            plotSolution(S,std_u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Ux',formats,renderer);
        case {'statichori3','statichori4','fatigue3','fatigue4'}
            plotSolution(S,mean_u,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uy',formats,renderer);
            
            plotSolution(S,std_u,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uy',formats,renderer);
    end
    
    % plotSolution(S,mean_u,'rotation',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'mean_Rx',formats,renderer);
    %
    % plotSolution(S,std_u,'rotation',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'std_Rx',formats,renderer);
    %
    % plotSolution(S,mean_u,'rotation',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'mean_Ry',formats,renderer);
    %
    % plotSolution(S,std_u,'rotation',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'std_Ry',formats,renderer);
end

%% Display convergence Monte-Carlo
if displayCv
    N = size(u,2);
    switch lower(test)
        case {'stability1','stability2','stability3','stability4',...
                'staticvert','impact','drop'}
            means_u = arrayfun(@(x) eval_sol(S,mean(u(:,1:x),2),P,'UZ'),1:N);
            stds_u = arrayfun(@(x) eval_sol(S,std(u(:,1:x),0,2),P,'UZ'),1:N);
            lowercis_u = arrayfun(@(x) eval_sol(S,quantile(u(:,1:x),probs(1),2),P,'UZ'),1:N);
            uppercis_u = arrayfun(@(x) eval_sol(S,quantile(u(:,1:x),probs(2),2),P,'UZ'),1:N);
            
        case {'statichori1','statichori2','fatigue1','fatigue2'}
            means_u = arrayfun(@(x) eval_sol(S,mean(u(:,1:x),2),P,'UY'),1:N);
            stds_u = arrayfun(@(x) eval_sol(S,std(u(:,1:x),0,2),P,'UY'),1:N);
            lowercis_u = arrayfun(@(x) eval_sol(S,quantile(u(:,1:x),probs(1),2),P,'UY'),1:N);
            uppercis_u = arrayfun(@(x) eval_sol(S,quantile(u(:,1:x),probs(2),2),P,'UY'),1:N);
        case {'statichori3','statichori4','fatigue3','fatigue4'}
            means_u = arrayfun(@(x) eval_sol(S,mean(u(:,1:x),2),P,'UX'),1:N);
            stds_u = arrayfun(@(x) eval_sol(S,std(u(:,1:x),0,2),P,'UX'),1:N);
            lowercis_u = arrayfun(@(x) eval_sol(S,quantile(u(:,1:x),probs(1),2),P,'UX'),1:N);
            uppercis_u = arrayfun(@(x) eval_sol(S,quantile(u(:,1:x),probs(2),2),P,'UX'),1:N);
    end
    
    figure('Name','Convergence solution')
    clf
    ciplot(lowercis_u,uppercis_u,1:N,'b');
    hold on
    alpha(0.2)
    plot(1:N,means_u,'-b','LineWidth',1)
    if strcmpi(test,'staticvert')
        plot(1:N,repmat(uz_exp,1,N),'-r','LineWidth',1)
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',16)
    % xlabel('Nombre de r\''ealisations','Interpreter','latex')
    xlabel('Number of samples','Interpreter','latex')
    ylabel('Solution','Interpreter','latex')
    if strcmpi(test,'staticvert')
        legend({[num2str((probs(2)-probs(1))*100) '% confidence interval'],'mean value','experimental value'})
    else
        legend({[num2str((probs(2)-probs(1))*100) '% confidence interval'],'mean value'})
    end
    mysaveas(pathname,'convergence_solution','fig');
    mymatlab2tikz(pathname,'convergence_solution.tex');
    
    figure('Name','Convergence mean')
    clf
    plot(1:N,means_u,'-b','LineWidth',1)
    grid on
    box on
    set(gca,'FontSize',16)
    % xlabel('Nombre de r\''ealisations','Interpreter','latex')
    % ylabel('Moyenne','Interpreter','latex')
    xlabel('Number of samples','Interpreter','latex')
    ylabel('Mean','Interpreter','latex')
    mysaveas(pathname,'convergence_mean','fig');
    mymatlab2tikz(pathname,'convergence_mean.tex');
    
    figure('Name','Convergence standard deviation')
    clf
    plot(1:N,stds_u,'-r','LineWidth',1)
    grid on
    box on
    set(gca,'FontSize',16)
    % xlabel('Nombre de r\''ealisations','Interpreter','latex')
    % ylabel('Ecart-type','Interpreter','latex')
    xlabel('Number of samples','Interpreter','latex')
    ylabel('Standard deviation','Interpreter','latex')
    mysaveas(pathname,'convergence_std','fig');
    mymatlab2tikz(pathname,'convergence_std.tex');
end

end

myparallel('stop');
