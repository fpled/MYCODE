%% FCBA table circular deterministic linear elasticity %%
%%-----------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');
myparallel('start');

%% Input data
solveProblem = true;
displaySolution = true;
displayCv = true;

% test = 'Stability1'; % stability test under vertical load 1
% test = 'Stability2'; % stability test under vertical load 2
% test = 'Stability3'; % stability test under vertical load 3
% test = 'Stability4'; % stability test under vertical load 4
% test = 'StaticHori1'; % test under static horizontal load 1
% test = 'StaticHori2'; % test under static horizontal load 2
% test = 'StaticHori3'; % test under static horizontal load 3
% test = 'StaticHori4'; % test under static horizontal load 4
test = 'StaticVert'; % test under static vertical load
% test = 'Fatigue1'; % fatigue test under horizontal load 1
% test = 'Fatigue2'; % fatigue test under horizontal load 2
% test = 'Fatigue3'; % fatigue test under horizontal load 3
% test = 'Fatigue4'; % fatigue test under horizontal load 4
% test = 'Impact'; % vertical impact test
% test = 'Drop'; % drop test

belt = 1; % belt modelisation

formats = {'fig','epsc2'};
renderer = 'OpenGL';

filename = ['FCBATableCircStoLinElas' test];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,...
    'results',filesep,'plate',filesep,filename,filesep);
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
        S_plate = gmshFCBAtablecirc(C,Q_belt,C_masse,Pb,Pe,[],Pi,cl_plate,cl_belt,cl_plate,cl_plate,cl_plate,cl_plate,cl_plate,[pathname 'gmsh_plate_circ_' elemtype  '_cl_' num2str(cl_plate)],3);
    else
        S_plate = gmshFCBAtablecirc(C,Q_belt,C_masse,Pb,Pe,[],Pi,cl_plate,cl_belt,cl_plate,cl_plate,cl_plate,cl_plate,cl_plate,[pathname 'gmsh_plate_circ_' elemtype  '_cl_' num2str(cl_plate)],3,'recombine');
    end
    S_plate = convertelem(S_plate,elemtype);
    
    % Beams meshes
    nbelem_beam = 80;
    S_beam = cellfun(@(L) build_model(L,'nbelem',nbelem_beam,'elemtype','BEAM'),L_beam,'UniformOutput',false);
    % cl_beam = l/80;
    % S_beam = cellfun(@(L,n) build_model(L,'cl',cl_beam,'elemtype','BEAM','filename',[pathname 'gmsh_beam_' num2str(n) '_cl_' num2str(cl_beam)]),L_beam,num2cell(1:length(L_beam)),'UniformOutput',false);
    
    % Belt mesh
    S_belt = cellfun(@(L,n) build_model(L,'cl',cl_belt,'elemtype','BEAM','filename',[pathname 'gmsh_belt_' num2str(n) '_cl_' num2str(cl_belt)]),L_belt,num2cell(1:length(L_belt)),'UniformOutput',false);
    
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
    mass_plate = 18.54;
    Sec_plate = pi*r^2;
    RHO = mass_plate/(Sec_plate*h);
    % Extensional stiffness (or Membrane rigidity)
    A_rig = E*h/(1-NU^2);
    % Bending stiffness (or Flexural rigidity)
    D_rig = E*h^3/(12*(1-NU^2));
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
    p_plate = RHO*g*h;
    p_beam = RHO_beam*g*Sec_beam;
    p_belt = RHO_beam*g*Sec_belt;
    switch lower(test)
        case {'stability1','stability2','stability3','stability4'}
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
    x_support = cellfun(@(L) getvertex(L,2)',L_beam,'UniformOutput',false);
    x_support = [x_support{:}]';
    P_support = POINT(x_support);
    
    S = final(S);
    switch lower(test)
        case 'stability1'
            S = addcl(S,P_support(4)); % addcl(S,P_support(4),{'U','R'},0);
            S = addcl(S,P_support(1),'U'); % addcl(S,P_support(1),'U',0);
            pmin = (p_plate*(-integral(@(x) 2*sqrt(r^2-x.^2).*(x-a/2),a/2,r,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(a/2-x),0,a/2,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(a/2+x),0,r,'RelTol',eps,'AbsTol',eps))...
                +p_belt*a*(a+b)+2*p_beam*a*l)/(-x_load_stab{1}(1)-a/2);
            if p < pmin % no lifting
                S = addcl(S,P_support([2 3]),'UZ'); % addcl(S,P_support([2 3]),'UZ',0);
            end
        case 'stability2'
            S = addcl(S,P_support(1)); % addcl(S,P_support(4),{'U','R'},0);
            S = addcl(S,P_support(2),'U'); % addcl(S,P_support(1),'U',0);
            pmin = (p_plate*(-integral(@(x) 2*sqrt(r^2-x.^2).*(x-b/2),b/2,r,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(b/2-x),0,b/2,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(b/2+x),0,r,'RelTol',eps,'AbsTol',eps))...
                +p_belt*b*(b+a)+2*p_beam*b*l)/(-x_load_stab{2}(2)-b/2);
            if p < pmin % no lifting
                S = addcl(S,P_support([3 4]),'UZ'); % addcl(S,P_support([3 4]),'UZ',0);
            end
        case 'stability3'
            S = addcl(S,P_support(2)); % addcl(S,P_support(4),{'U','R'},0);
            S = addcl(S,P_support(3),'U'); % addcl(S,P_support(1),'U',0);
            pmin = (p_plate*(-integral(@(x) 2*sqrt(r^2-x.^2).*(x-a/2),a/2,r,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(a/2-x),0,a/2,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(a/2+x),0,r,'RelTol',eps,'AbsTol',eps))...
                +p_belt*a*(a+b)+2*p_beam*a*l)/(x_load_stab{3}(1)-a/2);
            if p < pmin % no lifting
                S = addcl(S,P_support([1 4]),'UZ'); % addcl(S,P_support([1 4]),'UZ',0);
            end
        case 'stability4'
            S = addcl(S,P_support(3)); % addcl(S,P_support(4),{'U','R'},0);
            S = addcl(S,P_support(4),'U'); % addcl(S,P_support(1),'U',0);
            pmin = (p_plate*(-integral(@(x) 2*sqrt(r^2-x.^2).*(x-b/2),b/2,r,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(b/2-x),0,b/2,'RelTol',eps,'AbsTol',eps)...
                +integral(@(x) 2*sqrt(r^2-x.^2).*(b/2+x),0,r,'RelTol',eps,'AbsTol',eps))...
                +p_belt*b*(b+a)+2*p_beam*b*l)/(x_load_stab{4}(2)-b/2);
            if p < pmin % no lifting
                S = addcl(S,P_support([1 2]),'UZ'); % addcl(S,P_support([1 2]),'UZ',0);
            end
        case {'statichori1','statichori2'}
            S = addcl(S,P_support([3 4])); % addcl(S,P_support([3 4]),{'U','R'},0);
            S = addcl(S,P_support([1 2]),'UZ'); % addcl(S,P_support([1 2]),'UZ',0);
        case {'statichori3','statichori4'}
            S = addcl(S,P_support([4 1])); % addcl(S,P_support([4 1]),{'U','R'},0);
            S = addcl(S,P_support([2 3]),'UZ'); % addcl(S,P_support([2 3]),'UZ',0);
        case 'staticvert'
            S = addcl(S,P_support,'U'); % addcl(S,P_support,'U',0);
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            S = addcl(S,P_support); % addcl(S,P_support,{'U','R'},0);
        case {'impact','drop'}
            S = addcl(S,P_support); % addcl(S,P_support,{'U','R'},0);
    end
    
    %% Stiffness matrices
    A = cell(N,1);
    for i=1:N
        % Young modulus
        Ei = e(i);
        % Material
        mat_platei = setparam(mat_plate,'E',Ei);
        Si = setmaterial(S,mat_platei,[1,2,3]);
        % Stiffness matrix
        A{i} = calc_rigi(Si);
    end
    
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
            f = nodalload(S,P_load_stab{4},'FZ',-p*cosd(slope));
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
    
    %% Solution
    t = tic;
    parfor i=1:N
        u(:,i) = A{i}\f;
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
    
    mean_U = mean_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    mean_Ux = mean_u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,mean_u,S.node,'UX')),:);
    mean_Uy = mean_u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,mean_u,S.node,'UY')),:);
    mean_Uz = mean_u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,mean_u,S.node,'UZ')),:);
    mean_Ur = funr(mean_Ux,mean_Uy,t);
    mean_Ut = funt(mean_Ux,mean_Uy,t);
    
    mean_R = mean_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    mean_Rx = mean_u(findddl(S,'RX'),:); % Rx = double(squeeze(eval_sol(S,mean_u,S.node,'RX')),:);
    mean_Ry = mean_u(findddl(S,'RY'),:); % Ry = double(squeeze(eval_sol(S,mean_u,S.node,'RY')),:);
    mean_Rz = mean_u(findddl(S,'RZ'),:); % Rz = double(squeeze(eval_sol(S,mean_u,S.node,'RZ')),:);
    mean_Rr = funr(mean_Rx,mean_Ry,t);
    mean_Rt = funt(mean_Rx,mean_Ry,t);
    
    std_U = std_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    std_Ux = std_u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,std_u,S.node,'UX')),:);
    std_Uy = std_u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,std_u,S.node,'UY')),:);
    std_Uz = std_u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,std_u,S.node,'UZ')),:);
    std_Ur = funr(std_Ux,std_Uy,t);
    std_Ut = funt(std_Ux,std_Uy,t);
    
    std_R = std_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    std_Rx = std_u(findddl(S,'RX'),:); % Rx = double(squeeze(eval_sol(S,std_u,S.node,'RX')),:);
    std_Ry = std_u(findddl(S,'RY'),:); % Ry = double(squeeze(eval_sol(S,std_u,S.node,'RY')),:);
    std_Rz = std_u(findddl(S,'RZ'),:); % Rz = double(squeeze(eval_sol(S,std_u,S.node,'RZ')),:);
    std_Rr = funr(std_Rx,std_Ry,t);
    std_Rt = funt(std_Rx,std_Ry,t);
    
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
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','S_plate','elemtype','C','L_beam','L_belt','r','h','f','N');
    save(fullfile(pathname,'solution.mat'),'u','mean_u','std_u','time',...
        'mean_U','mean_Ux','mean_Uy','mean_Uz','mean_Ur','mean_Ut',...
        'mean_R','mean_Rx','mean_Ry','mean_Rz','mean_Rr','mean_Rt',...
        'std_U','std_Ux','std_Uy','std_Uz','std_Ur','std_Ut',...
        'std_R','std_Rx','std_Ry','std_Rz','std_Rr','std_Rt');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'mean_ux','mean_uy','mean_uz','mean_ur','mean_ut',...
        'mean_rx','mean_ry','mean_rz','mean_rr','mean_rt',...
        'std_ux','std_uy','std_uz','std_ur','std_ut',...
        'std_rx','std_ry','std_rz','std_rr','std_rt');
else
    load(fullfile(pathname,'problem.mat'),'S','S_plate','elemtype','C','L_beam','L_belt','r','h','f','N');
    load(fullfile(pathname,'solution.mat'),'u','mean_u','std_u','time',...
        'mean_U','mean_Ux','mean_Uy','mean_Uz','mean_Ur','mean_Ut',...
        'mean_R','mean_Rx','mean_Ry','mean_Rz','mean_Rr','mean_Rt',...
        'std_U','std_Ux','std_Uy','std_Uz','std_Ur','std_Ut',...
        'std_R','std_Rx','std_Ry','std_Rz','std_Rr','std_Rt');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'mean_ux','mean_uy','mean_uz','mean_ur','mean_ut',...
        'mean_rx','mean_ry','mean_rz','mean_rr','mean_rt',...
        'std_ux','std_uy','std_uz','std_ur','std_ut',...
        'std_rx','std_ry','std_rz','std_rr','std_rt');
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
fprintf('mean(ux) = %g, std(ux) = %g\n',mean_ux,std_ux);
fprintf('mean(uy) = %g, std(uy) = %g\n',mean_uy,std_uy);
fprintf('mean(uz) = %g, std(uz) = %g\n',mean_uz,std_uz);
if strcmpi(test,'staticvert')
    uz_exp = -2.35e-3;
    err_uz = norm(mean_uz-uz_exp)/norm(uz_exp);
    fprintf('uz_exp   = %g, error    = %.3e\n',uz_exp,err_uz);
end
fprintf('mean(ur) = %g, std(ur) = %g\n',mean_ur,std_ur);
fprintf('mean(ut) = %g, std(ut) = %g\n',mean_ut,std_ut);
fprintf('\n');

disp('Rotation r at point'); disp(P);
fprintf('mean(rx) = %g, std(rx) = %g\n',mean_rx,std_rx);
fprintf('mean(ry) = %g, std(ry) = %g\n',mean_ry,std_ry);
fprintf('mean(rz) = %g, std(rz) = %g\n',mean_rz,std_rz);
fprintf('mean(rr) = %g, std(rr) = %g\n',mean_rr,std_rr);
fprintf('mean(rt) = %g, std(rt) = %g\n',mean_rt,std_rt);
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
    
    plotSolution(S,mean_u,'displ',3,'ampl',ampl,options{:});
    mysaveas(pathname,'mean_Uz',formats,renderer);
    
    plotSolution(S,std_u,'displ',3,'ampl',ampl,options{:});
    mysaveas(pathname,'std_Uz',formats,renderer);
    
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
