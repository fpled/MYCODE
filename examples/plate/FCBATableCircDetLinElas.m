%% FCBA table circular deterministic linear elasticity %%
%%-----------------------------------------------------%%

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

test = 'Stability1'; % stability test under vertical load 1
% test = 'Stability2'; % stability test under vertical load 2
% test = 'Stability3'; % stability test under vertical load 3
% test = 'Stability4'; % stability test under vertical load 4
% test = 'StaticHori1'; % test under static horizontal load 1
% test = 'StaticHori2'; % test under static horizontal load 2
% test = 'StaticHori3'; % test under static horizontal load 3
% test = 'StaticHori4'; % test under static horizontal load 4
% test = 'StaticVert'; % test under static vertical load
% test = 'Fatigue1'; % fatigue test under horizontal load 1
% test = 'Fatigue2'; % fatigue test under horizontal load 2
% test = 'Fatigue3'; % fatigue test under horizontal load 3
% test = 'Fatigue4'; % fatigue test under horizontal load 4
% test = 'Impact'; % vertical impact test
% test = 'Drop'; % drop test

belt = 1; % belt modelisation

formats = {'fig','epsc2'};
renderer = 'OpenGL';

filename = ['FCBATableCircDetLinElas' test];
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
    
    %% Materials
    % Gravitational acceleration
    g = 9.81;
    
    % Plate
    % Young modulus
    E = 2.9914e9;
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
    
    % Beams
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
    
    %% Stiffness matrices and sollicitation vectors
    A = calc_rigi(S);
    
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
    u = A\f;
    time = toc(t);
    
    x = getcoord(S.node);
    t = cart2pol(x(:,1),x(:,2),x(:,3));
    funr = @(x,y,theta) dot([cos(theta),sin(theta)],[x,y],2);
    funt = @(x,y,theta) dot([-sin(theta),cos(theta)],[x,y],2);
    funx = @(r,t,theta) dot([cos(theta),-sin(theta)],[r,t],2);
    funy = @(r,t,theta) dot([sin(theta),cos(theta)],[r,t],2);
    
    u = unfreevector(S,u);
    
    U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    Ux = u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,u,S.node,'UX')));
    Uy = u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,u,S.node,'UY')));
    Uz = u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,u,S.node,'UZ')));
    Ur = funr(Ux,Uy,t);
    Ut = funt(Ux,Uy,t);
    
    R = u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    Rx = u(findddl(S,'RX'),:); % Rx = double(squeeze(eval_sol(S,u,S.node,'RX')));
    Ry = u(findddl(S,'RY'),:); % Ry = double(squeeze(eval_sol(S,u,S.node,'RY')));
    Rz = u(findddl(S,'RZ'),:); % Rz = double(squeeze(eval_sol(S,u,S.node,'RZ')));
    Rr = funr(Rx,Ry,t);
    Rt = funt(Rx,Ry,t);
    
    %% Test solution
    P = getcenter(C);
    xP = double(getcoord(P));
    tP = cart2pol(xP(:,1),xP(:,2),xP(:,3));
    
    ux = eval_sol(S,u,P,'UX');
    uy = eval_sol(S,u,P,'UY');
    uz = eval_sol(S,u,P,'UZ');
    ur = funr(ux,uy,tP);
    ut = funt(ux,uy,tP);
    
    rx = eval_sol(S,u,P,'RX');
    ry = eval_sol(S,u,P,'RY');
    rz = eval_sol(S,u,P,'RZ');
    rr = funr(rx,ry,tP);
    rt = funt(rx,ry,tP);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','S_plate','elemtype','C','L_beam','L_belt','r','h','f');
    save(fullfile(pathname,'solution.mat'),'u','time',...
        'U','Ux','Uy','Uz','Ur','Ut',...
        'R','Rx','Ry','Rz','Rr','Rt');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz','ur','ut',...
        'rx','ry','rz','rr','rt');
else
    load(fullfile(pathname,'problem.mat'),'S','S_plate','elemtype','C','L_beam','L_belt','r','h','f');
    load(fullfile(pathname,'solution.mat'),'u','time',...
        'U','Ux','Uy','Uz','Ur','Ut',...
        'R','Rx','Ry','Rz','Rr','Rt');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz','ur','ut',...
        'rx','ry','rz','rr','rt');
end

%% Outputs
fprintf('\nCircular table\n');
fprintf(['test : ' test '\n']);
fprintf(['mesh : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(S_plate));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S_plate));
fprintf('span-to-thickness ratio = %g\n',r/h);
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

disp('Displacement u at point'); disp(P);
fprintf('ux    = %g\n',ux);
fprintf('uy    = %g\n',uy);
fprintf('uz    = %g\n',uz);
if strcmpi(test,'staticvert')
    uz_exp = -2.35e-3;
    err_uz = norm(uz-uz_exp)/norm(uz_exp);
    fprintf('uz_exp= %g, error = %.3e\n',uz_exp,err_uz);
end
fprintf('ur    = %g\n',ur);
fprintf('ut    = %g\n',ut);
fprintf('\n');

disp('Rotation r at point'); disp(P);
fprintf('rx    = %g\n',rx);
fprintf('ry    = %g\n',ry);
fprintf('rz    = %g\n',rz);
fprintf('rr    = %g\n',rr);
fprintf('rt    = %g\n',rt);
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
    % legend([hD,hN],[legD,legN],'Location','northeastoutside')
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    ampl = getsize(S)/max(abs(u))/10;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*u,'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    %% Display solution
    % ampl = 0;
    ampl = getsize(S)/max(abs(u))/10;
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
            plotSolution(S,u,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'Uy',formats,renderer);
        case {'fatigue3','fatigue4'}
            plotSolution(S,u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'Ux',formats,renderer);
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
