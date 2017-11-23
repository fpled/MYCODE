%% Table circular deterministic linear elasticity %%
%%------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');

%% Input data
solveProblem = true;
displaySolution = true;

% loadings = {'Uniform'};
% loadings = {'Concentrated'};
loadings = {'Uniform','Concentrated'};
% elemtypes = {'DKT'};
% elemtypes = {'DKQ'};
% elemtypes = {'DST'};
% elemtypes = {'DSQ'};
% elemtypes = {'COQ4'};
% elemtypes = {'DKT','DKQ'}; % Kirchhoff-Love (classical) plate theory
% elemtypes = {'DST','DSQ','COQ4'}; % Reissner-Mindlin (first-order shear) plate theory
elemtypes = {'DKT','DKQ','DST','DSQ','COQ4'}; % Both plate theories

formats = {'fig','epsc2'};
renderer = 'OpenGL';

for il=1:length(loadings)
    loading = loadings{il};
    filename = ['tableCircDetLinElas' loading];
    close all
    
for ie=1:length(elemtypes)
    elemtype = elemtypes{ie};
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','plate',filename,...
        elemtype);
    if ~exist(pathname,'dir')
        mkdir(pathname);
    end

%% Problem
if solveProblem
    %% Domains and meshes
    % Plate
    r = 1;
    C = CIRCLE(0.0,0.0,0.0,r);
    
    % Points
    P_load = POINT([-r/2,0.0,0.0]);
    P_beam = getcenter(C);
    x_load = double(getcoord(P_load));
    x_beam = double(getcoord(P_beam));
    
    % Beam
    l = 1;
    L_beam = LIGNE(P_beam,P_beam+POINT([0.0,0.0,-l]));
    
    % Plate mesh
    cl_plate = r/10;
    switch lower(loading)
        case 'uniform'
            points = x_beam;
        case 'concentrated'
            points = {x_beam,x_load};
    end
    S_plate = build_model(C,'cl',cl_plate,'elemtype',elemtype,'filename',fullfile(pathname,['gmsh_plate_circ_' elemtype '_cl_' num2str(cl_plate)]),'points',points);
    
    % Beam mesh
    nbelem_beam = 10;
    S_beam = build_model(L_beam,'nbelem',nbelem_beam,'elemtype','BEAM');
    % cl_beam = 0.1;
    % S_beam = build_model(L_beam,'cl',cl_beam,'elemtype','BEAM','filename',fullfile(pathname,['gmsh_beam_cl_' num2str(cl_beam)]));
    
    %% Materials
    % Gravitational acceleration
    g = 9.81;
    
    % Plate
    % Young modulus
    E = 1;
    % Poisson ratio
    NU = 0.3;
    % Density
    RHO = 1;
    % Thickness
    h = 0.1;
    % Material
    mat_plate = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',h,'k',5/6);
    mat_plate = setnumber(mat_plate,1);
    S_plate = setmaterial(S_plate,mat_plate);
    
    % Beam
    % Young modulus
    E_beam = 1;
    % Poisson ratio
    NU_beam = 0.3;
    % Density
    RHO_beam = 1;
    % Radius
    r_beam = 0.1;
    % Cross-section area
    Sec_beam = pi*r_beam^2;
    % Planar second moment of area (or Planar area moment of inertia)
    IY = pi*r_beam^4/4;
    IZ = IY;
    % Polar second moment of area (or Polar area moment of inertia)
    IX = IY+IZ;
    % Material
    mat_beam = ELAS_BEAM('E',E_beam,'NU',NU_beam,'S',Sec_beam,'IZ',IZ,'IY',IY,'IX',IX,'RHO',RHO_beam);
    mat_beam = setnumber(mat_beam,2);
    S_beam = setmaterial(S_beam,mat_beam);
    
    S = union(S_plate,S_beam);
    
    %% Dirichlet boundary conditions
    x_support = getvertex(L_beam,2);
    P_support = POINT(x_support);
    
    S = final(S);
    S = addcl(S,P_support);
    
    %% Stiffness matrices and sollicitation vectors
    % Uniform or Concentrated load
    switch lower(loading)
        case 'uniform'
            p = RHO*g*h;
        case 'concentrated'
            Sec = pi*r^2;
            p = RHO*g*h*Sec;
    end
    
    A = calc_rigi(S);
    switch lower(loading)
        case 'uniform'
            f = bodyload(keepgroupelem(S,1),[],'FZ',-p);
        case 'concentrated'
            f = nodalload(S,P_load,'FZ',-p);
            if isempty(ispointin(P_load,POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
    end
    
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
    save(fullfile(pathname,'problem.mat'),'S','S_plate','C','L_beam','r','h','f');
    save(fullfile(pathname,'solution.mat'),'u','time',...
        'U','Ux','Uy','Uz','Ur','Ut',...
        'R','Rx','Ry','Rz','Rr','Rt');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz','ur','ut',...
        'rx','ry','rz','rr','rt');
else
    load(fullfile(pathname,'problem.mat'),'S','S_plate','C','L_beam','r','h','f');
    load(fullfile(pathname,'solution.mat'),'u','time',...
        'U','Ux','Uy','Uz','Ur','Ut',...
        'R','Rx','Ry','Rz','Rr','Rt');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz','ur','ut',...
        'rx','ry','rz','rr','rt');
end

%% Outputs
fprintf('\nCircular table\n');
fprintf(['load : ' loading '\n']);
fprintf(['mesh : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(S_plate));
fprintf('nb nodes    = %g\n',getnbnode(S_plate));
fprintf('nb dofs     = %g\n',getnbddl(S_plate));
fprintf('span-to-thickness ratio = %g\n',r/h);
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

disp('Displacement u at point'); disp(P);
fprintf('ux    = %g\n',ux);
fprintf('uy    = %g\n',uy);
fprintf('uz    = %g\n',uz);
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
    plotDomain(C,L_beam,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    switch lower(loading)
        case 'uniform'
            ampl = 2;
        case 'concentrated'
            ampl = 0.2;
    end
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
    % legend([hD,hN],[legD,legN],'Location','northeastoutside')
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
    % ampl = 0;
    ampl = getsize(S)/max(abs(u))/10;
    options = {'solid',true};
    % options = {};
    
    plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
    mysaveas(pathname,'Uz',formats,renderer);
    
    % plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'Rx',formats,renderer);
    %
    % plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'Ry',formats,renderer);
end

end
end
