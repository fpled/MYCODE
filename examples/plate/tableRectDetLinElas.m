%% Table rectangular deterministic linear elasticity %%
%%---------------------------------------------------%%

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
% meshtypes = {'Structured'};
% meshtypes = {'Unstructured'};
meshtypes = {'Structured','Unstructured'};

formats = {'fig','epsc2'};
renderer = 'OpenGL';

for il=1:length(loadings)
    loading = loadings{il};
    filename = ['tableRectDetLinElas' loading];
    close all
    
for ie=1:length(elemtypes)
    elemtype = elemtypes{ie};
    
for im=1:length(meshtypes)
    meshtype = meshtypes{im};
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','plate',filename,[elemtype meshtype]);
    if ~exist(pathname,'dir')
        mkdir(pathname);
    end

%% Problem
if solveProblem
    %% Domains and meshes
    % Plate
    a = 1;
    b = 1;
    Q = QUADRANGLE([0.0,0.0,0.0],[a,0.0,0.0],[a,b,0.0],[0.0,b,0.0]);
    
    % Beams
    l = 1;
    L_beam{1} = LIGNE([1/5*a,1/5*b,0.0],[1/5*a,1/5*b,-l]);
    L_beam{2} = LIGNE([4/5*a,1/5*b,0.0],[4/5*a,1/5*b,-l]);
    L_beam{3} = LIGNE([4/5*a,4/5*b,0.0],[4/5*a,4/5*b,-l]);
    L_beam{4} = LIGNE([1/5*a,4/5*b,0.0],[1/5*a,4/5*b,-l]);
    
    % Points
    x_beam = cellfun(@(L) getvertex(L,1),L_beam,'UniformOutput',false);
    P_beam = cellfun(@(x) POINT(x),x_beam,'UniformOutput',false);
    P_load = getcenter(Q);
    x_load = double(getcoord(P_load));
    
    % Plate mesh
    switch lower(meshtype)
        case 'structured'
            nbelem_plate = [20,20];
            S_plate = build_model(Q,'nbelem',nbelem_plate,'elemtype',elemtype);
        case 'unstructured'
            cl_plate = min(a,b)/20;
            switch lower(loading)
                case 'uniform'
                    points = x_beam;
                case 'concentrated'
                    points = [x_beam,{x_load}];
            end
            S_plate = build_model(Q,'cl',cl_plate,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_rect_' elemtype '_cl_' num2str(cl_plate)],'points',points);
    end
    
    % Beams meshes
    nbelem_beam = 10;
    S_beam = cellfun(@(L) build_model(L,'nbelem',nbelem_beam,'elemtype','BEAM'),L_beam,'UniformOutput',false);
    % cl_beam = 0.1;
    % S_beam = cellfun(@(L,n) build_model(L,'cl',cl_beam,'elemtype','BEAM','filename',[pathname 'gmsh_beam_' num2str(n) '_cl_' num2str(cl_beam)]),L_beam,num2cell(1:length(L_beam)),'UniformOutput',false);
    
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
    
    % Beams
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
    S_beam = cellfun(@(S) setmaterial(S,mat_beam),S_beam,'UniformOutput',false);
    
    S = union(S_plate,S_beam{:});
    
    %% Dirichlet boundary conditions
    x_support = cellfun(@(L) getvertex(L,2)',L_beam,'UniformOutput',false);
    x_support = [x_support{:}]';
    P_support = POINT(x_support);
    
    S = final(S);
    S = addcl(S,P_support); % S = addcl(S,P_support,{'U','R'},0);
    
    %% Stiffness matrices and sollicitation vectors
    % Uniform or Concentrated load
    switch lower(loading)
        case 'uniform'
            p = RHO*g*h;
        case 'concentrated'
            Sec = a*b;
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
    P = getcenter(Q);
    
    ux = eval_sol(S,u,P,'UX');
    uy = eval_sol(S,u,P,'UY');
    uz = eval_sol(S,u,P,'UZ');
    
    rx = eval_sol(S,u,P,'RX');
    ry = eval_sol(S,u,P,'RY');
    rz = eval_sol(S,u,P,'RZ');
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','S_plate','Q','L_beam','a','b','h','f');
    save(fullfile(pathname,'solution.mat'),'u','time',...
        'U','Ux','Uy','Uz',...
        'R','Rx','Ry','Rz');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz',...
        'rx','ry','rz');
else
    load(fullfile(pathname,'problem.mat'),'S','S_plate','Q','L_beam','a','b','h','f');
    load(fullfile(pathname,'solution.mat'),'u','time',...
        'U','Ux','Uy','Uz',...
        'R','Rx','Ry','Rz');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'ux','uy','uz',...
        'rx','ry','rz');
end

%% Outputs
fprintf('\nRectangular table\n');
fprintf(['load : ' loading '\n']);
fprintf(['mesh : ' elemtype ' ' meshtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(S_plate));
fprintf('nb nodes    = %g\n',getnbnode(S_plate));
fprintf('nb dofs     = %g\n',getnbddl(S_plate));
fprintf('span-to-thickness ratio = %g\n',max(a,b)/h);
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

disp('Displacement u at point'); disp(P);
fprintf('ux    = %g\n',ux);
fprintf('uy    = %g\n',uy);
fprintf('uz    = %g\n',uz);
fprintf('\n');

disp('Rotation r at point'); disp(P);
fprintf('rx    = %g\n',rx);
fprintf('ry    = %g\n',ry);
fprintf('rz    = %g\n',rz);
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain(Q,L_beam,'Color','w','legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    switch lower(loading)
        case 'uniform'
            ampl = 2;
        case 'concentrated'
            ampl = 0.5;
    end
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
    % legend([hD,hN],[legD,legN])
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'node',true,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    ampl = getsize(S)/max(abs(u))/5;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'node',true,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'node',true);
    plot(S+ampl*u,'Color','b','FaceColor','b','FaceAlpha',0.1,'node',true);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    %% Display solution
    % ampl = 0;
    ampl = getsize(S)/max(abs(u))/5;
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
end
