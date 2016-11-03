%% Table rectangular deterministic linear elasticity %%
%%---------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% myparallel('start');

%% Input data

% loadings = {'uniform'};
% loadings = {'concentrated'};
loadings = {'uniform','concentrated'};
% elemtypes = {'DKT'};
% elemtypes = {'DKQ'};
% elemtypes = {'COQ4'};
% elemtypes = {'DKT','DKQ'};
elemtypes = {'DKT','DKQ','COQ4'};
% meshtypes = {'structured'};
% meshtypes = {'unstructured'};
meshtypes = {'structured','unstructured'};

formats = {'fig','epsc2'};
renderer = 'OpenGL';

for il=1:length(loadings)
    loading = loadings{il};
    filename = ['table_rect_det_lin_elas_' loading];
    close all
    
for ie=1:length(elemtypes)
    elemtype = elemtypes{ie};
    
for im=1:length(meshtypes)
    meshtype = meshtypes{im};
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep,[elemtype '_' meshtype],filesep);
    if ~exist(pathname,'dir')
        mkdir(pathname);
    end

%% Domains and meshes

% Plate
a = 1;
b = 10;
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
switch meshtype
    case 'structured'
        nbelem_plate = [20,20];
        S_plate = build_model(Q,'nbelem',nbelem_plate,'elemtype',elemtype);
    case 'unstructured'
        cl_plate = min(a,b)/20;
        switch loading
            case 'uniform'
                points = x_beam;
            case 'concentrated'
                points = [x_beam,{x_load}];
        end
        S_plate = build_model(Q,'cl',cl_plate,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_rect_' elemtype  '_cl_' num2str(cl_plate)],'points',points);
end

% Beam meshes
nbelem_beam = 10;
S_beam = cellfun(@(L) build_model(L,'nbelem',nbelem_beam,'elemtype','BEAM'),L_beam,'UniformOutput',false);
% cl_beam = 0.1;
% S_beam = cellfun(@(L,n) build_model(L,'cl',cl_plate,'elemtype','BEAM','filename',[pathname 'gmsh_beam_' num2str(n) '_cl_' num2str(cl_beam)]),L_beam,num2cell(1:length(L_beam)),'UniformOutput',false);

%% Materials

% Gravitational acceleration
g = 10;

% Plate
% Young modulus
E = 1;
% Poisson ratio
NU = 0.3;
% Density
RHO = 1;
% Thickness
h = 0.1;
% Extensional stiffness (or Membrane rigidity)
A = E*h/(1-NU^2);
% Bending stiffness (or Flexural rigidity)
D = E*h^3/(12*(1-NU^2));
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
S_beam = cellfun(@(S) setmaterial(S,mat_beam),S_beam,'UniformOutput',false);

problem.S = union(S_plate,S_beam{:});

%% Dirichlet boundary conditions

x_support = cellfun(@(L) getvertex(L,2)',L_beam,'UniformOutput',false);
x_support = [x_support{:}]';
P_support = POINT(x_support);

problem.S = final(problem.S);
problem.S = addcl(problem.S,P_support); % addcl(problem.S,P_support,{'U','R'},0);

%% Stiffness matrices and sollicitation vectors

% Uniform or Concentrated load
switch loading
    case 'uniform'
        p = RHO*g*h;
    case 'concentrated'
        p = RHO*g*h*a*b;
end

problem.A = calc_rigi(problem.S);
switch loading
    case 'uniform'
        problem.b = bodyload(keepgroupelem(problem.S,1),[],'FZ',-p);
    case 'concentrated'
        problem.b = nodalload(problem.S,P_load,'FZ',-p);
        if isempty(ispointin(P_load,POINT(problem.S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
end

%% Resolution

t = tic;
u = solveSystem(problem);
time = toc(t);

%% Outputs

u = unfreevector(problem.S,u);

U = u(findddl(problem.S,DDL(DDLVECT('U',problem.S.syscoord,'TRANS'))),:);
Ux = u(findddl(problem.S,'UX'),:); % Ux = double(squeeze(eval_sol(problem.S,u,problem.S.node,'UX')));
Uy = u(findddl(problem.S,'UY'),:); % Uy = double(squeeze(eval_sol(problem.S,u,problem.S.node,'UY')));
Uz = u(findddl(problem.S,'UZ'),:); % Uz = double(squeeze(eval_sol(problem.S,u,problem.S.node,'UZ')));

R = u(findddl(problem.S,DDL(DDLVECT('R',problem.S.syscoord,'ROTA'))),:);
Rx = u(findddl(problem.S,'RX'),:); % Rx = double(squeeze(eval_sol(problem.S,u,problem.S.node,'RX'))));
Ry = u(findddl(problem.S,'RY'),:); % Ry = double(squeeze(eval_sol(problem.S,u,problem.S.node,'RY'))));
Rz = u(findddl(problem.S,'RZ'),:); % Rz = double(squeeze(eval_sol(problem.S,u,problem.S.node,'RZ'))));

P = getcenter(Q);

ux = eval_sol(problem.S,u,P,'UX');
uy = eval_sol(problem.S,u,P,'UY');
uz = eval_sol(problem.S,u,P,'UZ');

rx = eval_sol(problem.S,u,P,'RX');
ry = eval_sol(problem.S,u,P,'RY');
rz = eval_sol(problem.S,u,P,'RZ');

fprintf('\nRectangular table\n');
fprintf(['Load : ' loading '\n']);
fprintf(['Mesh : ' elemtype ' ' meshtype ' elements\n']);
fprintf('Nb elements = %g\n',getnbelem(S_plate));
fprintf('Span-to-thickness ratio = %g\n',max(a,b)/h);
fprintf('Elapsed time = %f s\n',time);
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

%% Save variables

save(fullfile(pathname,'solution.mat'),'u','U','R');

%% Display domains, boundary conditions and meshes

plotDomain(Q,L_beam,'Color','w','legend',false);
mysaveas(pathname,'domain',formats,renderer);
mymatlab2tikz(pathname,'domain.tex');

[hD,legD] = plotBoundaryConditions(problem.S,'legend',false);
switch loading
    case 'uniform'
        ampl = 2;
    case 'concentrated'
        ampl = 0.5;
end
[hN,legN] = vectorplot(problem.S,'F',problem.b,ampl,'r');
% legend([hD,hN],'Dirichlet','Neumann')
% legend([hD,hN],[legD,legN])
mysaveas(pathname,'boundary_conditions',formats,renderer);

plotModel(problem.S,'Color','k','FaceColor','k','FaceAlpha',0.1,'node',true,'legend',false);
mysaveas(pathname,'mesh',formats,renderer);

ampl = getsize(problem.S)/max(abs(u))/5;
plotModelDeflection(problem.S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'node',true,'legend',false);
mysaveas(pathname,'mesh_deflected',formats,renderer);

figure('Name','Meshes')
clf
plot(problem.S,'Color','k','FaceColor','k','FaceAlpha',0.1,'node',true);
plot(problem.S+ampl*u,'Color','b','FaceColor','b','FaceAlpha',0.1,'node',true);
mysaveas(pathname,'meshes_deflected',formats,renderer);

% plotFacets(problem.S);
% plotRidges(problem.S);

%% Display solution

% ampl = 0;
ampl = getsize(problem.S)/max(abs(u))/5;
options = {'solid',true};
% options = {};

plotSolution(problem.S,u,'displ',3,'ampl',ampl,options{:});
mysaveas(pathname,'Uz',formats,renderer);

% plotSolution(problem.S,u,'rotation',1,'ampl',ampl,options{:});
% mysaveas(pathname,'Rx',formats,renderer);

% plotSolution(problem.S,u,'rotation',2,'ampl',ampl,options{:});
% mysaveas(pathname,'Ry',formats,renderer);

end
end
end

% myparallel('stop');
