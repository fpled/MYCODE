%% Table circular deterministic linear elasticity %%
%%---------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% myparallel('start');

%% Input data

% loadings = {'uniform'};
% loadings = {'concentrated'};
loadings={'uniform','concentrated'};
% elemtypes = {'DKT'};
% elemtypes = {'DKQ'};
% elemtypes = {'COQ4'};
% elemtypes = {'DKT','DKQ'};
elemtypes = {'DKT','DKQ','COQ4'};

formats = {'fig','epsc2'};
renderer = 'OpenGL';

for il=1:length(loadings)
    loading = loadings{il};
    filename = ['table_circ_det_lin_elas_' loading];
    
for ie=1:length(elemtypes)
    elemtype = elemtypes{ie};
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep,elemtype,filesep);
    if ~exist(pathname,'dir')
        mkdir(pathname);
    end

%% Domains and meshes

r = 1;
C = CIRCLE(0.0,0.0,0.0,r);

P_load = POINT([-r/2,0.0,0.0]);
x_load = double(getcoord(P_load));

P_beam = getcenter(C);
% P_beam = POINT([r/2,0.0,0.0]);
x_beam = double(getcoord(P_beam));

l = 1;
L_beam = LIGNE(P_beam,P_beam+POINT([0.0,0.0,-l]));

cl_plate = 0.1;
switch loading
    case 'uniform'
        points = x_beam;
    case 'concentrated'
        points = {x_beam,x_load};
end
S_plate = build_model(C,'cl',cl_plate,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_circ_' elemtype  '_cl_' num2str(cl_plate)],'points',points);

nbelem_beam = 10;
S_beam = build_model(L_beam,'nbelem',nbelem_beam,'elemtype','BEAM');
% cl_beam = 0.1;
% S_beam = build_model(L_beam,'cl',cl_beam,'elemtype','BEAM','filename',[pathname 'gmsh_beam_cl_' num2str(cl_beam)]);

%% Materials

% Gravitational acceleration
g = 10;
% Young modulus
E = 1;
% Poisson ratio
NU = 0.3;
% Thickness
h = 0.1;
% Density
RHO = 1;
% Extensional stiffness (or Membrane rigidity)
A = E*h/(1-NU^2);
% Bending stiffness (or Flexural rigidity)
D = E*h^3/(12*(1-NU^2));

% Plate
mat_plate = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',h,'k',5/6);
mat_plate = setnumber(mat_plate,1);
S_plate = setmaterial(S_plate,mat_plate);

% Radius
r_beam = 0.1;
% Section
Sec_beam = pi*r_beam^2;
% Planar second moment of area (or Planar area moment of inertia)
IY = pi*r_beam^4/2;
IZ = IY;
% Polar second moment of area (or Polar area moment of inertia)
IX = IY+IZ;

% Beam
mat_beam = ELAS_BEAM('E',E,'NU',NU,'S',Sec_beam,'IZ',IZ,'IY',IY,'IX',IX,'RHO',RHO);
mat_beam = setnumber(mat_beam,2);
S_beam = setmaterial(S_beam,mat_beam);

problem.S = union(S_plate,S_beam);

%% Dirichlet boundary conditions

x_support = getvertex(L_beam,2);
P_support = POINT(x_support);

problem.S = final(problem.S);
problem.S = addcl(problem.S,P_support); % addcl(problem.S,P_support,{'U','R'},0);

%% Stiffness matrices and sollicitation vectors

% Uniform or Concentrated load
switch loading
    case 'uniform'
        p = RHO*g*h;
    case 'concentrated'
        p = RHO*g*h*r^2;
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

U = u(findddl(problem.S,DDL(DDLVECT('U',problem.S.syscoord,'TRANS'))));
Ux = u(findddl(problem.S,'UX'),:); % Ux = double(squeeze(eval_sol(problem.S,u,problem.S.node,'UX')));
Uy = u(findddl(problem.S,'UY'),:); % Uy = double(squeeze(eval_sol(problem.S,u,problem.S.node,'UY')));
Uz = u(findddl(problem.S,'UZ'),:); % Uz = double(squeeze(eval_sol(problem.S,u,problem.S.node,'UZ')));

R = u(findddl(problem.S,DDL(DDLVECT('R',problem.S.syscoord,'ROTA'))));
Rx = u(findddl(problem.S,'RX'),:); % Rx = double(squeeze(eval_sol(problem.S,u,problem.S.node,'RX'))));
Ry = u(findddl(problem.S,'RY'),:); % Ry = double(squeeze(eval_sol(problem.S,u,problem.S.node,'RY'))));
Rz = u(findddl(problem.S,'RZ'),:); % Rz = double(squeeze(eval_sol(problem.S,u,problem.S.node,'RZ'))));

P = getcenter(C);

ux = eval_sol(problem.S,u,P,'UX');
uy = eval_sol(problem.S,u,P,'UY');
uz = eval_sol(problem.S,u,P,'UZ');

rx = eval_sol(problem.S,u,P,'RX');
ry = eval_sol(problem.S,u,P,'RY');
rz = eval_sol(problem.S,u,P,'RZ');

fprintf('\nCircular table\n');
fprintf(['Load : ' loading '\n']);
fprintf(['Mesh : ' elemtype ' elements\n']);
fprintf('Nb elements = %g\n',getnbelem(problem.S));
fprintf('Span-to-thickness ratio = %g\n',r/h);
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
save(fullfile(pathname,'all.mat'));

%% Display domains, boundary conditions and meshes

plotDomain(C,L_beam,'legend',false);
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
axis image
mysaveas(pathname,'boundary_conditions',formats,renderer);

plotModel(problem.S,'Color','k','FaceColor','k','FaceAlpha',0.1,'node',true,'legend',false);
mysaveas(pathname,'mesh',formats,renderer);

ampl = max(getsize(problem.S))/max(abs(u))/10;
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
ampl = max(getsize(problem.S))/max(abs(u))/10;
% options = {'solid',true};
options = {};

plotSolution(problem.S,u,'displ',3,'ampl',ampl,options{:});
mysaveas(pathname,'Uz',formats,renderer);

% plotSolution(problem.S,u,'rotation',1,'ampl',ampl,options{:});
% mysaveas(pathname,'Rx',formats,renderer);

% plotSolution(problem.S,u,'rotation',2,'ampl',ampl,options{:});
% mysaveas(pathname,'Ry',formats,renderer);

end
end

% myparallel('stop');
