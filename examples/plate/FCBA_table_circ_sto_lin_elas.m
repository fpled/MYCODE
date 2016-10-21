%% FCBA table circular stochastic linear elasticity %%
%%--------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% myparallel('start');

%% Input data

formats = {'fig','epsc2'};
renderer = 'OpenGL';

filename = 'FCBA_table_circ_sto_lin_elas';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Domains and meshes

r = 600e-3;
C = CIRCLE(0.0,0.0,0.0,r);

x_load_hori = getvertices(C); % essai de charge statique horizontale
x_load_vert = double(getcenter(C)); % essai de charge statique verticale
x_load_stab = [0.0,-r+50e-3,0.0]; % essai de stabilité sous charge verticale
x_load = [x_load_hori,{x_load_vert},{x_load_stab}];
P_load = cellfun(@(x) POINT(x),x_load,'UniformOutput',false);

l = (750-25/2)*1e-3;
a = (800-43)*1e-3;
L_beam{1} = LIGNE([-a/2,-a/2,0.0],[-a/2,-a/2,-l]);
L_beam{2} = LIGNE([a/2,-a/2,0.0],[a/2,-a/2,-l]);
L_beam{3} = LIGNE([a/2,a/2,0.0],[a/2,a/2,-l]);
L_beam{4} = LIGNE([-a/2,a/2,0.0],[-a/2,a/2,-l]);

x_beam = cellfun(@(L) getvertex(L,1),L_beam,'UniformOutput',false);
P_beam = cellfun(@(x) POINT(x),x_beam,'UniformOutput',false);

cl_plate = 0.1;
elemtype = 'DKT';
points = [x_beam,x_load];
S_plate = build_model(C,'cl',cl_plate,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_circ_' elemtype  '_cl_' num2str(cl_plate)],'points',points);

nbelem_beam = 10;
S_beam = build_model(L_beam,'nbelem',nbelem_beam,'elemtype','BEAM');
% cl_beam = 0.1;
% S_beam = build_model(L_beam,'cl',cl_beam,'elemtype','BEAM','filename',[pathname 'gmsh_beam_cl_' num2str(cl_beam)]);

%% Materials

% Gravitational acceleration
g = 10;

% Plate
% Young modulus
E = 4e9;
% Poisson ratio
NU = 0.3;
% Density
RHO = 1;
% Thickness
h = 25e-3;
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
E_beam = 4e9;
% Poisson ratio
NU_beam = 0.3;
% Density
RHO_beam = 1;
% Base (along y)
b_beam = 43e-3;
% Height (along z)
h_beam = 43e-3;
% Section
Sec_beam = b_beam*h_beam;
% Planar second moment of area (or Planar area moment of inertia)
IY = b_beam*h_beam^3/12;
IZ = h_beam*b_beam^3/12;
% Polar second moment of area (or Polar area moment of inertia)
IX = IY+IZ;
% Material
mat_beam = ELAS_BEAM('E',E_beam,'NU',NU_beam,'S',Sec_beam,'IZ',IZ,'IY',IY,'IX',IX,'RHO',RHO_beam);
mat_beam = setnumber(mat_beam,2);
S_beam = setmaterial(S_beam,mat_beam);

problem.S = union(S_plate,S_beam);

%% Dirichlet boundary conditions

x_support = cellfun(@(L) getvertex(L,2),L_beam,'UniformOutput',false);
P_support = cellfun(@(x) POINT(x),x_support,'UniformOutput',false);

problem.S = final(problem.S);
for k=1:length(P_support)
    problem.S = addcl(problem.S,P_support{k}); % addcl(problem.S,P_support{k},{'U','R'},0);
end

%% Stiffness matrices and sollicitation vectors

% Essai de charge verticale
p = 400;
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
fprintf('Nb elements = %g\n',getnbelem(S_plate));
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

ampl = max(getsize(problem.S))/max(abs(u))/2;
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
options = {'solid',true};
% options = {};

plotSolution(problem.S,u,'displ',3,'ampl',ampl,options{:});
mysaveas(pathname,'Uz',formats,renderer);

% plotSolution(problem.S,u,'rotation',1,'ampl',ampl,options{:});
% mysaveas(pathname,'Rx',formats,renderer);

% plotSolution(problem.S,u,'rotation',2,'ampl',ampl,options{:});
% mysaveas(pathname,'Ry',formats,renderer);

% myparallel('stop');
