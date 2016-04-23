%% Table rectangular deterministic linear elasticity %%
%%---------------------------------------------------%%

% clc
clear all
close all

%% Input data
% loadings = {'uniform'};
% loadings = {'concentrated'};
loadings={'uniform','concentrated'};
% elemtypes = {'DKT'};
% elemtypes = {'DKQ'};
% elemtypes = {'COQ4'};
elemtypes = {'DKT','DKQ'};
% elemtypes = {'DKT','DKQ','COQ4'};
% meshtypes = 'structured';
% meshtypes = {'unstructured'};
meshtypes = {'structured','unstructured'};

for indexl=1:length(loadings)
    loading = loadings{indexl};
for indexe=1:length(elemtypes)
    elemtype = elemtypes{indexe};
for indexm=1:length(meshtypes)
    meshtype = meshtypes{indexm};

close all

filename = strcat('table_rect_det_lin_elas_',loading,'_',elemtype,'_',meshtype);
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'RESULTS',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
% set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
renderer = 'OpenGL';

% Parallel computing
% myparallel('start');

%% Domains and meshes

a = 1;
b = 1;
Q = QUADRANGLE([0.0,0.0,0.0],[a,0.0,0.0],[a,b,0.0],[0.0,b,0.0]);

P_load = getcenter(Q);
x_load = double(getcoord(P_load));

l = 1;
L_beam{1} = LIGNE([1/5*a,1/5*b,0.0],[1/5*a,1/5*b,-l]);
L_beam{2} = LIGNE([4/5*a,1/5*b,0.0],[4/5*a,1/5*b,-l]);
L_beam{3} = LIGNE([4/5*a,4/5*b,0.0],[4/5*a,4/5*b,-l]);
L_beam{4} = LIGNE([1/5*a,4/5*b,0.0],[1/5*a,4/5*b,-l]);

x_beam = cellfun(@(L) getvertex(L,1),L_beam,'UniformOutput',false);
P_beam = cellfun(@(x) POINT(x),x_beam,'UniformOutput',false);

switch meshtype
    case 'structured'
        nbelem_plate = [20,20];
        S_plate = build_model(Q,'nbelem',nbelem_plate,'elemtype',elemtype);
    case 'unstructured'
        cl_plate = 0.05;
        switch loading
            case 'uniform'
                points = x_beam;
            case 'concentrated'
                points = [x_beam(:)',{x_load}];
        end
        S_plate = build_model(Q,'cl',cl_plate,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_rect_' elemtype  '_cl_' num2str(cl_plate)],'points',points);
end

nbelem_beam = [10,10];
S_beam = cellfun(@(L) build_model(L,'nbelem',nbelem_beam,'elemtype','BEAM'),L_beam,'UniformOutput',false);
% cl_beam = 0.1;
% S_beam = cellfun(@(L,n) build_model(L,'cl',cl_plate,'elemtype','BEAM','filename',[pathname 'gmsh_beam_' num2str(n) '_cl_' num2str(cl_beam)]),L_beam,num2cell(1:length(L_beam)),'UniformOutput',false);

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
S_beam = cellfun(@(S) setmaterial(S,mat_beam),S_beam,'UniformOutput',false);

S_beams = union(S_beam{:});
S_beams = concatgroupelem(S_beams);
system.S = union(S_plate,S_beams);

%% Dirichlet boundary conditions

x_support = cellfun(@(L) getvertex(L,2),L_beam,'UniformOutput',false);
P_support = cellfun(@(x) POINT(x),x_support,'UniformOutput',false);

system.S = final(system.S);
for k=1:length(P_support)
    system.S = addcl(system.S,P_support{k}); % addcl(system.S,P_support{k},{'U','R'},0);
end

%% Stiffness matrices and sollicitation vectors

% Uniform or Concentrated load
switch loading
    case 'uniform'
        p = RHO*g*h;
    case 'concentrated'
        p = RHO*g*h*a*b;
end

system.A = calc_rigi(system.S);
switch loading
    case 'uniform'
        system.b = bodyload(keepgroupelem(system.S,1),[],'FZ',-p);
    case 'concentrated'
        system.b = nodalload(system.S,P_load,'FZ',-p);
        if isempty(ispointin(P_load,POINT(system.S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
end

%% Resolution

t = tic;
u = solve_system(system);
time = toc(t);
fprintf('\nRectangular table\n');
fprintf(['Load : ' loading '\n']);
fprintf(['Mesh : ' meshtype ' with ' elemtype ' elements\n']);
fprintf('Span-to-thickness ratio = %g\n',max(a,b)/h);
fprintf('Elapsed time = %f s\n',time);
fprintf('\n');

%% Outputs

u = unfreevector(system.S,u);

U = u(findddl(system.S,DDL(DDLVECT('U',system.S.syscoord,'TRANS'))));
Ux = u(findddl(system.S,'UX'),:); % Ux = double(squeeze(eval_sol(system.S,u,system.S.node,'UX')));
Uy = u(findddl(system.S,'UY'),:); % Uy = double(squeeze(eval_sol(system.S,u,system.S.node,'UY')));
Uz = u(findddl(system.S,'UZ'),:); % Uz = double(squeeze(eval_sol(system.S,u,system.S.node,'UZ')));

R = u(findddl(system.S,DDL(DDLVECT('R',system.S.syscoord,'ROTA'))));
Rx = u(findddl(system.S,'RX'),:); % Rx = double(squeeze(eval_sol(system.S,u,system.S.node,'RX'))));
Ry = u(findddl(system.S,'RY'),:); % Ry = double(squeeze(eval_sol(system.S,u,system.S.node,'RY'))));
Rz = u(findddl(system.S,'RZ'),:); % Rz = double(squeeze(eval_sol(system.S,u,system.S.node,'RZ'))));

P = getcenter(Q);

ux = eval_sol(system.S,u,P,'UX');
uy = eval_sol(system.S,u,P,'UY');
uz = eval_sol(system.S,u,P,'UZ');

rx = eval_sol(system.S,u,P,'RX');
ry = eval_sol(system.S,u,P,'RY');
rz = eval_sol(system.S,u,P,'RZ');

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

plot_domain(Q,L_beam,'color','w','nolegend');
mysaveas(pathname,'domain',{'fig','epsc2'},renderer);
mymatlab2tikz(pathname,'domain.tex');

[hD,legD] = plot_boundary_conditions(system.S,'nolegend');
switch loading
    case 'uniform'
        ampl = 2;
    case 'concentrated'
        ampl = 0.5;
end
[hN,legN] = vectorplot(system.S,'F',system.b,ampl,'r');
% legend([hD,hN],'Dirichlet','Neumann')
legend([hD,hN],[legD,legN])
axis image
mysaveas(pathname,'boundary_conditions',{'fig','epsc2'},renderer);

plot_model(system.S,'color','k','facecolor','k','facealpha',0.1,'node','nolegend');
mysaveas(pathname,'mesh',{'fig','epsc2'},renderer);

ampl = max(getsize(system.S))/max(abs(u))/2;
plot_model_deflection(system.S,u,'ampl',ampl,'color','b','facecolor','b','facealpha',0.1,'node','nolegend');
mysaveas(pathname,'mesh_deflected',{'fig','epsc2'},renderer);

figure('Name','Meshes')
clf
plot(system.S,'color','k','facecolor','k','facealpha',0.1,'node');
plot(system.S+ampl*u,'color','b','facecolor','b','facealpha',0.1,'node');
mysaveas(pathname,'meshes_deflected',{'fig','epsc2'},renderer);

% plot_facets(system.S);
% plot_ridges(system.S);

%% Display solution

% ampl = 0;
ampl = max(getsize(system.S))/max(abs(u))/2;
% options = {'solid'};
options = {};

plot_solution(system.S,u,'displ',3,'ampl',ampl,options{:});
mysaveas(pathname,'Uz',{'fig','epsc2'},renderer);

% plot_solution(system.S,u,'rotation',1,'ampl',ampl,options{:});
% mysaveas(pathname,'Rx',{'fig','epsc2'},renderer);

% plot_solution(system.S,u,'rotation',2,'ampl',ampl,options{:});
% mysaveas(pathname,'Ry',{'fig','epsc2'},renderer);

end
end
end

% myparallel('stop');
