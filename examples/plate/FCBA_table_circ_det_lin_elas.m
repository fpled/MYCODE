%% FCBA table circular deterministic linear elasticity %%
%%-----------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% myparallel('start');

%% Input data

% test = 'stability_1'; % stability test under vertical load 1
% test = 'stability_2'; % stability test under vertical load 2
% test = 'stability_3'; % stability test under vertical load 3
% test = 'stability_4'; % stability test under vertical load 4
% test = 'static_hori_1'; % test under static horizontal load 1
% test = 'static_hori_2'; % test under static horizontal load 2
% test = 'static_hori_3'; % test under static horizontal load 3
% test = 'static_hori_4'; % test under static horizontal load 4
test = 'static_vert'; % test under static vertical load
% test = 'fatigue_1'; % fatigue test under horizontal load 1
% test = 'fatigue_2'; % fatigue test under horizontal load 2
% test = 'fatigue_3'; % fatigue test under horizontal load 3
% test = 'fatigue_4'; % fatigue test under horizontal load 4
% test = 'impact'; % vertical impact test
% test = 'drop'; % drop test

formats = {'fig','epsc2'};
renderer = 'OpenGL';

filename = ['FCBA_table_circ_det_lin_elas_' test];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Domains and meshes

% Plate
% Radius
r = 600e-3;
% Thickness
h = 25e-3;
C = CIRCLE(0.0,0.0,0.0,r);

% Beams
% Cross-section base
b_beam_top = 48e-3;
b_beam_bot = 38e-3;
b_beam = (b_beam_top+b_beam_bot)/2;
% Cross-section height
h_beam_top = 48e-3;
h_beam_bot = 38e-3;
h_beam = (h_beam_top+h_beam_bot)/2;
% Length
l = 710e-3+h/2;
% Distance between beams
a = 800e-3-h_beam_top;
b = 800e-3-b_beam_top;
L_beam{1} = LIGNE([-a/2,-b/2,0.0],[-a/2,-b/2,-l]);
L_beam{2} = LIGNE([a/2,-b/2,0.0],[a/2,-b/2,-l]);
L_beam{3} = LIGNE([a/2,b/2,0.0],[a/2,b/2,-l]);
L_beam{4} = LIGNE([-a/2,b/2,0.0],[-a/2,b/2,-l]);

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
cl_plate = r/15;
elemtype = 'DKT';
r_masse = 150e-3;
C_masse = CIRCLE(0.0,0.0,0.0,r_masse);
% points = [x_beam,x_load_stab,x_load_vert,x_load_fati{3:4}];
% S_plate = gmshcirclewithinclusionandpoints(C,C_masse,points,cl_plate,cl_plate,cl_plate,[pathname 'gmsh_plate_circ_' elemtype  '_cl_' num2str(cl_plate)],3);
Pi = [x_beam,x_load_stab,x_load_vert];
Pb = {getvertex(C,1),getvertex(C,2),getvertex(C,3),x_load_fati{4},getvertex(C,4),x_load_fati{3}};
S_plate = gmshFCBAtablecirc(C,C_masse,Pi,Pb,cl_plate,cl_plate,cl_plate,cl_plate,[pathname 'gmsh_plate_circ_' elemtype  '_cl_' num2str(cl_plate)],3);
S_plate = convertelem(S_plate,elemtype);

% Beam meshes
nbelem_beam = 80;
S_beam = cellfun(@(L) build_model(L,'nbelem',nbelem_beam,'elemtype','BEAM'),L_beam,'UniformOutput',false);
% cl_beam = 0.1;
% S_beam = cellfun(@(L,n) build_model(L,'cl',cl_plate,'elemtype','BEAM','filename',[pathname 'gmsh_beam_' num2str(n) '_cl_' num2str(cl_beam)]),L_beam,num2cell(1:length(L_beam)),'UniformOutput',false);

%% Materials

% Gravitational acceleration
g = 10;

% Plate
% Young modulus
E = 3.797e9;
% Poisson ratio
NU = 0.3;
% Density
mass_plate = 18.54;
RHO = mass_plate/(pi*r^2*h);
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
E_beam = 15e9;
% Poisson ratio
NU_beam = 0.3;
% Cross-section area
Sec_beam_top = b_beam_top*h_beam_top;
Sec_beam_bot = b_beam_bot*h_beam_bot;
Sec_beam = b_beam*h_beam;
% Density
Vol_beam = (l-h/2)*(Sec_beam_top+Sec_beam_bot+sqrt(Sec_beam_top*Sec_beam_bot))/3;
b_belt = 30e-3;
h_belt = 80e-3;
Vol_belt = 2*(a-h_beam_top)*b_belt*h_belt + 2*(b-b_beam_top)*b_belt*h_belt;
mass_beams = 8.48;
RHO_beam = mass_beams/(length(L_beam)*Vol_beam + Vol_belt);
% Planar second moment of area (or Planar area moment of inertia)
IY = h_beam*b_beam^3/12;
IZ = b_beam*h_beam^3/12;
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
switch test
    case 'stability_1'
        problem.S = addcl(problem.S,P_support(4)); % addcl(problem.S,P_support(4),{'U','R'},0);
        problem.S = addcl(problem.S,P_support(1),'U'); % addcl(problem.S,P_support(1),'U',0);
    case 'stability_2'
        problem.S = addcl(problem.S,P_support(1)); % addcl(problem.S,P_support(4),{'U','R'},0);
        problem.S = addcl(problem.S,P_support(2),'U'); % addcl(problem.S,P_support(1),'U',0);
    case 'stability_3'
        problem.S = addcl(problem.S,P_support(2)); % addcl(problem.S,P_support(4),{'U','R'},0);
        problem.S = addcl(problem.S,P_support(3),'U'); % addcl(problem.S,P_support(1),'U',0);
    case 'stability_4'
        problem.S = addcl(problem.S,P_support(3)); % addcl(problem.S,P_support(4),{'U','R'},0);
        problem.S = addcl(problem.S,P_support(4),'U'); % addcl(problem.S,P_support(1),'U',0);
    case {'static_hori_1','static_hori_2'}
        problem.S = addcl(problem.S,P_support([3 4])); % addcl(problem.S,P_support([3 4]),{'U','R'},0);
        problem.S = addcl(problem.S,P_support([1 2]),'UZ'); % addcl(problem.S,P_support([1 2]),'UZ',0);
    case {'static_hori_3','static_hori_4'}
        problem.S = addcl(problem.S,P_support([4 1])); % addcl(problem.S,P_support([4 1]),{'U','R'},0);
        problem.S = addcl(problem.S,P_support([2 3]),'UZ'); % addcl(problem.S,P_support([2 3]),'UZ',0);
    case 'static_vert'
        problem.S = addcl(problem.S,P_support,'U'); % addcl(problem.S,P_support,'U',0);
    case {'fatigue_1','fatigue_2','fatigue_3','fatigue_4',...
            'impact','drop'}
        problem.S = addcl(problem.S,P_support); % addcl(problem.S,P_support,{'U','R'},0);
end

%% Stiffness matrices and sollicitation vectors

p_plate = RHO*g*h;
p_beam = cellfun(@(S) RHO_beam*g*Sec_beam,S_beam,'UniformOutput',false);
switch test
    case {'stability_1','stability_2','stability_3','stability_4'}
        p = 400;
    case {'static_hori_1','static_hori_2','static_hori_3','static_hori_4'}
        masse = 50;
        Sec_masse = pi*r_masse^2;
        p_masse = masse*g/Sec_masse;
        p = 400;
        slope = 0;
    case 'static_vert'
        p = 1200;
    case {'fatigue_1','fatigue_2','fatigue_3','fatigue_4'}
        masse = 50;
        Sec_masse = pi*r_masse^2;
        p_masse = masse*g/Sec_masse;
        p = 300;
    case 'impact'
        H = 180e-3;
    case 'drop'
        H = 100e-3;
end

problem.A = calc_rigi(problem.S);
switch test
    case 'stability_1'
        problem.b = nodalload(problem.S,P_load_stab{1},'FZ',-p);
        if isempty(ispointin(P_load_stab{1},POINT(problem.S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
    case 'stability_2'
        problem.b = nodalload(problem.S,P_load_stab{2},'FZ',-p);
        if isempty(ispointin(P_load_stab{2},POINT(problem.S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
    case 'stability_3'
        problem.b = nodalload(problem.S,P_load_stab{3},'FZ',-p);
        if isempty(ispointin(P_load_stab{3},POINT(problem.S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
    case 'stability_4'
        problem.b = nodalload(problem.S,P_load_stab{4},'FZ',-p*cosd(slope));
        if isempty(ispointin(P_load_stab{4},POINT(problem.S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
    case {'static_hori_1','static_hori_2','static_hori_3','static_hori_4'}
        if strcmp(test,'static_hori_1')
            problem.b = nodalload(problem.S,P_load_hori{1},{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
            if isempty(ispointin(P_load_hori{1},POINT(problem.S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        elseif strcmp(test,'static_hori_2')
            problem.b = nodalload(problem.S,P_load_hori{2},{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
            if isempty(ispointin(P_load_hori{2},POINT(problem.S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        elseif strcmp(test,'static_hori_3')
            problem.b = nodalload(problem.S,P_load_hori{3},{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
            if isempty(ispointin(P_load_hori{3},POINT(problem.S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        elseif strcmp(test,'static_hori_4')
            problem.b = nodalload(problem.S,P_load_hori{4},{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
            if isempty(ispointin(P_load_hori{4},POINT(problem.S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        end
        problem.b = problem.b + bodyload(keepgroupelem(problem.S,2),[],'FZ',-p_masse);
    case 'static_vert'
        problem.b = nodalload(problem.S,P_load_vert,'FZ',-p);
        if isempty(ispointin(P_load_vert,POINT(problem.S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
    case {'fatigue_1','fatigue_2','fatigue_3','fatigue_4'}
        if strcmp(test,'fatigue_1')
            problem.b = nodalload(problem.S,P_load_fati{1},'FY',-p);
            if isempty(ispointin(P_load_fati{1},POINT(problem.S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        elseif strcmp(test,'fatigue_2')
            problem.b = nodalload(problem.S,P_load_fati{2},'FY',p);
            if isempty(ispointin(P_load_fati{2},POINT(problem.S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        elseif strcmp(test,'fatigue_3')
            problem.b = nodalload(problem.S,P_load_fati{3},'FX',p);
            if isempty(ispointin(P_load_fati{3},POINT(problem.S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        elseif strcmp(test,'fatigue_4')
            problem.b = nodalload(problem.S,P_load_fati{4},'FX',-p);
            if isempty(ispointin(P_load_fati{4},POINT(problem.S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        end
        problem.b = problem.b + bodyload(keepgroupelem(problem.S,2),[],'FZ',-p_masse);
    case {'impact','drop'}
        error('Not implemented')
end
problem.b = problem.b + bodyload(keepgroupelem(problem.S,[1,2]),[],'FZ',-p_plate);
for k=1:length(P_beam)
    problem.b = problem.b + bodyload(keepgroupelem(problem.S,2+k),[],'FX',p_beam{k});
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

P = getcenter(C);

ux = eval_sol(problem.S,u,P,'UX');
uy = eval_sol(problem.S,u,P,'UY');
uz = eval_sol(problem.S,u,P,'UZ');

rx = eval_sol(problem.S,u,P,'RX');
ry = eval_sol(problem.S,u,P,'RY');
rz = eval_sol(problem.S,u,P,'RZ');

fprintf('\nCircular table\n');
fprintf(['Test : ' test '\n']);
fprintf(['Mesh : ' elemtype ' elements\n']);
fprintf('Nb elements = %g\n',getnbelem(S_plate));
fprintf('Span-to-thickness ratio = %g\n',r/h);
fprintf('Elapsed time = %f s\n',time);
fprintf('\n');

disp('Displacement u at point'); disp(P);
fprintf('ux    = %g\n',ux);
fprintf('uy    = %g\n',uy);
fprintf('uz    = %g\n',uz);
if strcmp(test,'static_vert')
    uz_exp = -2.35e-3;
    err_uz = norm(uz-uz_exp)/norm(uz_exp);
    fprintf('uz_exp= %g\n',uz_exp);
    fprintf('error = %g\n',err_uz);
end
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
switch test
    case {'stability_1','stability_2','stability_3','stability_4',...
            'static_hori_1','static_hori_2','static_hori_3','static_hori_4',...
            'static_vert',...
            'fatigue_1','fatigue_2','fatigue_3','fatigue_4',...
            'impact','drop'}
        ampl = 5;
end
[hN,legN] = vectorplot(problem.S,'F',problem.b,ampl,'r');
% legend([hD,hN],'Dirichlet','Neumann')
% legend([hD,hN],[legD,legN])
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
options = {'solid',true};
% options = {};

plotSolution(problem.S,u,'displ',3,'ampl',ampl,options{:});
mysaveas(pathname,'Uz',formats,renderer);

% plotSolution(problem.S,u,'rotation',1,'ampl',ampl,options{:});
% mysaveas(pathname,'Rx',formats,renderer);

% plotSolution(problem.S,u,'rotation',2,'ampl',ampl,options{:});
% mysaveas(pathname,'Ry',formats,renderer);

% myparallel('stop');