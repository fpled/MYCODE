%% FCBA table circular deterministic linear elasticity %%
%%-----------------------------------------------------%%

% clc
% clear all
close all
% set(0,'DefaultFigureVisible','off');

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

belt = 1; % belt modelisation

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
RHO = mass_plate/(pi*r^2*h);
% Extensional stiffness (or Membrane rigidity)
A_rig = E*h/(1-NU^2);
% Bending stiffness (or Flexural rigidity)
D_rig = E*h^3/(12*(1-NU^2));
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

%% Dirichlet boundary conditions

x_support = cellfun(@(L) getvertex(L,2)',L_beam,'UniformOutput',false);
x_support = [x_support{:}]';
P_support = POINT(x_support);

S = final(S);
switch test
    case 'stability_1'
        S = addcl(S,P_support(4)); % addcl(S,P_support(4),{'U','R'},0);
        S = addcl(S,P_support(1),'U'); % addcl(S,P_support(1),'U',0);
    case 'stability_2'
        S = addcl(S,P_support(1)); % addcl(S,P_support(4),{'U','R'},0);
        S = addcl(S,P_support(2),'U'); % addcl(S,P_support(1),'U',0);
    case 'stability_3'
        S = addcl(S,P_support(2)); % addcl(S,P_support(4),{'U','R'},0);
        S = addcl(S,P_support(3),'U'); % addcl(S,P_support(1),'U',0);
    case 'stability_4'
        S = addcl(S,P_support(3)); % addcl(S,P_support(4),{'U','R'},0);
        S = addcl(S,P_support(4),'U'); % addcl(S,P_support(1),'U',0);
    case {'static_hori_1','static_hori_2'}
        S = addcl(S,P_support([3 4])); % addcl(S,P_support([3 4]),{'U','R'},0);
        S = addcl(S,P_support([1 2]),'UZ'); % addcl(S,P_support([1 2]),'UZ',0);
    case {'static_hori_3','static_hori_4'}
        S = addcl(S,P_support([4 1])); % addcl(S,P_support([4 1]),{'U','R'},0);
        S = addcl(S,P_support([2 3]),'UZ'); % addcl(S,P_support([2 3]),'UZ',0);
    case 'static_vert'
        S = addcl(S,P_support,'U'); % addcl(S,P_support,'U',0);
    case {'fatigue_1','fatigue_2','fatigue_3','fatigue_4',...
            'impact','drop'}
        S = addcl(S,P_support); % addcl(S,P_support,{'U','R'},0);
end

%% Stiffness matrices and sollicitation vectors

p_plate = RHO*g*h;
p_beam = RHO_beam*g*Sec_beam;
p_belt = RHO_beam*g*Sec_belt;
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

A = calc_rigi(S);
switch test
    case 'stability_1'
        f = nodalload(S,P_load_stab{1},'FZ',-p);
        if isempty(ispointin(P_load_stab{1},POINT(S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
    case 'stability_2'
        f = nodalload(S,P_load_stab{2},'FZ',-p);
        if isempty(ispointin(P_load_stab{2},POINT(S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
    case 'stability_3'
        f = nodalload(S,P_load_stab{3},'FZ',-p);
        if isempty(ispointin(P_load_stab{3},POINT(S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
    case 'stability_4'
        f = nodalload(S,P_load_stab{4},'FZ',-p*cosd(slope));
        if isempty(ispointin(P_load_stab{4},POINT(S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
    case {'static_hori_1','static_hori_2','static_hori_3','static_hori_4'}
        if strcmp(test,'static_hori_1')
            f = nodalload(S,P_load_hori{1},{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
            if isempty(ispointin(P_load_hori{1},POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        elseif strcmp(test,'static_hori_2')
            f = nodalload(S,P_load_hori{2},{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
            if isempty(ispointin(P_load_hori{2},POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        elseif strcmp(test,'static_hori_3')
            f = nodalload(S,P_load_hori{3},{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
            if isempty(ispointin(P_load_hori{3},POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        elseif strcmp(test,'static_hori_4')
            f = nodalload(S,P_load_hori{4},{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
            if isempty(ispointin(P_load_hori{4},POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        end
        f = f + bodyload(keepgroupelem(S,3),[],'FZ',-p_masse);
    case 'static_vert'
        f = nodalload(S,P_load_vert,'FZ',-p);
        if isempty(ispointin(P_load_vert,POINT(S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
    case {'fatigue_1','fatigue_2','fatigue_3','fatigue_4'}
        if strcmp(test,'fatigue_1')
            f = nodalload(S,P_load_fati{1},'FY',-p);
            if isempty(ispointin(P_load_fati{1},POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        elseif strcmp(test,'fatigue_2')
            f = nodalload(S,P_load_fati{2},'FY',p);
            if isempty(ispointin(P_load_fati{2},POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        elseif strcmp(test,'fatigue_3')
            f = nodalload(S,P_load_fati{3},'FX',p);
            if isempty(ispointin(P_load_fati{3},POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
        elseif strcmp(test,'fatigue_4')
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
    f = f + bodyload(keepgroupelem(S,3+k),[],'FX',p_beam);
end
f = f + bodyload(keepgroupelem(S,3+length(L_beam)+1),[],'FY',p_belt);
f = f + bodyload(keepgroupelem(S,3+length(L_beam)+2),[],'FZ',p_belt);
f = f + bodyload(keepgroupelem(S,3+length(L_beam)+3),[],'FY',-p_belt);
f = f + bodyload(keepgroupelem(S,3+length(L_beam)+4),[],'FZ',-p_belt);

%% Resolution

t = tic;
u = A\f;
time = toc(t);

%% Outputs

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

fprintf('\nCircular table\n');
fprintf(['Test : ' test '\n']);
fprintf(['Mesh : ' elemtype ' elements\n']);
fprintf('Nb elements = %g\n',getnbelem(S_plate));
fprintf('Nb dofs     = %g\n',getnbddl(S_plate));
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

%% Save variables

save(fullfile(pathname,'solution.mat'),'u','U','R');

%% Display domains, boundary conditions and meshes

if ~belt
    plotDomain(C,L_beam,'legend',false);
else
    plotDomain(C,[L_beam,L_belt],'legend',false);
end
mysaveas(pathname,'domain',formats,renderer);
mymatlab2tikz(pathname,'domain.tex');

[hD,legD] = plotBoundaryConditions(S,'legend',false);
switch test
    case {'stability_1','stability_2','stability_3','stability_4',...
            'static_hori_1','static_hori_2','static_hori_3','static_hori_4',...
            'static_vert',...
            'fatigue_1','fatigue_2','fatigue_3','fatigue_4',...
            'impact','drop'}
        ampl = 5;
end
[hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
% legend([hD,hN],'Dirichlet','Neumann')
% legend([hD,hN],[legD,legN])
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
