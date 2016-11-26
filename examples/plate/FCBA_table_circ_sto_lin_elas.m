%% FCBA table circular deterministic linear elasticity %%
%%-----------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
rng('default');
myparallel('start');

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

filename = ['FCBA_table_circ_sto_lin_elas_' test];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

t = tic;

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

%% Random variables

% Experimental data
samples_E = [4.211 4.057 3.685 3.921 3.839 3.845 3.795...
    3.406 3.389 3.299 3.485 3.319 3.267 3.349 3.307...
    4.684 4.245 4.076 4.407 4.283 4.054 4.226 4.041...
    4.104 4.075 3.556 3.319 3.848 3.707 3.664 3.493 3.550]*1e9;
% Parameters for Gamma distribution
phat = gamfit(samples_E);
% Number of samples
N = 1e3;
% Sample set
e = gamrnd(phat(1),phat(2),1,N);

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
RHO = mass_plate/(pi*r^2*h);
% Extensional stiffness (or Membrane rigidity)
A_rig = E*h/(1-NU^2);
% Bending stiffness (or Flexural rigidity)
D_rig = E*h^3/(12*(1-NU^2));

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

S = union(S_plate,S_beam{:});

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

A = cell(1,N);
for i=1:N
    % Young modulus
    Ei = e(i);
    
    % Material
    mat_platei = setparam(mat_plate,'E',Ei);
    Si = setmaterial(S,mat_platei,[1,2]);
    
    % Stiffness matrix
    A{i} = calc_rigi(Si);
end

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
        f = f + bodyload(keepgroupelem(S,2),[],'FZ',-p_masse);
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
        f = f + bodyload(keepgroupelem(S,2),[],'FZ',-p_masse);
    case {'impact','drop'}
        error('Not implemented')
end
f = f + bodyload(keepgroupelem(S,[1,2]),[],'FZ',-p_plate);
for k=1:length(P_beam)
    f = f + bodyload(keepgroupelem(S,2+k),[],'FX',p_beam);
end

%% Resolution

parfor i=1:N
    u(:,i) = A{i}\f;
end

time = toc(t);

%% Convergence Monte-Carlo

fontsize = 16;
interpreter = 'latex';

figure('Name','Convergence empirical mean')
clf
mean_u = arrayfun(@(x) norm(mean(u(:,1:x),2)),1:N);
plot(1:N,mean_u,'-b','LineWidth',1)
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
ylabel('Moyenne empirique','Interpreter',interpreter)
mysaveas(pathname,'convergence_empirical_mean','fig');
mymatlab2tikz(pathname,'convergence_empirical_mean.tex');

figure('Name','Convergence empirical standard deviation')
clf
std_u = arrayfun(@(x) norm(std(u(:,1:x),0,2)),1:N);
plot(1:N,std_u,'-r','LineWidth',1)
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
ylabel('Ecart-type empirique','Interpreter',interpreter)
mysaveas(pathname,'convergence_empirical_std','fig');
mymatlab2tikz(pathname,'convergence_empirical_std.tex');

%% Outputs

mean_u = mean(u,2);
std_u = std(u,0,2);

mean_u = unfreevector(S,mean_u);
std_u = unfreevector(S,std_u);

mean_U = mean_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
mean_Ux = mean_u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,mean_u,S.node,'UX')),:);
mean_Uy = mean_u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,mean_u,S.node,'UY')),:);
mean_Uz = mean_u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,mean_u,S.node,'UZ')),:);

mean_R = mean_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
mean_Rx = mean_u(findddl(S,'RX'),:); % Rx = double(squeeze(eval_sol(S,mean_u,S.node,'RX'))),:);
mean_Ry = mean_u(findddl(S,'RY'),:); % Ry = double(squeeze(eval_sol(S,mean_u,S.node,'RY'))),:);
mean_Rz = mean_u(findddl(S,'RZ'),:); % Rz = double(squeeze(eval_sol(S,mean_u,S.node,'RZ'))),:);

std_U = std_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
std_Ux = std_u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,std_u,S.node,'UX')),:);
std_Uy = std_u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,std_u,S.node,'UY')),:);
std_Uz = std_u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,std_u,S.node,'UZ')),:);

std_R = std_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
std_Rx = std_u(findddl(S,'RX'),:); % Rx = double(squeeze(eval_sol(S,std_u,S.node,'RX'))),:);
std_Ry = std_u(findddl(S,'RY'),:); % Ry = double(squeeze(eval_sol(S,std_u,S.node,'RY'))),:);
std_Rz = std_u(findddl(S,'RZ'),:); % Rz = double(squeeze(eval_sol(S,std_u,S.node,'RZ'))),:);

P = getcenter(C);

mean_ux = eval_sol(S,mean_u,P,'UX');
mean_uy = eval_sol(S,mean_u,P,'UY');
mean_uz = eval_sol(S,mean_u,P,'UZ');

mean_rx = eval_sol(S,mean_u,P,'RX');
mean_ry = eval_sol(S,mean_u,P,'RY');
mean_rz = eval_sol(S,mean_u,P,'RZ');

std_ux = eval_sol(S,std_u,P,'UX');
std_uy = eval_sol(S,std_u,P,'UY');
std_uz = eval_sol(S,std_u,P,'UZ');

std_rx = eval_sol(S,std_u,P,'RX');
std_ry = eval_sol(S,std_u,P,'RY');
std_rz = eval_sol(S,std_u,P,'RZ');

fprintf('\nCircular table\n');
fprintf(['Test : ' test '\n']);
fprintf(['Mesh : ' elemtype ' elements\n']);
fprintf('Nb elements = %g\n',getnbelem(S_plate));
fprintf('Span-to-thickness ratio = %g\n',r/h);
fprintf('Nb samples = %g\n',N);
fprintf('Elapsed time = %f s\n',time);
fprintf('\n');

disp('Displacement u at point'); disp(P);
fprintf('mean(ux) = %g\n',mean_ux);
fprintf('std(ux)  = %g\n',std_ux);
fprintf('mean(uy) = %g\n',mean_uy);
fprintf('std(uy)  = %g\n',std_uy);
fprintf('mean(uz) = %g\n',mean_uz);
fprintf('std(uz)  = %g\n',std_uz);
if strcmp(test,'static_vert')
    uz_exp = -2.35e-3;
    err_uz = norm(mean_uz-uz_exp)/norm(uz_exp);
    fprintf('uz_exp   = %g\n',uz_exp);
    fprintf('error    = %g\n',err_uz);
end
fprintf('\n');

disp('Rotation r at point'); disp(P);
fprintf('mean(rx) = %g\n',mean_rx);
fprintf('std(rx)  = %g\n',std_rx);
fprintf('mean(ry) = %g\n',mean_ry);
fprintf('std(ry)  = %g\n',std_ry);
fprintf('mean(rz) = %g\n',mean_rz);
fprintf('std(rz)  = %g\n',std_rz);
fprintf('\n');

%% Save variables

save(fullfile(pathname,'solution.mat'),...
    'mean_u','mean_U','mean_R',...
    'std_u','std_U','std_R');

%% Display domains, boundary conditions and meshes

plotDomain(C,L_beam,'legend',false);
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

ampl = max(getsize(S))/max(abs(mean_u))/10;
plotModelDeflection(S,mean_u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'node',true,'legend',false);
mysaveas(pathname,'mesh_deflected',formats,renderer);

figure('Name','Meshes')
clf
plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'node',true);
plot(S+ampl*mean_u,'Color','b','FaceColor','b','FaceAlpha',0.1,'node',true);
mysaveas(pathname,'meshes_deflected',formats,renderer);

% plotFacets(S);
% plotRidges(S);

%% Display solution

% ampl = 0;
ampl = max(getsize(S))/max(abs(mean_u))/10;
options = {'solid',true};
% options = {};

plotSolution(S,mean_u,'displ',3,'ampl',ampl,options{:});
mysaveas(pathname,'mean_Uz',formats,renderer);

plotSolution(S,std_u,'displ',3,'ampl',ampl,options{:});
mysaveas(pathname,'std_Uz',formats,renderer);

% myparallel('stop');
