%% Portique plan hyperstatique de degré 3 %%
%%----------------------------------------%%

% clc
clearvars
close all

filename = 'portalFrame';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','plate',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Domains and meshes
% Beam 1 and 2 (rounded dimensions mm)
L1 = 750;
b1 = 400;
L2 = 1000;
b2 = 500;
h = 15;

% Points
P1_beam1_left = POINT([0.0,0.0]);
P2_beam1_left = POINT([0.0,L1]); % the same point of P1_beam2
P_load = POINT([L2/2,L1]);
P1_beam1_right = POINT([L2,0.0]);
P2_beam1_right = POINT([L2,L1]); % the same point of P2_beam2

L_beam{1} = LIGNE(P1_beam1_left, P2_beam1_left); % L_beam1_left
L_beam{2} = LIGNE(P2_beam1_left, P_load);
L_beam{3} = LIGNE(P_load, P2_beam1_right); % L_beam2
L_beam{4} = LIGNE(P1_beam1_right, P2_beam1_right); % L_beam1_right

cl_beam = h/2;
S_beam = cellfun(@(L,n) build_model(L,'cl',cl_beam,'elemtype','BEAM','filename',fullfile(pathname,['gmsh_beam_' num2str(n) '_cl_' num2str(cl_beam)])),L_beam,num2cell(1:length(L_beam)),'UniformOutput',false);

%% Materials
% Gravitational acceleration
g = 9.81;
% Young modulus
E_beam = 1;
% Poisson ratio
NU_beam = 0.3;
% Density
RHO_beam = 1;
% Cross-section area
Sec_beam1 = b1*h;
Sec_beam2 = b2*h;
% Planar second moment of area (or Planar area moment of inertia)
IY1 = h*b1^3/12;
IY2 = h*b1^3/12;
IZ1 = b1*h^3/12;
IZ2 = b2*h^3/12;
IX1 = IY1+IZ1;
IX2 = IY2+IZ2;

mat1_beam = ELAS_BEAM('E',E_beam,'NU',NU_beam,'S',Sec_beam1,'IZ',IZ1,'IY',IY1,'IX',IX1,'RHO',RHO_beam);
mat1_beam = setnumber(mat1_beam,1);
mat2_beam = ELAS_BEAM('E',E_beam,'NU',NU_beam,'S',Sec_beam2,'IZ',IZ2,'IY',IY2,'IX',IX2,'RHO',RHO_beam);
mat2_beam = setnumber(mat2_beam,2);
S_beam([1,4]) = cellfun(@(S) setmaterial(S,mat1_beam),S_beam([1,4]),'UniformOutput',false);
S_beam([2,3]) = cellfun(@(S) setmaterial(S,mat2_beam),S_beam([2,3]),'UniformOutput',false);

S = union(S_beam{:});

%% Neumann boundary conditions
p = 1;

%% Dirichlet boundary conditions
P_support = [P1_beam1_left P1_beam1_right];
S = final(S);
S = addcl(S,P_support);

%% Stiffness matrix and sollicitation vector
A = calc_rigi(S);
f = nodalload(S,P_load,'FY',-p);

%% Solution
t = tic;
u = A\f;
time = toc(t);

u = unfreevector(S,u);

e = calc_epsilon(S,u,'smooth');
s = calc_sigma(S,u,'smooth');

%% Outputs
fprintf('\nCrossbar\n');
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

%% Display domains, boundary conditions and meshes
plotDomain(S,'legend',false);
mysaveas(pathname,'domain',formats,renderer);
mymatlab2tikz(pathname,'domain.tex');

[hD,legD] = plotBoundaryConditions(S,'legend',false);
ampl = 300;
[hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
hP = plot(P_load,'g+');
legend([hD,hN,hP],[legD,legN,'measure'],'Location','NorthEastOutside')
mysaveas(pathname,'boundary_conditions',formats,renderer);

plotModel(S,'Color','k','FaceColor','k','node',true,'legend',false);
mysaveas(pathname,'mesh',formats,renderer);

ampl = getsize(S)/max(abs(u))/5;
plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','node',true,'legend',false);
mysaveas(pathname,'mesh_deflected',formats,renderer);

figure('Name','Meshes')
clf
plot(S,'Color','k','FaceColor','k','node',true);
plot(S+ampl*u,'Color','b','FaceColor','b','node',true);
mysaveas(pathname,'meshes_deflected',formats,renderer);

ampl = 0;
% ampl = getsize(S)/max(abs(u))/5;

plotSolution(S,u,'displ',1,'ampl',ampl);
mysaveas(pathname,'u_x',formats,renderer);

plotSolution(S,u,'displ',2,'ampl',ampl);
mysaveas(pathname,'u_y',formats,renderer);

plotSolution(S,u,'displ',3,'ampl',ampl);
mysaveas(pathname,'r_z',formats,renderer)

plotSolution(S,u,'epsilon',1,'ampl',ampl);
mysaveas(pathname,'eps_x',formats,renderer);

plotSolution(S,u,'epsilon',2,'ampl',ampl);
mysaveas(pathname,'gam_z',formats,renderer);

plotSolution(S,u,'sigma',1,'ampl',ampl);
mysaveas(pathname,'F_x',formats,renderer);

plotSolution(S,u,'sigma',2,'ampl',ampl);
mysaveas(pathname,'M_z',formats,renderer);

u = unfreevector(S,u);

figure('Name','Solution eps_x')
clf
plot(e,S+ampl*u,'compo','EPSX')
colorbar
set(gca,'FontSize',fontsize)
mysaveas(pathname,'eps_x',formats,renderer);

figure('Name','Solution gam_z')
clf
plot(e,S+ampl*u,'compo','GAMZ')
colorbar
set(gca,'FontSize',fontsize)
mysaveas(pathname,'gam_z',formats,renderer);

figure('Name','Solution eff_x')
clf
plot(s,S+ampl*u,'compo','EFFX')
colorbar
set(gca,'FontSize',fontsize)
mysaveas(pathname,'eff_x',formats,renderer);

figure('Name','Solution mom_z')
clf
plot(s,S+ampl*u,'compo','MOMZ')
colorbar
set(gca,'FontSize',fontsize)
mysaveas(pathname,'mom_z',formats,renderer);
