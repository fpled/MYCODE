%% Portal frame %%
%%--------------%%

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
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Domains and meshes
% Beams 1 and 2 (rounded dimensions in mm)
L1 = 750;
b1 = 400;
L2 = 1000;
b2 = 500;
h = 15;

% Points
P1 = POINT([0.0,0.0]);
P2 = POINT([0.0,L1]);
P_load = POINT([L2/2,L1]);
P3 = POINT([L2,L1]);
P4 = POINT([L2,0.0]);

L_beam{1} = LIGNE(P1,P2);
L_beam{2} = LIGNE(P2,P_load);
L_beam{3} = LIGNE(P_load,P3);
L_beam{4} = LIGNE(P4,P3);

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
S1 = b1*h;
S2 = b2*h;
% Planar second moment of area (or Planar area moment of inertia)
IY1 = h*b1^3/12;
IY2 = h*b1^3/12;
IZ1 = b1*h^3/12;
IZ2 = b2*h^3/12;
IX1 = IY1+IZ1;
IX2 = IY2+IZ2;

mat_beam{1} = ELAS_BEAM('E',E_beam,'NU',NU_beam,'S',S1,'IZ',IZ1,'IY',IY1,'IX',IX1,'RHO',RHO_beam);
mat_beam{1} = setnumber(mat_beam{1},1);
mat_beam{2} = ELAS_BEAM('E',E_beam,'NU',NU_beam,'S',S2,'IZ',IZ2,'IY',IY2,'IX',IX2,'RHO',RHO_beam);
mat_beam{2} = setnumber(mat_beam{2},2);
S_beam([1,4]) = cellfun(@(S) setmaterial(S,mat_beam{1}),S_beam([1,4]),'UniformOutput',false);
S_beam([2,3]) = cellfun(@(S) setmaterial(S,mat_beam{2}),S_beam([2,3]),'UniformOutput',false);

S = union(S_beam{:});

%% Neumann boundary conditions
p = 1;

%% Dirichlet boundary conditions
P_support = [P1 P4];
S = final(S);
S = addcl(S,P_support);

%% Stiffness matrix and sollicitation vector
A = calc_rigi(S);
f = nodalload(S,P_load,'FY',-p);

a_add = 0; % additonal junction rigidity
% [~,numnode2,~] = intersect(S,P2,'strict',false);
% [~,numnode3,~] = intersect(S,P3,'strict',false);
numddl2 = findddl(S,'RZ',P2,'free');
numddl3 = findddl(S,'RZ',P3,'free');
A(numddl2,numddl2) = A(numddl2,numddl2) + a_add;
A(numddl3,numddl3) = A(numddl3,numddl3) + a_add;

%% Solution
t = tic;
u = A\f;
time = toc(t);

e = calc_epsilon(S,u,'smooth');
s = calc_sigma(S,u,'smooth');

%% Test solution
ux = eval_sol(S,u,P2,'UX');
uy = eval_sol(S,u,P2,'UY');
rz = eval_sol(S,u,P2,'RZ');

N = s(1);
Mz = s(2);
N_max = max(abs(N));
Mz_max = max(abs(Mz));

[~,numnode,~] = intersect(S,P2,'strict',false);
N_corner = 0;
Mz_corner = 0;
for i=1:getnbgroupelem(S)
    Ni = reshape(abs(N{i}),[getnbnode(S),1]);
    Mzi = reshape(abs(Mz{i}),[getnbnode(S),1]);
    Ni_corner = double(Ni(numnode));
    Mzi_corner = double(Mzi(numnode));
    N_corner = max(N_corner,Ni_corner);
    Mz_corner = max(Mz_corner,Mzi_corner);
end
N = N_corner;
Mz = Mz_corner;

XA = 3/8*p*L2^2*IZ1/(L1*(L1*IZ2+2*L2*IZ1) + 3/S2*L2/L1*IZ1*(2*IZ2+L2/L1*IZ1));
YA = p/2;
MA = XA*(-L1/3+L2/L1^2*IZ1/S2);
% MA = p/8*L2^2*IZ1*(-L1+3*L2/L1^2*IZ1/S2)/(L1*(L1*IZ2+2*L2*IZ1) + 3/S2*L2/L1*IZ1*(2*IZ2+L2/L1*IZ1));

ux1_ex = @(x) -p/(2*E_beam*S1)*x;
ux2_ex = @(x) XA/(E_beam*S2)*(-x+L2/2);
rz1_ex = @(x) -1/(E_beam*IZ1)*(XA*x.^2/2+MA*x);
rz2_ex = @(x) 1/(E_beam*IZ2)*(p/4*(x+L2/2) - (XA*L1+MA)).*(x-L2/2);
% rz2_ex = @(x) 1/(E_beam*IZ2)*(p/4*x.^2 - (XA*L1+MA)*x) - 1/(E_beam*IZ1)*(XA*L1^2/2+MA*L1);
uy1_ex = @(x) -1/(E_beam*IZ1)*(XA*x.^3/6+MA*x.^2/2);
uy2_ex = @(x) 1/(E_beam*IZ2)*(p/2*x.^3/6 - (XA*L1+MA)*x.^2/2) - L1/(E_beam*IZ1)*(XA*L1/2+MA)*x - p/(2*E_beam*S1)*L1;

ux_ex = ux2_ex(0); % ux_ex = -uy1_ex(L1);
uy_ex = uy2_ex(0); % uy_ex = ux1_ex(L1);
rz_ex = rz2_ex(0); % rz_ex = rz1_ex(L1);
err_ux = norm(ux-ux_ex)/norm(ux_ex);
err_uy = norm(uy-uy_ex)/norm(uy_ex);
err_rz = norm(rz-rz_ex)/norm(rz_ex);

N1_ex = @(x) -YA;
N2_ex = @(x) -XA;
Ty1_ex = @(x) XA;
Ty2_ex = @(x) -YA;
Mz1_ex = @(x) -XA*x-MA;
Mz2_ex = @(x) p/2*x-XA*L1-MA;

[~,N1_max_ex] = fminbnd(@(x) -abs(N1_ex(x)),0,L1);
[~,N2_max_ex] = fminbnd(@(x) -abs(N2_ex(x)),0,L2/2);
[~,Mz1_max_ex] = fminbnd(@(x) -abs(Mz1_ex(x)),0,L1);
[~,Mz2_max_ex] = fminbnd(@(x) -abs(Mz2_ex(x)),0,L2/2);
N1_max_ex = -N1_max_ex;
N2_max_ex = -N2_max_ex;
Mz1_max_ex = -Mz1_max_ex;
Mz2_max_ex = -Mz2_max_ex;

N_ex = max(abs(N1_ex(L1)),abs(N2_ex(0)));
Mz_ex = max(abs(Mz1_ex(L1)),abs(Mz2_ex(0)));
err_N = norm(N-N_ex)/norm(N_ex);
err_Mz = norm(Mz-Mz_ex)/norm(Mz_ex);

N_max_ex = max(N1_max_ex,N2_max_ex);
Mz_max_ex = max(Mz1_max_ex,Mz2_max_ex);
err_N_max = norm(N_max-N_max_ex)/norm(N_max_ex);
err_Mz_max = norm(Mz_max-Mz_max_ex)/norm(Mz_max_ex);

%% Outputs
fprintf('\nCrossbar\n');
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

disp('Displacement u and rotation r at point'); disp(P2);
fprintf('ux    = %g mm\n',ux);
fprintf('ux_ex = %g mm, error = %g\n',ux_ex,err_ux);
fprintf('uy    = %g mm\n',uy);
fprintf('uy_ex = %g mm, error = %g\n',uy_ex,err_uy);
fprintf('rz    = %g rad = %g deg\n',rz,rad2deg(rz));
fprintf('rz_ex = %g rad = %g deg, error = %g\n',rz_ex,rad2deg(rz_ex),err_rz);
fprintf('\n');

disp('Force N and moment Mz at point'); disp(P2);
fprintf('N     = %g N\n',N);
fprintf('N_ex  = %g N, error = %g\n',N_ex,err_N);
fprintf('Mz    = %g N.mm\n',Mz);
fprintf('Mz_ex = %g N.mm, error = %g\n',Mz_ex,err_Mz);
fprintf('\n');

disp('Maximum force N and moment Mz');
fprintf('N_max     = %g N\n',N_max);
fprintf('N_max_ex  = %g N, error = %g\n',N_max_ex,err_N_max);
fprintf('Mz_max    = %g N.mm\n',Mz_max);
fprintf('Mz_max_ex = %g N.mm, error = %g\n',Mz_max_ex,err_Mz_max);
fprintf('\n');

%% Display domains, boundary conditions and meshes
plotDomain(S,'legend',false);
mysaveas(pathname,'domain',formats,renderer);
mymatlab2tikz(pathname,'domain.tex');

[hD,legD] = plotBoundaryConditions(S,'legend',false);
ampl = 300;
[hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
hP = plot(P2,'g+');
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

% DO NOT WORK WITH BEAM ELEMENTS
% plotSolution(S,u,'epsilon',1,'ampl',ampl);
% mysaveas(pathname,'eps_x',formats,renderer);
% 
% plotSolution(S,u,'epsilon',2,'ampl',ampl);
% mysaveas(pathname,'gam_z',formats,renderer);
% 
% plotSolution(S,u,'sigma',1,'ampl',ampl);
% mysaveas(pathname,'eff_x',formats,renderer);
% 
% plotSolution(S,u,'sigma',2,'ampl',ampl);
% mysaveas(pathname,'mom_z',formats,renderer);

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


