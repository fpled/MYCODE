%% FCBA table 3D deterministic linear elasticity %%
%%-----------------------------------------------%%

% clc
clearvars
close all

formats = {'fig','epsc'};
renderer = 'OpenGL';

filename = 'Table';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','plate',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Mesh
S = gmsh2femobject(3,'Table.msh');

%% Materials
% Material symmetry
materialSym = 'isotTrans';

switch lower(materialSym)
    case 'isot'
        % Young modulus
        E = 1e9; % Pa
        % Poisson ratio
        NU = 0.3;
        % Shear modulus
        G = E/(2*(1+NU)); % Pa
        % Material
        mat = ELAS_ISOT('E',E,'NU',NU);
        mat = setnumber(mat,1);
        S = setmaterial(S,mat);
    case 'isottrans'
        % Transverse Young modulus
        ET = 1e9; % Pa
        % Longitudinal shear modulus
        GL = 1e6; % Pa
        % Longitudinal Young modulus
        EL = 1e6; % Pa
        % Longitudinal Poisson ratio
        NUL = 0.3;
        % Transverse Poisson ratio
        NUT = 0.2;
        % Material
        mat1 = ELAS_ISOT_TRANS('AXISL',[1;0;0],'EL',EL,'ET',ET,'NUL',NUL,'NUT',NUT,'GL',GL);
        mat1 = setnumber(mat1,1);
        mat2 = ELAS_ISOT_TRANS('AXISL',[0;0;1],'EL',EL,'ET',ET,'NUL',NUL,'NUT',NUT,'GL',GL);
        mat2 = setnumber(mat2,2);
        S = setmaterial(S,mat1,1:4);
        S = setmaterial(S,mat2,5);
    otherwise
        error('Wrong material symmetry !')
end

%% Neumann boundary conditions
RHO = 1;
g = 9.81;
p_vol = RHO*g; % body load

r = 40e-3;
C = CIRCLE(500e-3,250e-3,620e-3,r);
Sec = pi*r^2;
p = 200;
p_surf = p/Sec; % surf load

%% Dirichlet boundary conditions
S = final(S);
P1 = POINT([0.0,0.0,0.0]);
P2 = POINT([0.0,1.0,0.0]);
P3 = POINT([1.0,0.0,0.0]);
P = PLAN(P1,P2,P3);
% S = addcl(S,P,{'UX','UY','UZ'},0);
S = addcl(S,P);
% [~,numnode] = intersect(S,P);
% S = addcl(S,numnode,{'UX','UY','UZ'},0);
% S = addcl(S,numnode,'U',0);
% S = addcl(S,numnode,'U');
% S = addcl(S,numnode);

%% Stiffness matrix and sollicitation vector
A = calc_rigi(S);

[S_surf,numnode_surf,numelem_surf] = intersect(S,C,'strict',0);
f = surfload(S,S_surf,'FZ',-p_surf);
f = f + bodyload(S,[],'FZ',-p_vol);

%% Solution
t = tic;
u = A\f;
time = toc(t);

%% Outputs
fprintf('\nTable\n');
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
ampl = 5;
[hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
mysaveas(pathname,'boundary_conditions',formats,renderer);

plotModel(S,'legend',false);
mysaveas(pathname,'mesh',formats,renderer);

ampl = getsize(S)/max(abs(u))/10;
plotModelDeflection(S,u,'ampl',ampl,'FaceColor','b','legend',false);
mysaveas(pathname,'mesh_deflected',formats,renderer);

figure('Name','Meshes')
clf
plot(S);
plot(S+ampl*u,'FaceColor','b');
mysaveas(pathname,'meshes_deflected',formats,renderer);

%% Display solution
% ampl = 0;
ampl = getsize(S)/max(abs(u))/10;
options = {'solid',true};
% options = {};

% Displacements
plotSolution(S,u,'displ',1,'ampl',ampl,options{:});
mysaveas(pathname,'Ux',formats,renderer);

plotSolution(S,u,'displ',2,'ampl',ampl,options{:});
mysaveas(pathname,'Uy',formats,renderer);

plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
mysaveas(pathname,'Uz',formats,renderer);

% Stresses
plotSolution(S,u,'sigma',1,'ampl',ampl,options{:});
mysaveas(pathname,'Sigxx',formats,renderer);

plotSolution(S,u,'sigma',2,'ampl',ampl,options{:});
mysaveas(pathname,'Sigyy',formats,renderer);

plotSolution(S,u,'sigma',3,'ampl',ampl,options{:});
mysaveas(pathname,'Sigzz',formats,renderer);

plotSolution(S,u,'sigma',4,'ampl',ampl,options{:});
mysaveas(pathname,'Sigyz',formats,renderer);

plotSolution(S,u,'sigma',5,'ampl',ampl,options{:});
mysaveas(pathname,'Sigxz',formats,renderer);

plotSolution(S,u,'sigma',6,'ampl',ampl,options{:});
mysaveas(pathname,'Sigxy',formats,renderer);

plotSolution(S,u,'sigma','mises','ampl',ampl,options{:});
mysaveas(pathname,'SigVM',formats,renderer);

