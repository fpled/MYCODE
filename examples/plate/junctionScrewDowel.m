%% Modeling of junction for screw and dowel %%
%%------------------------------------------%%

% clc
clearvars
close all

formats = {'fig','epsc'};
renderer = 'OpenGL';

filename = 'junctionScrewDowel';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','plate',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Domains and meshes
% Plates
a1 = 142.5e-3; % m
b = 113e-3;
a2 = 67.5e-3;
h = 15e-3;

Q1 = QUADRANGLE([0,0,0],[0,b,0],...
    [0,b,a1],[0,0,a1]);
Q2 = QUADRANGLE([0,b,a1],[0,0,a1],...
    [a2,0,a1],[a2,b,a1]);

elemtype = 'DKT';
cl = h;

S1 = build_model(Q1,'cl',cl,'elemtype',elemtype,...
    'filename',fullfile(pathname,['gmsh_junction_1_' elemtype '_cl_' num2str(cl)]));
%
S2 = build_model(Q2,'cl',cl,'elemtype',elemtype,...
    'filename',fullfile(pathname,['gmsh_junction_2_' elemtype '_cl_' num2str(cl)]));

%% Materials
% Gravitational acceleration
g = 9.81;

% Density
RHO = 707.1384;
Vol_total = h*(150e-3*113e-3+60e-3*113e-3);
Mass_total = Vol_total*RHO; % kg

% Data
filenameAna = 'data_ET_GL.mat';
filenameNum = 'data_EL_NUL.mat';
pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathnameIdentification,filenameAna));
load(fullfile(pathnameIdentification,filenameNum));

% Sample number
sample = 'B';
numSample = 13;

% Material symmetry
materialSym = 'isotTrans';

switch lower(materialSym)
    case 'isot'
        % Young modulus
        E = mean_ET_data(numSample)*1e6; % Pa
        %E = 1.7e9; % Pa
        % Shear modulus
        G = mean_GL_data(numSample)*1e6*13; % Pa
        % Poisson ratio
        NU = E./(2*G)-1;
        %NU = 0.25;
        % Material
        mat = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',h,'k',5/6);
    case 'isottrans'
        % Transverse Young modulus
        ET = mean_ET_data(numSample)*1e6; % Pa
        % Longitudinal shear modulus
        GL = mean_GL_data(numSample)*1e6; % Pa
        % Longitudinal Young modulus
        % EL = mean_EL_data(numSample)*1e6; % Pa
        % Longitudinal Poisson ratio
        % NUL = mean_NUL_data(numSample);
        % Transverse Poisson ratio
        NUT = 0.25;
        % Material
        mat = ELAS_SHELL_ISOT_TRANS('ET',ET,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',h,'k',5/6);
    otherwise
        error('Wrong material symmetry !')
end

mat = setnumber(mat,1);
S = union(S1,S2);
S = setmaterial(S,mat);

%% Neumann boundary conditions
p_plate = RHO*g*h; % surface load (body load for plates)
L_load = b;

junction_type = 'S1';

switch lower(junction_type)
    case 's1'
        p = [24 65 114 163 200 252 291];
    case 's2'
        p = [14 54 113 159 199 249 299];
    case 's3'
        p = [15 65 111 153 203 243 295 341];
    case 's4'
        p = [23 34 91 141 180 229 290];
    case 'd1'
        P = [40 46 101 146];
    case 'd2'
        P = [6 63 110 175];
end

%% Dirichlet boundary conditions
S = final(S);
L1 = getedge(Q1,1);
[~,numnode1] = intersect(S,L1);

S = addcl(S,numnode1);

%% Stiffness matrix and sollicitation vector
A = calc_rigi(S);

p = p(1)/L_load; % line load (surface load for plates)

Lf = LIGNE([a2,0,a1],[a2,b,a1]);

f = surfload(S,Lf,'FZ',-p);

f = f + bodyload(S,[],'FZ',-p_plate);

%% Solution
t = tic;
u = A\f;
time = toc(t);

u = unfreevector(S,u);

%% Display domains, boundary conditions and meshes
plotDomain(S,'legend',false);
mysaveas(pathname,'domain',formats,renderer);
mymatlab2tikz(pathname,'domain.tex');

[hD,legD] = plotBoundaryConditions(S,'legend',false);
ampl = 3;
[hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
mysaveas(pathname,'boundary_conditions',formats,renderer);

plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh',formats,renderer);

U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
ampl = getsize(S)/max(abs(U))/10;
plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh_deflected',formats,renderer);

figure('Name','Meshes')
clf
plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
plot(S+ampl*u,'Color','b','FaceColor','b','FaceAlpha',0.1);
mysaveas(pathname,'meshes_deflected',formats,renderer);

%% Display solution
% ampl = 0;
ampl = getsize(S)/max(abs(U))/10;
options = {'solid',true};
% options = {};

% plotSolution(S,u,'displ',1,'ampl',ampl,options{:});
% mysaveas(pathname,'Ux',formats,renderer);
% 
% plotSolution(S,u,'displ',2,'ampl',ampl,options{:});
% mysaveas(pathname,'Uy',formats,renderer);
% 
plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
mysaveas(pathname,'Uz',formats,renderer);

% plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
% mysaveas(pathname,'Rx',formats,renderer);
%
% plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
% mysaveas(pathname,'Ry',formats,renderer);

% plotSolution(S,u,'epsilon','mises','ampl',ampl,options{:});
% mysaveas(pathname,'EpsVM',formats,renderer);
% 
% plotSolution(S,u,'sigma','mises','ampl',ampl,options{:});
% mysaveas(pathname,'SigVM',formats,renderer);
