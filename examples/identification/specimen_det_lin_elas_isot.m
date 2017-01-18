%% Specimen determinsitic linear elasticity isotropy %%
%%---------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');

%% Input data

filename = 'specimen_det_lin_elas_isot';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Domains and meshes

L = 1;
D = DOMAIN(2,[0.0,0.0],[L,L]);

% elemtype = 'TRI3';
elemtype = 'QUA4';
option = 'DEFO'; % plane strain
% option = 'CONT'; % plane stress
nbelem = [20,20];
S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'option',option);
% cl = 0.05;
% S = build_model(D,'cl',cl,'elemtype',elemtype,'option',option,'filename',[pathname 'gmsh_domain']);

%% Materials

% Poisson ratio
NU = 0.3;
% Thickness
DIM3 = 1;
% Density
RHO = 1;
% Young modulus
E = 1;

% Material
mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
mat = setnumber(mat,1);
S = setmaterial(S,mat);

%% Dirichlet boundary conditions

LU = LIGNE([0.0,L],[L,L]);
LL = LIGNE([0.0,0.0],[L,0.0]);

S = final(S);
S = addcl(S,LL);

loading = 'Dirichlet'; % Imposed displacement
% loading = 'Neumann'; % Traction force density
if strcmp(loading,'Dirichlet')
    udmax = 1e-2;
    S = addcl(S,LU,'UY',@(x) 4*udmax/L^2*x(:,1).*(x(:,1)-L));
end

%% Stiffness matrices and sollicitation vectors

switch loading
    case 'Neumann'
        A = calc_rigi(S);
        f = -1;
        b = surfload(S,LU,'FY',f);
    case 'Dirichlet' % Imposed displacement
        [A,b] = calc_rigi(S);
        b = -b;
end

%% Resolution

t = tic;
u = A\b;
time = toc(t);

e = calc_epsilon(S,u);
s = calc_sigma(S,u);

%% Outputs

fprintf('\nSquare specimen\n');
fprintf(['Load     : ' loading '\n']);
fprintf(['Mesh     : ' elemtype ' elements\n']);
fprintf('Nb elements = %g\n',getnbelem(S));
fprintf('Nb dofs     = %g\n',getnbddl(S));
fprintf('Elapsed time = %f s\n',time);
fprintf('\n');

%% Save variables

save(fullfile(pathname,'solution.mat'),'u');

%% Display domains, boundary conditions and meshes

plotDomain(D,'legend',false);
mysaveas(pathname,'domain',formats,renderer);
mymatlab2tikz(pathname,'domain.tex');

[hD,legD] = plotBoundaryConditions(S,'legend',false);
ampl = 0.5;
switch loading
    case 'Neumann'
        [hN,legN] = vectorplot(S,'F',b,ampl,'r','LineWidth',1);
    case 'Dirichlet'
        v = calc_init_dirichlet(S);
        [hN,legN] = vectorplot(S,'U',v,ampl,'r','LineWidth',1);
end
% legend([hD,hN],[legD,legN])
mysaveas(pathname,'boundary_conditions',formats,renderer);

% plotModel(S,'legend',false);
% mysaveas(pathname,'mesh',formats,renderer);

plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh',formats,renderer);

ampl = getsize(S)/max(abs(u))/5;
plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh_deflected',formats,renderer);

figure('Name','Meshes')
clf
plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
plot(S+ampl*unfreevector(S,u),'Color','b','FaceColor','b','FaceAlpha',0.1);
mysaveas(pathname,'meshes_deflected',formats,renderer);

%% Display solution

% ampl = 0;
ampl = getsize(S)/max(abs(u))/5;

% for i=1:2
%     plotSolution(S,u,'displ',i,'ampl',ampl);
%     mysaveas(pathname,['u_' num2str(i)],formats,renderer);
% end
% 
% for i=1:3
%     plotSolution(S,u,'epsilon',i,'ampl',ampl);
%     mysaveas(pathname,['eps_' num2str(i)],formats,renderer);
%     
%     plotSolution(S,u,'sigma',i,'ampl',ampl);
%     mysaveas(pathname,['sig_' num2str(i)],formats,renderer);
% end

% figure('Name','Solution eps_xx')
% clf
% plot(e,S+ampl*u,'compo','EPXX')
% colorbar
% set(gca,'FontSize',16)
% mysaveas(pathname,'eps_xx',formats,renderer);
% 
% figure('Name','Solution eps_yy')
% clf
% plot(e,S+ampl*u,'compo','EPYY')
% colorbar
% set(gca,'FontSize',16)
% mysaveas(pathname,'eps_yy',formats,renderer);
% 
% figure('Name','Solution eps_xy')
% clf
% plot(e,S+ampl*u,'compo','EPXY')
% colorbar
% set(gca,'FontSize',16)
% mysaveas(pathname,'eps_xy',formats,renderer);

% figure('Name','Solution sig_xx')
% clf
% plot(s,S+ampl*u,'compo','SMXX')
% colorbar
% set(gca,'FontSize',16)
% mysaveas(pathname,'sig_xx',formats,renderer);
% 
% figure('Name','Solution sig_yy')
% clf
% plot(s,S+ampl*u,'compo','SMYY')
% colorbar
% set(gca,'FontSize',16)
% mysaveas(pathname,'sig_yy',formats,renderer);
% 
% figure('Name','Solution sig_xy')
% clf
% plot(s,S+ampl*u,'compo','SMXY')
% colorbar
% set(gca,'FontSize',16)
% mysaveas(pathname,'sig_xy',formats,renderer);
