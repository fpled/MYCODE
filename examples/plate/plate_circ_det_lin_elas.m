%% Plate circular deterministic linear elasticity %%
%%------------------------------------------------%%
% Code_Aster v3.03.100.pdf
% SSLS100 - Plaque circulaire encastrée soumise à une pression uniforme
% Code_Aster v3.03.101.pdf
% SSLS101 - Plaque circulaire posée soumise à une pression uniforme

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% myparallel('start');

%% Input data

% boundaries = {'simply_supported'};
% boundaries = {'clamped'};
boundaries = {'simply_supported','clamped'};
% loadings = {'uniform'};
% loadings = {'concentrated'};
loadings = {'uniform','concentrated'};
% elemtypes = {'DKT'};
% elemtypes = {'DKQ'};
% elemtypes = {'COQ4'};
% elemtypes = {'DKT','DKQ'};
elemtypes = {'DKT','DKQ','COQ4'};

formats = {'fig','epsc2'};
renderer = 'OpenGL';

for ib=1:length(boundaries)
    boundary = boundaries{ib};
    
for il=1:length(loadings)
    loading = loadings{il};
    filename = ['plate_circ_det_lin_elas_' boundary '_' loading];
    close all
    
for ie=1:length(elemtypes)
    elemtype = elemtypes{ie};
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep,elemtype,filesep);
    if ~exist(pathname,'dir')
        mkdir(pathname);
    end

%% Domains and meshes

r = 1;
C = CIRCLE(0.0,0.0,0.0,r);

P_load = getcenter(C);
x_load = double(getcoord(P_load));

cl = r/10;
switch loading
    case 'uniform'
        S = build_model(C,'cl',cl,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_circ_' elemtype '_cl_' num2str(cl)]);
    case 'concentrated'
        S = build_model(C,'cl',cl,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_circ_' elemtype '_cl_' num2str(cl)],'points',x_load);
end

%% Materials

% Gravitational acceleration
g = 10;
% Young modulus
E = 1;
% Poisson ratio
NU = 0.3;
% Density
RHO = 1;
% Thickness
h = 0.1;
% Extensional stiffness (or Membrane rigidity)
A_rig = E*h/(1-NU^2);
% Bending stiffness (or Flexural rigidity)
D_rig = E*h^3/(12*(1-NU^2));

% Material
mat = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',h,'k',5/6);
S = setmaterial(S,mat);

%% Dirichlet boundary conditions

S = final(S);
switch boundary
    case 'clamped'
        S = addcl(S,[]); % addcl(S,[],{'U','R'},0);
    case 'simply_supported'
        S = addcl(S,[],'U'); % S = addcl(S,[],{'UX','UY','UZ'},0);
end
% S = addcl(S,[],'R'); % S = addcl(S,[],{'RX','RY','RZ'},0);

%% Stiffness matrices and sollicitation vectors

% Uniform or Concentrated load
switch loading
    case 'uniform'
        p = RHO*g*h;
    case 'concentrated'
        Sec = pi*r^2;
        p = RHO*g*h*Sec;
end
% Moment per unit length
c = 0;

A = calc_rigi(S);
switch loading
    case 'uniform'
        f = bodyload(S,[],'FZ',-p);
    case 'concentrated'
        f = nodalload(S,P_load,'FZ',-p);
        if isempty(ispointin(P_load,POINT(S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
end
if strcmp(boundary,'simply_supported')
    f = f + surfload(S,[],{'MX','MY'},-c*[1;1]);
end

%% Resolution

t = tic;
u = A\f;
time = toc(t);

%% Outputs

u = unfreevector(S,u);

U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
Ux = u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,u,S.node,'UX')));
Uy = u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,u,S.node,'UY')));
Uz = u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,u,S.node,'UZ')));

R = u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
Rx = u(findddl(S,'RX'),:); % Rx = double(squeeze(eval_sol(S,u,S.node,'RX'))));
Ry = u(findddl(S,'RY'),:); % Ry = double(squeeze(eval_sol(S,u,S.node,'RY'))));
Rz = u(findddl(S,'RZ'),:); % Rz = double(squeeze(eval_sol(S,u,S.node,'RZ'))));

switch loading
    case 'uniform'
        switch elemtype
            case {'DKT','DKQ'} % Kirchhoff-Love
                phi = 0;
            case {'COQ4'} % Reissner-Mindlin
                phi = 16/5*(h/r)^2/(1-NU);
        end
        switch boundary
            case 'clamped'
                w = @(x) -p/(64*D_rig) * (r^2 - (x(:,1).^2+x(:,2).^2)).*(r^2 - (x(:,1).^2+x(:,2).^2) + r^2*phi);
            case 'simply_supported'
                w = @(x) -1/(2*D_rig*(1+NU)) * (r^2 - (x(:,1).^2+x(:,2).^2)) .* (p/32*((5+NU)*r^2 - (1+NU)*(x(:,1).^2+x(:,2).^2) + (1+NU)*r^2*phi) + c);
        end
    case 'concentrated'
        switch boundary
            case 'clamped'
                w = @(x) -p/(16*pi*D_rig) * (r^2 - (x(:,1).^2+x(:,2).^2) - 2*(x(:,1).^2+x(:,2).^2).*log(r./sqrt(x(:,1).^2+x(:,2).^2)));
            case 'simply_supported'
                w = @(x) -p/(16*pi*D_rig) * ((3+NU)/(1+NU)*(r^2 - (x(:,1).^2+x(:,2).^2)) - 2*(x(:,1).^2+x(:,2).^2).*log(r./sqrt(x(:,1).^2+x(:,2).^2))) - c/(2*D_rig*(1+NU))*(r^2 - (x(:,1).^2+x(:,2).^2));
        end
end
x = getcoord(S.node);
Uz_ex = w(x);

ind = find(~isnan(Uz) & ~isnan(Uz_ex));
err = norm(Uz(ind)-Uz_ex(ind))/norm(Uz_ex(ind));

P = getcenter(C);

ux = eval_sol(S,u,P,'UX');
uy = eval_sol(S,u,P,'UY');
uz = eval_sol(S,u,P,'UZ');
uz_ex = w(double(P));
err_uz = norm(uz-uz_ex)/norm(uz_ex);

rx = eval_sol(S,u,P,'RX');
ry = eval_sol(S,u,P,'RY');
rz = eval_sol(S,u,P,'RZ');

fprintf('\nCircular plate\n');
fprintf(['Boundary : ' boundary '\n']);
fprintf(['Load     : ' loading '\n']);
fprintf(['Mesh     : ' elemtype ' elements\n']);
fprintf('Nb elements = %g\n',getnbelem(S));
fprintf('Span-to-thickness ratio = %g\n',r/h);
fprintf('Error = %g\n',err);
fprintf('Elapsed time = %f s\n',time);
fprintf('\n');

disp('Displacement u at point'); disp(P);
fprintf('ux    = %g\n',ux);
fprintf('uy    = %g\n',uy);
fprintf('uz    = %g\n',uz);
fprintf('uz_ex = %g\n',uz_ex);
fprintf('error = %g\n',err_uz);
fprintf('\n');

disp('Rotation r at point'); disp(P);
fprintf('rx    = %g\n',rx);
fprintf('ry    = %g\n',ry);
fprintf('rz    = %g\n',rz);
fprintf('\n');

%% Save variables

save(fullfile(pathname,'solution.mat'),'u','U','R');

%% Display domains, boundary conditions and meshes

plotDomain(C,'solid',true,'legend',false);
mysaveas(pathname,'domain',formats,renderer);
mymatlab2tikz(pathname,'domain.tex');

[hD,legD] = plotBoundaryConditions(S,'legend',false);
switch loading
    case 'uniform'
        ampl = 2;
    case 'concentrated'
        ampl = 0.2;
end
[hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
% legend([hD,hN],'Dirichlet','Neumann')
% legend([hD,hN],[legD,legN])
mysaveas(pathname,'boundary_conditions',formats,renderer);

plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh',formats,renderer);

ampl = getsize(S)/max(abs(u))/5;
plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh_deflected',formats,renderer);

figure('Name','Meshes')
clf
plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
plot(S+ampl*u,'Color','b','FaceColor','b','FaceAlpha',0.1);
mysaveas(pathname,'meshes_deflected',formats,renderer);

% plotFacets(S);
% plotRidges(S);

%% Display solution

% ampl = 0;
ampl = getsize(S)/max(abs(u))/5;
options = {'solid',true};
% options = {};

plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
mysaveas(pathname,'Uz',formats,renderer);

figure('Name','Solution u_3_ex')
clf
plot(FENODEFIELD(w(x)),S+ampl*u,options{:});
colorbar
set(gca,'FontSize',16)
mysaveas(pathname,'Uz_ex',formats,renderer);

% plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
% mysaveas(pathname,'Rx',formats,renderer);

% plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
% mysaveas(pathname,'Ry',formats,renderer);

end
end
end

% myparallel('stop');
