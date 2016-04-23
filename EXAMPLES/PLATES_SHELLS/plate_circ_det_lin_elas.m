%% Plate circular deterministic linear elasticity %%
%%------------------------------------------------%%
% Code_Aster v3.03.100.pdf
% SSLS100 - Plaque circulaire encastrée soumise à une pression uniforme
% Code_Aster v3.03.101.pdf
% SSLS101 - Plaque circulaire posée soumise à une pression uniforme

% clc
clear all
close all

%% Input data
% boundaries = {'simply_supported'};
% boundaries = {'clamped'};
boundaries = {'simply_supported','clamped'};
% loadings = {'uniform'};
% loadings = {'concentrated'};
loadings={'uniform','concentrated'};
% elemtypes = {'DKT'};
% elemtypes = {'DKQ'};
% elemtypes = {'COQ4'};
elemtypes = {'DKT','DKQ'};
% elemtypes = {'DKT','DKQ','COQ4'};

for indexb=1:length(boundaries)
    boundary = boundaries{indexb};
for indexl=1:length(loadings)
    loading = loadings{indexl};
for indexe=1:length(elemtypes)
    elemtype = elemtypes{indexe};

close all

filename = strcat('plate_circ_det_lin_elas_',boundary,'_',loading,'_',elemtype);
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'RESULTS',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
% set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
renderer = 'OpenGL';

% Parallel computing
% myparallel('start');

%% Domains and meshes

r = 1;
C = CIRCLE(0.0,0.0,0.0,r);

P_load = getcenter(C);
x_load = double(getcoord(P_load));

cl = 0.1;
switch loading
    case 'uniform'
        system.S = build_model(C,'cl',cl,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_circ_' elemtype '_cl_' num2str(cl)]);
    case 'concentrated'
        system.S = build_model(C,'cl',cl,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_circ_' elemtype '_cl_' num2str(cl)],'points',x_load);
end

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

% Material
mat = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',h,'k',5/6);
system.S = setmaterial(system.S,mat);

%% Dirichlet boundary conditions

system.S = final(system.S);
switch boundary
    case 'clamped'
        system.S = addcl(system.S,[]); % addcl(system.S,[],{'U','R'},0);
    case 'simply_supported'
        system.S = addcl(system.S,[],'U'); % system.S = addcl(system.S,[],{'UX','UY','UZ'},0);
end
% system.S = addcl(system.S,[],'R'); % system.S = addcl(system.S,[],{'RX','RY','RZ'},0);

%% Stiffness matrices and sollicitation vectors

% Uniform or Concentrated load
switch loading
    case 'uniform'
        p = RHO*g*h;
    case 'concentrated'
        p = RHO*g*h*r^2;
end
% Moment per unit length
c = 0;

system.A = calc_rigi(system.S);
switch loading
    case 'uniform'
        system.b = bodyload(system.S,[],'FZ',-p);
    case 'concentrated'
        system.b = nodalload(system.S,P_load,'FZ',-p);
        if isempty(ispointin(P_load,POINT(system.S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
end
if strcmp(boundary,'simply_supported')
    system.b = system.b + surfload(system.S,[],{'MX','MY'},-c*[1;1]);
end

%% Resolution

t = tic;
u = solve_system(system);
time = toc(t);
fprintf('\nCircular plate\n');
fprintf(['Boundary : ' boundary '\n']);
fprintf(['Load : ' loading '\n']);
fprintf(['Mesh : unstructured with ' elemtype ' elements\n']);
fprintf('Span-to-thickness ratio = %g\n',r/h);
fprintf('Elapsed time = %f s\n',time);

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

switch elemtype
    case {'DKT','DKQ'} % Kirchhoff-Love
        phi = 0;
    case {'COQ4'} % Reissner-Mindlin
        phi = 16/5*(h/r)^2/(1-NU);
end
switch loading
    case 'uniform'
        switch boundary
            case 'clamped'
                w = @(x) -p/(64*D) * (r^2 - (x(:,1).^2+x(:,2).^2)).*(r^2 - (x(:,1).^2+x(:,2).^2) + phi);
            case 'simply_supported'
                w = @(x) -1/(2*D*(1+NU)) * (r^2 - (x(:,1).^2+x(:,2).^2)) .* (p/32*((5+NU)*r^2 - (1+NU)*(x(:,1).^2+x(:,2).^2) + phi*(1+NU)) + c);
        end
    case 'concentrated'
        switch boundary
            case 'clamped'
                w = @(x) -p/(16*pi*D) * (r^2 - (x(:,1).^2+x(:,2).^2) - 2*(x(:,1).^2+x(:,2).^2).*log(r./sqrt(x(:,1).^2+x(:,2).^2)));
            case 'simply_supported'
                w = @(x) -p/(16*pi*D) * ((3+NU)/(1+NU)*(r^2 - (x(:,1).^2+x(:,2).^2)) - 2*(x(:,1).^2+x(:,2).^2).*log(r./sqrt(x(:,1).^2+x(:,2).^2))) - c/(2*D*(1+NU))*(r^2 - (x(:,1).^2+x(:,2).^2));
        end
end
x = getcoord(system.S.node);
Uz_ex = w(x);
error_Uz = norm(Uz-Uz_ex)/norm(Uz_ex);
fprintf('\n');
fprintf('error = %g\n',error_Uz);
fprintf('\n');

P = getcenter(C);

ux = eval_sol(system.S,u,P,'UX');
uy = eval_sol(system.S,u,P,'UY');
uz = eval_sol(system.S,u,P,'UZ');
uz_ex = w(double(P));
error_uz = norm(uz-uz_ex)/norm(uz_ex);

rx = eval_sol(system.S,u,P,'RX');
ry = eval_sol(system.S,u,P,'RY');
rz = eval_sol(system.S,u,P,'RZ');

disp('Displacement u at point'); disp(P);
fprintf('ux    = %g\n',ux);
fprintf('uy    = %g\n',uy);
fprintf('uz    = %g\n',uz);
fprintf('uz_ex = %g\n',uz_ex);
fprintf('error = %g\n',error_uz);
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

plot_domain(C,'solid','nolegend');
mysaveas(pathname,'domain',{'fig','epsc2'},renderer);
mymatlab2tikz(pathname,'domain.tex');

[hD,legD] = plot_boundary_conditions(system.S,'nolegend');
switch loading
    case 'uniform'
        ampl = 2;
    case 'concentrated'
        ampl = 1;
end
[hN,legN] = vectorplot(system.S,'F',system.b,ampl,'r');
% legend([hD,hN],'Dirichlet','Neumann')
legend([hD,hN],[legD,legN])
axis image
mysaveas(pathname,'boundary_conditions',{'fig','epsc2'},renderer);

plot_model(system.S,'color','k','facecolor','k','facealpha',0.1,'nolegend');
mysaveas(pathname,'mesh',{'fig','epsc2'},renderer);

ampl = max(getsize(system.S))/max(abs(u));
plot_model_deflection(system.S,u,'ampl',ampl,'color','b','facecolor','b','facealpha',0.1,'nolegend');
mysaveas(pathname,'mesh_deflected',{'fig','epsc2'},renderer);

figure('Name','Meshes')
clf
plot(system.S,'color','k','facecolor','k','facealpha',0.1);
plot(system.S+ampl*u,'color','b','facecolor','b','facealpha',0.1);
mysaveas(pathname,'meshes_deflected',{'fig','epsc2'},renderer);

% plot_facets(system.S);
% plot_ridges(system.S);

%% Display solution

% ampl = 0;
ampl = max(getsize(system.S))/max(abs(u));
options = {'solid'};
% options = {};

plot_solution(system.S,u,'displ',3,'ampl',ampl,options{:});
mysaveas(pathname,'Uz',{'fig','epsc2'},renderer);

figure('Name','Solution u_3_ex')
clf
plot(FENODEFIELD(w(x)),system.S+ampl*u,options{:});
colorbar
set(gca,'FontSize',16)
mysaveas(pathname,'Uz_ex',{'fig','epsc2'},renderer);

% plot_solution(system.S,u,'rotation',1,'ampl',ampl,options{:});
% mysaveas(pathname,'Rx',{'fig','epsc2'},renderer);

% plot_solution(system.S,u,'rotation',2,'ampl',ampl,options{:});
% mysaveas(pathname,'Ry',{'fig','epsc2'},renderer);

end
end
end

% myparallel('stop');
