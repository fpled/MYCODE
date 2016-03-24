%% Plate/Shell deterministic linear elasticity disk %%
%%--------------------------------------------------%%

% clc
clear all
close all

%% Input data
filename = 'plate_shell_det_lin_elas_disk';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'RESULTS',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
% set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
renderer = 'OpenGL';

% Parallel computing
% myparallel('start');

%% Domain and mesh definition

r = 1;
C = CIRCLE(0.0,0.0,0.0,r);

elemtype = 'DKT'; % DKT, DKQ, COQ4
cl = 0.1;
system.S = build_model(C,'cl',cl,'elemtype',elemtype,'filename',[pathname 'gmsh_disk_' elemtype]);

%% Materials

% Gravitational acceleration
g = 10;
% Young modulus
E = 1;
% Poisson ratio
NU = 0.3;
% Thickness
H = 0.1;
% Density
RHO = 1;
% Extensional stiffness (or Membrane rigidity)
A = E*H/(1-NU^2);
% Bending stiffness (or Flexural rigidity)
D = E*H^3/(12*(1-NU^2));

% Material
% a(u,v) = int( epsilon(u) : K : epsilon(v) )
mat = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',H,'k',5/6);
mat = setnumber(mat,1);
system.S = setmaterial(system.S,mat);

%% Dirichlet boundary conditions

% bctype = 'clamped';
bctype = 'simply supported';

system.S = final(system.S);
switch bctype
    case 'clamped'
        system.S = addcl(system.S,[]); % addcl(system.S,[],{'U','R'});
    case 'simply supported'
        system.S = addcl(system.S,[],'U'); % system.S = addcl(system.S,[],{'UX','UY','UZ'});
end
% system.S = addcl(system.S,[],'R'); % system.S = addcl(system.S,[],{'RX','RY','RZ'});

%% Stiffness matrices and sollicitation vectors

% Uniform or Concentrated load p
forceload = 'uniform';
% forceload = 'concentrated';
switch forceload
    case 'uniform'
        p = RHO*g*H;
    case 'concentrated'
        p = RHO*g*H*r^2;
end
Pload = getcenter(C);
xload = double(getcoord(Pload));
% Moment per unit length c
c = 0;

% Stiffness matrix system.A and sollicitation vector system.b associated to mesh system.S
system.A = calc_rigi(system.S);
switch forceload
    case 'uniform'
        system.b = bodyload(system.S,[],'FZ',-p);
    case 'concentrated'
        system.b = nodalload(system.S,Pload,'FZ',-p);
        if isempty(ispointin(Pload,POINT(system.S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
end
if strcmp(bctype,'simply supported')
    system.b = system.b + surfload(system.S,[],{'MX','MY'},-c*[1;1]);
end

%% Resolution

u = solve_system(system);

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
        phi = 16/5*(H/r)^2/(1-NU);
end
switch forceload
    case 'uniform'
        switch bctype
            case 'clamped'
                w = @(x) -p/(64*D) * (r^2 - (x(:,1).^2+x(:,2).^2)).*(r^2 - (x(:,1).^2+x(:,2).^2) + phi);
            case 'simply supported'
                w = @(x) -1/(2*D*(1+NU)) * (r^2 - (x(:,1).^2+x(:,2).^2)) .* (p/32*((5+NU)*r^2 - (1+NU)*(x(:,1).^2+x(:,2).^2) + phi*(1+NU)) + c);
        end
    case 'concentrated'
        switch bctype
            case 'clamped'
                w = @(x) -p/(16*pi*D) * (r^2 - (x(:,1).^2+x(:,2).^2) - 2*(x(:,1).^2+x(:,2).^2).*log(r/sqrt(x(:,1).^2+x(:,2).^2)));
            case 'simply supported'
                w = @(x) -p/(16*pi*D) * ((3+NU)/(1+NU)*(r^2 - (x(:,1).^2+x(:,2).^2)) - 2*(x(:,1).^2+x(:,2).^2).*log(r/sqrt(x(:,1).^2+x(:,2).^2))) - c/(2*D*(1+NU))*(r^2 - (x(:,1).^2+x(:,2).^2));
        end
end
x = getcoord(system.S.node);
Uz_ex = w(x);
error_Uz = norm(Uz-Uz_ex)/norm(Uz_ex);
fprintf('\n');
fprintf('error = %.4e\n',error_Uz);
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
fprintf('ux    = %.4e\n',ux);
fprintf('uy    = %.4e\n',uy);
fprintf('uz    = %.4e\n',uz);
fprintf('uz_ex = %.4e\n',uz_ex);
fprintf('error = %.4e\n',error_uz);
fprintf('\n');

disp('Rotation r at point'); disp(P);
fprintf('rx    = %.4e\n',rx);
fprintf('ry    = %.4e\n',ry);
fprintf('rz    = %.4e\n',rz);
fprintf('\n');

%% Save variables

save(fullfile(pathname,'solution.mat'),'u','U','R');
save(fullfile(pathname,'all.mat'));

%% Display domain, partition and mesh

% Display domain
plot_domain(C,'solid');
mysaveas(pathname,'domain',{'fig','epsc2'},renderer);
mymatlab2tikz(pathname,'domain.tex');

% Display mesh system.S
plot_model(system.S,'color','k','facecolor','k','facealpha',0.1,'nolegend');
mysaveas(pathname,'mesh',{'fig','epsc2'},renderer);

% Display deflection (or deflected shape) of mesh system.S
ampl = max(getsize(system.S))/max(abs(u));
plot_model_deflection(system.S,u,'ampl',ampl,'color','b','facecolor','b','facealpha',0.3,'nolegend');
mysaveas(pathname,'mesh_deflected',{'fig','epsc2'},renderer);

figure('Name','Meshes')
clf
plot(system.S,'color','k','facecolor','k','facealpha',0.1);
plot(system.S+ampl*u,'color','b','facecolor','b','facealpha',0.3);
mysaveas(pathname,'meshes_deflected',{'fig','epsc2'},renderer);

% Display boundary conditions
[hD,legD] = plot_boundary_conditions(system.S,'nolegend');
switch forceload
    case 'uniform'
        ampl = 2;
    case 'concentrated'
        ampl = 1/(2*system.S.nbnode);
end
[hN,legN] = vectorplot(system.S,'F',system.b,ampl,'r');
% legend([hD,hN],'Dirichlet','Neumann')
legend([hD,hN],[legD,legN])
axis image
mysaveas(pathname,'boundary_conditions',{'fig','epsc2'},renderer);

% Display facets of mesh system.S
% plot_facets(system.S);
% mysaveas(pathname,'facets',{'fig','epsc2'},renderer);

% Display ridges of mesh system.S
% plot_ridges(system.S);
% mysaveas(pathname,'ridges',{'fig','epsc2'},renderer);

%% Display solution u=(Ux,Uy,Uz,Rx,Ry,Rz)

ampl = max(getsize(system.S))/max(abs(u));

% Uz
% plot_solution(system.S,u,'displ',3,'solid');
% mysaveas(pathname,'Uz',{'fig','epsc2'},renderer);

plot_solution(system.S,u,'displ',3,'ampl',ampl,'solid');
mysaveas(pathname,'Uz_deflected',{'fig','epsc2'},renderer);

% Uz_ex
% figure('Name','Solution u_3_ex')
% clf
% plot(FENODEFIELD(w(x)),system.S,'solid');
% colorbar
% set(gca,'FontSize',16)
% mysaveas(pathname,'Uz_ex',{'fig','epsc2'},renderer);

figure('Name','Solution u_3_ex')
clf
plot(FENODEFIELD(w(x)),system.S+ampl*u,'solid');
colorbar
set(gca,'FontSize',16)
mysaveas(pathname,'Uz_ex_deflected',{'fig','epsc2'},renderer);

% Uz Uz_ex
% figure('Name','Nodal values')
% clf
% plot3(x(:,1),x(:,2),Uz,'b.');
% plot3(x(:,1),x(:,2),w(x),'r.');

if exist('momentdof','var')
    % Rx
%     plot_solution(system.S,u,'rotation',1,'solid');
%     mysaveas(pathname,'Rx',{'fig','epsc2'},renderer);
    
    plot_solution(system.S,u,'rotation',1,'ampl',ampl,'solid');
    mysaveas(pathname,'Rx_deflected',{'fig','epsc2'},renderer);
    
    % Ry
%     plot_solution(system.S,u,'rotation',2,'solid');
%     mysaveas(pathname,'Ry',{'fig','epsc2'},renderer);
    
    plot_solution(system.S,u,'rotation',2,'ampl',ampl,'solid');
    mysaveas(pathname,'Ry_deflected',{'fig','epsc2'},renderer);
end

% myparallel('stop');
