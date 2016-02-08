%% Shell deterministic linear elasticity disk %%
%%--------------------------------------------%%

% clc
% clear all
close all

%% Input data
filename = 'shell_det_lin_elas_disk';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'RESULTS',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
set(0,'DefaultFigureVisible','on'); % change the default figure properties of the MATLAB root object
renderer = 'OpenGL';

% Parallel computing
% myparallel('start');

%% Domain and mesh definition

r = 500e-3;
C = CIRCLE(0.0,0.0,r);
cl = 0.05;
elemtype = 'DKT'; % DKT, DKQ, COQ4, TRI3, QUA4
if strcmp(elemtype,'DKT') || strcmp(elemtype,'DKQ') || strcmp(elemtype,'COQ4')
    shell = 1;
    system.S = build_model(C,'cl',cl,'elemtype',elemtype,'indim',3,'filename',[pathname 'gmsh_disk_' elemtype]);
else
    shell = 0;
    system.S = build_model(C,'cl',cl,'elemtype',elemtype,'indim',2,'filename',[pathname 'gmsh_disk_' elemtype]);
end

%% Materials

% Gravitational acceleration
g = 9.81;
% Young modulus
E = 11.5e9;
% Poisson ratio
NU = 0.23;
% Thickness
H = 40e-3;
% Density
q = -500;
RHO = -q/(H*g);
% Constants
D = E*H^3/(12*(1-NU^2));
A = E*H/(1-NU^2);

% Material mat_out associated to outside subdomain
% a(u,v) = int( epsilon(u) : K : epsilon(v) )
if shell
    mat = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',H); % 3D shell or plate model
else
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',H); % 2D linear elasticity model
end
mat = setnumber(mat,1);
system.S = setmaterial(system.S,mat);

%% Dirichlet boundary conditions

system.S = final(system.S);
% system.S = addcl(system.S,[]); % addcl(system.S,[],{'U','R'},0);
system.S = addcl(system.S,[],'U',0); % system.S = addcl(system.S,[],{'UX','UY','UZ'},0);
% system.S = addcl(system.S,[],'R',0); % system.S = addcl(system.S,[],{'RX','RY','RZ'},0);

%% Stiffness matrices and sollicitation vectors

% Body force field q
q = -RHO*H*g;
% Moment per unit length c
c = -100;

if shell
    forcedof = 'FZ'; % FX, FY, FZ
    momentdof = {'MX','MY'}; % MX, MY, MZ
else
    forcedof = 'FX'; % FX, FY
end

% Stiffness matrix system.A and sollicitation vector system.b associated to mesh system.S
system.A = calc_rigi(system.S);
system.b = bodyload(system.S,[],forcedof,q);
if exist('momentdof','var')
    fun = @(x,c) c*[-x(:,2) x(:,1)]'/norm(x);
    system.b = system.b + surfload(system.S,[],momentdof,fun,c);
end

%% Resolution

u = solve_system(calc_system(system));

%% Outputs

u = unfreevector(system.S,u);

U = u(findddl(system.S,DDL(DDLVECT('U',system.S.syscoord,'TRANS'))));
Ux = u(findddl(system.S,'UX'),:); % Ux = double(squeeze(eval_sol(system.S,u,system.S.node,'UX')));
Uy = u(findddl(system.S,'UY'),:); % Uy = double(squeeze(eval_sol(system.S,u,system.S.node,'UY')));
if shell
    Uz = u(findddl(system.S,'UZ'),:); % Uz = double(squeeze(eval_sol(system.S,u,system.S.node,'UZ')));
    if ~exist('momentdof','var')
        w = @(x) q/(64*D) * ((x(:,1).^2+x(:,2).^2) - r^2) .* (x(:,1).^2+x(:,2).^2 - (5+NU)/(1+NU)*r^2);
    else
        w = @(x) q/(64*D) * ((x(:,1).^2+x(:,2).^2).^2 - r^4) - ((x(:,1).^2+x(:,2).^2) - r^2) .* (c + q*r^2*(3+NU)/8) / (2*D*(1+NU));
    end
    x = getcoord(system.S.node);
    Uz_ex = w(x);
    error_Uz = norm(Uz-Uz_ex)/norm(Uz_ex);
    fprintf('error = %.4e\n',error_Uz);
    fprintf('\n');
    
    R = u(findddl(system.S,DDL(DDLVECT('R',system.S.syscoord,'ROTA'))));
    Rx = u(findddl(system.S,'RX'),:); % Rx = double(squeeze(eval_sol(system.S,u,system.S.node,'RX'))));
    Ry = u(findddl(system.S,'RY'),:); % Ry = double(squeeze(eval_sol(system.S,u,system.S.node,'RY'))));
    Rz = u(findddl(system.S,'RZ'),:); % Rz = double(squeeze(eval_sol(system.S,u,system.S.node,'RZ'))));
end

if shell % getindim(system.S) == 3
    P = POINT([0.0,0.0,0.0]);
else % getindim(system.S) == 2
    P = POINT([0.0,0.0]);
end

disp('Displacement u at point');
disp(P);
ux = eval_sol(system.S,u,P,'UX');
uy = eval_sol(system.S,u,P,'UY');
fprintf('ux    = %.4e\n',ux);
fprintf('uy    = %.4e\n',uy);
if shell
    uz = eval_sol(system.S,u,P,'UZ');
    uz_ex = w(double(P));
    error_uz = norm(uz-uz_ex)/norm(uz_ex);
    fprintf('uz    = %.4e\n',uz);
    fprintf('uz_ex = %.4e\n',uz_ex);
    fprintf('error = %.4e\n',error_uz);
    fprintf('\n');
    
    disp('Rotation r at point');
    disp(P);
    rx = eval_sol(system.S,u,P,'RX');
    ry = eval_sol(system.S,u,P,'RY');
    rz = eval_sol(system.S,u,P,'RZ');
    fprintf('rx    = %.4e\n',rx);
    fprintf('ry    = %.4e\n',ry);
    fprintf('rz    = %.4e\n',rz);
    fprintf('\n');
end

%% Save variables

if shell
    save(fullfile(pathname,'solution.mat'),'u','U','R');
else
    save(fullfile(pathname,'solution.mat'),'u','U');
end
save(fullfile(pathname,'all.mat'));

%% Display domain, partition and mesh

% Display domain
% plot_domain(C);
% mysaveas(pathname,'domain',{'fig','epsc2'},renderer);
% mymatlab2tikz(pathname,'domain.tex');

% Display mesh system.S
plot_model(system.S,'color','k','facecolor','k','facealpha',0.3,'nolegend');
mysaveas(pathname,'mesh',{'fig','epsc2'},renderer);

% Display deflection (or deflected shape) of mesh system.S
ampl = max(getsize(system.S))/max(abs(u));

plot_model_deflection(system.S,u,'ampl',ampl,'color','b','facecolor','b','facealpha',0.3,'nolegend');
mysaveas(pathname,'mesh_deflected',{'fig','epsc2'},renderer);

figure('Name','Meshes')
clf
plot(system.S,'color','k','facecolor','k','facealpha',0.1);
plot(system.S+ampl*u,'color','b','facecolor','b','facealpha',0.3);
if shell
    view(3)
end
mysaveas(pathname,'meshes_deflected',{'fig','epsc2'},renderer);

% Display boundary conditions
plot_model(system.S,'color','k','facecolor','w','facealpha',0.3,'nolegend');
hold on
% hD = plot(getnode(create_boundary(system.S)),'b*');
[hD,leg] = plotbcond(system.S);
hN = vectorplot(system.S,'F',system.b,'r');
hold off
% legend([hD,hN],'Dirichlet','Neumann')
legend([hD,hN],leg{:},'Neumann')
mysaveas(pathname,'boundary_conditions',{'fig','epsc2'},renderer);

% Display facets of mesh system.S
plot_facets(system.S);
mysaveas(pathname,'facets',{'fig','epsc2'},renderer);

% Display ridges of mesh system.S
plot_ridges(system.S);
mysaveas(pathname,'ridges',{'fig','epsc2'},renderer);

%% Display solution u

if shell % u=(Ux,Uy,Uz,Rx,Ry,Rz)
    % Uz
    plot_solution(system.S,u,'displ',3,'solid');
    mysaveas(pathname,['U_' num2str(3)],{'fig','epsc2'},renderer);
    
    plot_solution(system.S,u,'displ',3,'ampl',ampl,'solid');
    mysaveas(pathname,['U_' num2str(3) '_deflected'],{'fig','epsc2'},renderer);
    
    figure('Name',['Solution u_ex_' num2str(3)])
    clf
    plot(FENODEFIELD(w(x)),system.S,'solid');
    mysaveas(pathname,['U_ex_' num2str(3)],{'fig','epsc2'},renderer);
    
    figure('Name',['Solution u_ex_' num2str(3)])
    clf
    plot(FENODEFIELD(w(x)),system.S+ampl*u,'solid');
    mysaveas(pathname,['U_ex_' num2str(3) '_deflected'],{'fig','epsc2'},renderer);
    
    figure
    clf
    plot3(x(:,1),x(:,2),Uz,'b.');
    plot3(x(:,1),x(:,2),w(x),'r.');
    
    
    if exist('momentdof','var')
        % Rx
        plot_solution(system.S,u,'rotation',1,'solid');
        mysaveas(pathname,['R_' num2str(1)],{'fig','epsc2'},renderer);
        
        plot_solution(system.S,u,'rotation',1,'ampl',ampl,'solid');
        mysaveas(pathname,['R_' num2str(1) '_deflected'],{'fig','epsc2'},renderer);
        
        % Ry
        plot_solution(system.S,u,'rotation',2,'solid');
        mysaveas(pathname,['R_' num2str(2)],{'fig','epsc2'},renderer);
        
        plot_solution(system.S,u,'rotation',2,'ampl',ampl,'solid');
        mysaveas(pathname,['R_' num2str(2) '_deflected'],{'fig','epsc2'},renderer);
    end
else % u=(Ux,Uy)
    % Ux
    plot_solution(system.S,u,'displ',1);
    mysaveas(pathname,['u_' num2str(1)],{'fig','epsc2'},renderer);
    
    plot_solution(system.S,u,'displ',1,'ampl',ampl);
    mysaveas(pathname,['u_' num2str(1) '_deflected'],{'fig','epsc2'},renderer);
end

% myparallel('stop');
