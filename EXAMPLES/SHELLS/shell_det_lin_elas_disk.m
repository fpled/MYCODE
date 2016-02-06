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
cl = 0.1;
elemtype = 'DKT'; % DKT, DKQ, COQ4, TRI3, QUA4
if strcmp(elemtype,'DKT') || strcmp(elemtype,'DKQ') || strcmp(elemtype,'COQ4')
    shell = 1;
    view3 = 'view3';
    system.S = build_model(C,'cl',cl,'elemtype',elemtype,'indim',3,'filename',[pathname 'gmsh_disk_' elemtype]);
else
    shell = 0;
    view3 = [];
    system.S = build_model(C,'cl',cl,'elemtype',elemtype,'indim',2,'filename',[pathname 'gmsh_disk_' elemtype]);
end

%% Materials

% Young modulus E
E = 11.5e9;
% Poisson ratio
NU = 0.23;
% Thickness
H = 40e-3;
% Density
RHO = 1;
% Constant
D = E*H^3/(12*(1-NU^2));

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
q = -500;
% Moment per unit length c
c = -10;

if shell
    forcedof = 'FZ'; % FX, FY, FZ
%     momentdof = {'MX','MY'}; % MX, MY, MZ
else
    forcedof = 'FX'; % FX, FY
end

% Stiffness matrix system.A and sollicitation vector system.b associated to mesh system.S
system.A = calc_rigi(system.S);
system.b = bodyload(system.S,[],forcedof,q);
if exist('momentdof','var')
    fun = @(x) c*[x(:,1);x(:,2)]/norm(x);
    system.b = system.b + surfload(system.S,[],momentdof,fun,c);
end

%% Resolution

u = solve_system(calc_system(system));

u = unfreevector(system.S,u);

syscoord = getsyscoord(system.S.node);
U = u(findddl(system.S,DDL(DDLVECT('U',syscoord))));
U = reshape(U,getindim(syscoord),size(U,1)/getindim(syscoord));
U = U(:);

Ux = u(findddl(system.S,'UX'),:);
Uy = u(findddl(system.S,'UY'),:);
if shell
    Uz = u(findddl(system.S,'UZ'),:);
    if ~exist('momentdof','var')
        w = @(x) -q/(64*D) * (r^2 - (x(:,1).^2+x(:,2).^2)) .* (x(:,1).^2+x(:,2).^2 - (5+NU)/(1+NU)*r^2);
    else
        w = @(x) -q/(64*D) * ((x(:,1).^2+x(:,2).^2).^2 - r^4) + ((x(:,1).^2+x(:,2).^2) - r^2) .* (q*r^2*(3+NU)/8) / (2*D*(1+NU));
    end
    P = POINT(system.S.node);
    p = double(squeeze(P))';
    Uz_ex = w(p);
    error_Uz = norm(Uz-Uz_ex)/norm(Uz_ex)
    
    R = u(findddl(system.S,DDL(DDLVECT('R',syscoord))));
    R = reshape(R,getindim(syscoord),size(R,1)/getindim(syscoord));
    R = R(:);
    
    Rx = u(findddl(system.S,'RX'),:);
    Ry = u(findddl(system.S,'RY'),:);
    Rz = u(findddl(system.S,'RZ'),:);
end

%% Outputs

if shell % getindim(system.S) == 3
    P = POINT([0.0,0.0,0.0]);
else % getindim(system.S) == 2
    P = POINT([0.0,0.0]);
end

ux = eval_sol(system.S,u,P,'UX')
uy = eval_sol(system.S,u,P,'UY')
if shell
    uz = eval_sol(system.S,u,P,'UZ')
    uz_ex = w(double(P))
    error_uz = norm(uz-uz_ex)/norm(uz_ex)
    
    rx = eval_sol(system.S,u,P,'RX')
    ry = eval_sol(system.S,u,P,'RY')
    rz = eval_sol(system.S,u,P,'RZ')
end

%% Save variables

save(fullfile(pathname,'solution.mat'),'u','U','R');
save(fullfile(pathname,'all.mat'));

%% Display domain, partition and mesh

% Display domain
% plot_domain(C);
% mysaveas(pathname,'domain',{'fig','epsc2'},renderer);
% mymatlab2tikz(pathname,'domain.tex');

% Display mesh system.S
plot_model(system.S,'color','k','facecolor','k','facealpha',0.3,'facelighting','gouraud',view3,'nolegend');
mysaveas(pathname,'mesh',{'fig','epsc2'},renderer);

% Display deflected shape of mesh system.S
ampl = max(getsize(system.S))/max(abs(u))/4;

plot_model_deflected(system.S,u,'ampl',ampl,'color','b','facecolor','b','facealpha',0.3,'facelighting','gouraud',view3,'nolegend');
mysaveas(pathname,'mesh_deflected',{'fig','epsc2'},renderer);

figure
clf
plot(system.S,'color','k','facecolor','k','facealpha',0.1,'facelighting','gouraud');
plot(system.S+ampl*u,'color','b','facecolor','b','facealpha',0.3,'facelighting','gouraud');
if shell
    view(3)
end
mysaveas(pathname,'meshes_deflected',{'fig','epsc2'},renderer);

% Display boundary conditions
plot_model(system.S,'color','k','facecolor','w','facealpha',0.3,'facelighting','gouraud',view3,'nolegend');
hold on
h1 = plot(getnode(create_boundary(system.S)),'b*');
h2 = vectorplot(system.S,'F',system.b,'r');
hold off
legend([h1(1),h2(1)],'Dirichlet','Neumann')
mysaveas(pathname,'boundary_conditions',{'fig','epsc2'},renderer);

% Display facets of mesh system.S
plot_facets(system.S,view3);
mysaveas(pathname,'facets',{'fig','epsc2'},renderer);

% Display ridges of mesh system.S
plot_ridges(system.S,view3);
mysaveas(pathname,'ridges',{'fig','epsc2'},renderer);

%% Display solution u

if shell % u=(Ux,Uy,Uz,Rx,Ry,Rz)
    % Uz
    plot_solution(system.S,u,'displ',3,'facealpha',1,view3);
    mysaveas(pathname,['U_' num2str(3)],{'fig','epsc2'},renderer);
    
    plot_solution(system.S,u,'displ',3,'ampl',ampl,'facealpha',1,view3);
    mysaveas(pathname,['U_' num2str(3) '_deflected'],{'fig','epsc2'},renderer);
    
    if exist('momentdof','var')
        % Rx
        plot_solution(system.S,u,'rotation',1,'facealpha',1,view3);
        mysaveas(pathname,['R_' num2str(1)],{'fig','epsc2'},renderer);
        
        plot_solution(system.S,u,'rotation',1,'ampl',ampl,'facealpha',1,view3);
        mysaveas(pathname,['R_' num2str(1) '_deflected'],{'fig','epsc2'},renderer);
        
        % Ry
        plot_solution(system.S,u,'rotation',2,'facealpha',1,view3);
        mysaveas(pathname,['R_' num2str(2)],{'fig','epsc2'},renderer);
        
        plot_solution(system.S,u,'rotation',2,'ampl',ampl,'facealpha',1,view3);
        mysaveas(pathname,['R_' num2str(2) '_deflected'],{'fig','epsc2'},renderer);
    end
else % u=(Ux,Uy)
    % Ux
    plot_solution(system.S,u,'displ',1,'facealpha',1,view3);
    mysaveas(pathname,['u_' num2str(1)],{'fig','epsc2'},renderer);
    
    plot_solution(system.S,u,'displ',1,'ampl',ampl,'facealpha',1,view3);
    mysaveas(pathname,['u_' num2str(1) '_deflected'],{'fig','epsc2'},renderer);
end

% myparallel('stop');
