%% Shell deterministic linear elasticity rect %%
%%--------------------------------------------%%

% clc
% clear all
close all

%% Input data
filename = 'shell_det_lin_elas_rect';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'RESULTS',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
% set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
renderer = 'OpenGL';

% Parallel computing
% myparallel('start');

%% Domain and mesh definition

a = 1;
b = 1;
Q = QUADRANGLE([0.0,0.0,0.0],[a,0.0,0.0],[a,b,0.0],[0.0,b,0.0]);

elemtype = 'DKT'; % DKT, DKQ, COQ4
nbelem = [30,30];
system.S = build_model(Q,'nbelem',nbelem,'elemtype',elemtype);
% cl = 0.05;
% system.S = build_model(Q,'cl',cl,'elemtype',elemtype,'filename',[pathname 'gmsh_rect_' elemtype]);

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
RHO = 500/(g*H);
% Extensional stiffness (or In-plane plate rigidity)
A = E*H/(1-NU^2);
% Bending stiffness (or Flexural rigidity)
D = E*H^3/(12*(1-NU^2));

% Material
% a(u,v) = int( epsilon(u) : K : epsilon(v) )
mat = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',H);
mat = setnumber(mat,1);
system.S = setmaterial(system.S,mat);

%% Dirichlet boundary conditions

L1 = LIGNE([0.0,0.0,0.0],[a,0.0,0.0]);
L2 = LIGNE([a,0.0,0.0],[a,b,0.0]);
L3 = LIGNE([a,b,0.0],[0.0,b,0.0]);
L4 = LIGNE([0.0,b,0.0],[0.0,0.0,0.0]);

bctype = 'clamped';
% bctype = 'simply supported';

system.S = final(system.S);
switch bctype
    case 'clamped'
        system.S = addcl(system.S,[]); % addcl(system.S,[],{'U','R'},0);
    case 'simply supported'
        system.S = addcl(system.S,[],'U',0); % system.S = addcl(system.S,[],{'UX','UY','UZ'},0);
%         system.S = addcl(system.S,L1,'RY',0);
%         system.S = addcl(system.S,L3,'RY',0);
%         system.S = addcl(system.S,L2,'RX',0);
%         system.S = addcl(system.S,L4,'RX',0);
end
% system.S = addcl(system.S,[],'R',0); % system.S = addcl(system.S,[],{'RX','RY','RZ'},0);

%% Stiffness matrices and sollicitation vectors

% Uniform or Concentrated load p
% forceload = 'uniform';
forceload = 'concentrated';
switch forceload
    case 'uniform'
        p = RHO*g*H;
    case 'concentrated'
        p = RHO*g*H*a*b;
end
Pload = getcenter(Q);
xload = double(getcoord(Pload));

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

%% Resolution

u = solve_system(calc_system(system));

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

w = @(x) 0;
m_max = 10;
n_max = 10;
switch forceload
    case 'uniform'
        switch bctype
            case 'clamped'
%                 for n=1:n_max
%                     w = @(x) w(x) - p/(4*D*pi^4*n^4*(3*(1/a^4+1/b^4)+2/(a^2*b^2))) .* (1-cos(2*n*pi*x(:,1)/a)) .* (1-cos(2*n*pi*x(:,2)/b));
%                 end
                w = @(x) -p/(4*D*pi^4*(3*(1/a^4+1/b^4)+2/(a^2*b^2))) .* (1-cos(2*pi*x(:,1)/a)) .* (1-cos(2*pi*x(:,2)/b));
            case 'simply supported'
                for m=1:m_max
                    for n=1:n_max
                        w = @(x) w(x) - 16*p/(D*pi^6*m*n*(m^2/a^2+n^2/b^2)^2) * sin(m*pi/2)^2 * sin(n*pi/2)^2 .* sin(m*pi*x(:,1)/a) .* sin(n*pi*x(:,2)/b);
                    end
                end
        end
    case 'concentrated'
        switch bctype
            case 'clamped'
%                 for n=1:n_max
%                     w = @(x) w(x) - p/(D*pi^4*n^4*a*b*(3*(1/a^4+1/b^4)+2/(a^2*b^2))) .* (1-cos(2*n*pi*x(:,1)/a)) .* (1-cos(2*n*pi*x(:,2)/b));
%                 end
                w = @(x) -p/(D*pi^4*a*b*(3*(1/a^4+1/b^4)+2/(a^2*b^2))) .* (1-cos(2*pi*x(:,1)/a)) .* (1-cos(2*pi*x(:,2)/b));
            case 'simply supported'
                for m=1:m_max
                    for n=1:n_max
                        w = @(x) w(x) - 4*p/(D*pi^4*a*b*(m^2/a^2+n^2/b^2)^2) * sin(m*pi*xload(1)/a) * sin(n*pi*xload(2)/b) .* sin(m*pi*x(:,1)/a) .* sin(n*pi*x(:,2)/b);
                    end
                end
        end
end
x = getcoord(system.S.node);
Uz_ex = w(x);
error_Uz = norm(Uz-Uz_ex)/norm(Uz_ex);
fprintf('error = %.4e\n',error_Uz);
fprintf('\n');

P = getcenter(Q);

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

% P1 = POINT([0.0,0.0,0.0]);
% P2 = POINT([a,0.0,0.0]);
% P3 = POINT([a,b,0.0]);
% P4 = POINT([0.0,b,0.0]);
% P5 = POINT([a/2,0.0,0.0]);
% P6 = POINT([a,b/2,0.0]);
% P7 = POINT([a/2,b,0.0]);
% P8 = POINT([0.0,b/2,0.0]);
% 
% r1x = eval_sol(system.S,u,P1,'RX');
% r1y = eval_sol(system.S,u,P1,'RY');
% r1z = eval_sol(system.S,u,P1,'RZ');
% 
% r2x = eval_sol(system.S,u,P2,'RX');
% r2y = eval_sol(system.S,u,P2,'RY');
% r2z = eval_sol(system.S,u,P2,'RZ');
% 
% r3x = eval_sol(system.S,u,P3,'RX');
% r3y = eval_sol(system.S,u,P3,'RY');
% r3z = eval_sol(system.S,u,P3,'RZ');
% 
% r4x = eval_sol(system.S,u,P4,'RX');
% r4y = eval_sol(system.S,u,P4,'RY');
% r4z = eval_sol(system.S,u,P4,'RZ');
% 
% r5x = eval_sol(system.S,u,P5,'RX');
% r5y = eval_sol(system.S,u,P5,'RY');
% r5z = eval_sol(system.S,u,P5,'RZ');
% 
% r6x = eval_sol(system.S,u,P6,'RX');
% r6y = eval_sol(system.S,u,P6,'RY');
% r6z = eval_sol(system.S,u,P6,'RZ');
% 
% r7x = eval_sol(system.S,u,P7,'RX');
% r7y = eval_sol(system.S,u,P7,'RY');
% r7z = eval_sol(system.S,u,P7,'RZ');
% 
% r8x = eval_sol(system.S,u,P8,'RX');
% r8y = eval_sol(system.S,u,P8,'RY');
% r8z = eval_sol(system.S,u,P8,'RZ');
% 
% disp('Rotation r at point'); disp(P1);
% fprintf('r1x    = %.4e\n',r1x);
% fprintf('r1y    = %.4e\n',r1y);
% fprintf('r1z    = %.4e\n',r1z);
% fprintf('\n');
% 
% disp('Rotation r at point'); disp(P2);
% fprintf('r2x    = %.4e\n',r2x);
% fprintf('r2y    = %.4e\n',r2y);
% fprintf('r2z    = %.4e\n',r2z);
% fprintf('\n');
% 
% disp('Rotation r at point'); disp(P3);
% fprintf('r3x    = %.4e\n',r3x);
% fprintf('r3y    = %.4e\n',r3y);
% fprintf('r3z    = %.4e\n',r3z);
% fprintf('\n');
% 
% disp('Rotation r at point'); disp(P4);
% fprintf('r4x    = %.4e\n',r4x);
% fprintf('r4y    = %.4e\n',r4y);
% fprintf('r4z    = %.4e\n',r4z);
% fprintf('\n');
% 
% disp('Rotation r at point'); disp(P5);
% fprintf('r5x    = %.4e\n',r5x);
% fprintf('r5y    = %.4e\n',r5y);
% fprintf('r5z    = %.4e\n',r5z);
% fprintf('\n');
% 
% disp('Rotation r at point'); disp(P6);
% fprintf('r6x    = %.4e\n',r6x);
% fprintf('r6y    = %.4e\n',r6y);
% fprintf('r6z    = %.4e\n',r6z);
% fprintf('\n');
% 
% disp('Rotation r at point'); disp(P7);
% fprintf('r7x    = %.4e\n',r7x);
% fprintf('r7y    = %.4e\n',r7y);
% fprintf('r7z    = %.4e\n',r7z);
% fprintf('\n');
% 
% disp('Rotation r at point'); disp(P8);
% fprintf('r8x    = %.4e\n',r8x);
% fprintf('r8y    = %.4e\n',r8y);
% fprintf('r8z    = %.4e\n',r8z);
% fprintf('\n');

%% Save variables

save(fullfile(pathname,'solution.mat'),'u','U','R');
save(fullfile(pathname,'all.mat'));

%% Display domain, partition and mesh

% Display domain
plot_domain(Q,'solid');
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
