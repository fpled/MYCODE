%% Table deterministic linear elasticity rectangular %%
%%---------------------------------------------------%%

% clc
clear all
close all

%% Input data
loading = 'uniform';
% loading = 'concentrated';

filename = ['table_rect_det_lin_elas_' loading];
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

Pload = getcenter(Q);
xload = double(getcoord(Pload));

Pbeam{1} = POINT([1/4*a,1/4*b]);
Pbeam{2} = POINT([3/4*a,1/4*b]);
Pbeam{3} = POINT([3/4*a,3/4*b]);
Pbeam{4} = POINT([1/4*a,3/4*b]);
xbeam = cellfun(@(P) double(getcoord(P)),Pbeam,'UniformOutput',false);

elemtype = 'DKT'; % DKT, DKQ, COQ4
% nbelem = [30,30];
% S_plate = build_model(Q,'nbelem',nbelem,'elemtype',elemtype);
cl = 0.05;
switch loading
    case 'uniform'
        S_plate = build_model(Q,'cl',cl,'elemtype',elemtype,'filename',[pathname 'gmsh_rect_' elemtype  '_cl_' num2str(cl)]);
    case 'concentrated'
        S_plate = build_model(Q,'points',xload,'cl',cl,'elemtype',elemtype,'filename',[pathname 'gmsh_rect_' elemtype  '_cl_' num2str(cl)]);
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

% Material mat_plate associated to plate
mat_plate = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',h,'k',5/6);
mat_plate = setnumber(mat_plate,1);
S_plate = setmaterial(S_plate,mat_plate);



d = 0.01;
IY = pi*d^4/2;
IZ = IY;

% Material mat_beam associated to beams
mat_beam = ELAS_BEAM('E',E,'NU',NU,'S',d*l,'IZ',IZ,'IY',IY,'IX',IY+IZ,'RHO',RHO);
mat_beam = setnumber(mat_beam,2);
S_beam = setmaterial(S_beam,mat_beam);

mat = MATERIALS();
mat{1} = mat_plate;
mat{2} = mat_beam;

system.S = union(S_plate,S_beam);

%% Dirichlet boundary conditions

L1 = LIGNE([0.0,0.0,0.0],[a,0.0,0.0]);
L2 = LIGNE([a,0.0,0.0],[a,b,0.0]);
L3 = LIGNE([a,b,0.0],[0.0,b,0.0]);
L4 = LIGNE([0.0,b,0.0],[0.0,0.0,0.0]);

system.S = final(system.S);
switch boundary
    case 'clamped'
        system.S = addcl(system.S,[]); % addcl(system.S,[],{'U','R'},0);
    case 'simply_supported'
        system.S = addcl(system.S,[],'U'); % system.S = addcl(system.S,[],{'UX','UY','UZ'});
%         system.S = addcl(system.S,L1,'RY');
%         system.S = addcl(system.S,L3,'RY');
%         system.S = addcl(system.S,L2,'RX');
%         system.S = addcl(system.S,L4,'RX');
end
% system.S = addcl(system.S,[],'R'); % system.S = addcl(system.S,[],{'RX','RY','RZ'});

%% Stiffness matrices and sollicitation vectors

% Uniform or Concentrated load p
switch loading
    case 'uniform'
        p = RHO*g*h;
    case 'concentrated'
        p = RHO*g*h*a*b;
end

% Stiffness matrix system.A and sollicitation vector system.b associated to mesh system.S
system.A = calc_rigi(system.S);
switch loading
    case 'uniform'
        system.b = bodyload(system.S,[],'FZ',-p);
    case 'concentrated'
        system.b = nodalload(system.S,Pload,'FZ',-p);
        if isempty(ispointin(Pload,POINT(system.S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
end

%% Resolution

t = tic;
u = solve_system(system);
time = toc(t);
fprintf(['\nRectangular ' boundary ' plate under ' loading ' load\n']);
fprintf('Span-to-thickness ratio = %.3e\n',max(a,b)/h);
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

w = @(x) 0;
m_max = 10;
n_max = 10;
switch loading
    case 'uniform'
        switch boundary
            case 'clamped'
%                 for n=1:n_max
%                     w = @(x) w(x) - p/(4*D*pi^4*n^4*(3*(1/a^4+1/b^4)+2/(a^2*b^2))) .* (1-cos(2*n*pi*x(:,1)/a)) .* (1-cos(2*n*pi*x(:,2)/b));
%                 end
                w = @(x) -p/(4*D*pi^4*(3*(1/a^4+1/b^4)+2/(a^2*b^2))) .* (1-cos(2*pi*x(:,1)/a)) .* (1-cos(2*pi*x(:,2)/b));
            case 'simply_supported'
                for m=1:m_max
                    for n=1:n_max
                        w = @(x) w(x) - 16*p/(D*pi^6*m*n*(m^2/a^2+n^2/b^2)^2) * sin(m*pi/2)^2 * sin(n*pi/2)^2 .* sin(m*pi*x(:,1)/a) .* sin(n*pi*x(:,2)/b);
                    end
                end
        end
    case 'concentrated'
        switch boundary
            case 'clamped'
%                 for n=1:n_max
%                     w = @(x) w(x) - p/(D*pi^4*n^4*a*b*(3*(1/a^4+1/b^4)+2/(a^2*b^2))) .* (1-cos(2*n*pi*x(:,1)/a)) .* (1-cos(2*n*pi*x(:,2)/b));
%                 end
                w = @(x) -p/(D*pi^4*a*b*(3*(1/a^4+1/b^4)+2/(a^2*b^2))) .* (1-cos(2*pi*x(:,1)/a)) .* (1-cos(2*pi*x(:,2)/b));
            case 'simply_supported'
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
fprintf('\n');
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
switch loading
    case 'uniform'
        ampl = 2;
    case 'concentrated'
        ampl = 0.5;
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

% Rx
% plot_solution(system.S,u,'rotation',1,'solid');
% mysaveas(pathname,'Rx',{'fig','epsc2'},renderer);

% plot_solution(system.S,u,'rotation',1,'ampl',ampl,'solid');
% mysaveas(pathname,'Rx_deflected',{'fig','epsc2'},renderer);

% Ry
% plot_solution(system.S,u,'rotation',2,'solid');
% mysaveas(pathname,'Ry',{'fig','epsc2'},renderer);

% plot_solution(system.S,u,'rotation',2,'ampl',ampl,'solid');
% mysaveas(pathname,'Ry_deflected',{'fig','epsc2'},renderer);

% myparallel('stop');
