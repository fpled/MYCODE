%% Plate rectangular deterministic linear elasticity %%
%%---------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');

%% Input data

boundaries = {'simply_supported'};
% boundaries = {'clamped'};
% boundaries = {'simply_supported','clamped'};
loadings = {'uniform'};
% loadings = {'concentrated'};
% loadings = {'uniform','concentrated'};
% elemtypes = {'DKT'};
% elemtypes = {'DKQ'};
elemtypes = {'DKT','DKQ'}; % Kirchhoff-Love (classical) plate theory
% meshtypes = {'structured'};
% meshtypes = {'unstructured'};
meshtypes = {'structured','unstructured'};
nbelems = 2.^(1:5);

formats = {'fig','epsc2'};
renderer = 'OpenGL';

for ib=1:length(boundaries)
    boundary = boundaries{ib};
    
for il=1:length(loadings)
    loading = loadings{il};
    filename = ['plate_rect_det_lin_elas_' boundary '_' loading];
    hcvUz = figure('Name','Evolution of error indicator for Uz w.r.t number of elements');
    clf
    hcvRx = figure('Name','Evolution of error indicator for Rx w.r.t number of elements');
    clf
    hcvRy = figure('Name','Evolution of error indicator for Ry w.r.t number of elements');
    clf
    htime = figure('Name','Evolution of CPU time w.r.t number of elements');
    clf
    leg = cell(1,length(elemtypes)*length(meshtypes));
    
for ie=1:length(elemtypes)
    elemtype = elemtypes{ie};
    
for im=1:length(meshtypes)
    meshtype = meshtypes{im};
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep,[elemtype '_' meshtype],filesep);
    if ~exist(pathname,'dir')
        mkdir(pathname);
    end
    iem = (ie-1)*length(meshtypes)+im;
    leg{iem} = [elemtype ' ' meshtype];

%% Domains and meshes

a = 1;
b = 1;
Q = QUADRANGLE([0.0,0.0,0.0],[a,0.0,0.0],[a,b,0.0],[0.0,b,0.0]);

P_load = getcenter(Q);
x_load = double(getcoord(P_load));

err_Uz = zeros(1,length(nbelems));
err_Rx = zeros(1,length(nbelems));
err_Ry = zeros(1,length(nbelems));
time = zeros(1,length(nbelems));
Nbelem = zeros(1,length(nbelems));
for i=1:length(nbelems)

switch meshtype
    case 'structured'
        nbelem = [nbelems(i) nbelems(i)];
        S = build_model(Q,'nbelem',nbelem,'elemtype',elemtype);
    case 'unstructured'
        cl = min(a,b)./nbelems(i);
        switch loading
            case 'uniform'
                S = build_model(Q,'cl',cl,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_rect_' elemtype  '_cl_' num2str(cl)]);
            case 'concentrated'
                S = build_model(Q,'cl',cl,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_rect_' elemtype  '_cl_' num2str(cl)],'points',x_load);
        end
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

% L{1} = LIGNE([0.0,0.0,0.0],[a,0.0,0.0]);
% L{2} = LIGNE([a,0.0,0.0],[a,b,0.0]);
% L{3} = LIGNE([a,b,0.0],[0.0,b,0.0]);
% L{4} = LIGNE([0.0,b,0.0],[0.0,0.0,0.0]);
L = getedges(Q);

S = final(S);
switch boundary
    case 'clamped'
        S = addcl(S,[]); % addcl(S,[],{'U','R'},0);
    case 'simply_supported'
        % Soft support
        S = addcl(S,[],'U'); % S = addcl(S,[],{'UX','UY','UZ'},0);
        % Hard support
%         S = addcl(S,L{1},{'U','RY'});
%         S = addcl(S,L{2},{'U','RX'});
%         S = addcl(S,L{3},{'U','RY'});
%         S = addcl(S,L{4},{'U','RX'});
end
% S = addcl(S,[],'R'); % S = addcl(S,[],{'RX','RY','RZ'},0);

%% Stiffness matrices and sollicitation vectors

switch loading
    case 'uniform' % Uniform transverse load per unit area applied on the plate surface
        p = RHO*g*h;
    case 'concentrated' % Concentrated transverse load applied at point P_load
        p = RHO*g*h*a*b;
end

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

%% Resolution

t = tic;
u = A\f;
time(i) = toc(t);

%% Outputs

u = unfreevector(S,u);

U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
Ux = u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,u,S.node,'UX')));
Uy = u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,u,S.node,'UY')));
Uz = u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,u,S.node,'UZ')));

R = u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
Rx = u(findddl(S,'RX'),:); % Rx = double(squeeze(eval_sol(S,u,S.node,'RX')));
Ry = u(findddl(S,'RY'),:); % Ry = double(squeeze(eval_sol(S,u,S.node,'RY')));
Rz = u(findddl(S,'RZ'),:); % Rz = double(squeeze(eval_sol(S,u,S.node,'RZ')));

switch boundary
    case 'clamped'
        % Galerkin approximation (based on weak formulation)
        n = 1:10;
        switch loading
            case 'uniform'
                A = - p/(4*D_rig*pi^4*(3*(1/a^4+1/b^4)+2/(a^2*b^2))) ./ n.^4;
            case 'concentrated'
                A = - p/(D_rig*pi^4*a*b*(3*(1/a^4+1/b^4)+2/(a^2*b^2))) ./ n.^4;
        end
        fun_Uz = @(x) dot(repmat(A,size(x,1),1) .* (1-cos(2*pi/a*x(:,1)*n)),1-cos(2*pi/b*x(:,2)*n),2);
        fun_Rx = @(x) dot(repmat(A,size(x,1),1) .* (1-cos(2*pi/a*x(:,1)*n)),2*pi/b * repmat(n,size(x,1),1) .* sin(2*pi/b*x(:,2)*n),2);
        fun_Ry = @(x) -dot(repmat(A,size(x,1),1) .* (1-cos(2*pi/b*x(:,2)*n)),2*pi/a * repmat(n,size(x,1),1) .* sin(2*pi/a*x(:,1)*n),2);
%         A = zeros(length(n),1);
%         fun_Uz = @(x) 0;
%         fun_Rx = @(x) 0;
%         fun_Ry = @(x) 0;
%         for k=1:length(n)
%             switch loading
%                 case 'uniform'
%                     A(k) = - p/(4*D_rig*pi^4*n(k)^4*(3*(1/a^4+1/b^4)+2/(a^2*b^2)));
%                 case 'concentrated'
%                     A(k) = - p/(D_rig*pi^4*n(k)^4*a*b*(3*(1/a^4+1/b^4)+2/(a^2*b^2)));
%             end
%             fun_Uz = @(x) fun_Uz(x) + A(k) .* (1-cos(2*n(k)*pi*x(:,1)/a)) .* (1-cos(2*n(k)*pi*x(:,2)/b));
%             fun_Rx = @(x) fun_Rx(x) + A(k) * 2*n(k)*pi/b .* (1-cos(2*n(k)*pi*x(:,1)/a)) .* sin(2*n(k)*pi*x(:,2)/b);
%             fun_Ry = @(x) fun_Ry(x) - A(k) * 2*n(k)*pi/a .* sin(2*n(k)*pi*x(:,1)/a) .* (1-cos(2*n(k)*pi*x(:,2)/b));
%         end
    case 'simply_supported'
        % Fourier series representation (based on strong formulation)
        m = 1:10;
        n = 1:10;
        switch loading
            case 'uniform'
                A = - 16*p/(D_rig*pi^6) ./ (m'*n .* (repmat(m,length(n),1)'.^2/a^2+repmat(n,length(m),1).^2/b^2).^2) .* (sin(m*pi/2)'.^2 * sin(n*pi/2).^2);
            case 'concentrated'
                A = - 4*p/(D_rig*pi^4*a*b) ./ (repmat(m,length(n),1)'.^2/a^2+repmat(n,length(m),1).^2/b^2).^2 .* (sin(m*pi/a*x_load(1))' * sin(n*pi/b*x_load(2)));
        end
        fun_Uz = @(x) dot(sin(pi/a*x(:,1)*m) * A,sin(pi/b*x(:,2)*n),2);
        fun_Rx = @(x) dot(sin(pi/a*x(:,1)*m) * A,pi/b * repmat(n,size(x,1),1) .* cos(pi/b*x(:,2)*n),2);
        fun_Ry = @(x) -dot(sin(pi/b*x(:,2)*n) * A',pi/a * repmat(m,size(x,1),1) .* cos(pi/a*x(:,1)*m),2);
%         A = zeros(length(m),length(n));
%         fun_Uz = @(x) 0;
%         fun_Rx = @(x) 0;
%         fun_Ry = @(x) 0;
%         for j=1:length(m)
%             for k=1:length(n)
%                 switch loading
%                     case 'uniform'
%                         A(j,k) = - 16*p/(D_rig*pi^6*m(j)*n(k)*(m(j)^2/a^2+n(k)^2/b^2)^2) * sin(m(j)*pi/2)^2 * sin(n(k)*pi/2)^2;
%                     case 'concentrated'
%                         A(j,k) = - 4*p/(D_rig*pi^4*a*b*(m(j)^2/a^2+n(k)^2/b^2)^2) * sin(m(j)*pi*x_load(1)/a) * sin(n(k)*pi*x_load(2)/b);
%                 end
%                 fun_Uz = @(x) fun_Uz(x) + A(j,k) .* sin(m(j)*pi*x(:,1)/a) .* sin(n(k)*pi*x(:,2)/b);
%                 fun_Rx = @(x) fun_Rx(x) + A(j,k) .* n(k)*pi/b .* sin(m(j)*pi*x(:,1)/a) .* cos(n(k)*pi*x(:,2)/b);
%                 fun_Ry = @(x) fun_Ry(x) - A(j,k) .* m(j)*pi/a .* cos(m(j)*pi*x(:,1)/a) .* sin(n(k)*pi*x(:,2)/b);
%             end
%         end
end
fun_Uz = MultiVariateFunction(fun_Uz,3);
fun_Rx = MultiVariateFunction(fun_Rx,3);
fun_Ry = MultiVariateFunction(fun_Ry,3);
fun_Uz.evaluationAtMultiplePoints = true;
fun_Rx.evaluationAtMultiplePoints = true;
fun_Ry.evaluationAtMultiplePoints = true;

x = getcoord(S.node);
Uz_ex = fun_Uz.functionEval(x);
Rx_ex = fun_Rx.functionEval(x);
Ry_ex = fun_Ry.functionEval(x);

err_Uz(i) = norm(Uz-Uz_ex)/norm(Uz_ex);
err_Rx(i) = norm(Rx-Rx_ex)/norm(Rx_ex);
err_Ry(i) = norm(Ry-Ry_ex)/norm(Ry_ex);

P = getcenter(Q);
xP = double(getcoord(P));

ux = eval_sol(S,u,P,'UX');
uy = eval_sol(S,u,P,'UY');
uz = eval_sol(S,u,P,'UZ');

rx = eval_sol(S,u,P,'RX');
ry = eval_sol(S,u,P,'RY');
rz = eval_sol(S,u,P,'RZ');

uz_ex = fun_Uz.functionEval(xP);
rx_ex = fun_Rx.functionEval(xP);
ry_ex = fun_Ry.functionEval(xP);

err_uz = norm(uz-uz_ex)/norm(uz_ex);
err_rx = norm(rx-rx_ex)/norm(rx_ex);
err_ry = norm(ry-ry_ex)/norm(ry_ex);

fprintf('\nRectangular plate\n');
fprintf(['Boundary : ' boundary '\n']);
fprintf(['Load     : ' loading '\n']);
fprintf(['Mesh     : ' elemtype ' ' meshtype ' elements\n']);
Nbelem(i) = getnbelem(S);
fprintf('Nb elements = %g\n',Nbelem(i));
fprintf('Span-to-thickness ratio = %g\n',max(a,b)/h);
fprintf('Error = %.3e for Uz\n',err_Uz(i));
fprintf('      = %.3e for Rx\n',err_Rx(i));
fprintf('      = %.3e for Ry\n',err_Ry(i));
fprintf('Elapsed time = %f s\n',time(i));
fprintf('\n');

disp('Displacement u at point'); disp(P);
fprintf('ux    = %g\n',ux);
fprintf('uy    = %g\n',uy);
fprintf('uz    = %g\n',uz);
fprintf('uz_ex = %g, error = %.3e\n',uz_ex,err_uz);
fprintf('\n');

disp('Rotation r at point'); disp(P);
fprintf('rx    = %g\n',rx);
fprintf('rx_ex = %g, error = %.3e\n',rx_ex,err_rx);
fprintf('ry    = %g\n',ry);
fprintf('ry_ex = %g, error = %.3e\n',ry_ex,err_ry);
fprintf('rz    = %g\n',rz);
fprintf('\n');

%% Save variables

save(fullfile(pathname,['solution_' num2str(i) '.mat']),'u','U','R');

%% Display domains, boundary conditions and meshes

% plotDomain(Q,'solid',true,'legend',false);
% mysaveas(pathname,['domain_' num2str(i)],formats,renderer);
% mymatlab2tikz(pathname,['domain_' num2str(i) '.tex']);
% 
% [hD,legD] = plotBoundaryConditions(S,'legend',false);
% switch loading
%     case 'uniform'
%         ampl = 2;
%     case 'concentrated'
%         ampl = 0.5;
% end
% [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
% % legend([hD,hN],'Dirichlet','Neumann')
% % legend([hD,hN],[legD,legN])
% mysaveas(pathname,['boundary_conditions_' num2str(i)],formats,renderer);
% 
% plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
% mysaveas(pathname,['mesh_' num2str(i)],formats,renderer);
% 
% ampl = getsize(S)/max(abs(u))/2;
% plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
% mysaveas(pathname,['mesh_deflected_' num2str(i)],formats,renderer);
% 
% figure('Name','Meshes')
% clf
% plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
% plot(S+ampl*u,'Color','b','FaceColor','b','FaceAlpha',0.1);
% mysaveas(pathname,['meshes_deflected_' num2str(i)],formats,renderer);

%% Display solution

% % ampl = 0;
% ampl = getsize(S)/max(abs(u))/2;
% options = {'solid',true};
% % options = {};
% 
% plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
% mysaveas(pathname,['Uz_' num2str(i)],formats,renderer);
% 
% figure('Name','Solution u_3_ex')
% clf
% plot(FENODEFIELD(Uz_ex),S+ampl*u,options{:});
% colorbar
% set(gca,'FontSize',16)
% mysaveas(pathname,['Uz_ex_' num2str(i)],formats,renderer);
% 
% plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
% mysaveas(pathname,['Rx_' num2str(i)],formats,renderer);
% 
% figure('Name','Solution r_1_ex')
% clf
% plot(FENODEFIELD(Rx_ex),S+ampl*u,options{:});
% colorbar
% set(gca,'FontSize',16)
% mysaveas(pathname,['Rx_ex_' num2str(i)],formats,renderer);
% 
% plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
% mysaveas(pathname,['Ry_' num2str(i)],formats,renderer);
% 
% figure('Name','Solution r_2_ex')
% clf
% plot(FENODEFIELD(Ry_ex),S+ampl*u,options{:});
% colorbar
% set(gca,'FontSize',16)
% mysaveas(pathname,['Ry_ex_' num2str(i)],formats,renderer);

end

figure(hcvUz)
loglog(Nbelem,err_Uz,'-','Color',getfacecolor(iem+1),'LineWidth',1);
hold on

figure(hcvRx)
loglog(Nbelem,err_Rx,'-','Color',getfacecolor(iem+1),'LineWidth',1);
hold on

figure(hcvRy)
loglog(Nbelem,err_Ry,'-','Color',getfacecolor(iem+1),'LineWidth',1);
hold on

figure(htime)
loglog(Nbelem,time,'-','Color',getfacecolor(iem+1),'LineWidth',1);
hold on

end
end

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);

figure(hcvUz)
grid on
box on
set(gca,'FontSize',16)
xlabel('Number of elements')
ylabel('Error')
legend(leg{:})
mysaveas(pathname,'error_Uz','fig');
mymatlab2tikz(pathname,'error_Uz.tex');

figure(hcvRx)
grid on
box on
set(gca,'FontSize',16)
xlabel('Number of elements')
ylabel('Error')
legend(leg{:})
mysaveas(pathname,'error_Rx','fig');
mymatlab2tikz(pathname,'error_Rx.tex');

figure(hcvRy)
grid on
box on
set(gca,'FontSize',16)
xlabel('Number of elements')
ylabel('Error')
legend(leg{:})
mysaveas(pathname,'error_Ry','fig');
mymatlab2tikz(pathname,'error_Ry.tex');

figure(htime)
grid on
box on
set(gca,'FontSize',16)
xlabel('Number of elements')
ylabel('CPU time (s)')
legend(leg{:})
mysaveas(pathname,'cputime','fig');
mymatlab2tikz(pathname,'cputime.tex');

end
end
