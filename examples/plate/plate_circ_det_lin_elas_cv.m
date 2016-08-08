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
loadings={'uniform','concentrated'};
% elemtypes = {'DKT'};
% elemtypes = {'DKQ'};
% elemtypes = {'COQ4'};
% elemtypes = {'DKT','DKQ'};
elemtypes = {'DKT','DKQ','COQ4'};
nbelems = 2.^(1:6);

formats = {'fig','epsc2'};
renderer = 'OpenGL';

for ib=1:length(boundaries)
    boundary = boundaries{ib};
    
for il=1:length(loadings)
    loading = loadings{il};
    filename = ['plate_circ_det_lin_elas_' boundary '_' loading];
    hcv = figure('Name','Evolution of error indicator w.r.t number of elements');
    clf
    htime = figure('Name','Evolution of CPU time w.r.t number of elements');
    clf
    leg = cell(1,length(elemtypes));
    
for ie=1:length(elemtypes)
    elemtype = elemtypes{ie};
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep,elemtype,filesep);
    if ~exist(pathname,'dir')
        mkdir(pathname);
    end
    leg{ie} = elemtype;

%% Domains and meshes

r = 1;
C = CIRCLE(0.0,0.0,0.0,r);

P_load = getcenter(C);
x_load = double(getcoord(P_load));

err = zeros(1,length(nbelems));
time = zeros(1,length(nbelems));
Nbelem = zeros(1,length(nbelems));
for i=1:length(nbelems)

cl = r./nbelems(i);
switch loading
    case 'uniform'
        problem.S = build_model(C,'cl',cl,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_circ_' elemtype '_cl_' num2str(cl)]);
    case 'concentrated'
        problem.S = build_model(C,'cl',cl,'elemtype',elemtype,'filename',[pathname 'gmsh_plate_circ_' elemtype '_cl_' num2str(cl)],'points',x_load);
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
problem.S = setmaterial(problem.S,mat);

%% Dirichlet boundary conditions

problem.S = final(problem.S);
switch boundary
    case 'clamped'
        problem.S = addcl(problem.S,[]); % addcl(problem.S,[],{'U','R'},0);
    case 'simply_supported'
        problem.S = addcl(problem.S,[],'U'); % problem.S = addcl(problem.S,[],{'UX','UY','UZ'},0);
end
% problem.S = addcl(problem.S,[],'R'); % problem.S = addcl(problem.S,[],{'RX','RY','RZ'},0);

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

problem.A = calc_rigi(problem.S);
switch loading
    case 'uniform'
        problem.b = bodyload(problem.S,[],'FZ',-p);
    case 'concentrated'
        problem.b = nodalload(problem.S,P_load,'FZ',-p);
        if isempty(ispointin(P_load,POINT(problem.S.node)))
            error('Pointwise load must be applied to a node of the mesh')
        end
end
if strcmp(boundary,'simply_supported')
    problem.b = problem.b + surfload(problem.S,[],{'MX','MY'},-c*[1;1]);
end

%% Resolution

t = tic;
u = solveSystem(problem);
time(i) = toc(t);

%% Outputs

u = unfreevector(problem.S,u);

U = u(findddl(problem.S,DDL(DDLVECT('U',problem.S.syscoord,'TRANS'))));
Ux = u(findddl(problem.S,'UX'),:); % Ux = double(squeeze(eval_sol(problem.S,u,problem.S.node,'UX')));
Uy = u(findddl(problem.S,'UY'),:); % Uy = double(squeeze(eval_sol(problem.S,u,problem.S.node,'UY')));
Uz = u(findddl(problem.S,'UZ'),:); % Uz = double(squeeze(eval_sol(problem.S,u,problem.S.node,'UZ')));

R = u(findddl(problem.S,DDL(DDLVECT('R',problem.S.syscoord,'ROTA'))));
Rx = u(findddl(problem.S,'RX'),:); % Rx = double(squeeze(eval_sol(problem.S,u,problem.S.node,'RX'))));
Ry = u(findddl(problem.S,'RY'),:); % Ry = double(squeeze(eval_sol(problem.S,u,problem.S.node,'RY'))));
Rz = u(findddl(problem.S,'RZ'),:); % Rz = double(squeeze(eval_sol(problem.S,u,problem.S.node,'RZ'))));

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
x = getcoord(problem.S.node);
Uz_ex = w(x);

ind = find(~isnan(Uz) & ~isnan(Uz_ex));
err(i) = norm(Uz(ind)-Uz_ex(ind))/norm(Uz_ex(ind));

P = getcenter(C);

ux = eval_sol(problem.S,u,P,'UX');
uy = eval_sol(problem.S,u,P,'UY');
uz = eval_sol(problem.S,u,P,'UZ');
uz_ex = w(double(P));
err_uz = norm(uz-uz_ex)/norm(uz_ex);

rx = eval_sol(problem.S,u,P,'RX');
ry = eval_sol(problem.S,u,P,'RY');
rz = eval_sol(problem.S,u,P,'RZ');

fprintf('\nCircular plate\n');
fprintf(['Boundary : ' boundary '\n']);
fprintf(['Load     : ' loading '\n']);
fprintf(['Mesh     : ' elemtype ' elements\n']);
Nbelem(i) = getnbelem(problem.S);
fprintf('Nb elements = %g\n',Nbelem(i));
fprintf('Span-to-thickness ratio = %g\n',r/h);
fprintf('Error = %g\n',err(i));
fprintf('Elapsed time = %f s\n',time(i));
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

plotModel(problem.S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,['mesh_' num2str(i)],formats,renderer);

end

figure(hcv)
loglog(Nbelem,err,'-','Color',getfacecolor(ie+1),'LineWidth',1);
hold on

figure(htime)
loglog(Nbelem,time,'-','Color',getfacecolor(ie+1),'LineWidth',1);
hold on

%% Save variables

save(fullfile(pathname,'solution.mat'),'u','U','R');

%% Display domains, boundary conditions and meshes

% plotDomain(C,'solid',true,'legend',false);
% mysaveas(pathname,'domain',formats,renderer);
% mymatlab2tikz(pathname,'domain.tex');
% 
% [hD,legD] = plotBoundaryConditions(problem.S,'legend',false);
% switch loading
%     case 'uniform'
%         ampl = 2;
%     case 'concentrated'
%         ampl = 1;
% end
% [hN,legN] = vectorplot(problem.S,'F',problem.b,ampl,'r');
% % legend([hD,hN],'Dirichlet','Neumann')
% % legend([hD,hN],[legD,legN])
% axis image
% mysaveas(pathname,'boundary_conditions',formats,renderer);
% 
% plotModel(problem.S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
% mysaveas(pathname,'mesh',formats,renderer);
% 
% ampl = max(getsize(problem.S))/max(abs(u));
% plotModelDeflection(problem.S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
% mysaveas(pathname,'mesh_deflected',formats,renderer);
% 
% figure('Name','Meshes')
% clf
% plot(problem.S,'Color','k','FaceColor','k','FaceAlpha',0.1);
% plot(problem.S+ampl*u,'Color','b','FaceColor','b','FaceAlpha',0.1);
% mysaveas(pathname,'meshes_deflected',formats,renderer);
% 
% % plotFacets(problem.S);
% % plotRidges(problem.S);
% 
% %% Display solution
% 
% % ampl = 0;
% ampl = max(getsize(problem.S))/max(abs(u));
% options = {'solid',true};
% % options = {};
% 
% plotSolution(problem.S,u,'displ',3,'ampl',ampl,options{:});
% mysaveas(pathname,'Uz',formats,renderer);
% 
% figure('Name','Solution u_3_ex')
% clf
% plot(FENODEFIELD(w(x)),problem.S+ampl*u,options{:});
% colorbar
% set(gca,'FontSize',16)
% mysaveas(pathname,'Uz_ex',formats,renderer);
% 
% % plotSolution(problem.S,u,'rotation',1,'ampl',ampl,options{:});
% % mysaveas(pathname,'Rx',formats,renderer);
% 
% % plotSolution(problem.S,u,'rotation',2,'ampl',ampl,options{:});
% % mysaveas(pathname,'Ry',formats,renderer);

end

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);

figure(hcv)
grid on
box on
set(gca,'FontSize',16)
xlabel('Number of elements')
ylabel('Error')
legend(leg{:})
mysaveas(pathname,'error','fig');
mymatlab2tikz(pathname,'error.tex');

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

% myparallel('stop');
