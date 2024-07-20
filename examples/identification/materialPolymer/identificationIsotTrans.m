%% Identification of transversely isotropic elastic properties %%
%%-------------------------------------------------------------%%

% clc
clearvars
close all
% rng('default');

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

% structure = 'plate';
structure = 'plate_hole';

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialPolymer',structure);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Domains and meshes
D = DOMAIN(2,[36.3462,66.6026],[89.6795,83.0128]); % [mm]
C = CIRCLE(36.3462,66.6026,10); % [mm]

elemtype = 'TRI3';
% elemtype = 'QUA4';
% option = 'DEFO'; % plane strain
option = 'CONT'; % plane stress
nbelem = [14,4];
cl = min(getsize(D)./nbelem);

if strcmp(structure,'plate')
    S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'option',option);
    % S = build_model(D,'cl',cl,'elemtype',elemtype,'option',option,'filename',fullfile(pathname,'gmsh_plate'));
elseif strcmp(structure,'plate_hole')
    G = GMSHFILE();
    G = setfile(G,fullfile(pathname,'gmsh_plate_hole'));
    PD = getvertices(D);
    PC = getvertices(C);
    G = createpoint(G,getcenter(C),cl,1);
    G = createpoints(G,PC(3:4),cl,2:3);
    G = createpoints(G,PD(2:4),cl,4:6);
    G = createcircle(G,1,2:3,1);
    G = createlines(G,[[2 4];[4 5];[5 6];[6 3]],2:5);
    G = createlineloop(G,[-1 2:5],1);
    G = createplanesurface(G,1,1);
    S = gmsh2femobject(2,G,2);
    S = convertelem(S,elemtype,'elemtype');
end

%% Materials
% Young modulus
EL_exp = 2e3; % [MPa]
ET_exp = 4e3; % [MPa]
% Poisson ratio
NUL_exp = 0.2;
% Shear modulus
GL_exp = 1e3; % [MPa]
% Thickness
DIM3 = 1; % [m]
% Density
RHO = 1; % [kg/m3]

% Material
mat = ELAS_ISOT_TRANS('AXISL',[0;1],'AXIST',[1;0],'EL',EL_exp,'ET',ET_exp,'NUL',NUL_exp,'GL',GL_exp,'RHO',RHO,'DIM3',DIM3);
mat = setnumber(mat,1);
S = setmaterial(S,mat);

%% Dirichlet boundary conditions for experimental solution
L = getedges(D);
S = final(S);

S_exp = addcl(S,L{4},'UX');
S_exp = addcl(S_exp,L{1},'UY');

%% Stiffness matrices and sollicitation vectors
A = calc_rigi(S_exp);
fx = 3.55; % [N/mm2]
fy = 1.775; % [N/mm2]
% b = surfload(S_exp,L{2},'FX',fx);
% b = b + surfload(S_exp,L{3},'FY',fy);
b = surfload(S_exp,L{3},'FY',fy);

%% Experimental solution
u_exp_in = A\b;
u_exp = unfreevector(S_exp,u_exp_in);

%% Dirichlet boundary conditions for numerical solution
I = create_boundary(S);
I = final(I);
P = calc_P(S_exp,I);
u_exp_b = P*u_exp;
S = addcl(S,[],'U',u_exp_b);
u_exp_in = freevector(S,u_exp);

%% Display experimental solution
[hD,legD] = plotBoundaryConditions(S_exp,'legend',false);
ampl = 0.5;
[hN,legN] = vectorplot(S_exp,'F',b,ampl,'r','LineWidth',1);
legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
mysaveas(pathname,'boundary_conditions_exp',formats,renderer);

plotModel(S_exp,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh_exp',formats,renderer);

ampl = getsize(S_exp)/max(abs(u_exp))/10;
plotModelDeflection(S_exp,u_exp,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh_deflected_exp',formats,renderer);

figure('Name','Meshes')
clf
plot(S_exp,'Color','k','FaceColor','k','FaceAlpha',0.1);
plot(S_exp+ampl*u_exp,'Color','b','FaceColor','b','FaceAlpha',0.1);
mysaveas(pathname,'meshes_deflected_exp',formats,renderer);

ampl = 0;
plotSolution(S_exp,u_exp,'displ',1,'ampl',ampl);
mysaveas(pathname,'Ux_exp',formats,renderer);

plotSolution(S_exp,u_exp,'displ',2,'ampl',ampl);
mysaveas(pathname,'Uy_exp',formats,renderer);

%% Identification
% initial guess
EL0 = 1e3; % longitudinal Young modulus [MPa]
NUL0 = 0.1; % longitudinal Poisson ratio
GL0 = 5e2; % longitudinal shear modulus [MPa]

disp('Initial parameters');
disp('------------------');
fprintf('EL  = %g MPa\n',EL0);
fprintf('NUL = %g\n',NUL0);
fprintf('GL  = %g MPa\n',GL0);

x0 = [EL0 NUL0 GL0];
lb = [0 0 0];
ub = [Inf 0.5 Inf];

optimFun = 'lsqnonlin';  % optimization function
% optimFun = 'fminsearch';
% optimFun = 'fminunc';
% optimFun = 'fmincon';

display = 'off';
% display = 'iter';
% display = 'iter-detailed';
% display = 'final';
% display = 'final-detailed';

tolX = 1e-14; % tolerance on the parameter value
tolFun = 1e-14; % tolerance on the function value

switch optimFun
    case {'lsqnonlin','fminunc','fmincon'}
        options  = optimoptions(optimFun,'Display',display,'TolX',tolX,'TolFun',tolFun);
    case 'fminsearch'
        options = optimset('Display',display,'TolX',tolX,'TolFun',tolFun);
    otherwise
        error(['Wrong optimization function' optimFun])
end

t = tic;
switch optimFun
    case 'lsqnonlin'
        fun = @(x) funlsqnonlinIsotTrans(x,@solveTractionIsotTrans,u_exp_in,S);
        [x,err,~,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
    case 'fminsearch'
        fun = @(x) funoptimIsotTrans(x,@solveTractionIsotTrans,u_exp_in,S);
        [x,err,exitflag,output] = fminsearch(fun,x0,options);
    case 'fminunc'
        fun = @(x) funoptimIsotTrans(x,@solveTractionIsotTrans,u_exp_in,S);
        [x,err,exitflag,output] = fminunc(fun,x0,options);
    case 'fmincon'
        fun = @(x) funoptimIsotTrans(x,@solveTractionIsotTrans,u_exp_in,S);
        [x,err,exitflag,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
end
toc(t)

EL = x(1); % [MPa]
NUL = x(2);
GL = x(3); % [MPa]
err = sqrt(err)./norm(u_exp_in);

disp('Optimal parameters');
disp('------------------');
fprintf('EL  = %g MPa\n',EL);
fprintf('NUL = %g\n',NUL);
fprintf('GL  = %g MPa\n',GL);
fprintf('err = %g\n',err);
% fprintf('exitflag = %g\n',exitflag);
% disp(output);

%% Numerical solution
[u_in,S] = solveTractionIsotTrans(x,S);
u = unfreevector(S,u_in);

%% Display numerical solution
ampl = 0.5;
v_exp = calc_init_dirichlet(S);
figure('Name','Imposed experimental displacement')
clf
h = plot(S,'FaceColor','w','LineWidth',0.5);
hold on
[hD,legD] = vectorplot(S,'U',v_exp,ampl,'r','LineWidth',1);
hold off
set(gca,'FontSize',fontsize)
hg = hggroup;
set([h(:),hD],'Parent',hg);
axis image
l = legend(hD,'{\boldmath$u$}$^{\mathrm{exp}}$','Location','NorthEastOutside');
% l = legend(hD,legD,'Location','NorthEastOutside');
set(l,'Interpreter',interpreter);
mysaveas(pathname,'boundary_conditions',formats,renderer);

plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh',formats,renderer);

ampl = getsize(S)/max(abs(u))/10;
plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh_deflected',formats,renderer);

figure('Name','Meshes')
clf
plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
plot(S+ampl*u,'Color','b','FaceColor','b','FaceAlpha',0.1);
mysaveas(pathname,'meshes_deflected',formats,renderer);

ampl = 0;
plotSolution(S,u,'displ',1,'ampl',ampl);
mysaveas(pathname,'Ux',formats,renderer);

plotSolution(S,u,'displ',2,'ampl',ampl);
mysaveas(pathname,'Uy',formats,renderer);

%% Test numerical solution
EL_series = linspace(EL*0.5,EL*1.5,50); % [MPa]
NUL_series = linspace(NUL*0.5,NUL*1.5,50);
GL_series = linspace(GL*0.5,GL*1.5,50); % [MPa]

% Plot error EL GL
err = zeros(length(EL_series),length(GL_series));
for m=1:length(EL_series)
    for n=1:length(GL_series)
        x = [EL_series(m) NUL GL_series(n)];
        u_in = solveTractionIsotTrans(x,S);
        err(m,n) = norm(u_exp_in - u_in);
    end
end
err = err./norm(u_exp_in);
[errmin,I] = min(err);
[errmin,c] = min(errmin);
r = I(c);

figure
surfc(GL_series,EL_series,err,'EdgeColor','none');
colorbar
hold on
scatter3(GL_series(c),EL_series(r),errmin,'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off
set(gca,'FontSize',fontsize)
% set(gca,'ZScale','log')
xlabel('$G_L$ [MPa]','Interpreter',interpreter)
ylabel('$E_L$ [MPa]','Interpreter',interpreter)
zlabel('Error','Interpreter',interpreter)
% zlabel('Erreur','Interpreter',interpreter)
mysaveas(pathname,'error_EL_GL_3D',formats,renderer);

figure
contourf(GL_series,EL_series,err,30);
colorbar
hold on
scatter(GL_series(c),EL_series(r),'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off
set(gca,'FontSize',fontsize)
% set(gca,'ZScale','log')
xlabel('$G_L$ [MPa]','Interpreter',interpreter)
ylabel('$E_L$ [MPa]','Interpreter',interpreter)
mysaveas(pathname,'error_EL_GL_2D',formats,renderer);

% Plot error EL NUL
err = zeros(length(EL_series),length(NUL_series));
for m=1:length(EL_series)
    for n=1:length(NUL_series)
        x = [EL_series(m) NUL_series(n) GL];
        u_in = solveTractionIsotTrans(x,S);
        err(m,n) = norm(u_exp_in - u_in);
    end
end
err = err./norm(u_exp_in);
[errmin,I] = min(err);
[errmin,c] = min(errmin);
r = I(c);

figure
surfc(NUL_series,EL_series,err,'EdgeColor','none');
colorbar
hold on
scatter3(NUL_series(c),EL_series(r),errmin,'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off
set(gca,'FontSize',fontsize)
% set(gca,'ZScale','log')
xlabel('$\nu_L$','Interpreter',interpreter)
ylabel('$E_L$ [MPa]','Interpreter',interpreter)
zlabel('Error','Interpreter',interpreter)
% zlabel('Erreur','Interpreter',interpreter)
mysaveas(pathname,'error_EL_NUL_3D',formats,renderer);

figure
contourf(NUL_series,EL_series,err,30);
colorbar
hold on
scatter(NUL_series(c),EL_series(r),'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off
set(gca,'FontSize',fontsize)
% set(gca,'ZScale','log')
xlabel('$\nu_L$','Interpreter',interpreter)
ylabel('$E_L$ [MPa]','Interpreter',interpreter)
mysaveas(pathname,'error_EL_NUL_2D',formats,renderer);

% Plot error NUL GL
err = zeros(length(NUL_series),length(GL_series));
for m=1:length(NUL_series)
    for n=1:length(GL_series)
        x = [EL NUL_series(m) GL_series(n)];
        u_in = solveTractionIsotTrans(x,S);
        err(m,n) = norm(u_exp_in - u_in);
    end
end
err = err./norm(u_exp_in);
[errmin,I] = min(err);
[errmin,c] = min(errmin);
r = I(c);

figure
surfc(GL_series,NUL_series,err,'EdgeColor','none');
colorbar
hold on
scatter3(GL_series(c),NUL_series(r),errmin,'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off
set(gca,'FontSize',fontsize)
% set(gca,'ZScale','log')
xlabel('$G_L$ [MPa]','Interpreter',interpreter)
ylabel('$\nu_L$','Interpreter',interpreter)
zlabel('Error','Interpreter',interpreter)
% zlabel('Erreur','Interpreter',interpreter)
mysaveas(pathname,'error_NUL_GL_3D',formats,renderer);

figure
contourf(GL_series,NUL_series,err,30);
colorbar
hold on
scatter(GL_series(c),NUL_series(r),'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off
set(gca,'FontSize',fontsize)
% set(gca,'ZScale','log')
xlabel('$G_L$ [MPa]','Interpreter',interpreter)
ylabel('$\nu_L$','Interpreter',interpreter)
mysaveas(pathname,'error_NUL_GL_2D',formats,renderer);
