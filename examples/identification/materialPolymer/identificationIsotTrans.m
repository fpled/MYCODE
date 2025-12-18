%% Identification of transversely isotropic elastic properties %%
%%-------------------------------------------------------------%%

% clc
clearvars
close all
% rng('default');
myparallel('start');

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};

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
    G = createpoints(G,PC(1:2),cl,2:3);
    G = createpoints(G,PD(2:4),cl,4:6);
    G = createcirclearc(G,1,2:3,1);
    G = createlines(G,[[2 4];[4 5];[5 6];[6 3]],2:5);
    G = createcurveloop(G,[-1 2:5],1);
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
mysaveas(pathname,'boundary_conditions_exp',formats);

plotModel(S_exp,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh_exp',formats);

ampl = getsize(S_exp)/max(abs(u_exp))/10;
plotModelDeflection(S_exp,u_exp,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh_deflected_exp',formats);

figure('Name','Meshes')
clf
plot(S_exp,'Color','k','FaceColor','k','FaceAlpha',0.1);
plot(S_exp+ampl*u_exp,'Color','b','FaceColor','b','FaceAlpha',0.1);
mysaveas(pathname,'meshes_deflected_exp',formats);

ampl = 0;
plotSolution(S_exp,u_exp,'displ',1,'ampl',ampl);
mysaveas(pathname,'Ux_exp',formats);

plotSolution(S_exp,u_exp,'displ',2,'ampl',ampl);
mysaveas(pathname,'Uy_exp',formats);

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

optimFun = 'lsqnonlin'; % optimization function
% optimFun = 'fminsearch';
% optimFun = 'fminunc';
% optimFun = 'fmincon';

display = 'off';
% display = 'iter';
% display = 'iter-detailed';
% display = 'notify'; % only for fmincon and fminunc
% display = 'notify-detailed'; % only for fmincon and fminunc
% display = 'final';
% display = 'final-detailed';

algo = 'trust-region-reflective'; % default for lsqnonlin
% algo = 'interior-point'; % default for fmincon
% algo = 'active-set'; % only for fmincon
% algo = 'sqp'; % only for fmincon
% algo = 'sqp-legacy'; % only for fmincon
% algo = 'levenberg-marquardt'; % only for lsqnonlin
% algo = 'quasi-newton'; % default for fminunc
% algo = 'trust-region'; % only for fminunc

tolX = 1e-14; % tolerance on the parameter value
tolFun = 1e-14; % tolerance on the function value
maxIters = Inf; % maximum number of iterations
maxFunEvals = Inf; % maximum number of function evaluations

switch optimFun
    case {'lsqnonlin','fminunc','fmincon'}
        % options = optimoptions(optimFun,'Display',display,'Algorithm',algo,...
        %     'TolX',tolX,'TolFun',tolFun,'MaxIter',maxIters,'MaxFunEvals',maxFunEvals);
        options = optimoptions(optimFun,'Display',display,'Algorithm',algo,...
            'StepTolerance',tolX,'FunctionTolerance',tolFun,'OptimalityTolerance',tolFun,...
            'MaxIterations',maxIters,'MaxFunctionEvaluations',maxFunEvals);
    case 'fminsearch'
        options = optimset('Display',display,'TolX',tolX,'TolFun',tolFun,...
            'MaxIter',maxIters,'MaxFunEvals',maxFunEvals);
    otherwise
        error(['Wrong optimization function' optimFun])
end

t = tic;
switch optimFun
    case 'lsqnonlin'
        fun = @(x) funlsqnonlinIsotTrans(x,@solveTractionIsotTrans,u_exp_in,S);
        [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
    case 'fminsearch'
        fun = @(x) funoptimIsotTrans(x,@solveTractionIsotTrans,u_exp_in,S);
        [x,resnorm,exitflag,output] = fminsearch(fun,x0,options);
    case 'fminunc'
        fun = @(x) funoptimIsotTrans(x,@solveTractionIsotTrans,u_exp_in,S);
        [x,resnorm,exitflag,output] = fminunc(fun,x0,options);
    case 'fmincon'
        fun = @(x) funoptimIsotTrans(x,@solveTractionIsotTrans,u_exp_in,S);
        [x,resnorm,exitflag,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
end
toc(t)

EL = x(1); % [MPa]
NUL = x(2);
GL = x(3); % [MPa]
err = sqrt(resnorm)./norm(u_exp_in);

fprintf('\n');
disp('Optimal parameters');
disp('------------------');
fprintf('EL  = %g MPa\n',EL);
fprintf('NUL = %g\n',NUL);
fprintf('GL  = %g MPa\n',GL);
fprintf('resnorm = %g\n',resnorm);
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
legend(hD,'{\boldmath$u$}$^{\mathrm{exp}}$','Location','NorthEastOutside','Interpreter',interpreter)
% legend(hD,legD,'Location','NorthEastOutside','Interpreter',interpreter)
mysaveas(pathname,'boundary_conditions',formats);

plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh',formats);

ampl = getsize(S)/max(abs(u))/10;
plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
mysaveas(pathname,'mesh_deflected',formats);

figure('Name','Meshes')
clf
plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
plot(S+ampl*u,'Color','b','FaceColor','b','FaceAlpha',0.1);
mysaveas(pathname,'meshes_deflected',formats);

ampl = 0;
plotSolution(S,u,'displ',1,'ampl',ampl);
mysaveas(pathname,'Ux',formats);

plotSolution(S,u,'displ',2,'ampl',ampl);
mysaveas(pathname,'Uy',formats);

%% Test numerical solution
EL_series = linspace(EL*0.5,EL*1.5,50); % [MPa]
NUL_series = linspace(NUL*0.5,NUL*1.5,50);
GL_series = linspace(GL*0.5,GL*1.5,50); % [MPa]

% Plot relative error with respect to EL and GL
err_series = zeros(length(EL_series),length(GL_series));
for m=1:length(EL_series)
    EL_m = EL_series(m);
    parfor n=1:length(GL_series)
        GL_n = GL_series(n);
        x = [EL_m NUL GL_n];
        u_in = solveTractionIsotTrans(x,S);
        err_series(m,n) = norm(u_exp_in - u_in);
    end
end
err_series = err_series./norm(u_exp_in);
% [err_min,I] = min(err_series);
% [err_min,c] = min(err_min);
% r = I(c);
% GL_min = GL_series(c);
% EL_min = EL_series(r);

figure('Name','Surface plot: Relative error with respect to EL and GL')
clf
surfc(GL_series,EL_series,err_series,'EdgeColor','none');
colorbar
set(gca,'ColorScale','log')
view(-37.5,30) % default view
hold on
% scatter3(GL_min,EL_min,err_min,'MarkerEdgeColor','k','MarkerFaceColor','r');
scatter3(GL,EL,err,'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off
set(gca,'FontSize',fontsize)
set(gca,'ZScale','log')
xlabel('$G_L$ [MPa]','Interpreter',interpreter)
ylabel('$E_L$ [MPa]','Interpreter',interpreter)
zlabel('Relative error','Interpreter',interpreter)
%zlabel('Erreur relative','Interpreter',interpreter)
mysaveas(pathname,'error_EL_GL_3D',formats);

figure('Name','Contour plot: Relative error with respect to EL and GL')
clf
contourf(GL_series,EL_series,err_series,50);
colorbar
set(gca,'ColorScale','log')
hold on
% scatter(GL_min,EL_min,'MarkerEdgeColor','k','MarkerFaceColor','r');
scatter(GL,EL,'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off
set(gca,'FontSize',fontsize)
xlabel('$G_L$ [MPa]','Interpreter',interpreter)
ylabel('$E_L$ [MPa]','Interpreter',interpreter)
mysaveas(pathname,'error_EL_GL_2D',formats);

% Plot relative error with respect to EL and NUL
err_series = zeros(length(EL_series),length(NUL_series));
for m=1:length(EL_series)
    EL_m = EL_series(m);
    parfor n=1:length(NUL_series)
        NUL_n = NUL_series(n);
        x = [EL_m NUL_n GL];
        u_in = solveTractionIsotTrans(x,S);
        err_series(m,n) = norm(u_exp_in - u_in);
    end
end
err_series = err_series./norm(u_exp_in);
% [err_min,I] = min(err_series);
% [err_min,c] = min(err_min);
% r = I(c);
% NUL_min = NUL_series(c);
% EL_min = EL_series(r);

figure('Name','Surface plot: Relative error with respect to EL and NUL')
clf
surfc(NUL_series,EL_series,err_series,'EdgeColor','none');
colorbar
set(gca,'ColorScale','log')
view(-37.5,30) % default view
hold on
% scatter3(NUL_min,EL_min,err_min,'MarkerEdgeColor','k','MarkerFaceColor','r');
scatter3(NUL,EL,err,'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off
set(gca,'FontSize',fontsize)
set(gca,'ZScale','log')
xlabel('$\nu_L$','Interpreter',interpreter)
ylabel('$E_L$ [MPa]','Interpreter',interpreter)
zlabel('Relative error','Interpreter',interpreter)
%zlabel('Erreur relative','Interpreter',interpreter)
mysaveas(pathname,'error_EL_NUL_3D',formats);

figure('Name','Contour plot: Relative error with respect to EL and NUL')
clf
contourf(NUL_series,EL_series,err_series,50);
colorbar
set(gca,'ColorScale','log')
hold on
% scatter(NUL_min,EL_min,'MarkerEdgeColor','k','MarkerFaceColor','r');
scatter(NUL,EL,'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off
set(gca,'FontSize',fontsize)
xlabel('$\nu_L$','Interpreter',interpreter)
ylabel('$E_L$ [MPa]','Interpreter',interpreter)
mysaveas(pathname,'error_EL_NUL_2D',formats);

% Plot relative error with respect to NUL and GL
err_series = zeros(length(NUL_series),length(GL_series));
for m=1:length(NUL_series)
    NUL_m = NUL_series(m);
    parfor n=1:length(GL_series)
        GL_n = GL_series(n);
        x = [EL NUL_m GL_n];
        u_in = solveTractionIsotTrans(x,S);
        err_series(m,n) = norm(u_exp_in - u_in);
    end
end
err_series = err_series./norm(u_exp_in);
% [err_min,I] = min(err_series);
% [err_min,c] = min(err_min);
% r = I(c);
% GL_min = GL_series(c);
% NUL_min = NUL_series(r);

figure('Name','Surface plot: Relative error with respect to NUL and GL')
clf
surfc(GL_series,NUL_series,err_series,'EdgeColor','none');
colorbar
set(gca,'ColorScale','log')
view(-37.5,30) % default view
hold on
% scatter3(GL_min,NUL_min,err_min,'MarkerEdgeColor','k','MarkerFaceColor','r');
scatter3(GL,NUL,err,'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off
set(gca,'FontSize',fontsize)
set(gca,'ZScale','log')
xlabel('$G_L$ [MPa]','Interpreter',interpreter)
ylabel('$\nu_L$','Interpreter',interpreter)
zlabel('Relative error','Interpreter',interpreter)
%zlabel('Erreur relative','Interpreter',interpreter)
mysaveas(pathname,'error_NUL_GL_3D',formats);

figure('Name','Contour plot: Relative error with respect to NUL and GL')
clf
contourf(GL_series,NUL_series,err_series,50);
colorbar
set(gca,'ColorScale','log')
hold on
% scatter(GL_min,NUL_min,'MarkerEdgeColor','k','MarkerFaceColor','r');
scatter(GL,NUL,'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off
set(gca,'FontSize',fontsize)
xlabel('$G_L$ [MPa]','Interpreter',interpreter)
ylabel('$\nu_L$','Interpreter',interpreter)
mysaveas(pathname,'error_NUL_GL_2D',formats);

myparallel('stop');
