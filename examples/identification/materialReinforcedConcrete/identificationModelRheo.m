%% Identification of an equivalent rheological model                     %%
% Reinforced Concrete (RC) wall under cyclic loading                     %%
% 2 springs in series with stiffnesses kc and kc*z(d) (damaged concrete) %%
% working in parallel with 2 springs in series with stiffness ks (steel) %%
% with a sliding pad representing concrete-steel debonding               %%
%%-----------------------------------------------------------------------%%

% clc
clearvars
close all

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

damageFun = 2; % choice for damage function

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard',['damageFunction' num2str(damageFun)]);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Experimental data
pathnameExp = fileparts(mfilename('fullpath'));
filename = 'resu_exp.11';
[displ,force_exp] = readData(fullfile(pathnameExp,filename)); % displacement [m], force [N]
force_exp = force_exp*1e-6; % force [MN]
L = 3; % wall length [m]
h = 1.2; % wall height [m]
t = 0.2; % wall thickness [m]
l = sqrt(L^2+h^2); % rod length [m]
angle = atan(h/L); % angle between rod and axis x [rad]
epsilon = displ/l; % strain

%% Identification
% Initial point
S = sqrt(L*h)*t; % equivalent rod section [m^2]
Es = 200e3; % elastic modulus of steel [MPa]
Ec = 28.6e3; % elastic modulus of concrete [MPa]
ft = 3.32; % concrete failure under tension [MPa]
fc = 46.4; % concrete failure under compression [MPa]
sig0t = 40/100*ft; % first damage stress threshold for concrete under tension [MPa]
sig0c = 10/100*fc; % first damage stress threshold for concrete under compression [MPa]
Dsig0 = 1; % first debonding stress threshold for concrete-steel bond [MPa]

ps = 0.8/100; % proportion of steel in RC section
pc = 1-ps; % proportion of concrete in RC section
ks = Es*ps; % elastic modulus of steel times proportion of steel in RC section [MPa]
kc = Ec*pc; % elastic modulus of concrete times proportion of concrete in RC section [MPa]
F0t = sig0t*S; % first damage force threshold for concrete under tension [MN]
F0c = sig0c*S; % first damage force threshold for concrete under compression [MN]
DF0 = Dsig0*S; % first debonding force threshold for concrete-steel bond [MN]

switch damageFun
    case 1
        % damage parameters for concrete under tension
        alphaGt = 1;
        alphadt = 1.5;
        % damage parameters for concrete under compression
        alphaGc = 0.5;
        alphadc = 1.5;
    case 2
        gammat = -0.01; % damage parameter for concrete under tension
        gammac = 0.8; % damage parameter for concrete under compression
end

disp(['Damage function = ' num2str(damageFun)]);
disp('-------------------');

disp('Initial parameters');
disp('------------------');
fprintf('S  = %g m^2\n',S);
fprintf('Es = %g GPa\n',Es*1e-3);
fprintf('Ec = %g GPa\n',Ec*1e-3);
fprintf('ks = %g GPa\n',ks*1e-3);
fprintf('kc = %g GPa\n',kc*1e-3);
fprintf('sig0t = %g MPa\n',sig0t);
fprintf('sig0c = %g MPa\n',sig0c);
fprintf('Dsig0 = %g MPa\n',Dsig0);
fprintf('F0t = %g MN\n',F0t);
fprintf('F0c = %g MN\n',F0c);
fprintf('DF0 = %g MN\n',DF0);
switch damageFun
    case 1
        fprintf('alphaGt = %g\n',alphaGt);
        fprintf('alphadt = %g\n',alphadt);
        fprintf('alphaGc = %g\n',alphaGc);
        fprintf('alphadc = %g\n',alphadc);
    case 2
        fprintf('gammat = %g\n',gammat);
        fprintf('gammac = %g\n',gammac);
end

x0 = [S ks kc sig0t sig0c Dsig0];
% lb = [0 0 0 0 0 0];
% ub = [Inf Inf Inf Inf Inf Inf];
lb = [0 160e3*ps 15e3*pc 0 0 0];
ub = [Inf 240e3*ps 45e3*pc Inf Inf Inf];
switch damageFun
    case 1
        x0 = [x0 alphaGt alphadt alphaGc alphadc];
        lb = [lb 1e-2 1e-2 1e-2 1e-2];
        ub = [ub Inf Inf Inf Inf];
    case 2
        x0 = [x0 gammat gammac];
        lb = [lb -Inf -Inf];
        ub = [ub 1 1];
end
optimFun = 'lsqnonlin'; % optimization function
% optimFun = 'fminsearch';
% optimFun = 'fminunc';
% optimFun = 'fmincon';

% display = 'off';
% display = 'iter';
display = 'iter-detailed';
% display = 'final';
% display = 'final-detailed';

tolX = 1e-6; % tolerance on the parameter value in optimization algorithm
tolFun = 1e-6; % tolerance on the function value in optimization algorithm
tol = 1e-10; % tolerance for fixed-point algorithm
maxFunEvals = 1e3; % maximum number of function evaluations

switch optimFun
    case {'lsqnonlin','fminunc','fmincon'}
        options  = optimoptions(optimFun,'Display',display,'TolX',tolX,'TolFun',tolFun,'MaxFunEvals',maxFunEvals);
        % options  = optimoptions(optimFun,'Display',display,'StepTolerance',tolX,'FunctionTolerance',tolFun,'OptimalityTolerance',tolFun,'MaxFunctionEvaluations',maxFunEvals);
    case 'fminsearch'
        options = optimset('Display',display,'TolX',tolX,'TolFun',tolFun,'MaxFunEvals',maxFunEvals);
    otherwise
        error(['Wrong optimization function' optimFun])
end

time = tic;
switch optimFun
    case 'lsqnonlin'
        fun = @(x) funlsqnonlin(x,force_exp,epsilon,angle,damageFun,tol);
        [x,err,~,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
    case 'fminsearch'
        fun = @(x) funoptim(x,force_exp,epsilon,angle,damageFun,tol);
        [x,err,exitflag,output] = fminsearch(fun,x0,options);
    case 'fminunc'
        fun = @(x) funoptim(x,force_exp,epsilon,angle,damageFun,tol);
        [x,err,exitflag,output] = fminunc(fun,x0,options);
    case 'fmincon'
        fun = @(x) funoptim(x,force_exp,epsilon,angle,damageFun,tol);
        [x,err,exitflag,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
end
% x = x0;
% err = 0;

S = x(1);
ks = x(2);
kc = x(3);
sig0t = x(4);
sig0c = x(5);
Dsig0 = x(6);
switch damageFun
    case 1
        alphaGt = x(7);
        alphadt = x(8);
        alphaGc = x(9);
        alphadc = x(10);
    case 2
        gammat = x(7);
        gammac = x(8);
end
err = sqrt(err)./norm(force_exp);

Es = ks/ps;
Ec = kc/pc;
F0t = sig0t*S;
F0c = sig0c*S;
DF0 = Dsig0*S;

disp('Optimal parameters');
disp('------------------');
fprintf('S  = %g m^2\n',S);
fprintf('Es = %g GPa\n',Es*1e-3);
fprintf('Ec = %g GPa\n',Ec*1e-3);
fprintf('ks = %g GPa\n',ks*1e-3);
fprintf('kc = %g GPa\n',kc*1e-3);
fprintf('sig0t = %g MPa\n',sig0t);
fprintf('sig0c = %g MPa\n',sig0c);
fprintf('Dsig0 = %g MPa\n',Dsig0);
fprintf('F0t = %g MN\n',F0t);
fprintf('F0c = %g MN\n',F0c);
fprintf('DF0 = %g MN\n',DF0);
switch damageFun
    case 1
        fprintf('alphaGt = %g\n',alphaGt);
        fprintf('alphadt = %g\n',alphadt);
        fprintf('alphaGc = %g\n',alphaGc);
        fprintf('alphadc = %g\n',alphadc);
    case 2
        fprintf('gammat = %g\n',gammat);
        fprintf('gammac = %g\n',gammac);
end
fprintf('err = %g\n',err);
% fprintf('exitflag = %g\n',exitflag);
% disp(output);
toc(time)

%% Solution
strain = epsilon(:)*cos(angle)*[1 -1]; 
stress = zeros(size(strain));
stressDebonding = zeros(size(strain));
strainDebonding = zeros(size(strain));
damage = zeros(size(strain));
dissipatedEnergy = zeros(size(strain));
dissipatedEnergyDamage = zeros(size(strain));
dissipatedEnergyDebonding = zeros(size(strain));
parfor i=1:2
    [stress(:,i),stressDebonding(:,i),strainDebonding(:,i),damage(:,i),...
        dissipatedEnergy(:,i),dissipatedEnergyDamage(:,i),dissipatedEnergyDebonding(:,i)]...
        = solveModelRheo(x,strain(:,i),damageFun,tol,false);
end
sigma = (stress(:,1)-stress(:,2))/(2*cos(angle)); % [MPa]
force = sigma*S; % [MN]
forceDebonding = stressDebonding*S; % [MN]
displDebonding = strainDebonding*l; % [m]

%% Display solution
% Force-displacement curves
figure
clf
plot(displ*1e3,force_exp,'-r','Linewidth',linewidth)
hold on
plot(displ*1e3,force,'-b','Linewidth',linewidth)
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Deplacement [mm]','Interpreter',interpreter)
ylabel('Force [MN]','Interpreter',interpreter)
leg = legend('$F_{\mathrm{exp}}$','$F$','Location','NorthWest');
set(leg,'Interpreter',interpreter);
mysaveas(pathname,'force_displacement',formats);

% Evolution of displacement
figure
clf
plot(displ*1e3,'-b','Linewidth',linewidth)
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Nombre d''iterations','Interpreter',interpreter)
ylabel('Deplacement [mm]','Interpreter',interpreter)
mysaveas(pathname,'displacement',formats);

% Evolution of force
figure
clf
plot(force_exp,'-r','Linewidth',linewidth)
hold on
plot(force,'-b','Linewidth',linewidth)
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Nombre d''iterations','Interpreter',interpreter)
ylabel('Force [MN]','Interpreter',interpreter)
leg = legend('$F_{\mathrm{exp}}$','$F$','Location','NorthWest');
set(leg,'Interpreter',interpreter);
mysaveas(pathname,'force',formats);

% Evolution of damage
figure
clf
plot(damage(:,1),'-b','Linewidth',linewidth)
hold on
plot(damage(:,2),'-r','Linewidth',linewidth)
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Nombre d''iterations','Interpreter',interpreter)
ylabel('Endommagement','Interpreter',interpreter)
leg = legend('barre 1','barre 2','Location','NorthWest');
% set(leg,'Interpreter',interpreter);
mysaveas(pathname,'damage',formats);

% Evolution of debonding displacement
figure
clf
plot(displDebonding(:,1)*1e3,'-b','Linewidth',linewidth)
hold on
plot(displDebonding(:,2)*1e3,'-r','Linewidth',linewidth)
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Nombre d''iterations','Interpreter',interpreter)
ylabel('Glissement [mm]','Interpreter',interpreter)
leg = legend('barre 1','barre 2','Location','NorthWest');
% set(leg,'Interpreter',interpreter);
mysaveas(pathname,'displacement_debonding',formats);

% Evolution of debonding force
figure
clf
plot(forceDebonding(:,1),'-b','Linewidth',linewidth)
hold on
plot(forceDebonding(:,2),'-r','Linewidth',linewidth)
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Nombre d''iterations','Interpreter',interpreter)
ylabel('Force de glissement [MN]','Interpreter',interpreter)
leg = legend('barre 1','barre 2','Location','NorthWest');
% set(leg,'Interpreter',interpreter);
mysaveas(pathname,'force_debonding',formats);

% Evolution of dissipated energy
figure
clf
plot(dissipatedEnergy(:,1),'-b','Linewidth',linewidth)
hold on
plot(dissipatedEnergyDamage(:,1),'--b','Linewidth',linewidth)
plot(dissipatedEnergyDebonding(:,1),'-.b','Linewidth',linewidth)
plot(dissipatedEnergy(:,2),'-r','Linewidth',linewidth)
plot(dissipatedEnergyDamage(:,2),'--r','Linewidth',linewidth)
plot(dissipatedEnergyDebonding(:,2),'-.r','Linewidth',linewidth)
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Nombre d''iterations','Interpreter',interpreter)
ylabel('Energie dissipee [J/m$^3$]','Interpreter',interpreter)
leg = legend('barre 1 - total','barre 1 - endommagement','barre 1 - glissement',...
    'barre 2 - total','barre 2 - endommagement','barre 2 - glissement','Location','NorthWest');
% set(leg,'Interpreter',interpreter);
mysaveas(pathname,'disspated_energy',formats);
