%% Identification of cyclic compression behavior of wood in the transverse direction %%
%%-----------------------------------------------------------------------------------%%

% clc
clearvars
close all

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialWoodCyclicCompression');
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Input data
%F = [120 240 360 480 600 720];
%F = [120 360 600];
%F = [240 480 720];
F = 300; % [N]
h = 50; % [mm]

%% Experimental data
pathnameExp = fileparts(mfilename('fullpath'));
filename = 'BG1-A-10k-300N-4000.txt';
% filename  ='BG1-A-10k-300N-propre.txt';
[time,displ,force] = readDataCyclicCompression(fullfile(pathnameExp,filename));

flag = true;
displ_max = [];
displ_min = [];
N = 0;
for i=1:length(displ)-1
    if flag && force(i+1)<=force(i)
        N = N+1;
        displ_max = [displ_max displ(i)];
        flag = false;
    elseif ~flag && force(i+1)>=force(i)
        displ_min = [displ_min displ(i)];
        flag = true;
    end
end
N = N-1;
displ_max(end) = [];

%% Identification
delta_exp = displ_min;
N = 1:length(delta_exp);

% Initial point
delta0 = 0.08; % [mm]
Kinf = 1e-5; % [mm]
Nref = 30;

disp('Initial parameters');
disp('------------------');
fprintf('delta0 = %g mm\n',delta0);
fprintf('Kinf   = %g mm\n',Kinf);
fprintf('Nref   = %g\n',Nref);

x0 = [delta0 Kinf Nref];
lb = [0 0 0];
ub = [h Inf length(N)];

optimFun = 'lsqnonlin'; % optimization function
% optimFun = 'fminsearch';
% optimFun = 'fminunc';
% optimFun = 'fmincon';

% display = 'off';
% display = 'iter';
display = 'iter-detailed';
% display = 'final';
% display = 'final-detailed';

tolX = eps; % tolerance on the parameter value
tolFun = eps; % tolerance on the function value
tolOpt = eps; % tolerance on the first-order optimality
maxIters = Inf; % maximum number of iterations
maxFunEvals = Inf; % maximum number of function evaluations

switch optimFun
    case {'lsqnonlin','fminunc','fmincon'}
        % options = optimoptions(optimFun,'Display',display,'TolX',tolX,'TolFun',tolFun,'MaxFunEvals',maxFunEvals);
        options = optimoptions(optimFun,'Display',display,...
            'StepTolerance',tolX,'FunctionTolerance',tolFun,'OptimalityTolerance',tolOpt,...
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
        fun = @(x) funlsqnonlinCyclicCompression(x,delta_exp,N);
        [x,resnorm,~,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
    case 'fminsearch'
        fun = @(x) funoptimCyclicCompression(x,delta_exp,N);
        [x,resnorm,exitflag,output] = fminsearch(fun,x0,options);
    case 'fminunc'
        fun = @(x) funoptimCyclicCompression(x,delta_exp,N);
        [x,resnorm,exitflag,output] = fminunc(fun,x0,options);
    case 'fmincon'
        fun = @(x) funoptimCyclicCompression(x,delta_exp,N);
        [x,resnorm,exitflag,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
end
toc(t)

delta0 = x(1);
Kinf = x(2);
Nref = x(3);
err = sqrt(resnorm)./norm(delta_exp);

fprintf('\n');
disp('Optimal parameters');
disp('------------------');
fprintf('delta0 = %g mm\n',delta0);
fprintf('Kinf   = %g mm\n',Kinf);
fprintf('Nref   = %g\n',Nref);
fprintf('resnorm = %g\n',resnorm);
fprintf('err = %g\n',err);
% fprintf('exitflag = %g\n',exitflag);
% disp(output);

%% Solution
delta = solveCyclicCompression(x,N);

%% Display solution
% Evolution of experimental displacement
figure
plot(time,displ,'-r','LineWidth',linewidth)
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Time (s)')
ylabel('Displacement [mm]')
mysaveas(pathname,'displacement',formats);

% Evolution of force
figure
plot(time,force,'-b','LineWidth',linewidth)
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Time [s]')
ylabel('Load [N]')
mysaveas(pathname,'force',formats);

% Evolution of experimental mininmal and maximal displacements
figure
plot(N,displ_max,'-r','LineWidth',linewidth)
hold on
plot(N,displ_min,'-b','LineWidth',linewidth)
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Number of cycles')
ylabel('Displacement [mm]')
legend('max','min');
mysaveas(pathname,'displacement_min_max',formats);

% Evolution of experimental and model damage
figure
plot(N,delta_exp,'.b')
hold on
plot(N,delta,'-r','LineWidth',linewidth)
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Number of cycles')
ylabel('Damage [mm]')
legend('experimental','model');
mysaveas(pathname,'solution',formats);
