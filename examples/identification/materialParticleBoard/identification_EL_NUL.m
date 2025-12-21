%% Identification of EL and NUL - Finite Element Model Updating (FEMU) method %%
%%----------------------------------------------------------------------------%%
% Call after Digital Image Correlation (DIC) RT3 and Identification of ET and GL

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

filenameAna = 'data_ET_GL.mat';
filenameNum = 'data_EL_NUL.mat';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathname,filenameAna));
pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','materialParticleBoard','resultsDIC');

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};

%% Identification
if solveProblem
% Initial parameter values
EL0 = 1e2; % longitudinal Young modulus [MPa]
NUL0 = 1e-2; % longitudinal Poisson ratio

fprintf('\n');
fprintf('Initial parameter values\n');
fprintf('EL  = %g MPa\n',EL0);
fprintf('NUL = %g\n',NUL0);

x0 = [EL0 NUL0];
lb = [0 0];
ub = [Inf 0.5];

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
tolOpt = 1e-14; % tolerance on the first-order optimality
maxIters = Inf; % maximum number of iterations
maxFunEvals = Inf; % maximum number of function evaluations

switch optimFun
    case {'lsqnonlin','fminunc','fmincon'}
        % options = optimoptions(optimFun,'Display',display,'Algorithm',algo,...
        %     'TolX',tolX,'TolFun',tolFun,'MaxIter',maxIters,'MaxFunEvals',maxFunEvals);
        options = optimoptions(optimFun,'Display',display,'Algorithm',algo,...
            'StepTolerance',tolX,'FunctionTolerance',tolFun,'OptimalityTolerance',tolOpt,...
            'MaxIterations',maxIters,'MaxFunctionEvaluations',maxFunEvals);
    case 'fminsearch'
        options = optimset('Display',display,'TolX',tolX,'TolFun',tolFun,...
            'MaxIter',maxIters,'MaxFunEvals',maxFunEvals);
    otherwise
        error(['Wrong optimization function' optimFun])
end

% Geometric dimensions
b = 50; % sample width [mm]
h = 15; % sample thickness [mm]
Iz = b*h^3/12; % planar second moment of area (or planar area moment of inertia) [mm^4]
d = 20; % distance between the support and the region of interest (ROI) [mm]

numSamples = 27;
EL_data = cell(numSamples,1);
NUL_data = cell(numSamples,1);
resnorm_num = cell(numSamples,1);
err_num = cell(numSamples,1);

mean_EL_data = zeros(numSamples,1);
mean_NUL_data = zeros(numSamples,1);
std_EL_data = zeros(numSamples,1);
std_NUL_data = zeros(numSamples,1);

for j=1:numSamples
    
    numSample = ['B' num2str(j)];
    F = appliedLoad(numSample); % applied load [N]
    
    t = tic;
    
    numImages = length(F);
    EL = zeros(numImages,1);
    NUL = zeros(numImages,1);
    resnorm = zeros(numImages,1);
    err = zeros(numImages,1);
    
    for k=1:numImages
        
        numImage = num2str(k,'%02d');
        filenameDIC = [numSample '_00-' numImage '-Mesh'];
        load(fullfile(pathnameDIC,filenameDIC));
        
        [u_exp,coord] = extractCorreliElas(Job,Mesh,U,h,d); % [mm]
        
        node = NODE(coord,1:size(coord,1));
        elem = Mesh.TRI;
        elemtype = 'TRI3';
        % option = 'DEFO'; % plane strain
        option = 'CONT'; % plane stress
        S = MODEL('PLANE');
        S = addnode(S,node);
        S = addelem(S,elemtype,elem,'option',option);
        
        ET = ET_data{j}(k); % [MPa]
        GL = GL_data{j}(k); % [MPa]
        mat = ELAS_ISOT_TRANS('AXISL',[0;1],'AXIST',[1;0],'EL',[],'ET',ET,'NUL',[],'GL',GL,'DIM3',h);
        mat = setnumber(mat,1);
        S = setmaterial(S,mat);
        S = final(S);
        I = create_boundary(S);
        I = final(I);
        P = calc_P(S,I);
        u_exp_b = P*u_exp;
        S = addcl(S,[],'U',u_exp_b);
        u_exp_in = freevector(S,u_exp);
        
        switch optimFun
            case 'lsqnonlin'
                fun = @(x) funlsqnonlinThreePointBending(x,@solveThreePointBendingNum,u_exp_in,S);
                [x,resnorm(k),residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
            case 'fminsearch'
                fun = @(x) funoptimThreePointBending(x,@solveThreePointBendingNum,u_exp_in,S);
                [x,resnorm(k),exitflag,output] = fminsearch(fun,x0,options);
            case 'fminunc'
                fun = @(x) funoptimThreePointBending(x,@solveThreePointBendingNum,u_exp_in,S);
                [x,resnorm(k),exitflag,output] = fminunc(fun,x0,options);
            case 'fmincon'
                fun = @(x) funoptimThreePointBending(x,@solveThreePointBendingNum,u_exp_in,S);
                [x,resnorm(k),exitflag,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
        end
        
        % Optimal parameter values
        EL(k) = x(1); % [MPa]
        NUL(k) = x(2);
        err(k) = sqrt(resnorm(k))./norm(u_exp_in);
    end
    
    %% Outputs
    fprintf('\n')
    fprintf('Sample B%d\n',j)
    disp('+--------+-----------------+------------+---------------+----------------+')
    disp('| Image  | Young''s modulus | Poisson''s  | Objective fun | Relative error |')
    disp('| number |     EL [MPa]    | ratio  NUL |     value     | U_num vs U_exp |')
    disp('+--------+-----------------+------------+---------------+----------------+')
    for k=1:numImages
        fprintf('| %6d | %15.4f | %10.4f | %13.4e | %14.4e |\n',k,EL(k),NUL(k),resnorm(k),err(k))
    end
    disp('+--------+-----------------+------------+---------------+----------------+')
    
    toc(t)
    
    EL_data{j} = EL;
    NUL_data{j} = NUL;
    resnorm_num{j} = resnorm;
    err_num{j} = err;
    mean_EL_data(j) = mean(EL(initImage:end));
    mean_NUL_data(j) = mean(NUL(initImage:end));
    std_EL_data(j) = std(EL(initImage:end));
    std_NUL_data(j) = std(NUL(initImage:end));
end
close all

%% Save variables
save(fullfile(pathname,filenameNum),'EL_data','NUL_data','resnorm_num','err_num',...
    'mean_EL_data','mean_NUL_data','std_EL_data','std_NUL_data');
else
%% Load variables
load(fullfile(pathname,filenameNum),'EL_data','NUL_data','resnorm_num','err_num',...
    'mean_EL_data','mean_NUL_data','std_EL_data','std_NUL_data');
end

%% Statistics
fprintf('\n');
fprintf('Longitudinal Young''s modulus EL\n');
fprintf('mean(EL) = %g GPa\n',mean(mean_EL_data));
fprintf('var(EL)  = %g (GPa)^2\n',var(mean_EL_data));
fprintf('std(EL)  = %g GPa\n',std(mean_EL_data));
fprintf('cv(EL)   = %g\n',std(mean_EL_data)/mean(mean_EL_data));

fprintf('\n');
fprintf('Longitudinal Poisson''s ratio NUL\n');
fprintf('mean(NUL) = %g MPa\n',mean(mean_NUL_data));
fprintf('var(NUL)  = %g (MPa)^2\n',var(mean_NUL_data));
fprintf('std(NUL)  = %g MPa\n',std(mean_NUL_data));
fprintf('cv(NUL)   = %g\n',std(mean_NUL_data)/mean(mean_NUL_data));

%% Plot data
if displaySolution
    numSamples = length(EL_data);
    for j=1:numSamples
        numSample = ['B' num2str(j)];
        
        figure('Name',['Sample ' numSample ': Longitudinal Young''s modulus'])
        clf
        bar(EL_data{j});
        %bar(EL_data{j}(initImage:end));
        grid on
        set(gca,'FontSize',fontsize)
        %legend(numSample,'Location','NorthEastOutside');
        xlabel('Image number','Interpreter',interpreter);
        ylabel('Longitudinal Young''s modulus $E_L$ [MPa]','Interpreter',interpreter);
        %xlabel('Num\''ero d''image','Interpreter',interpreter);
        %ylabel('Module d''Young longitudinal $E_L$ [MPa]','Interpreter',interpreter);
        mysaveas(pathname,['data_EL_' numSample],formats);
        mymatlab2tikz(pathname,['data_EL_' numSample '.tex']);
        
        figure('Name',['Sample ' numSample ': Longitudinal Poisson''s ratio'])
        clf
        bar(NUL_data{j});
        %bar(NUL_data{j}(initImage:end));
        grid on
        set(gca,'FontSize',fontsize)
        %legend(numSample,'Location','NorthEastOutside');
        xlabel('Image number','Interpreter',interpreter);
        ylabel('Longitudinal Poisson''s ratio $\nu_L$','Interpreter',interpreter);
        %xlabel('Num\''ero d''image','Interpreter',interpreter);
        %ylabel('Coefficient de Poisson longitudinal $\nu_L$','Interpreter',interpreter);
        mysaveas(pathname,['data_NUL_' numSample],formats);
        mymatlab2tikz(pathname,['data_NUL_' numSample '.tex']);
    end
    close all
    
    figure('Name','Longitudinal Young''s modulus')
    clf
    bar(mean_EL_data);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Longitudinal Young''s modulus $E_L$ [MPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module d''Young longitudinal $E_L$ [MPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_EL',formats);
    mymatlab2tikz(pathname,'data_EL.tex');
    
    figure('Name','Longitudinal Poisson''s ratio')
    clf
    bar(mean_NUL_data);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Longitudinal Poisson''s ratio $\nu_L$','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Coefficient de Poisson longitudinal $\nu_L$','Interpreter',interpreter);
    mysaveas(pathname,'data_NUL',formats);
    mymatlab2tikz(pathname,'data_NUL.tex');
end
