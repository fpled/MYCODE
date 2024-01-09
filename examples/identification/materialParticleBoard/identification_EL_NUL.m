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
renderer = 'OpenGL';

%% Identification
if solveProblem
% initial guess
EL0 = 1e2; % longitudinal Young modulus [MPa]
NUL0 = 1e-2; % longitudinal Poisson ratio

disp('Initial parameters');
disp('------------------');
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

numSamples = 27;
EL_data = cell(numSamples,1);
NUL_data = cell(numSamples,1);
err_num_data = cell(numSamples,1);

mean_EL_data = zeros(numSamples,1);
mean_NUL_data = zeros(numSamples,1);
std_EL_data = zeros(numSamples,1);
std_NUL_data = zeros(numSamples,1);

for j=1:numSamples
% for j=14
    
    numSample = ['B' num2str(j)];
    F = appliedLoad(numSample);
    [b,h,d,Iz] = dimSample(numSample);
    
    t = tic;
    
    numImages = length(F);
    EL = zeros(numImages,1);
    NUL = zeros(numImages,1);
    err = zeros(numImages,1);
    
    for k=1:numImages
    % for k=6
        
        numImage = num2str(k,'%02d');
        filenameDIC = [numSample '_00-' numImage '-Mesh'];
        load(fullfile(pathnameDIC,filenameDIC));
        
        [u_exp,coord] = extractCorreli(Job,Mesh,U,h,d); % [mm]
        
        node = NODE(coord,1:size(coord,1));
        elem = Mesh.TRI;
        elemtype = 'TRI3';
        % option = 'DEFO'; % plane strain
        option = 'CONT'; % plane stress
        S = MODEL('PLAN');
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
                [x,err(k),~,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
            case 'fminsearch'
                fun = @(x) funoptimThreePointBending(x,@solveThreePointBendingNum,u_exp_in,S);
                [x,err(k),exitflag,output] = fminsearch(fun,x0,options);
            case 'fminunc'
                fun = @(x) funoptimThreePointBending(x,@solveThreePointBendingNum,u_exp_in,S);
                [x,err(k),exitflag,output] = fminunc(fun,x0,options);
            case 'fmincon'
                fun = @(x) funoptimThreePointBending(x,@solveThreePointBendingNum,u_exp_in,S);
                [x,err(k),exitflag,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
        end
        
        EL(k) = x(1); % [MPa]
        NUL(k) = x(2);
        err(k) = sqrt(err(k))./norm(u_exp_in);
    end
    
    %% Outputs
    fprintf('\n')
    disp('+-----------------+')
    fprintf('| Sample #%2d      |\n',j)
    disp('+-----------------+-----------------+-----------------+')
    disp('| Young''s modulus | Poisson''s ratio |  Error between  |')
    disp('|     EL [MPa]    |       NUL       | U_num and U_exp |')
    disp('+-----------------+-----------------+-----------------+')
    for k=1:numImages
        fprintf('| %15.4f | %15.4f | %15.4e |\n',EL(k),NUL(k),err(k))
    end
    disp('+-----------------+-----------------+-----------------+')
    
    toc(t)
    
    EL_data{j} = EL;
    NUL_data{j} = NUL;
    err_num_data{j} = err;
    mean_EL_data(j) = mean(EL(initImage:end));
    mean_NUL_data(j) = mean(NUL(initImage:end));
    std_EL_data(j) = std(EL(initImage:end));
    std_NUL_data(j) = std(NUL(initImage:end));
end

%% Save variables
save(fullfile(pathname,filenameNum),'EL_data','NUL_data','err_num_data',...
    'mean_EL_data','mean_NUL_data','std_EL_data','std_NUL_data');
else
%% Load variables
load(fullfile(pathname,filenameNum),'EL_data','NUL_data','err_num_data',...
    'mean_EL_data','mean_NUL_data','std_EL_data','std_NUL_data');
end

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
        ylabel('Longitudinal Young''s modulus $E^L$ [MPa]','Interpreter',interpreter);
        %xlabel('Num\''ero d''image','Interpreter',interpreter);
        %ylabel('Module d''Young longitudinal $E^L$ [MPa]','Interpreter',interpreter);
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
        ylabel('Longitudinal Poisson''s ratio $\nu^L$','Interpreter',interpreter);
        %xlabel('Num\''ero d''image','Interpreter',interpreter);
        %ylabel('Coefficient de Poisson longitudinal $\nu^L$','Interpreter',interpreter);
        mysaveas(pathname,['data_NUL_' numSample],formats);
        mymatlab2tikz(pathname,['data_NUL_' numSample '.tex']);
    end
    
    figure('Name','Longitudinal Young''s modulus')
    clf
    bar(mean_EL_data);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Longitudinal Young''s modulus $E^L$ [MPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module d''Young longitudinal $E^L$ [MPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_EL',formats);
    mymatlab2tikz(pathname,'data_EL.tex');
    
    figure('Name','Longitudinal Poisson''s ratio')
    clf
    bar(mean_NUL_data);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Longitudinal Poisson''s ratio $\nu^L$','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Coefficient de Poisson longitudinal $\nu^L$','Interpreter',interpreter);
    mysaveas(pathname,'data_NUL',formats);
    mymatlab2tikz(pathname,'data_NUL.tex');
end
