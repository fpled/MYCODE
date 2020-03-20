%% Identification of ET and GL - Analytical beam model %%
%%-----------------------------------------------------%%
% Call after Digital Image Correlation (DIC) RT3
% Coordinate system
% CORRELI:    right     down
% MATLAB:     right     up

% clc
clearvars
close all

%% Input data
displaySolution = true;

filename = 'data_ET_GL.mat';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
if ~exist(pathname,'dir')
    mkdir(pathname);
end
pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','materialParticleBoard','resultsDIC');

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};

%% Identification
% initial guess
ET0 = 1e3; % transverse Young modulus [MPa]
GL0 = 1e2; % longitudinal shear modulus [MPa]
Phi0 = 0; % rigid body rotation around z direction [rad]
U0 = 0; % rigid body displacement along x direction [mm]
V0 = 0; % rigid body displacement along y direction [mm]

disp('Initial parameters');
disp('------------------');
fprintf('ET  = %g MPa\n',ET0);
fprintf('GL  = %g MPa\n',GL0);
fprintf('Phi = %g rad = %g deg\n',Phi0,rad2deg(Phi0));
fprintf('U   = %g mm\n',U0);
fprintf('V   = %g mm\n',V0);

x0 = [ET0 GL0 Phi0 U0 V0];
lb = [0 0 -Inf -Inf -Inf];
ub = [Inf Inf Inf Inf Inf];

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

sample = 'B';
numSamples = 27;
ET_data = cell(numSamples,1);
GL_data = cell(numSamples,1);
Phi_data = cell(numSamples,1);
U0_data = cell(numSamples,1);
V0_data = cell(numSamples,1);
err_ana_data = cell(numSamples,1);

mean_ET_data = zeros(numSamples,1);
mean_GL_data = zeros(numSamples,1);
std_ET_data = zeros(numSamples,1);
std_GL_data = zeros(numSamples,1);

initImage = 3;

for j=1:numSamples
    
    numSample = [sample num2str(j)];
    F = appliedLoad(numSample);
    [b,h,d,Iz] = dimSample(numSample);
    
    t = tic;
    
    numImages = length(F);
    ET = zeros(numImages,1);
    GL = zeros(numImages,1);
    Phi = zeros(numImages,1);
    U0 = zeros(numImages,1);
    V0 = zeros(numImages,1);
    err = zeros(numImages,1);
    
    for k=1:numImages
        
        numImage = num2str(k,'%02d');
        filenameDIC = [numSample '_00-' numImage '-Mesh'];
        load(fullfile(pathnameDIC,filenameDIC));
        
        [u_exp,coord] = extractCorreli(Job,Mesh,U,h,d); % [mm]
        
        switch optimFun
            case 'lsqnonlin'
                fun = @(x) funlsqnonlin(x,@solveThreePointBendingAna,u_exp,coord,F(k),Iz,h);
                [x,err(k),~,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
            case 'fminsearch'
                fun = @(x) funoptim(x,@solveThreePointBendingAna,u_exp,coord,F(k),Iz,h);
                [x,err(k),exitflag,output] = fminsearch(fun,x0,options);
            case 'fminunc'
                fun = @(x) funoptim(x,@solveThreePointBendingAna,u_exp,coord,F(k),Iz,h);
                [x,err(k),exitflag,output] = fminunc(fun,x0,options);
            case 'fmincon'
                fun = @(x) funoptim(x,@solveThreePointBendingAna,u_exp,coord,F(k),Iz,h);
                [x,err(k),exitflag,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
        end
        
        ET(k) = x(1); % [MPa]
        GL(k) = x(2); % [MPa]
        Phi(k) = x(3); % [rad]
        U0(k) = x(4); % [mm]
        V0(k) = x(5); % [mm]
        err(k) = sqrt(err(k))./norm(u_exp);
    end
    
    %% Outputs
    fprintf('\n')
    disp('+-----------------+')
    fprintf('| Sample %s%2d      |\n',sample,j)
    disp('+-----------------+---------------+-----------------+')
    disp('| Young''s modulus | Shear modulus |  Error between  |')
    disp('|     ET [GPa]    |    GL [MPa]   | U_ana and U_exp |')
    disp('+-----------------+---------------+-----------------+')
    for k=1:numImages
        fprintf('| %15.4f | %13.4f | %15.4e |\n',ET(k)*1e-3,GL(k),err(k))
    end
    disp('+-----------------+---------------+-----------------+')
    
    toc(t)
    
    ET_data{j} = ET;
    GL_data{j} = GL;
    Phi_data{j} = Phi;
    U0_data{j} = U0;
    V0_data{j} = V0;
    err_ana_data{j} = err;
    mean_ET_data(j) = mean(ET(initImage:end));
    mean_GL_data(j) = mean(GL(initImage:end));
    std_ET_data(j) = std(ET(initImage:end));
    std_GL_data(j) = std(GL(initImage:end));
end

%% Save variables
save(fullfile(pathname,filename),'ET_data','GL_data',...
    'Phi_data','U0_data','V0_data','err_ana_data','initImage',...
    'mean_ET_data','mean_GL_data','std_ET_data','std_GL_data');

%% Plot data
if displaySolution
%     for j=1:numSamples
%         numSample = ['B' num2str(j)];
%         
%         figure
%         clf
%         bar(ET_data{j}(initImage:end)*1e-3);
%         grid on
%         set(gca,'FontSize',fontsize)
%         legend(numSample,'Location','NorthEastOutside');
%         xlabel('Image number','Interpreter',interpreter);
%         ylabel('Young''s modulus $E^T$ [GPa]','Interpreter',interpreter);
%         %xlabel('Num\''ero d''image','Interpreter',interpreter);
%         %ylabel('Module d''Young $E^T$ [GPa]','Interpreter',interpreter);
%         mysaveas(pathname,['data_ET_' numSample],formats);
%         mymatlab2tikz(pathname,['data_ET_' numSample '.tex']);
%         
%         figure
%         clf
%         bar(GL_data{j}(initImage:end));
%         grid on
%         set(gca,'FontSize',fontsize)
%         legend(numSample,'Location','NorthEastOutside');
%         xlabel('Image number','Interpreter',interpreter);
%         ylabel('Shear modulus $G^L$ [MPa]','Interpreter',interpreter);
%         %xlabel('Num\''ero d''image','Interpreter',interpreter);
%         %ylabel('Module de cisaillement $G^L$ [MPa]','Interpreter',interpreter);
%         mysaveas(pathname,['data_GL_' numSample],formats);
%         mymatlab2tikz(pathname,['data_GL_' numSample '.tex']);
%     end
    
    figure
    clf
    bar(mean_ET_data*1e-3,'LineWidth',linewidth);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Young''s modulus $E^T$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module d''Young $E^T$ [GPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_ET',formats);
    mymatlab2tikz(pathname,'data_ET.tex');
    
    figure
    clf
    bar(mean_GL_data,'LineWidth',linewidth);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Shear modulus $G^L$ [MPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module de cisaillement $G^L$ [MPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_GL',formats);
    mymatlab2tikz(pathname,'data_GL.tex');
end
