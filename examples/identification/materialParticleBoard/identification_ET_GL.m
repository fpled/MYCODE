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
solveProblem = true;
displayMesh = true;
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
renderer = 'OpenGL';

Scal = 1;
Unitx = '[mm]';
UnitU = '[mm]';

%% Identification
if solveProblem
% initial guess
ET0 = 1e3; % transverse Young modulus [MPa]
GL0 = 1e2; % longitudinal shear modulus [MPa]
R0 = 0; % rigid body rotation around z direction [rad]
U0 = 0; % rigid body displacement along x direction [mm]
V0 = 0; % rigid body displacement along y direction [mm]

disp('Initial parameters');
disp('------------------');
fprintf('ET  = %g MPa\n',ET0);
fprintf('GL  = %g MPa\n',GL0);
fprintf('R   = %g rad = %g deg\n',R0,rad2deg(R0));
fprintf('U   = %g mm\n',U0);
fprintf('V   = %g mm\n',V0);

x0 = [ET0 GL0 R0 U0 V0];
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

% geometric dimensions
b = 50; % sample width [mm]
h = 15; % sample thickness [mm]
Iz = b*h^3/12; % planar second moment of area (or planar area moment of inertia) [mm^4]
d = 20; % distance between the support and the region of interest (ROI) [mm]
L = 85; % length of ROI [mm]

numSamples = 27;
ET_data = cell(numSamples,1);
GL_data = cell(numSamples,1);
R0_data = cell(numSamples,1);
U0_data = cell(numSamples,1);
V0_data = cell(numSamples,1);
err_ana_data = cell(numSamples,1);

mean_ET_data = zeros(numSamples,1);
mean_GL_data = zeros(numSamples,1);
mean_R0_data = zeros(numSamples,1);
mean_U0_data = zeros(numSamples,1);
mean_V0_data = zeros(numSamples,1);
std_ET_data = zeros(numSamples,1);
std_GL_data = zeros(numSamples,1);
std_R0_data = zeros(numSamples,1);
std_U0_data = zeros(numSamples,1);
std_V0_data = zeros(numSamples,1);

initImage = 3;

for j=1:numSamples
    
    numSample = ['B' num2str(j)];
    F = appliedLoad(numSample); % applied load [N]
    
    t = tic;
    
    numImages = length(F);
    ET = zeros(numImages,1);
    GL = zeros(numImages,1);
    R0 = zeros(numImages,1);
    U0 = zeros(numImages,1);
    V0 = zeros(numImages,1);
    err = zeros(numImages,1);
    
    for k=1:numImages
        
        numImage = num2str(k,'%02d');
        filenameDIC = [numSample '_00-' numImage '-Mesh'];
        load(fullfile(pathnameDIC,filenameDIC));
        
        [u_exp,coord] = extractCorreliElas(Job,Mesh,U,h,d); % [mm]
        
        %-------------------------------
        % Reference and deformed meshes
        %-------------------------------
        if displayMesh && (k==6 || k==numImages)
            figure('name',['Sample ' numSample ' - Image ' numImage ': Reference and deformed meshes'])
            triplot(Mesh.TRI,coord(:,1),coord(:,2),'k');
            hold on
            triplot(Mesh.TRI,coord(:,1)+Scal*u_exp(1:2:end),coord(:,2)+Scal*u_exp(2:2:end),'r');
            hold off
            axis image
            grid on
            set(gca,'XLim',[d-5,L+5])
            set(gca,'YLim',[-h,2*h/3])
            set(gca,'FontSize',fontsize)
            xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
            ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
            mysaveas(pathname,['meshes_' numSample '_' numImage],formats,renderer);
        end
        
        switch optimFun
            case 'lsqnonlin'
                fun = @(x) funlsqnonlinThreePointBending(x,@solveThreePointBendingAna,u_exp,coord,F(k),Iz,h);
                [x,err(k),~,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
            case 'fminsearch'
                fun = @(x) funoptimThreePointBending(x,@solveThreePointBendingAna,u_exp,coord,F(k),Iz,h);
                [x,err(k),exitflag,output] = fminsearch(fun,x0,options);
            case 'fminunc'
                fun = @(x) funoptimThreePointBending(x,@solveThreePointBendingAna,u_exp,coord,F(k),Iz,h);
                [x,err(k),exitflag,output] = fminunc(fun,x0,options);
            case 'fmincon'
                fun = @(x) funoptimThreePointBending(x,@solveThreePointBendingAna,u_exp,coord,F(k),Iz,h);
                [x,err(k),exitflag,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
        end
        
        ET(k) = x(1); % [MPa]
        GL(k) = x(2); % [MPa]
        R0(k) = x(3); % [rad]
        U0(k) = x(4); % [mm]
        V0(k) = x(5); % [mm]
        err(k) = sqrt(err(k))./norm(u_exp);
    end
    
    %% Outputs
    fprintf('\n')
    disp('+-----------------+')
    fprintf('| Sample #%2d      |\n',j)
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
    R0_data{j} = R0;
    U0_data{j} = U0;
    V0_data{j} = V0;
    err_ana_data{j} = err;
    mean_ET_data(j) = mean(ET(initImage:end));
    mean_GL_data(j) = mean(GL(initImage:end));
    mean_R0_data(j) = mean(R0(initImage:end));
    mean_U0_data(j) = mean(U0(initImage:end));
    mean_V0_data(j) = mean(V0(initImage:end));
    std_ET_data(j) = std(ET(initImage:end));
    std_GL_data(j) = std(GL(initImage:end));
    std_R0_data(j) = std(R0(initImage:end));
    std_U0_data(j) = std(U0(initImage:end));
    std_V0_data(j) = std(V0(initImage:end));
end
close all

%% Save variables
save(fullfile(pathname,filename),'ET_data','GL_data',...
    'R0_data','U0_data','V0_data','err_ana_data','initImage',...
    'mean_ET_data','mean_GL_data','std_ET_data','std_GL_data',...
    'mean_R0_data','mean_U0_data','mean_V0_data','std_R0_data','std_U0_data','std_V0_data');
else
%% Load variables
load(fullfile(pathname,filename),'ET_data','GL_data',...
    'R0_data','U0_data','V0_data','err_ana_data','initImage',...
    'mean_ET_data','mean_GL_data','std_ET_data','std_GL_data',...
    'mean_R0_data','mean_U0_data','mean_V0_data','std_R0_data','std_U0_data','std_V0_data');
end

%% Statistics
fprintf('\nTransverse Young''s modulus ET\n');
fprintf('mean(ET) = %g GPa\n',mean(mean_ET_data*1e-3));
fprintf('var(ET)  = %g (GPa)^2\n',var(mean_ET_data*1e-3));
fprintf('std(ET)  = %g GPa\n',std(mean_ET_data*1e-3));
fprintf('cv(ET)   = %g\n',std(mean_ET_data)/mean(mean_ET_data));

fprintf('\nLongitudinal shear modulus GL\n');
fprintf('mean(GL) = %g MPa\n',mean(mean_GL_data));
fprintf('var(GL) = %g (MPa)^2\n',var(mean_GL_data));
fprintf('std(GL)  = %g MPa\n',std(mean_GL_data));
fprintf('cv(GL)   = %g\n',std(mean_GL_data)/mean(mean_GL_data));

%% Plot data
if displaySolution
    numSamples = length(ET_data);
    for j=1:numSamples
        numSample = ['B' num2str(j)];
        
        figure('Name',['Sample ' numSample ': Transverse Young''s modulus'])
        clf
        bar(ET_data{j}*1e-3);
        %bar(ET_data{j}(initImage:end)*1e-3);
        grid on
        set(gca,'FontSize',fontsize)
        %legend(numSample,'Location','NorthEastOutside');
        xlabel('Image number','Interpreter',interpreter);
        ylabel('Transverse Young''s modulus $E_T$ [GPa]','Interpreter',interpreter);
        %xlabel('Num\''ero d''image','Interpreter',interpreter);
        %ylabel('Module d''Young transverse $E_T$ [GPa]','Interpreter',interpreter);
        mysaveas(pathname,['data_ET_' numSample],formats);
        mymatlab2tikz(pathname,['data_ET_' numSample '.tex']);
        
        figure('Name',['Sample ' numSample ': Longitudinal shear modulus'])
        clf
        bar(GL_data{j});
        %bar(GL_data{j}(initImage:end));
        grid on
        set(gca,'FontSize',fontsize)
        %legend(numSample,'Location','NorthEastOutside');
        xlabel('Image number','Interpreter',interpreter);
        ylabel('Longitudinal shear modulus $G_L$ [MPa]','Interpreter',interpreter);
        %xlabel('Num\''ero d''image','Interpreter',interpreter);
        %ylabel('Module de cisaillement longitudinal $G_L$ [MPa]','Interpreter',interpreter);
        mysaveas(pathname,['data_GL_' numSample],formats);
        mymatlab2tikz(pathname,['data_GL_' numSample '.tex']);
        
        figure('Name',['Sample ' numSample ': Rigid body rotation'])
        clf
        bar(rad2deg(R0_data{j}));
        %bar(rad2deg(R0_data{j}(initImage:end)));
        grid on
        set(gca,'FontSize',fontsize)
        %legend(numSample,'Location','NorthEastOutside');
        xlabel('Image number','Interpreter',interpreter);
        ylabel('Rigid body rotation $\psi_0$ [deg]','Interpreter',interpreter);
        %xlabel('Num\''ero d''image','Interpreter',interpreter);
        %ylabel('Rotation de corps rigide $\psi_0$ [deg]','Interpreter',interpreter);
        mysaveas(pathname,['data_R0_' numSample],formats);
        mymatlab2tikz(pathname,['data_R0_' numSample '.tex']);
        
        figure('Name',['Sample ' numSample ': Horizontal rigid body translation'])
        clf
        bar(U0_data{j});
        %bar(U0_data{j}(initImage:end)));
        grid on
        set(gca,'FontSize',fontsize)
        %legend(numSample,'Location','NorthEastOutside');
        xlabel('Image number','Interpreter',interpreter);
        ylabel('Horizontal rigid body translation $u_0$ [mm]','Interpreter',interpreter);
        %xlabel('Num\''ero d''image','Interpreter',interpreter);
        %ylabel('Translation horizontale de corps rigide $u_0$ [mm]','Interpreter',interpreter);
        mysaveas(pathname,['data_U0_' numSample],formats);
        mymatlab2tikz(pathname,['data_U0_' numSample '.tex']);
        
        figure('Name',['Sample ' numSample ': Vertical rigid body translation'])
        clf
        bar(V0_data{j});
        %bar(V0_data{j}(initImage:end)));
        grid on
        set(gca,'FontSize',fontsize)
        %legend(numSample,'Location','NorthEastOutside');
        xlabel('Image number','Interpreter',interpreter);
        ylabel('Vertical rigid body translation $v_0$ [mm]','Interpreter',interpreter);
        %xlabel('Num\''ero d''image','Interpreter',interpreter);
        %ylabel('Translation verticale de corps rigide $v_0$ [mm]','Interpreter',interpreter);
        mysaveas(pathname,['data_V0_' numSample],formats);
        mymatlab2tikz(pathname,['data_V0_' numSample '.tex']);
    end
    close all
    
    figure('Name','Transverse Young''s modulus')
    clf
    bar(mean_ET_data*1e-3);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Transverse Young''s modulus $E_T$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module d''Young transverse $E_T$ [GPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_ET',formats);
    mymatlab2tikz(pathname,'data_ET.tex');
    
    figure('Name','Longitudinal shear modulus')
    clf
    bar(mean_GL_data);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Longitudinal shear modulus $G_L$ [MPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module de cisaillement longitudinal $G_L$ [MPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_GL',formats);
    mymatlab2tikz(pathname,'data_GL.tex');
    
    figure('Name','Rigid body rotation')
    clf
    bar(rad2deg(mean_R0_data));
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Rigid body rotation $\psi_0$ [deg]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Rotation de corps rigide $\psi_0$ [deg]','Interpreter',interpreter);
    mysaveas(pathname,'data_R0',formats);
    mymatlab2tikz(pathname,'data_R0.tex');
    
    figure('Name','Horizontal rigid body translation')
    clf
    bar(mean_U0_data);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Horizontal rigid body translation $u_0$ [mm]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Translation horizontale de corps rigide $u_0$ [mm]','Interpreter',interpreter);
    mysaveas(pathname,'data_U0',formats);
    mymatlab2tikz(pathname,'data_U0.tex');
    
    figure('Name','Vertical rigid body translation')
    clf
    bar(mean_V0_data);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Vertical rigid body translation $v_0$ [mm]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Translation verticale de corps rigide $v_0$ [mm]','Interpreter',interpreter);
    mysaveas(pathname,'data_V0',formats);
    mymatlab2tikz(pathname,'data_V0.tex');
end
