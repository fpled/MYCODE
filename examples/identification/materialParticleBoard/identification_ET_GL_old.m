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
% for j=14
    
    numSample = ['B' num2str(j)];
    F = appliedLoad(numSample);
    [b,h,d,Iz] = dimSample(numSample);
    
    t = tic;
    
    numImages = length(F);
    ET = zeros(numImages,1);
    GL = zeros(numImages,1);
    R0 = zeros(numImages,1);
    U0 = zeros(numImages,1);
    V0 = zeros(numImages,1);
    err = zeros(numImages,1);
    
    for k=1:numImages
    % for k=[6,numImages]
        
        numImage = num2str(k,'%02d');
        filenameDIC = [numSample '_00-' numImage '-Mesh'];
        load(fullfile(pathnameDIC,filenameDIC));
        
        [u_exp,coord] = extractCorreli(Job,Mesh,U,h,d); % [mm]
        
        %-------------------------------
        % Reference and deformed meshes
        %-------------------------------
        % figure('name','Reference and deformed meshes')
        % triplot(Mesh.TRI,coord(:,1),coord(:,2),'r');
        % hold on
        % triplot(Mesh.TRI,coord(:,1)+Scal*u_exp(1:2:end),coord(:,2)+Scal*u_exp(2:2:end),'k');
        % hold off
        % axis image
        % grid on
        % set(gca,'YLim',[-15,10])
        % set(gca,'FontSize',fontsize)
        % xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
        % ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
        % mysaveas(pathname,['meshes_' numSample '_image_' numImage],formats,renderer);
        
        A = zeros(5);
        B = zeros(5,1);
        for i=1:Mesh.NNodeTot
            xi = coord(i,1);
            yi = coord(i,2);
            ui = u_exp(2*i-1);
            vi = u_exp(2*i);
            
            A_UE = -F(k)/(4*Iz)*yi*xi^2;
            A_UG = F(k)/(12*Iz)*yi^3;
            A_VE = F(k)/(12*Iz)*xi^3;
            A_VG = -F(k)*h^2/(16*Iz)*xi;
            
            A = A + [ A_UE^2+A_VE^2        A_UE*A_UG+A_VE*A_VG  A_UE*yi-A_VE*xi  A_UE  A_VE
                A_UE*A_UG+A_VE*A_VG  A_UG^2+A_VG^2        A_UG*yi-A_VG*xi  A_UG  A_VG
                A_UE*yi-A_VE*xi      A_UG*yi-A_VG*xi      yi^2+xi^2        yi    -xi
                A_UE                 A_UG                 yi               1     0
                A_VE                 A_VG                 -xi              0     1 ];
            
            B = B + [ A_UE*ui + A_VE*vi
                A_UG*ui + A_VG*vi
                yi*ui   - xi*vi
                ui
                vi ];
        end
        
        param = A\B;
        
        ET(k) = 1/param(1); % [MPa]
        GL(k) = 1/param(2); % [MPa]
        R0(k) = param(3); % [rad]
        U0(k) = param(4); % [mm]
        V0(k) = param(5); % [mm]
        
        x = [ET(k) GL(k) R0(k) U0(k) V0(k)];
        u = solveThreePointBendingAna(x,coord,F(k),Iz,h);
        err(k) = norm(u - u_exp);
        err(k) = err(k)./norm(u_exp);
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
        ylabel('Rigid body rotation $\psi_0$ [$^{\circ}$]','Interpreter',interpreter);
        %xlabel('Num\''ero d''image','Interpreter',interpreter);
        %ylabel('Rotation de corps rigide $\psi_0$ [$^{\circ}$]','Interpreter',interpreter);
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
    ylabel('Rigid body rotation $\psi_0$ [$^{\circ}$]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Rotation de corps rigide $\psi_0$ [$^{\circ}$]','Interpreter',interpreter);
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
