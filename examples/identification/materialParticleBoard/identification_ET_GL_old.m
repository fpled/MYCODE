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
        
        [u_exp,coord] = extractCorreli(Job,Mesh,U,h,d);
        
set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',20,'DefaultTextInterpreter','tex');
figure('name','Reference and deformed mesh')
Scal=10;
Unitx=' (mm.)';
UnitU=' (mm.)';
triplot(Mesh.TRI,coord(:,1),coord(:,2),'r');
hold on
triplot(Mesh.TRI,coord(:,1)+Scal*u_exp(1:2:end),coord(:,2)+Scal*u_exp(2:2:end),'k');
hold off
axis equal
xlabel(['$y$ ',Unitx],'Interpreter','latex')
ylabel(['$z$ ',Unitx],'Interpreter','latex')
        
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
        
        ET(k) = 1/param(1); % MPa
        GL(k) = 1/param(2); % MPa
        Phi(k) = param(3);
        U0(k) = param(4); % mm
        V0(k) = param(5); % mm
        
        x = [ET(k) GL(k) Phi(k) U0(k) V0(k)];
        u = solveThreePointBendingAna(x,coord,F(k),Iz,h);
        err(k) = norm(u - u_exp);
        err(k) = err(k)./norm(u_exp);
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
    bar(mean_GL_data*1e-3,'LineWidth',linewidth);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Shear modulus $G^L$ [MPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module de cisaillement $G^L$ [MPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_GL',formats);
    mymatlab2tikz(pathname,'data_GL.tex');
end
