%% Stochastic modeling of random linear elasticity tensor %%
%  with transversely isotropic symmetry                   %%
%%--------------------------------------------------------%%

% clc
clearvars
close all

%% Input data
displaySolution = true;

MCMC = 'CUM';% Markov-Chain Monte Carlo (MCMC) method = 'MH', 'BUM' or 'CUM'
filename = ['modelStoLinElasIsotTrans_ElasTensor_' MCMC];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
linewidth = 1;
markersize = 36;
interpreter = 'latex';
formats = {'fig','epsc'};

%% Data
filenameAna = 'data_ET_GL.mat';
filenameNum = 'data_EL_NUL.mat';
pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathnameIdentification,filenameAna));
load(fullfile(pathnameIdentification,filenameNum));

ET_data = mean_ET_data*1e-3; % GPa
GL_data = mean_GL_data*1e-3; % GPa
EL_data = mean_EL_data*1e-3; % GPa
NUL_data = mean_NUL_data;
NUT_data = 0.1+0.2*rand(length(mean_ET_data),1); % artificial data for NUT varying from 0.1 to 0.3
GT_data = ET_data./(2*(1+NUT_data)); % GPa
kT_data = (EL_data.*ET_data)./(2*(1-NUT_data).*EL_data-4*ET_data.*(NUL_data).^2); % GPa
C1_data = EL_data + 4*(NUL_data.^2).*kT_data; % GPa
C2_data = 2*kT_data; % GPa
C3_data = 2*sqrt(2)*kT_data.*NUL_data; % GPa
C4_data = 2*GT_data; % GPa
C5_data = 2*GL_data; % GPa
fprintf('\nnb data = %d',length(ET_data));

%% Sample generation
N = 5e3; % number of samples
switch lower(MCMC)
    case 'mh'
        [C_sample,lambda] = mhsampleStoLinElasTensorIsotTrans(C1_data,C2_data,C3_data,C4_data,C5_data,N);
    case 'bum'
        [C_sample,lambda] = mhsampleStoLinElasTensorIsotTrans_BUM(C1_data,C2_data,C3_data,C4_data,C5_data,N);
    case 'cum'
        [C_sample,lambda] = mhsampleStoLinElasTensorIsotTrans_CUM(C1_data,C2_data,C3_data,C4_data,C5_data,N);
end

kT_sample = C_sample(2,:)/2;
NUL_sample = (C_sample(3,:)./kT_sample)/(2*sqrt(2));
EL_sample = C_sample(1,:) - 4*(NUL_sample.^2).*kT_sample;
GT_sample = C_sample(4,:)/2;
GL_sample = C_sample(5,:)/2;
k1 = 1./kT_sample+4*(NUL_sample.^2)./EL_sample;
k2 = 1./GT_sample;
ET_sample = 4./(k1+k2);
NUT_sample = (ET_sample./GT_sample)/2-1;

fprintf('\nlambda_1 = %.4f',lambda(1));
fprintf('\nlambda_2 = %.4f',lambda(2));
fprintf('\nlambda_3 = %.4f',lambda(3));
fprintf('\nlambda_4 = %.4f',lambda(4));
fprintf('\nlambda_5 = %.4f',lambda(5));
fprintf('\nlambda   = %.4f',lambda(6));
fprintf('\n');

a4 = 1-2*lambda(6);
b4 = 1/lambda(4);
a5 = 1-2*lambda(6);
b5 = 1/lambda(5);
fprintf('\nalpha_4 = %.4f',a4);
fprintf('\nbeta_4  = %.4f',b4);
fprintf('\n');

fprintf('\nalpha_5 = %.4f',a5);
fprintf('\nbeta_5  = %.4f',b5);
fprintf('\n');

%% Pdfs and cdfs
pdf_C4 = @(c4) gampdf(c4,a4,b4); % Gamma probability density function of C4
pdf_C5 = @(c5) gampdf(c5,a5,b5); % Gamma probability density function of C5
cdf_C4 = @(c4) gamcdf(c4,a4,b4); % Gamma cumulative density function of C4
cdf_C5 = @(c5) gamcdf(c5,a5,b5); % Gamma cumulative density function of C5

%% Display
if displaySolution
    %% Plot pdfs and cdfs
    N = 100;
    x4_min = max(0,mean(C4_data)-8*std(C4_data));
    x4_max = mean(C4_data)+8*std(C4_data);
    x5_min = max(0,mean(C5_data)-8*std(C5_data));
    x5_max = mean(C5_data)+8*std(C5_data);
    
    % Plot pdf of C4
    figure('Name','Probability density function of C4')
    clf
    x4 = linspace(x4_min,x4_max,N);
    y4 = pdf_C4(x4);
    plot(x4,y4,'-b');
    hold on
    plot(C4_data,pdf_C4(C4_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(x4),max(x4)])
    xlabel('$c_4$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_4}(c_4)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C4',formats);
    mymatlab2tikz(pathname,'pdf_C4.tex');
    
    % Plot pdf of C5
    figure('Name','Probability density function of C5')
    clf
    x5 = linspace(x5_min,x5_max,N);
    y5 = pdf_C5(x5);
    plot(x5,y5,'-b');
    hold on
    plot(C5_data,pdf_C5(C5_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(x5),max(x5)])
    xlabel('$c_5$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_5}(c_5)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C5',formats);
    mymatlab2tikz(pathname,'pdf_C5.tex');
    
    % Plot cdf of C4
    figure('name','cumulative distribution function of C4')
    clf
    z4 = cdf_C4(x4);
    plot(x4,z4,'-b');
    hold on
    plot(C4_data,cdf_C4(C4_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(x4),max(x4)])
    xlabel('$c_4$ [GPa]','Interpreter',interpreter)
    ylabel('$F_{C_4}(c_4)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_C4',formats);
    mymatlab2tikz(pathname,'cdf_C4.tex');
    
    % Plot cdf of C5
    figure('name','cumulative distribution function of C5')
    clf
    z5 = cdf_C5(x5);
    plot(x5,z5,'-b');
    hold on
    plot(C5_data,cdf_C5(C5_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(x5),max(x5)])
    xlabel('$c_5$ [GPa]','Interpreter',interpreter)
    ylabel('$F_{C_5}(c_5)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_C5',formats);
    mymatlab2tikz(pathname,'cdf_C5.tex');
    
    %% Plot samples
    figure('name','Samples of (ET,GL,NUT)')
    clf
    scatter3(ET_sample,GL_sample,NUT_sample,'b.')
    hold on
    scatter3(ET_data,GL_data,NUT_data,'r+')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$E^T$ [GPa]','Interpreter',interpreter);
    ylabel('$G^L$ [GPa]','Interpreter',interpreter);
    zlabel('$\nu^T$','Interpreter',interpreter);
%     xlabel('Young''s modulus $E^T$ [GPa]','Interpreter',interpreter);
%     ylabel('Shear modulus $G^L$ [GPa]','Interpreter',interpreter);
%     zlabel('Poisson''s ratio $\nu^T$','Interpreter',interpreter);
%     legend('samples','data');
    mysaveas(pathname,'samples_ET_GL_NUT',formats);
    mymatlab2tikz(pathname,'samples_ET_GL_NUT.tex');
end
