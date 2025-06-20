%% Stochastic modeling of random linear elasticity tensor %%
%  with transversely isotropic symmetry                   %%
%%--------------------------------------------------------%%

% clc
clearvars
close all
rng('default');

%% Input data
displaySolution = true;

MCMCalg = 'MH'; % algorithm for Markov-Chain Monte Carlo (MCMC) method = 'MH', 'BUM', 'CUM' or 'SS'
filename = ['modelStoLinElasIsotTrans_ElasTensor_' MCMCalg];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};

%% Data
filenameAna = 'data_ET_GL.mat';
filenameNum = 'data_EL_NUL.mat';
pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathnameIdentification,filenameAna));
load(fullfile(pathnameIdentification,filenameNum));

ET_data = mean_ET_data*1e-3; % [GPa]
GL_data = mean_GL_data*1e-3; % [GPa]
EL_data = mean_EL_data*1e-3; % [GPa]
NUL_data = mean_NUL_data;
NUT_data = 0.1+0.2*rand(length(mean_ET_data),1); % artificial data for NUT varying from 0.1 to 0.3
GT_data = ET_data./(2*(1+NUT_data)); % [GPa]
kT_data = (EL_data.*ET_data)./(2*(1-NUT_data).*EL_data-4*ET_data.*(NUL_data).^2); % [GPa]
C1_data = EL_data + 4*(NUL_data.^2).*kT_data; % [GPa]
C2_data = 2*kT_data; % [GPa]
C3_data = 2*sqrt(2)*kT_data.*NUL_data; % [GPa]
C4_data = 2*GT_data; % [GPa]
C5_data = 2*GL_data; % [GPa]
C_data = [C1_data(:) C2_data(:) C3_data(:) C4_data(:) C5_data(:)];

mC_data = mean(C_data,1);
vC_data = var(C_data,0,1);
% vC_data = size(C_data,1)/(size(C_data,1)-1)*moment(C_data,2,1);
sC_data = sqrt(norm(vC_data));
dC_data = sC_data/norm(mC_data);
phiC_data = log((C_data(:,1).*C_data(:,2)-C_data(:,3).^2).*(C_data(:,4).^2).*(C_data(:,5).^2));
nuC_data = mean(phiC_data,1);

%% Least-Squares estimation for computing Lagrange multipliers
lambda = lseStoLinElasTensorIsotTrans(C_data,MCMCalg);

la1 = lambda(1);
la2 = lambda(2);
la3 = lambda(3);
la4 = lambda(4);
la5 = lambda(5);
la  = lambda(6);

a = 1-2*la;
b4 = 1/la4;
b5 = 1/la5;

%% Pdfs and cdfs
pdf_C4 = @(c4) gampdf(c4,a,b4); % Gamma probability density function of C4
pdf_C5 = @(c5) gampdf(c5,a,b5); % Gamma probability density function of C5
cdf_C4 = @(c4) gamcdf(c4,a,b4); % Gamma cumulative density function of C4
cdf_C5 = @(c5) gamcdf(c5,a,b5); % Gamma cumulative density function of C5

%% Sample generation
N = 1e4; % number of samples
switch lower(MCMCalg)
    case 'mh'
        C_sample(:,1:3) = mhsampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),N);
    case 'bum'
        C_sample(:,1:3) = mhsampleStoLinElasTensorIsotTrans_BUM(lambda,C_data(:,1:3),N);
    case 'cum'
        C_sample(:,1:3) = mhsampleStoLinElasTensorIsotTrans_CUM(lambda,C_data(:,1:3),N);
    case 'ss'
        C_sample = slicesampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),N);
    otherwise
        error(['MCMC algorithm ' MCMC ' not implemented'])
end
C_sample(:,4) = gamrnd(a,b4,N,1);
C_sample(:,5) = gamrnd(a,b5,N,1);

kT_sample = C_sample(:,2)/2; % [GPa]
NUL_sample = (C_sample(:,3)./kT_sample)/(2*sqrt(2));
EL_sample = C_sample(:,1) - 4*(NUL_sample.^2).*kT_sample; % [GPa]
GT_sample = C_sample(:,4)/2; % [GPa]
GL_sample = C_sample(:,5)/2; % [GPa]
ET_sample = 4./(1./kT_sample+1./GT_sample+4*(NUL_sample.^2)./EL_sample); % [GPa]
NUT_sample = (ET_sample./GT_sample)/2-1;

mC_sample = mean(C_sample,1);
vC_sample = var(C_sample,0,1);
% vC_sample = size(C_sample,1)/(size(C_sample,1)-1)*moment(C_sample,2,1);
sC_sample = sqrt(norm(vC_sample));
dC_sample = sC_sample/norm(mC_sample);
phiC_sample = log((C_sample(:,1).*C_sample(:,2)-C_sample(:,3).^2).*(C_sample(:,4).^2).*(C_sample(:,5).^2));
nuC_sample = mean(phiC_sample,1);

%% Ouputs
fprintf('\nnb data   = %g',size(C_data,1));
fprintf('\nnb sample = %g',N);
fprintf('\n');

fprintf('\nlambda_1 = %.4f',la1);
fprintf('\nlambda_2 = %.4f',la2);
fprintf('\nlambda_3 = %.4f',la3);
fprintf('\nlambda_4 = %.4f',la4);
fprintf('\nlambda_5 = %.4f',la5);
fprintf('\nlambda   = %.4f',la);
fprintf('\n');
fprintf('\nalpha_4 = %.4f',a);
fprintf('\nbeta_4  = %.4f',b4);
fprintf('\n');
fprintf('\nalpha_5 = %.4f',a);
fprintf('\nbeta_5  = %.4f',b5);
fprintf('\n');

for i=1:5
    fprintf('\nmean(C%u_sample) = %.4f GPa',i,mC_sample(i));
    fprintf('\nmean(C%u_data)   = %.4f GPa',i,mC_data(i));
    fprintf('\nvar(C%u_sample)  = %.4f (GPa)^2',i,vC_sample(i));
    fprintf('\nvar(C%u_data)    = %.4f (GPa)^2',i,vC_data(i));
    fprintf('\nstd(C%u_sample)  = %.4f GPa',i,sqrt(vC_sample(i)));
    fprintf('\nstd(C%u_data)    = %.4f GPa',i,sqrt(vC_data(i)));
    fprintf('\ndisp(C%u_sample) = %.4f',i,sqrt(vC_sample(i))/mC_sample(i));
    fprintf('\ndisp(C%u_data)   = %.4f',i,sqrt(vC_data(i))/mC_data(i));
    fprintf('\n');
end
fprintf('\nnu(C_sample) = mean(log(det([C_sample])) = %.4f',nuC_sample);
fprintf('\nnu(C_data)   = mean(log(det([C_data]))   = %.4f',nuC_data);

err_mean_sample = norm(mC_sample - mC_data)^2/norm(mC_data)^2;
err_nu_sample = (nuC_sample - nuC_data)^2/(nuC_data)^2;
fprintf('\nerror on mean(C_sample) = %.4e',err_mean_sample);
fprintf('\nerror on mean(log(det([C_sample])) = %.4e',err_nu_sample);
fprintf('\n');

%% Display
if displaySolution
    %% Plot data
    figure('Name','Data for ET')
    clf
    bar(ET_data);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Transverse Young''s modulus $E_T$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module d''Young transverse $E_T$ [GPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_ET',formats);
    mymatlab2tikz(pathname,'data_ET.tex');
    
    figure('Name','Data for EL')
    clf
    bar(EL_data*1e3);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Longitudinal Young''s modulus $E_L$ [MPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module d''Young longitudinal $E_L$ [MPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_EL',formats);
    mymatlab2tikz(pathname,'data_EL.tex');
    
    figure('Name','Data for GL')
    clf
    bar(GL_data*1e3);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Longitudinal shear modulus $G_L$ [MPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module de cisaillement longitudinal $G_L$ [MPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_G','fig');
    mymatlab2tikz(pathname,'data_G.tex');
    
    figure('Name','Data for NUL')
    clf
    bar(NUL_data);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Longitudinal Poisson''s ratio $\nu_L$','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Coefficient de Poisson longitudinal $\nu_L$','Interpreter',interpreter);
    mysaveas(pathname,'data_NUL',formats);
    mymatlab2tikz(pathname,'data_NUL.tex');
    
    figure('Name','Data for NUT')
    clf
    bar(NUT_data);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Transverse Poisson''s ratio $\nu_T$','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Coefficient de Poisson transverse $\nu_T$','Interpreter',interpreter);
    mysaveas(pathname,'data_NUL',formats);
    mymatlab2tikz(pathname,'data_NUL.tex');
    
    %% Plot pdfs and cdfs
    N = 1e3;
    [f1,xi1,bw1] = ksdensity(C_sample(:,1),'NumPoints',N);
    [f2,xi2,bw2] = ksdensity(C_sample(:,2),'NumPoints',N);
    [f3,xi3,bw3] = ksdensity(C_sample(:,3),'NumPoints',N);
    [f4,xi4,bw4] = ksdensity(C_sample(:,4),'NumPoints',N);
    [f5,xi5,bw5] = ksdensity(C_sample(:,5),'NumPoints',N);
    
    % Plot pdf estimate of C1
    figure('Name','Probability density estimate of C1')
    clf
    plot(xi1,f1,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(xi1),max(xi1)])
    xlabel('$c_1$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_1}(c_1)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C1_ksdensity',formats);
    mymatlab2tikz(pathname,'pdf_C1_ksdensity.tex');
    
    % Plot pdf estimate of C2
    figure('Name','Probability density estimate of C2')
    clf
    plot(xi2,f2,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(xi2),max(xi2)])
    xlabel('$c_2$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_2}(c_2)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C2_ksdensity',formats);
    mymatlab2tikz(pathname,'pdf_C2_ksdensity.tex');
    
    % Plot pdf estimate of C3
    figure('Name','Probability density estimate of C3')
    clf
    plot(xi3,f3,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(xi3),max(xi3)])
    xlabel('$c_3$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_3}(c_3)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C3_ksdensity',formats);
    mymatlab2tikz(pathname,'pdf_C3_ksdensity.tex');
    
    % Plot pdf estimate of C4
    figure('Name','Probability density estimate of C4')
    clf
    plot(xi4,f4,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(xi4),max(xi4)])
    xlabel('$c_4$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_4}(c_4)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C4_ksdensity',formats);
    mymatlab2tikz(pathname,'pdf_C4_ksdensity.tex');
    
    % Plot pdf estimate of C5
    figure('Name','Probability density estimate of C5')
    clf
    plot(xi5,f5,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(xi5),max(xi5)])
    xlabel('$c_5$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_5}(c_5)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C5_ksdensity',formats);
    mymatlab2tikz(pathname,'pdf_C5_ksdensity.tex');
    
    x4_min = max(0,mean(C_data(:,4))-5*std(C_data(:,4)));
    x4_max = mean(C_data(:,4))+5*std(C_data(:,4));
    x5_min = max(0,mean(C_data(:,5))-5*std(C_data(:,5)));
    x5_max = mean(C_data(:,5))+5*std(C_data(:,5));
    
    % Plot pdf of C4
    figure('Name','Probability density function of C4')
    clf
    x4 = linspace(x4_min,x4_max,N);
    y4 = pdf_C4(x4);
    plot(x4,y4,'-b','LineWidth',linewidth);
    % hold on
    % plot(C_data(:,4),pdf_C4(C_data(:,4)),'k+','LineWidth',linewidth);
    % hold off
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
    plot(x5,y5,'-b','LineWidth',linewidth);
    % hold on
    % plot(C_data(:,5),pdf_C5(C_data(:,5)),'k+','LineWidth',linewidth);
    % hold off
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
    plot(x4,z4,'-r','LineWidth',linewidth);
    % hold on
    % plot(C_data(:,4),cdf_C4(C_data(:,4)),'k+','LineWidth',linewidth);
    % hold off
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
    plot(x5,z5,'-r','LineWidth',linewidth);
    % hold on
    % plot(C_data(:,5),cdf_C5(C_data(:,5)),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(x5),max(x5)])
    xlabel('$c_5$ [GPa]','Interpreter',interpreter)
    ylabel('$F_{C_5}(c_5)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_C5',formats);
    mymatlab2tikz(pathname,'cdf_C5.tex');
    
    %% Plot samples
    % c1min = min(min(C_sample(:,1)),min(C_data(:,1)));
    c1min = 0;
    c1max = max(max(C_sample(:,1)),max(C_data(:,1)));
    c1 = linspace(c1min,c1max,1e2);
    
    % c2min = min(min(C_sample(:,2)),min(C_data(:,2)));
    c2min = 0;
    c2max = max(max(C_sample(:,2)),max(C_data(:,2)));
    c2 = linspace(c2min,c2max,1e2);
    
    [C1,C2] = meshgrid(c1,c2);
    C3 = sqrt(C1.*C2);
    
    figure('name','Samples of (C_1,C_2,C_3)')
    clf
    scatter3(C_sample(:,1),C_sample(:,2),C_sample(:,3),'b.')
    hold on
    scatter3(C_data(:,1),C_data(:,2),C_data(:,3),'r+','LineWidth',linewidth)
    surf(C1,C2,C3,'FaceAlpha',0.5)
    surf(C1,C2,-C3,'FaceAlpha',0.5)
    hold off
    grid on
    box on
    set(gca,'Xdir','reverse','Ydir','reverse','FontSize',fontsize)
    xlabel('$c_1$ [GPa]','Interpreter',interpreter);
    ylabel('$c_2$ [GPa]','Interpreter',interpreter);
    zlabel('$c_3$ [GPa]','Interpreter',interpreter);
    legend('realizations','data','support');
    %legend('réalisations','données','support');
    mysaveas(pathname,'samples_C1_C2_C3',formats);
    % mymatlab2tikz(pathname,'samples_C1_C2_C3.tex');
    
    figure('name','Samples of (C_4,C_5)')
    clf
    scatter(C_sample(:,4),C_sample(:,5),'b.')
    hold on
    scatter(C_data(:,4),C_data(:,5),'r+','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$c_4$ [GPa]','Interpreter',interpreter);
    ylabel('$c_5$ [GPa]','Interpreter',interpreter);
    legend('realizations','data');
    %legend('réalisations','données');
    mysaveas(pathname,'samples_C4_C5',formats);
    % mymatlab2tikz(pathname,'samples_C4_C5.tex');
    
    figure('name','Samples of (ET,GL,EL)')
    clf
    scatter3(ET_sample,GL_sample*1e3,EL_sample*1e3,'b.')
    hold on
    scatter3(ET_data,GL_data*1e3,EL_data*1e3,'r+','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$E_T$ [GPa]','Interpreter',interpreter);
    ylabel('$G_L$ [MPa]','Interpreter',interpreter);
    zlabel('$E_L$ [MPa]','Interpreter',interpreter);
    legend('realizations','data');
    %legend('réalisations','données');
    mysaveas(pathname,'samples_ET_GL_EL',formats);
    % mymatlab2tikz(pathname,'samples_ET_GL_EL.tex');

    figure('name','Samples of (NUL,NUT)')
    clf
    scatter(NUL_sample,NUT_sample,'b.')
    hold on
    scatter(NUL_data,NUT_data,'r+','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$\nu_L$','Interpreter',interpreter);
    ylabel('$\nu_T$','Interpreter',interpreter);
    legend('realizations','data');
    %legend('réalisations','données');
    mysaveas(pathname,'samples_NUL_NUT',formats);
    % mymatlab2tikz(pathname,'samples_NUL_NUT.tex');
    
end
