%% Stochastic modeling of junction's rigidity %%
%%----------------------------------------%%

% clc
clearvars
close all
rng('default');

%% Input data
displaySolution = true;

filename = 'modelStoLinElas_JunctionRigidity';
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
filenameScrew = 'data_KS.mat';
filenameDowel = 'data_KD.mat';
pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathnameIdentification,filenameScrew));
load(fullfile(pathnameIdentification,filenameDowel));

KS_data = mean_KS_data*1e-3; % [kN/rad]
KD_data = mean_KD_data*1e-3; % [kN/rad]
% mean_KS = mean(KS_data);
% std_KS = std(KS_data);
% mean_KD = mean(KD_data);
% std_KD = std(KD_data);
fprintf('\nnb data for Dowel junction = %d',length(KD_data));
fprintf('\nnb data for Screw junction = %d',length(KS_data));

%% Maximum likelihood estimation
phatS = gamfit(KS_data);
phatD = gamfit(KD_data);

aS = phatS(1);
bS = phatS(2);
aD = phatD(1);
bD = phatD(2);
fprintf('\nalpha for screw = %.4f',aS);
fprintf('\nbeta for screw = %.4f',bS);
fprintf('\nalpha for dowel = %.4f',aD);
fprintf('\nbeta for dowel = %.4f',bD);
fprintf('\n');

mKS = aS*bS;
vKS = aS*bS^2;
sKS = sqrt(vKS);
dKS = sKS/mKS; 

mKD = aD*bD;
vKD = aD*bD^2;
sKD = sqrt(vKD);
dKD = sKD/mKD; 

fprintf('\nmean(KS) = %.4f kN/rad',mKS);
fprintf('\nvar(KS)  = %.4f (kN/rad)^2',vKS);
fprintf('\nstd(KS)  = %.4f kN/rad',sKS);
fprintf('\ndisp(KS) = %.4f',dKS);

fprintf('\nmean(KD) = %.4f kN/rad',mKD);
fprintf('\nvar(KD)  = %.4f (kN/rad)^2',vKD);
fprintf('\nstd(KD)  = %.4f kN/rad',sKD);
fprintf('\ndisp(KD) = %.4f',dKD);
fprintf('\n');

%% Pdf and cdf
pdf_KS = @(x) gampdf(x,aS,bS);
cdf_KS = @(x) gamcdf(x,aS,bS);
pdf_KD = @(x) gampdf(x,aD,bD);
cdf_KD = @(x) gamcdf(x,aD,bD);

%% Sample generation
N = 1e4; % number of samples
kS = gamrnd(aS,bS,N,1); % [kN/rad]
kD = gamrnd(aD,bD,N,1); % [kN/rad]

%% Display
if displaySolution
    %% Plot data
    figure('Name','Data for KS')
    clf
    bar(KS_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Bending  stiffness per unit length $k^S$ [kN/rad]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Rigidit\''e lin\''eique en flexion $k^S$ [kN/rad]','Interpreter',interpreter);
    mysaveas(pathname,'data_KS',formats);
    mymatlab2tikz(pathname,'data_KS.tex');
    
    figure('Name','Data for KD')
    clf
    bar(KD_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Bending stiffness per unit length $k^D$ [kN/rad]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Rigidit\''e lin\''eique en flexion $k^D$ [kN/rad]','Interpreter',interpreter);
    mysaveas(pathname,'data_KD',formats);
    mymatlab2tikz(pathname,'data_KD.tex');
    
    %% Plot pdfs and cdfs
    xminD = max(0,mKD-5*sKD);
    xmaxD = mKD+5*sKD;
    xD = linspace(xminD,xmaxD,1e3);
    
    xminS = max(0,mKS-5*sKS);
    xmaxS = mKS+5*sKS;
    xS = linspace(xminS,xmaxS,1e3);
    
    % Plot pdf of KS
    figure('Name','Probability density function of KS')
    clf
    plot(xS,pdf_KS(xS),'-b','LineWidth',linewidth);
    %hold on
    %plot(KS_data,pdf_KS(KS_data),'k+','LineWidth',linewidth);
    %hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xminS,xmaxS])
    xlabel('$k$ [kN/rad]','Interpreter',interpreter)
    ylabel('$p_{K^S}(k)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_KS',formats);
    mymatlab2tikz(pathname,'pdf_KS.tex',...
        'extraAxisOptions',{'ylabel style={overlay}'});
    
    % Plot cdf of KS
    figure('Name','Cumulative distribution function of KS')
    clf
    plot(xS,cdf_KS(xS),'-r','LineWidth',linewidth);
    %hold on
    %plot(KS_data,cdf_KS(KS_data),'k+','LineWidth',linewidth);
    %hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xminS,xmaxS])
    set(gca,'YLim',[0,1])
    xlabel('$k$ [kN/rad]','Interpreter',interpreter)
    ylabel('$F_{K^S}(k)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_KS',formats);
    mymatlab2tikz(pathname,'cdf_KS.tex',...
        'extraAxisOptions',{'ylabel style={overlay}'});
    
    % Plot pdf of KD
    figure('Name','Probability density function of KD')
    clf
    plot(xD,pdf_KD(xD),'-b','LineWidth',linewidth);
    %hold on
    %plot(KD_data,pdf_KD(KD_data),'k+','LineWidth',linewidth);
    %hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xminD,xmaxD])
    xlabel('$k$ [kN/rad]','Interpreter',interpreter)
    ylabel('$p_{K^D}(k)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_KD',formats);
    mymatlab2tikz(pathname,'pdf_KD.tex',...
        'extraAxisOptions',{'ylabel style={overlay}'});
    
    % Plot cdf of KD
    figure('Name','Cumulative distribution function of KD')
    clf
    plot(xD,cdf_KD(xD),'-r','LineWidth',linewidth);
    %hold on
    %plot(KD_data,cdf_KD(KD_data),'k+','LineWidth',linewidth);
    %hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xminD,xmaxD])
    set(gca,'YLim',[0,1])
    xlabel('$k$ [kN/rad]','Interpreter',interpreter)
    ylabel('$F_{K^D}(k)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_KD',formats);
    mymatlab2tikz(pathname,'cdf_KD.tex',...
        'extraAxisOptions',{'ylabel style={overlay}'});
   
    %% Plot samples
    figure('Name','Samples of KS')
    clf
    scatter(1:N,kS,'b.')
    hold on
    plot([1 N],[mKS mKS],'-r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of realizations','Interpreter',interpreter)
    ylabel('Bending stiffness per unit length $k^S$ [kN/rad]','Interpreter',interpreter)
    legend('samples','mean value');
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Rigidit\''e lin\''eique en flexion $k^S$ [kN/rad]','Interpreter',interpreter)
    %legend('réalisations','valeur moyenne');
    mysaveas(pathname,'samples_KS',formats);
    mymatlab2tikz(pathname,'samples_KS.tex');
    
    figure('Name','Samples of KD')
    clf
    scatter(1:N,kD,'b.')
    hold on
    plot([1 N],[mKD mKD],'-r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of realizations','Interpreter',interpreter)
    ylabel('Bending stiffness per unit length $k^D$ [kN/rad]','Interpreter',interpreter)
    legend('samples','mean value');
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Rigidit\''e lin\''eique en flexion $k^D$ [kN/rad]','Interpreter',interpreter)
    %legend('réalisations','valeur moyenne');
    mysaveas(pathname,'samples_KD',formats);
    mymatlab2tikz(pathname,'samples_KD.tex');
    
end
