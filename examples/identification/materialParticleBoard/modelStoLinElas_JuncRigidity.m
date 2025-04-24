%% Stochastic modeling of junction rigidity (bending stiffness) %%
%%--------------------------------------------------------------%%

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
% mKS_data = mean(KS_data);
% vKS_data = var(KS_data);
% sKS_data = std(KS_data);
% dKS_data = sKS_data/mKS_data;
% mKD_data = mean(KD_data);
% vKD_data = var(KD_data);
% sKD_data = std(KD_data);
% dKD_data = sKD_data/mKD_data;
fprintf('\nnb data = %d for screw junction',length(KS_data));
fprintf('\nnb data = %d for dowel junction',length(KD_data));
fprintf('\n');

%% Maximum likelihood estimation
paramS = gamfit(KS_data);
paramD = gamfit(KD_data);

% paramS = mle(KS_data,'pdf',@gampdf,'Start',[3 1],'LowerBound',[2 0]);
% paramD = mle(KD_data,'pdf',@gampdf,'Start',[3 1],'LowerBound',[2 0]);

aS = paramS(1);
bS = paramS(2);
aD = paramD(1);
bD = paramD(2);
fprintf('\nalpha = %.4f for screw junction',aS);
fprintf('\nbeta = %.4f for screw junction',bS);
fprintf('\n');
fprintf('\nalpha = %.4f for screw junction',aD);
fprintf('\nbeta = %.4f for screw junction',bD);
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
fprintf('\n');

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
KS_sample = gamrnd(aS,bS,N,1); % [kN/rad]
KD_sample = gamrnd(aD,bD,N,1); % [kN/rad]

%% Display
if displaySolution
    %% Plot data
    figure('Name','Data for KS')
    clf
    bar(KS_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Bending  stiffness per unit length $k_S$ [kN/rad]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Rigidit\''e lin\''eique en flexion $k_S$ [kN/rad]','Interpreter',interpreter);
    mysaveas(pathname,'data_KS',formats);
    mymatlab2tikz(pathname,'data_KS.tex');
    
    figure('Name','Data for KD')
    clf
    bar(KD_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Bending stiffness per unit length $k_D$ [kN/rad]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Rigidit\''e lin\''eique en flexion $k_D$ [kN/rad]','Interpreter',interpreter);
    mysaveas(pathname,'data_KD',formats);
    mymatlab2tikz(pathname,'data_KD.tex');
    
    %% Plot log-likelihood functions
    aS_series = linspace(aS*0.5,aS*1.5,1e2);
    bS_series = linspace(bS*0.5,bS*1.5,1e2);
    loglfS = zeros(length(aS_series),length(bS_series));
    for i=1:length(aS_series)
        aS_i = aS_series(i);
        for j=1:length(bS_series)
            bS_j = bS_series(j);
            paramS_ij = [aS_i,bS_j];
            loglfS(j,i) = -gamlike(paramS_ij,KS_data);
        end
    end
    
    aD_series = linspace(aD*0.5,aD*1.5,1e2);
    bD_series = linspace(bD*0.5,bD*1.5,1e2);
    loglfD = zeros(length(aD_series),length(bD_series));
    for i=1:length(aD_series)
        aD_i = aD_series(i);
        for j=1:length(bD_series)
            bD_j = bD_series(j);
            paramD_ij = [aD_i,bD_j];
            loglfD(j,i) = -gamlike(paramD_ij,KD_data);
        end
    end
    
    % Plot log-likelihood function loglf for KS
    figure('Name','Surface plot: Log-likelihood function for KS')
    clf
    surfc(aS_series,bS_series,loglfS,'EdgeColor','none');
    colorbar
    hold on
    scatter3(aS,bS,-gamlike(paramS,KS_data),'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    xlabel('$\alpha$','Interpreter',interpreter)
    ylabel('$\beta$ [kN/rad]','Interpreter',interpreter)
    zlabel('$\mathcal{L}_S(\alpha,\beta)$','Interpreter',interpreter)
    mysaveas(pathname,'loglf_KS_3D',formats);
    % mymatlab2tikz(pathname,'loglf_KS_3D.tex');
    
    figure('Name','Contour plot: Log-likelihood function for KS')
    clf
    contourf(aS_series,bS_series,loglfS,30);
    colorbar
    hold on
    scatter(aS,bS,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    xlabel('$\alpha$','Interpreter',interpreter)
    ylabel('$\beta$ [kN/rad]','Interpreter',interpreter)
    zlabel('$\mathcal{L}_S(\alpha,\beta)$','Interpreter',interpreter)
    mysaveas(pathname,'loglf_KS_2D',formats);
    % mymatlab2tikz(pathname,'loglf_KS_2D.tex');
    
    % Plot log-likelihood function loglf for KD
    figure('Name','Surface plot: Log-likelihood function for KD')
    clf
    surfc(aD_series,bD_series,loglfD,'EdgeColor','none');
    colorbar
    hold on
    scatter3(aD,bD,-gamlike(paramD,KD_data),'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    xlabel('$\alpha$','Interpreter',interpreter)
    ylabel('$\beta$ [kN/rad]','Interpreter',interpreter)
    zlabel('$\mathcal{L}_D(\alpha,\beta)$','Interpreter',interpreter)
    mysaveas(pathname,'loglf_KD_3D',formats);
    % mymatlab2tikz(pathname,'loglf_KD_3D.tex');
    
    figure('Name','Contour plot: Log-likelihood function for KD')
    clf
    contourf(aD_series,bD_series,loglfD,30);
    colorbar
    hold on
    scatter(aD,bD,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    xlabel('$\alpha$','Interpreter',interpreter)
    ylabel('$\beta$ [kN/rad]','Interpreter',interpreter)
    zlabel('$\mathcal{L}_D(\alpha,\beta)$','Interpreter',interpreter)
    mysaveas(pathname,'loglf_KD_2D',formats);
    % mymatlab2tikz(pathname,'loglf_KD_2D.tex');
    
    %% Plot pdfs and cdfs
    xminS = max(0,mKS-5*sKS);
    xmaxS = mKS+5*sKS;
    xS = linspace(xminS,xmaxS,1e3);
    
    xminD = max(0,mKD-5*sKD);
    xmaxD = mKD+5*sKD;
    xD = linspace(xminD,xmaxD,1e3);
    
    % Plot pdf of KS
    figure('Name','Probability density function of KS')
    clf
    plot(xS,pdf_KS(xS),'-b','LineWidth',linewidth);
    % hold on
    % plot(KS_data,pdf_KS(KS_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xminS,xmaxS])
    xlabel('$k$ [kN/rad]','Interpreter',interpreter)
    ylabel('$p_{K_S}(k)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_KS',formats);
    mymatlab2tikz(pathname,'pdf_KS.tex');
    
    % Plot cdf of KS
    figure('Name','Cumulative distribution function of KS')
    clf
    plot(xS,cdf_KS(xS),'-r','LineWidth',linewidth);
    % hold on
    % plot(KS_data,cdf_KS(KS_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xminS,xmaxS])
    set(gca,'YLim',[0,1])
    xlabel('$k$ [kN/rad]','Interpreter',interpreter)
    ylabel('$F_{K_S}(k)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_KS',formats);
    mymatlab2tikz(pathname,'cdf_KS.tex');
    
    % Plot pdf of KD
    figure('Name','Probability density function of KD')
    clf
    plot(xD,pdf_KD(xD),'-b','LineWidth',linewidth);
    % hold on
    % plot(KD_data,pdf_KD(KD_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xminD,xmaxD])
    xlabel('$k$ [kN/rad]','Interpreter',interpreter)
    ylabel('$p_{K_D}(k)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_KD',formats);
    mymatlab2tikz(pathname,'pdf_KD.tex');
    
    % Plot cdf of KD
    figure('Name','Cumulative distribution function of KD')
    clf
    plot(xD,cdf_KD(xD),'-r','LineWidth',linewidth);
    % hold on
    % plot(KD_data,cdf_KD(KD_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xminD,xmaxD])
    set(gca,'YLim',[0,1])
    xlabel('$k$ [kN/rad]','Interpreter',interpreter)
    ylabel('$F_{K_D}(k)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_KD',formats);
    mymatlab2tikz(pathname,'cdf_KD.tex');
   
    %% Plot samples
    % Plot samples of KS
    figure('Name','Samples of KS')
    clf
    scatter(1:N,KS_sample,'b.')
    hold on
    plot([1 N],[mKS mKS],'-r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of realizations','Interpreter',interpreter)
    %ylabel('Bending stiffness per unit length $k_S$ [kN/rad]','Interpreter',interpreter)
    ylabel('$k_S$ [kN/rad]','Interpreter',interpreter)
    legend('realizations','mean value');
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Rigidit\''e lin\''eique en flexion $k_S$ [kN/rad]','Interpreter',interpreter)
    %legend('réalisations','valeur moyenne');
    mysaveas(pathname,'samples_KS',formats);
    % mymatlab2tikz(pathname,'samples_KS.tex');
    
    % Plot samples of KD
    figure('Name','Samples of KD')
    clf
    scatter(1:N,KD_sample,'b.')
    hold on
    plot([1 N],[mKD mKD],'-r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of realizations','Interpreter',interpreter)
    %ylabel('Bending stiffness per unit length $k_D$ [kN/rad]','Interpreter',interpreter)
    ylabel('$k_D$ [kN/rad]','Interpreter',interpreter)
    legend('realizations','mean value');
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Rigidit\''e lin\''eique en flexion $k_D$ [kN/rad]','Interpreter',interpreter)
    %legend('réalisations','valeur moyenne');
    mysaveas(pathname,'samples_KD',formats);
    % mymatlab2tikz(pathname,'samples_KD.tex');
    
    figure('name','Samples of (KS,KD)')
    clf
    scatter(KS_sample,KD_sample,'b.')
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$k_S$ [kN/rad]','Interpreter',interpreter);
    ylabel('$k_D$ [kN/rad]','Interpreter',interpreter);
    mysaveas(pathname,'samples_KS_KD',formats);
    % mymatlab2tikz(pathname,'samples_KS_KD.tex');
    
end
