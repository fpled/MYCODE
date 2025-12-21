%% Stochastic modeling of junction rigidity (bending stiffness) %%
%%--------------------------------------------------------------%%

% clc
clearvars
close all
rng('default');

%% Input data
displaySolution = true;

filename = 'modelStoLinElas_JunctionStiffness';
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

% Empirical estimates
NS_data = length(KS_data);
mKS_data = mean(KS_data);
vKS_data = var(KS_data);
sKS_data = std(KS_data);
dKS_data = sKS_data/mKS_data;

ND_data = length(KD_data);
mKD_data = mean(KD_data);
vKD_data = var(KD_data);
sKD_data = std(KD_data);
dKD_data = sKD_data/mKD_data;

%% Maximum likelihood estimation
% display = 'off'; % default for gamfit and mle
display = 'iter';
% display = 'final';

tolX = 1e-8; % tolerance on the parameter value (1e-8 by default for gamfit and 1e-6 for mle)
tolFun = 1e-8; % tolerance on the function value (1e-8 by default for gamfit and 1e-6 for mle)
tolBnd = 1e-6; % tolerance on the parameter bound (1e-6 by default for gamfit and mle)
maxIters = Inf; % maximum number of iterations
maxFunEvals = Inf; % maximum number of function evaluations

options_gamfit = statset(statset('gamfit'),'Display',display,...
    'TolX',tolX,'TolFun',tolFun,'TolBnd',tolBnd,...
    'MaxIter',maxIters,'MaxFunEvals',maxFunEvals);

options_mlecustom = statset(statset('mlecustom'),'Display',display,...
    'TolX',tolX,'TolFun',tolFun,'TolBnd',tolBnd,...
    'MaxIter',maxIters,'MaxFunEvals',maxFunEvals);

% paramS = gamfit(KS_data,[],[],[],options_gamfit);
% paramD = gamfit(KD_data,[],[],[],options_gamfit);
% paramS = mle(KS_data,'Distribution','gam','Options',options_gamfit);
% paramD = mle(KD_data,'Distribution','gam','Options',options_gamfit);

% Initial parameter values using Method of moments estimates
a0S = mKS_data^2/vKS_data; % aS > 2
b0S = mKS_data/a0S;        % bS > 0
a0D = mKD_data^2/vKD_data; % aD > 2
b0D = mKD_data/a0D;        % bD > 0

param0S = [a0S b0S]; % initial parameter vector for screw junction
param0D = [a0D b0D]; % initial parameter vector for dowel junction
lb      = [2 0];     % lower bounds

loglfval0S = -gamlike(param0S,KS_data); % initial log-likelihood for screw junction
loglfval0D = -gamlike(param0D,KD_data); % initial log-likelihood for dowel junction

paramS = mle(KS_data,'pdf',@gampdf,'Start',param0S,'LowerBound',lb,'Options',options_mlecustom);
paramD = mle(KD_data,'pdf',@gampdf,'Start',param0D,'LowerBound',lb,'Options',options_mlecustom);
% paramS = mle(KS_data,'nloglf',@gamlike,'Start',param0S,'LowerBound',lb,'Options',options_mlecustom);
% paramD = mle(KD_data,'nloglf',@gamlike,'Start',param0D,'LowerBound',lb,'Options',options_mlecustom);

% Optimal parameter values
aS = paramS(1); % aS > 2
bS = paramS(2); % bS > 0
aD = paramD(1); % aD > 2
bD = paramD(2); % bD > 0

[mKS,vKS] = gamstat(aS,bS);
% mKS = aS*bS;
% vKS = aS*bS^2;
sKS = sqrt(vKS); % sKS = sqrt(aS)*bS;
dKS = sKS/mKS; % dKS = 1/sqrt(aS);

[mKD,vKD] = gamstat(aD,bD);
% mKD = aD*bD;
% vKD = aD*bD^2;
sKD = sqrt(vKD); % sKD = sqrt(aD)*bD;
dKD = sKD/mKD; % dKD = 1/sqrt(aD);

loglfvalS = -gamlike(paramS,KS_data); % optimal (maximal) log-likelihood for screw junction
loglfvalD = -gamlike(paramD,KD_data); % optimal (maximal) log-likelihood for dowel junction

%% Pdf and cdf
pdf_KS = @(x) gampdf(x,aS,bS);
cdf_KS = @(x) gamcdf(x,aS,bS);
pdf_KD = @(x) gampdf(x,aD,bD);
cdf_KD = @(x) gamcdf(x,aD,bD);

%% Sample generation
N = 1e4; % number of samples
KS_sample = gamrnd(aS,bS,N,1); % [kN/rad]
KD_sample = gamrnd(aD,bD,N,1); % [kN/rad]

mKS_sample = mean(KS_sample);
vKS_sample = var(KS_sample);
sKS_sample = std(KS_sample);
dKS_sample = sKS_sample/mKS_sample;

mKD_sample = mean(KD_sample);
vKD_sample = var(KD_sample);
sKD_sample = std(KD_sample);
dKD_sample = sKD_sample/mKD_sample;

%% Outputs
filenameResults = fullfile(pathname,'results.txt');
fid = fopen(filenameResults,'w');
fprintf(fid,'nb data    = %g for screw junction\n',NS_data);
fprintf(fid,'nb data    = %g for dowel junction\n',ND_data);
fprintf(fid,'nb samples = %g\n',N);

fprintf(fid,'\n');
fprintf(fid,'Initial parameter values for screw junction\n');
fprintf(fid,'alphaS = %.4f\n',a0S);
fprintf(fid,'betaS  = %.6f\n',b0S);
fprintf(fid,'loglfS = %.4f\n',loglfval0S);

fprintf(fid,'\n');
fprintf(fid,'Optimal parameter values for screw junction\n');
fprintf(fid,'alphaS = %.4f\n',aS);
fprintf(fid,'betaS  = %.6f\n',bS);
fprintf(fid,'loglfS = %.4f\n',loglfvalS);

fprintf(fid,'\n');
fprintf(fid,'Initial parameter values for dowel junction\n');
fprintf(fid,'alphaD = %.4f\n',a0D);
fprintf(fid,'betaD  = %.6f\n',b0D);
fprintf(fid,'loglfD = %.4f\n',loglfval0D);

fprintf(fid,'\n');
fprintf(fid,'Optimal parameter values for dowel junction\n');
fprintf(fid,'alphaD = %.4f\n',aD);
fprintf(fid,'betaD  = %.6f\n',bD);
fprintf(fid,'loglfD = %.4f\n',loglfvalD);

% Screw junction
fprintf(fid,'\n');
fprintf(fid,'mean(KS_data)   = %.4f kN/rad\n',mKS_data);
fprintf(fid,'mean(KS)        = %.4f kN/rad\n',mKS);
fprintf(fid,'mean(KS_sample) = %.4f kN/rad\n',mKS_sample);
fprintf(fid,'var(KS_data)    = %.4f (kN/rad)^2\n',vKS_data);
fprintf(fid,'var(KS)         = %.4f (kN/rad)^2\n',vKS);
fprintf(fid,'var(KS_sample)  = %.4f (kN/rad)^2\n',vKS_sample);
fprintf(fid,'std(KS_data)    = %.4f kN/rad\n',sKS_data);
fprintf(fid,'std(KS)         = %.4f kN/rad\n',sKS);
fprintf(fid,'std(KS_sample)  = %.4f kN/rad\n',sKS_sample);
fprintf(fid,'cv(KS_data)     = %.4f\n',dKS_data);
fprintf(fid,'cv(KS)          = %.4f\n',dKS);
fprintf(fid,'cv(KS_sample)   = %.4f\n',dKS_sample);

err_meanKS = (mKS - mKS_data)/mKS_data;
err_varKS = (vKS - vKS_data)/vKS_data;
err_stdKS = (sKS - sKS_data)/sKS_data;
err_cvKS = (dKS - dKS_data)/dKS_data;

err_meanKS_sample = (mKS_sample - mKS_data)/mKS_data;
err_varKS_sample = (vKS_sample - vKS_data)/vKS_data;
err_stdKS_sample = (sKS_sample - sKS_data)/sKS_data;
err_cvKS_sample = (dKS_sample - dKS_data)/dKS_data;

alpha = 1/2;
mseKS = alpha * (mKS - mKS_data)^2/(mKS_data)^2 + (1-alpha) * (dKS - dKS_data)^2/(dKS_data)^2;
mseKS_sample = alpha * (mKS_sample - mKS_data)^2/(mKS_data)^2 + (1-alpha) * (dKS_sample - dKS_data)^2/(dKS_data)^2;

fprintf(fid,'\n');
fprintf(fid,'relative error on mean(KS) = %.4e\n',err_meanKS);
fprintf(fid,'relative error on var(KS)  = %.4e\n',err_varKS);
fprintf(fid,'relative error on std(KS)  = %.4e\n',err_stdKS);
fprintf(fid,'relative error on cv(KS)   = %.4e\n',err_cvKS);
fprintf(fid,'mean-squared error mse(KS) = %.4e\n',mseKS);

fprintf(fid,'\n');
fprintf(fid,'relative error on mean(KS_sample) = %.4e\n',err_meanKS_sample);
fprintf(fid,'relative error on var(KS_sample)  = %.4e\n',err_varKS_sample);
fprintf(fid,'relative error on std(KS_sample)  = %.4e\n',err_stdKS_sample);
fprintf(fid,'relative error on cv(KS_sample)   = %.4e\n',err_cvKS_sample);
fprintf(fid,'mean-squared error mse(KS_sample) = %.4e\n',mseKS_sample);

% Dowel junction
fprintf(fid,'\n');
fprintf(fid,'mean(KD_data)   = %.4f kN/rad\n',mKD_data);
fprintf(fid,'mean(KD)        = %.4f kN/rad\n',mKD);
fprintf(fid,'mean(KD_sample) = %.4f kN/rad\n',mKD_sample);
fprintf(fid,'var(KD_data)    = %.4f (kN/rad)^2\n',vKD_data);
fprintf(fid,'var(KD)         = %.4f (kN/rad)^2\n',vKD);
fprintf(fid,'var(KD_sample)  = %.4f (kN/rad)^2\n',vKD_sample);
fprintf(fid,'std(KD_data)    = %.4f kN/rad\n',sKD_data);
fprintf(fid,'std(KD)         = %.4f kN/rad\n',sKD);
fprintf(fid,'std(KD_sample)  = %.4f kN/rad\n',sKD_sample);
fprintf(fid,'cv(KD_data)     = %.4f\n',dKD_data);
fprintf(fid,'cv(KD)          = %.4f\n',dKD);
fprintf(fid,'cv(KD_sample)   = %.4f\n',dKD_sample);

err_meanKD = (mKD - mKD_data)/mKD_data;
err_varKD = (vKD - vKD_data)/vKD_data;
err_stdKD = (sKD - sKD_data)/sKD_data;
err_cvKD = (dKD - dKD_data)/dKD_data;

err_meanKD_sample = (mKD_sample - mKD_data)/mKD_data;
err_varKD_sample = (vKD_sample - vKD_data)/vKD_data;
err_stdKD_sample = (sKD_sample - sKD_data)/sKD_data;
err_cvKD_sample = (dKD_sample - dKD_data)/dKD_data;

alpha = 1/2;
mseKD = alpha * (mKD - mKD_data)^2/(mKD_data)^2 + (1-alpha) * (dKD - dKD_data)^2/(dKD_data)^2;
mseKD_sample = alpha * (mKD_sample - mKD_data)^2/(mKD_data)^2 + (1-alpha) * (dKD_sample - dKD_data)^2/(dKD_data)^2;

fprintf(fid,'\n');
fprintf(fid,'relative error on mean(KD) = %.4e\n',err_meanKD);
fprintf(fid,'relative error on var(KD)  = %.4e\n',err_varKD);
fprintf(fid,'relative error on std(KD)  = %.4e\n',err_stdKD);
fprintf(fid,'relative error on cv(KD)   = %.4e\n',err_cvKD);
fprintf(fid,'mean-squared error mse(KD) = %.4e\n',mseKD);

fprintf(fid,'\n');
fprintf(fid,'relative error on mean(KD_sample) = %.4e\n',err_meanKD_sample);
fprintf(fid,'relative error on var(KD_sample)  = %.4e\n',err_varKD_sample);
fprintf(fid,'relative error on std(KD_sample)  = %.4e\n',err_stdKD_sample);
fprintf(fid,'relative error on cv(KD_sample)   = %.4e\n',err_cvKD_sample);
fprintf(fid,'mean-squared error mse(KD_sample) = %.4e\n',mseKD_sample);
fclose(fid);
type(filenameResults) % fprintf('%s', fileread(filenameResults))

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
    loglf = @(A,B,data) arrayfun(@(a,b) -gamlike([a,b],data), A, B);
    
    aS_series = linspace(aS*0.5,aS*1.5,1e2);
    bS_series = linspace(bS*0.5,bS*1.5,1e2);
    % logLS = zeros(length(aS_series),length(bS_series));
    % for i=1:length(aS_series)
    %     aS_i = aS_series(i);
    %     for j=1:length(bS_series)
    %         bS_j = bS_series(j);
    %         paramS_ij = [aS_i,bS_j];
    %         logLS(j,i) = -gamlike(paramS_ij,KS_data);
    %     end
    % end
    [AS,BS] = meshgrid(aS_series,bS_series);
    logLS = loglf(AS,BS,KS_data);
    
    aD_series = linspace(aD*0.5,aD*1.5,1e2);
    bD_series = linspace(bD*0.5,bD*1.5,1e2);
    % logLD = zeros(length(aD_series),length(bD_series));
    % for i=1:length(aD_series)
    %     aD_i = aD_series(i);
    %     for j=1:length(bD_series)
    %         bD_j = bD_series(j);
    %         paramD_ij = [aD_i,bD_j];
    %         logLD(j,i) = -gamlike(paramD_ij,KD_data);
    %     end
    % end
    [AD,BD] = meshgrid(aD_series,bD_series);
    logLD = loglf(AD,BD,KD_data);
    
    % Plot log-likelihood function loglf for KS
    figure('Name','Surface plot: Log-likelihood function for KS')
    clf
    surfc(aS_series,bS_series,logLS,'EdgeColor','none');
    colorbar
    hold on
    scatter3(aS,bS,loglfvalS,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    xlabel('$\alpha$','Interpreter',interpreter)
    ylabel('$\beta$ [kN/rad]','Interpreter',interpreter)
    zlabel('$\mathcal{L}_S(\alpha,\beta)$','Interpreter',interpreter)
    mysaveas(pathname,'loglf_KS_3D',formats);
    % mymatlab2tikz(pathname,'loglf_KS_3D.tex');
    
    figure('Name','Contour plot: Log-likelihood function for KS')
    clf
    contourf(aS_series,bS_series,logLS,50);
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
    surfc(aD_series,bD_series,logLD,'EdgeColor','none');
    colorbar
    hold on
    scatter3(aD,bD,loglfvalD,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    xlabel('$\alpha$','Interpreter',interpreter)
    ylabel('$\beta$ [kN/rad]','Interpreter',interpreter)
    zlabel('$\mathcal{L}_D(\alpha,\beta)$','Interpreter',interpreter)
    mysaveas(pathname,'loglf_KD_3D',formats);
    % mymatlab2tikz(pathname,'loglf_KD_3D.tex');
    
    figure('Name','Contour plot: Log-likelihood function for KD')
    clf
    contourf(aD_series,bD_series,logLD,50);
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
    hold on
    %scatter(KS_data,KD_data,'r+','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$k_S$ [kN/rad]','Interpreter',interpreter);
    ylabel('$k_D$ [kN/rad]','Interpreter',interpreter);
    %legend('realizations','data')
    %legend('réalisations','données')
    mysaveas(pathname,'samples_KS_KD',formats);
    % mymatlab2tikz(pathname,'samples_KS_KD.tex');
    
end
