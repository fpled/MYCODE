%% Stochastic modeling of Young's modulus %%
%%----------------------------------------%%

% clc
clearvars
close all
% rng('default');

%% Input data
displaySolution = true;

filename = 'modelStoLinElasIsot_YoungModulus';
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
pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathnameIdentification,filenameAna));

E_data = mean_ET_data*1e-3; % GPa

mE_data = mean(E_data,1);
vE_data = var(E_data,0,1);
% vE_data = size(E_data,1)/(size(E_data,1)-1)*moment(E_data,2,1);
sE_data = sqrt(vE_data);
% sE_data = std(E_data,0,1);
dE_data = sE_data/mE_data;

%% Maximum likelihood estimation
phat = gamfit(E_data);
% phat = mle(E_data,'distribution','gam');
% nloglf = @(phat,data,cens,freq) length(data)*gammaln(a)...
%     +length(data)*a*log(b)...
%     +(1-a)*sum(log(data))...
%     +1/b*sum(data);
% phat = mle(E_data,'nloglf',nloglf,'start',[2 0],'lowerbound',[2 0]);

a = phat(1);
b = phat(2);

mE = a*b;
vE = a*b^2;
sE = sqrt(vE);
dE = sE/mE; % dE = 1/sqrt(a);

%% Pdf and cdf
pdf_E = @(x) gampdf(x,a,b);
% pdf_E = @(x) pdf('gam',x,a,b);
% pdf_E = @(x) 1/(b^a*gamma(a))*(x.^(a-1)).*exp(-x./b);
cdf_E = @(x) gamcdf(x,a,b);

%% Sample generation
N = 1e3; % number of samples
E_sample = gamrnd(a,b,N,1);
% u = randn(N,1);
% E_sample = gaminv(normcdf(u),a,b);

mE_sample = mean(E_sample,1);
vE_sample = var(E_sample,0,1);
% vE_sample = size(E_sample,1)/(size(E_sample,1)-1)*moment(E_sample,2,1);
sE_sample = sqrt(vE_sample);
% sE_sample = std(E_sample,0,1);
dE_sample = sE_sample/mE_sample;

%% Outputs
fprintf('\nnb data = %g',length(E_data));
fprintf('\nnb sample = %g',N);
fprintf('\n');

fprintf('\nalpha = %.4f',a);
fprintf('\nbeta  = %.4f',b);
fprintf('\n');

fprintf('\nmean(E)        = %.4f GPa',mE);
fprintf('\nmean(E_sample) = %.4f GPa',mE_sample);
fprintf('\nmean(E_data)   = %.4f GPa',mE_data);
fprintf('\nvar(E)         = %.4f (GPa)^2',vE);
fprintf('\nvar(E_sample)  = %.4f (GPa)^2',vE_sample);
fprintf('\nvar(E_data)    = %.4f (GPa)^2',vE_data);
fprintf('\nstd(E)         = %.4f GPa',sE);
fprintf('\nstd(E_sample)  = %.4f GPa',sE_sample);
fprintf('\nstd(E_data)    = %.4f GPa',sE_data);
fprintf('\ndisp(E)        = %.4f',dE);
fprintf('\ndisp(E_sample) = %.4f',dE_sample);
fprintf('\ndisp(E_data)   = %.4f',dE_data);
fprintf('\n');
    
alpha = 1/2;
mse = alpha * (mE - mE_data)^2/(mE_data)^2 + (1-alpha) * (dE - dE_data)^2/(dE_data)^2;
mse_sample = alpha * (mE_sample - mE_data)^2/(mE_data)^2 + (1-alpha) * (dE_sample - dE_data)^2/(dE_data)^2;
fprintf('\nmean-squared error mse        = %.4e',mse);
fprintf('\nmean-squared error mse_sample = %.4e',mse_sample);
fprintf('\n');

%% Display
if displaySolution
    %% Plot data
    figure('Name','Data')
    clf
    bar(1:length(E_data),E_data)
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[0,length(E_data)+1])
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Young''s modulus $E$ [GPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_E',formats);
    mymatlab2tikz(pathname,'data_E.tex');
    
    %% Plot pdf and cdf
    xmin = max(0,mE-5*sE);
    xmax = mE+5*sE;
    x = linspace(xmin,xmax,1e3);
    
    % Plot pdf of E
    figure('Name','Probability density function')
    clf
    hold on
    plot(x,pdf_E(x),'-b','LineWidth',linewidth);
    plot(E_data,pdf_E(E_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xmin,xmax])
    xlabel('$e$ [GPa]','Interpreter',interpreter)
    ylabel('$p_E(e)$','Interpreter',interpreter)
    % l = legend('$p_E(e)$','$(e_i,p_E(e_i))_{i=1}^n$');
    % set(l,'Interpreter',interpreter,'Location','NorthWest');
    mysaveas(pathname,'pdf_E',formats);
    mymatlab2tikz(pathname,'pdf_E.tex',...
        'extraAxisOptions',{'ylabel style={overlay}'});
    
    % Plot cdf of E
    figure('Name','Cumulative distribution function')
    clf
    hold on
    plot(x,cdf_E(x),'-b','LineWidth',linewidth);
    plot(E_data,cdf_E(E_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xmin,xmax])
    set(gca,'YLim',[0,1])
    xlabel('$e$ [GPa]','Interpreter',interpreter)
    ylabel('$F_E(e)$','Interpreter',interpreter)
    % l = legend('$F_E(e)$','$(e_i,F_E(e_i))_{i=1}^n$');
    % set(l,'Interpreter',interpreter,'Location','NorthWest');
    mysaveas(pathname,'cdf_E',formats);
    mymatlab2tikz(pathname,'cdf_E.tex',...
        'extraAxisOptions',{'ylabel style={overlay}'});
    
    %% Plot samples
    figure('Name','Samples')
    clf
    scatter(1:N,E_sample,'b.')
    hold on
    plot([1 N],[mE mE],'-r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Young''s modulus $E$ [GPa]','Interpreter',interpreter)
    legend('samples','mean');
    mysaveas(pathname,'samples_E',formats);
    mymatlab2tikz(pathname,'samples_E.tex');
    
end
