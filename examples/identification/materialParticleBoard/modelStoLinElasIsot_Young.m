%% Stochastic modeling of Young modulus %%
%%--------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');

%% Input data
displaySolution = true;

filename = 'modelStoElasIsotYoung';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc2'};

%% Data
filenameAna = 'data_ET_GL.mat';
pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathnameIdentification,filenameAna));

sample = 'C';
E_data = zeros(1,20);
for j=1:20
    sampleNum = [sample num2str(j)];
    E_data(j) = eval(['mean_ET_' sampleNum '_data;']); % GPa
end
mean_E = mean(E_data);
std_E = std(E_data);

%% Plot data
if displaySolution
    figure('Name','Data')
    clf
    bar(1:length(E_data),E_data)
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[0,length(E_data)+1])
    xlabel('Sample','Interpreter',interpreter);
    ylabel('Young modulus $E$ (GPa)','Interpreter',interpreter);
    mysaveas(pathname,'data_E',formats);
    mymatlab2tikz(pathname,'data_E.tex');
end

%% Maximum likelihood estimation
phat = gamfit(E_data);
% phat = mle(E_data,'distribution','gam');
% nloglf = @(phat,data,cens,freq) length(data)*gammaln(phat(1))...
%     +length(data)*phat(1)*log(phat(2))...
%     +(1-phat(1))*sum(log(data))...
%     +1/phat(2)*sum(data);
% phat = mle(E_data,'nloglf',nloglf,'start',[2 0],'lowerbound',[2 0]);

fprintf('\nnb samples = %d',length(E_data));
fprintf('\nalpha = %.4f',phat(1));
fprintf('\nbeta  = %.4f',phat(2));

mE = phat(1)*phat(2);
vE = phat(1)*phat(2)^2;
sE = sqrt(vE);
dE = sE/mE; % dE = 1/sqrt(phat(1));
fprintf('\nmean(E) = %.4f',mE);
fprintf('\nvar(E)  = %.4f',vE);
fprintf('\nstd(E)  = %.4f',sE);
fprintf('\ndisp(E) = %.4f',dE);
fprintf('\n');

%% Plot pdf and cdf
if displaySolution
    pdf_E = @(x) gampdf(x,phat(1),phat(2));
    % pdf_E = @(x) pdf('gam',x,phat(1),phat(2));
    % pdf_E = @(x) 1/(phat(2)^phat(1)*gamma(phat(1)))*(x.^(phat(1)-1)).*exp(-x./phat(2));
    cdf_E = @(x) gamcdf(x,phat(1),phat(2));
    
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
    xlabel('$e$ (GPa)','Interpreter',interpreter)
    ylabel('$p_E(e)$','Interpreter',interpreter)
    % l = legend('$p_E(e)$','$(e_i,p_E(e_i))_{i=1}^n$');
    % set(l,'Interpreter',interpreter,'Location','northwest');
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
    xlabel('$e$ (GPa)','Interpreter',interpreter)
    ylabel('$F_E(e)$','Interpreter',interpreter)
    % l = legend('$F_E(e)$','$(e_i,F_E(e_i))_{i=1}^n$');
    % set(l,'Interpreter',interpreter,'Location','northwest');
    mysaveas(pathname,'cdf_E',formats);
    mymatlab2tikz(pathname,'cdf_E.tex',...
        'extraAxisOptions',{'ylabel style={overlay}'});
end

%% Sample generation
N = 1e4; % number of samples
e = gamrnd(phat(1),phat(2),N,1);
% u = randn(N,1);
% e = gaminv(normcdf(u),phat(1),phat(2));

%% Plot samples
if displaySolution
    figure('Name','Samples')
    clf
    scatter(1:N,e,'b.')
    hold on
    plot([1 N],[mE mE],'-r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Young modulus $E$ (GPa)','Interpreter',interpreter)
    mysaveas(pathname,'samples_E',formats);
    mymatlab2tikz(pathname,'samples_E.tex');
end
