%% Stochastic modeling of junction's rigidity %%
%%----------------------------------------%%

% clc
clearvars
close all
% rng('default');

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
filenameAna = 'data_Kjunction.mat';
pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathnameIdentification,filenameAna));

Kdowel_data = mean_Kdowel_data*1e-3; 
Kscrew_data = mean_Kscrew_data*1e-4; 
% mean_Kdowel = mean(Kdowel_data);
% std_Kdowel = std(Kdowel_data);
% mean_Kscrew = mean(Kscrew_data);
% std_Kscrew = std(Kscrew_data);
fprintf('\nnb data for Dowel junction = %d',length(Kdowel_data));
fprintf('\nnb data for Screw junction = %d',length(Kscrew_data));

%% Maximum likelihood estimation
phat_dowel = gamfit(Kdowel_data);
phat_screw = gamfit(Kscrew_data);

a_dowel = phat_dowel(1);
b_dowel = phat_dowel(2);
a_screw = phat_screw(1);
b_screw = phat_screw(2);
fprintf('\nalpha for dowel = %.4f',a_dowel);
fprintf('\nbeta for dowel = %.4f',b_dowel);
fprintf('\nalpha for screw = %.4f',a_screw);
fprintf('\nbeta for screw = %.4f',b_screw);
fprintf('\n');

mKdowel = a_dowel*b_dowel;
vKdowel = a_dowel*b_dowel^2;
sKdowel = sqrt(vKdowel);
dKdowel = sKdowel/mKdowel; 

mKscrew = a_screw*b_screw;
vKscrew = a_screw*b_screw^2;
sKscrew = sqrt(vKscrew);
dKscrew = sKscrew/mKscrew; 

fprintf('\nmean(Kdowel) = %.4f 10^3N/rad',mKdowel);
fprintf('\nvar(Kdowel)  = %.4f (10^6N/rad)^2',vKdowel);
fprintf('\nstd(Kdowel)  = %.4f 10^3N/rad',sKdowel);
fprintf('\ndisp(Kdowel) = %.4f',dKdowel);
fprintf('\nmean(Kscrew) = %.4f 10^4N/rad',mKscrew);
fprintf('\nvar(Kscrew)  = %.4f (10^8N/rad)^2',vKscrew);
fprintf('\nstd(Kscrew)  = %.4f 10^4N/rad',sKscrew);
fprintf('\ndisp(Kscrew) = %.4f',dKscrew);
fprintf('\n');

%% Pdf and cdf
pdf_Kdowel = @(x) gampdf(x,a_dowel,b_dowel);
cdf_Kdowel = @(x) gamcdf(x,a_dowel,b_dowel);
pdf_Kscrew = @(x) gampdf(x,a_screw,b_screw);
cdf_Kscrew = @(x) gamcdf(x,a_screw,b_screw);

%% Sample generation
N = 1e4; % number of samples
k_dowel = gamrnd(a_dowel,b_dowel,N,1);
k_screw = gamrnd(a_screw,b_screw,N,1);

%% Display
if displaySolution
    %% Plot data
    figure('Name','Dowel junction data')
    clf
    bar(1:length(Kdowel_data),Kdowel_data)
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[0,length(Kdowel_data)+1])
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Rigidity of dowel junction $k$ [$\times10^3$N/rad]','Interpreter',interpreter);
    mysaveas(pathname,'data_Kdowel',formats);
    mymatlab2tikz(pathname,'data_Kdowel.tex');

    figure('Name','Screw junction data')
    clf
    bar(1:length(Kscrew_data),Kscrew_data)
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[0,length(Kscrew_data)+1])
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Rigidity of screw junction $k$ [$\times10^4$N/rad]','Interpreter',interpreter);
    mysaveas(pathname,'data_Kscrew',formats);
    mymatlab2tikz(pathname,'data_Kscrew.tex');
    
  %% Plot pdf and cdf
    xmin = max(0,mKdowel-5*sKdowel);
    xmax = mKdowel+5*sKdowel;
    x = linspace(xmin,xmax,1e3);
    
    % Plot pdf of Kdowel
    figure('Name','Probability density function of Kdowel')
    clf
    hold on
    plot(x,pdf_Kdowel(x),'-b','LineWidth',linewidth);
    plot(Kdowel_data,pdf_Kdowel(Kdowel_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xmin,xmax])
    xlabel('$k$ [$\times10^3$N/rad]','Interpreter',interpreter)
    ylabel('$p_K(k)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_Kdowel',formats);
    mymatlab2tikz(pathname,'pdf_Kdowel.tex',...
        'extraAxisOptions',{'ylabel style={overlay}'});
    
    % Plot cdf of Kdowel
    figure('Name','Cumulative distribution function of Kdowel')
    clf
    hold on
    plot(x,cdf_Kdowel(x),'-b','LineWidth',linewidth);
    plot(Kdowel_data,cdf_Kdowel(Kdowel_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xmin,xmax])
    set(gca,'YLim',[0,1])
    xlabel('$k$ [$\times10^3$N/rad]','Interpreter',interpreter)
    ylabel('$F_K(k)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_Kdowel',formats);
    mymatlab2tikz(pathname,'cdf_Kdowel.tex',...
        'extraAxisOptions',{'ylabel style={overlay}'});  
    
    % Plot pdf of Kscrew
    figure('Name','Probability density function of Kscrew')
    clf
    hold on
    plot(x,pdf_Kscrew(x),'-b','LineWidth',linewidth);
    plot(Kscrew_data,pdf_Kscrew(Kscrew_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xmin,xmax])
    xlabel('$k$ [$\times10^4$N/rad]','Interpreter',interpreter)
    ylabel('$p_K(k)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_Kscrew',formats);
    mymatlab2tikz(pathname,'pdf_Kscrew.tex',...
        'extraAxisOptions',{'ylabel style={overlay}'});
    
    % Plot cdf of Kscrew
    figure('Name','Cumulative distribution function of Kscrew')
    clf
    hold on
    plot(x,cdf_Kscrew(x),'-b','LineWidth',linewidth);
    plot(Kscrew_data,cdf_Kscrew(Kscrew_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xmin,xmax])
    set(gca,'YLim',[0,1])
    xlabel('$k$ [$\times10^4$N/rad]','Interpreter',interpreter)
    ylabel('$F_K(k)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_Kscrew',formats);
    mymatlab2tikz(pathname,'cdf_Kscrew.tex',...
        'extraAxisOptions',{'ylabel style={overlay}'});     
   
    %% Plot samples
    figure('Name','Samples of Kdowel')
    clf
    scatter(1:N,k_dowel,'b.')
    hold on
    plot([1 N],[mKdowel mKdowel],'-r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Rigidity of dowel junction $k$ [$\times10^3$N/rad]','Interpreter',interpreter)
    legend('samples','mean');
    mysaveas(pathname,'samples_Kdowel',formats);
    mymatlab2tikz(pathname,'samples_Kdowel.tex');
    
    figure('Name','Samples of Kscrew')
    clf
    scatter(1:N,k_screw,'b.')
    hold on
    plot([1 N],[mKscrew mKscrew],'-r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Rigidity of screw junction $k$ [$\times10^4$N/rad]','Interpreter',interpreter)
    legend('samples','mean');
    mysaveas(pathname,'samples_Kscrew',formats);
    mymatlab2tikz(pathname,'samples_Kscrew.tex');
    
    
end
