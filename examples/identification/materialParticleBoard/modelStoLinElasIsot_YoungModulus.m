%% Stochastic modeling of Young's modulus %%
%%----------------------------------------%%

% clc
clearvars
close all
rng('default');

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
interpreter = 'latex';
formats = {'fig','epsc'};

%% Data
filenameAna = 'data_ET_GL.mat';
pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathnameIdentification,filenameAna));

E_data = mean_ET_data*1e-3; % [GPa]

mE_data = mean(E_data);
% mE_data = sum(E_data)/length(E_data);
vE_data = var(E_data);
% vE_data = length(E_data)/(length(E_data)-1)*moment(E_data,2);
sE_data = std(E_data);
% sE_data = sqrt(vE_data);
dE_data = sE_data/mE_data;

%% Maximum likelihood estimation
param = gamfit(E_data);
% param = mle(E_data,'Distribution','gam');
% param = mle(E_data,'pdf',@gampdf,'Start',[3 1],'LowerBound',[2 0]);
% param = mle(E_data,'nloglf',@gamlike,'Start',[3 1],'LowerBound',[2 0]);

custpdf = @(data,a,b) gampdf(data,a,b);
% custpdf = @(data,a,b) exp(-(a*log(b)+gammaln(a))) * data.^(a-1) .* exp(-data/b);
% param = mle(E_data,'pdf',custpdf,'Start',[3 1],'LowerBound',[2 0]);

nloglf = @(param,data,cens,freq) gamlike(param,data);
% nloglf = @(param,data,cens,freq) length(data)*param(1)*log(param(2))...
%     +length(data)*gammaln(param(1))...
%     +(1-param(1))*sum(log(data))...
%     +1/param(2)*sum(data);
% param = mle(E_data,'nloglf',nloglf,'Start',[3 1],'LowerBound',[2 0]);

a = param(1);
b = param(2);

% [mE,vE] = gamstat(a,b);
mE = a*b;
vE = a*b^2;
sE = sqrt(vE); % sE = sqrt(a)*b;
dE = sE/mE; % dE = 1/sqrt(a);

%% Pdf and cdf
pdf_E = @(x) gampdf(x,a,b);
% pdf_E = @(x) custpdf(x,a,b);
% pdf_E = @(x) pdf('gam',x,a,b);
cdf_E = @(x) gamcdf(x,a,b);
% cdf_E = @(x) cdf('gam',x,a,b);

%% Sample generation
N = 1e3; % number of samples
E_sample = gamrnd(a,b,N,1); % [GPa]
% u = randn(N,1);
% E_sample = gaminv(normcdf(u),a,b); % [GPa]

mE_sample = mean(E_sample);
% mE_sample = sum(E_sample)/length(E_sample);
vE_sample = var(E_sample);
% vE_sample = length(E_sample)/(length(E_sample)-1)*moment(E_sample,2);
sE_sample = std(E_sample);
% sE_sample = sqrt(vE_sample);
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
    bar(E_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Young''s modulus $E$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module d''Young $E$ [GPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_E',formats);
    mymatlab2tikz(pathname,'data_E.tex');
    
    %% Plot log-likelihood function
    a_series = linspace(a*0.5,a*1.5,1e2);
    b_series = linspace(b*0.5,b*1.5,1e2);
    loglf = zeros(length(a_series),length(b_series));
    for i=1:length(a_series)
        aS_i = a_series(i);
        for j=1:length(b_series)
            b_j = b_series(j);
            param_ij = [aS_i,b_j];
            loglf(i,j) = -gamlike(param_ij,E_data);
        end
    end
    
    % Plot log-likelihood function loglf for E
    figure('Name','Surface plot: Log-likelihood function for E')
    clf
    surfc(a_series,b_series,loglf,'EdgeColor','none');
    colorbar
    hold on
    scatter3(a,b,-gamlike(param,E_data),'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    xlabel('$\alpha$','Interpreter',interpreter)
    ylabel('$\beta$ [GPa]','Interpreter',interpreter)
    zlabel('$\mathcal{L}(\alpha,\beta)$','Interpreter',interpreter)
    mysaveas(pathname,'loglf_E_3D',formats);
    % mymatlab2tikz(pathname,'loglf_E_3D.tex');
    
    figure('Name','Contour plot: Log-likelihood function for KS')
    clf
    contourf(a_series,b_series,loglf,30);
    colorbar
    hold on
    scatter(a,b,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    xlabel('$\alpha$','Interpreter',interpreter)
    ylabel('$\beta$ [kN/rad]','Interpreter',interpreter)
    zlabel('$\mathcal{L}(\alpha,\beta)$','Interpreter',interpreter)
    mysaveas(pathname,'loglf_E_2D',formats);
    % mymatlab2tikz(pathname,'loglf_E_2D.tex');

    %% Plot pdf and cdf
    xmin = max(0,mE-5*sE);
    xmax = mE+5*sE;
    x = linspace(xmin,xmax,1e3);
    
    % Plot pdf of E
    figure('Name','Probability density function')
    clf
    plot(x,pdf_E(x),'-b','LineWidth',linewidth);
    % hold on
    % plot(E_data,pdf_E(E_data),'k+');
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[xmin,xmax])
    xlabel('$e$ [GPa]','Interpreter',interpreter)
    ylabel('$p_E(e)$','Interpreter',interpreter)
    % l = legend('$p_E(e)$','$(e_i,p_E(e_i))_{i=1}^n$');
    % set(l,'Interpreter',interpreter,'Location','NorthWest');
    mysaveas(pathname,'pdf_E',formats);
    mymatlab2tikz(pathname,'pdf_E.tex');
    
    % Plot cdf of E
    figure('Name','Cumulative distribution function')
    clf
    plot(x,cdf_E(x),'-r','LineWidth',linewidth);
    % hold on
    % plot(E_data,cdf_E(E_data),'k+');
    % hold off
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
    mymatlab2tikz(pathname,'cdf_E.tex');
    
    %% Plot samples
    figure('Name','Samples')
    clf
    scatter(1:length(E_sample),E_sample,'b.')
    hold on
    plot([1 length(E_sample)],[mE mE],'-r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of realizations','Interpreter',interpreter)
    %ylabel('Young''s modulus $E$ [GPa]','Interpreter',interpreter)
    ylabel('$E$ [GPa]','Interpreter',interpreter)
    legend('realizations','mean value');
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Module d''Young $E$ [GPa]','Interpreter',interpreter)
    %legend('réalisations','valeur moyenne');
    mysaveas(pathname,'samples_E',formats);
    % mymatlab2tikz(pathname,'samples_E.tex');
    
end
