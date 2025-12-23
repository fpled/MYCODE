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

% Empirical estimates
N_data = length(E_data);
mE_data = mean(E_data);
% mE_data = sum(E_data)/N_data;
vE_data = var(E_data);
% vE_data = N_data/(N_data-1)*moment(E_data,2);
sE_data = std(E_data);
% sE_data = sqrt(vE_data);
dE_data = sE_data/mE_data;

%% Maximum likelihood estimation
% display = 'off'; % default for gamfit and mle
display = 'iter';
% display = 'final';

tolX = 1e-8; % tolerance on the parameter value (default: 1e-8 for gamfit and 1e-6 for mle)
tolFun = 1e-8; % tolerance on the function value (default: 1e-8 for gamfit and 1e-6 for mle)
tolBnd = 1e-6; % tolerance on the parameter bound (default: 1e-6 for gamfit and mle)
maxIters = Inf; % maximum number of iterations
maxFunEvals = Inf; % maximum number of function evaluations

options_gamfit = statset(statset('gamfit'),'Display',display,...
    'TolX',tolX,'TolFun',tolFun,'TolBnd',tolBnd,...
    'MaxIter',maxIters,'MaxFunEvals',maxFunEvals);

options_mlecustom = statset(statset('mlecustom'),'Display',display,...
    'TolX',tolX,'TolFun',tolFun,'TolBnd',tolBnd,...
    'MaxIter',maxIters,'MaxFunEvals',maxFunEvals);

param = gamfit(E_data,[],[],[],options_gamfit);
% param = mle(E_data,'Distribution','gam','Options',options_gamfit);

% Initial parameter values using Method of moments estimates
a0 = mE_data^2/vE_data; % a > 2
b0 = mE_data/a0;        % b > 0

param0 = [a0 b0]; % initial parameter vector
lb     = [2 0];   % lower bounds

loglfval0 = -gamlike(param0,E_data); % initial log-likelihood

% param = mle(E_data,'pdf',@gampdf,'Start',param0,'LowerBound',lb,'Options',options_mlecustom);
% param = mle(E_data,'nloglf',@gamlike,'Start',param0,'LowerBound',lb,'Options',options_mlecustom);

custpdf = @(data,a,b) gampdf(data,a,b);
% custpdf = @(data,a,b) exp(-(a*log(b)+gammaln(a))) * data.^(a-1) .* exp(-data/b);
% param = mle(E_data,'pdf',custpdf,'Start',param0,'LowerBound',lb,'Options',options_mlecustom);

nloglf = @(param,data,cens,freq) gamlike(param,data);
% nloglf = @(param,data,cens,freq) length(data)*param(1)*log(param(2))...
%     +length(data)*gammaln(param(1))...
%     +(1-param(1))*sum(log(data))...
%     +1/param(2)*sum(data);
% param = mle(E_data,'nloglf',nloglf,'Start',param0,'LowerBound',lb,'Options',options_mlecustom);

% Optimal parameter values
a = param(1); % a > 2
b = param(2); % b > 0

[mE,vE] = gamstat(a,b);
% mE = a*b;
% vE = a*b^2;
sE = sqrt(vE); % sE = sqrt(a)*b;
dE = sE/mE; % dE = 1/sqrt(a);

loglfval = -gamlike(param,E_data); % optimal (maximal) log-likelihood = log-likelihood function evaluated at MLE param

%% Pdf and cdf
pdf_E = @(x) gampdf(x,a,b);
% pdf_E = @(x) custpdf(x,a,b);
% pdf_E = @(x) pdf('gam',x,a,b);
cdf_E = @(x) gamcdf(x,a,b);
% cdf_E = @(x) cdf('gam',x,a,b);

%% Sample generation
N = 1e4; % number of samples
E_sample = gamrnd(a,b,N,1); % [GPa]
% u = randn(N,1);
% E_sample = gaminv(normcdf(u),a,b); % [GPa]

mE_sample = mean(E_sample);
% mE_sample = sum(E_sample)/N;
vE_sample = var(E_sample);
% vE_sample = N/(N-1)*moment(E_sample,2);
sE_sample = std(E_sample);
% sE_sample = sqrt(vE_sample);
dE_sample = sE_sample/mE_sample;

%% Outputs
filenameResults = fullfile(pathname,'results.txt');
fid = fopen(filenameResults,'w');
fprintf(fid,'nb data    = %g\n',N_data);
fprintf(fid,'nb samples = %g\n',N);

fprintf(fid,'\n');
fprintf(fid,'Initial parameter values\n');
fprintf(fid,'alpha = %.4f\n',a0);
fprintf(fid,'beta  = %.6f\n',b0);
fprintf(fid,'loglf = %.4f\n',loglfval0);

fprintf(fid,'\n');
fprintf(fid,'Optimal parameter values\n');
fprintf(fid,'alpha = %.4f\n',a);
fprintf(fid,'beta  = %.6f\n',b);
fprintf(fid,'loglf = %.4f\n',loglfval);

fprintf(fid,'\n');
fprintf(fid,'mean(E_data)   = %.4f GPa\n',mE_data);
fprintf(fid,'mean(E)        = %.4f GPa\n',mE);
fprintf(fid,'mean(E_sample) = %.4f GPa\n',mE_sample);
fprintf(fid,'var(E_data)    = %.4f (GPa)^2\n',vE_data);
fprintf(fid,'var(E)         = %.4f (GPa)^2\n',vE);
fprintf(fid,'var(E_sample)  = %.4f (GPa)^2\n',vE_sample);
fprintf(fid,'std(E_data)    = %.4f GPa\n',sE_data);
fprintf(fid,'std(E)         = %.4f GPa\n',sE);
fprintf(fid,'std(E_sample)  = %.4f GPa\n',sE_sample);
fprintf(fid,'cv(E_data)     = %.4f\n',dE_data);
fprintf(fid,'cv(E)          = %.4f\n',dE);
fprintf(fid,'cv(E_sample)   = %.4f\n',dE_sample);

err_meanE = abs(mE - mE_data)/abs(mE_data);
err_varE = abs(vE - vE_data)/abs(vE_data);
err_stdE = abs(sE - sE_data)/abs(sE_data);
err_cvE = abs(dE - dE_data)/abs(dE_data);

err_meanE_sample = abs(mE_sample - mE_data)/abs(mE_data);
err_varE_sample = abs(vE_sample - vE_data)/abs(vE_data);
err_stdE_sample = abs(sE_sample - sE_data)/abs(sE_data);
err_cvE_sample = abs(dE_sample - dE_data)/abs(dE_data);

alpha = 1/2;
mse = alpha * (mE - mE_data)^2/(mE_data)^2 + (1-alpha) * (dE - dE_data)^2/(dE_data)^2;
mse_sample = alpha * (mE_sample - mE_data)^2/(mE_data)^2 + (1-alpha) * (dE_sample - dE_data)^2/(dE_data)^2;

fprintf(fid,'\n');
fprintf(fid,'relative error on mean(E) = %.4e\n',err_meanE);
fprintf(fid,'relative error on var(E)  = %.4e\n',err_varE);
fprintf(fid,'relative error on std(E)  = %.4e\n',err_stdE);
fprintf(fid,'relative error on cv(E)   = %.4e\n',err_cvE);
fprintf(fid,'mean-squared error mse(E) = %.4e\n',mse);

fprintf(fid,'\n');
fprintf(fid,'relative error on mean(E_sample) = %.4e\n',err_meanE_sample);
fprintf(fid,'relative error on var(E_sample)  = %.4e\n',err_varE_sample);
fprintf(fid,'relative error on std(E_sample)  = %.4e\n',err_stdE_sample);
fprintf(fid,'relative error on cv(E_sample)   = %.4e\n',err_cvE_sample);
fprintf(fid,'mean-squared error mse(E_sample) = %.4e\n',mse_sample);
fclose(fid);
type(filenameResults) % fprintf('%s', fileread(filenameResults))

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
    loglf = @(A,B,data) arrayfun(@(a,b) -gamlike([a,b],data), A, B);

    a_series = linspace(a*0.5,a*1.5,1e2);
    b_series = linspace(b*0.5,b*1.5,1e2);
    logL = zeros(length(a_series),length(b_series));
    for i=1:length(a_series)
        a_i = a_series(i);
        for j=1:length(b_series)
            b_j = b_series(j);
            param_ij = [a_i,b_j];
            logL(j,i) = -gamlike(param_ij,E_data);
        end
    end
    [A,B] = meshgrid(a_series,b_series);
    logL = loglf(A,B,E_data);
    
    % Plot log-likelihood function loglf for E
    figure('Name','Surface plot: Log-likelihood function for E')
    clf
    surfc(a_series,b_series,logL,'EdgeColor','none');
    colorbar
    hold on
    scatter3(a,b,loglfval,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    xlabel('$\alpha$','Interpreter',interpreter)
    ylabel('$\beta$ [GPa]','Interpreter',interpreter)
    zlabel('$\mathcal{L}(\alpha,\beta)$','Interpreter',interpreter)
    mysaveas(pathname,'loglf_E_3D',formats);
    % mymatlab2tikz(pathname,'loglf_E_3D.tex');
    
    figure('Name','Contour plot: Log-likelihood function for E')
    clf
    contourf(a_series,b_series,logL,50);
    colorbar
    hold on
    scatter(a,b,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    xlabel('$\alpha$','Interpreter',interpreter)
    ylabel('$\beta$ [GPa]','Interpreter',interpreter)
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
    % legend('$p_E(e)$','$(e_i,p_E(e_i))_{i=1}^n$','Location','NorthWest','Interpreter',interpreter)
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
    % legend('$F_E(e)$','$(e_i,F_E(e_i))_{i=1}^n$','Location','NorthWest','Interpreter',interpreter)
    mysaveas(pathname,'cdf_E',formats);
    mymatlab2tikz(pathname,'cdf_E.tex');
    
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
