%% Stochastic modeling of random linear elasticity tensor %%
%  with transversely isotropic symmetry                   %%
%%--------------------------------------------------------%%

% clc
clearvars
close all
rng('default');
s = rng;

%% Input data
displaySolution = true;
displayConvergenceInit = true;
displayConvergenceOpt = true;

useRedParam = true; % reduced parameterization
MCMCalgo = 'IMH'; % algorithm for Markov Chain Monte Carlo (MCMC) method = 'IMH', 'RWMH' or 'SS'
useParallel = true; % parallel computing to evaluate gradients of objective function
% plotFcn = []; % plot functions called at each iteration
plotFcn = {'optimplot',...
    'optimplotx','optimplotfunccount','optimplotfval','optimplotfvalconstr',...
    'optimplotconstrviolation','optimplotstepsize','optimplotfirstorderopt'}; % plot functions called at each iteration
globalSolver = []; % local solver
% globalSolver = 'PatternSearch'; % PatternSearch solver
% globalSolver = 'GlobalSearch'; % GlobalSearch solver
% globalSolver = 'MultiStart'; % MultiStart solver

filename = 'modelStoLinElasIsotTrans_ElasTensor_';
if useRedParam
    filename = [filename 'ReducedParam_'];
else
    filename = [filename 'FullParam_'];
end
filename = [filename MCMCalgo];
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

N_data = length(mean_ET_data);
ET_data = mean_ET_data*1e-3; % [GPa]
GL_data = mean_GL_data*1e-3; % [GPa]
EL_data = mean_EL_data*1e-3; % [GPa]
% NUL_data = mean_NUL_data;
NUL_data = 0.03+0.03*rand(N_data,1); % artificial data for NUL uniformly distributed from 0.03 to 0.06
NUT_data = 0.1+0.2*rand(N_data,1);   % artificial data for NUT uniformly distributed from 0.1 to 0.3
GT_data = ET_data./(2*(1+NUT_data)); % [GPa]
kT_data = (EL_data.*ET_data)./(2*(1-NUT_data).*EL_data-4*ET_data.*(NUL_data).^2); % [GPa]
C1_data = EL_data + 4*(NUL_data.^2).*kT_data; % [GPa]
C2_data = 2*kT_data; % [GPa]
C3_data = 2*sqrt(2)*kT_data.*NUL_data; % [GPa]
C4_data = 2*GT_data; % [GPa]
C5_data = 2*GL_data; % [GPa]
C_data = [C1_data(:) C2_data(:) C3_data(:) C4_data(:) C5_data(:)];

% Empirical estimates
mC_data = mean(C_data,1);
vC_data = var(C_data,0,1); % vC_data = size(C_data,1)/(size(C_data,1)-1)*moment(C_data,2,1);
stdC_data = std(C_data,0,1); % stdC_data = sqrt(vC_data);
cvC_data = stdC_data./mC_data;
% sC_data = sqrt(norm(vC_data));
mCnorm_data = norm(mC_data);
% dC_data = sC_data/mCnorm_data;
phiC_data = log((C1_data.*C2_data-C3_data.^2).*(C4_data.^2).*(C5_data.^2));
nuC_data = mean(phiC_data,1);
fC_data = [mC_data nuC_data];

%% Least-Squares estimation for computing Lagrange multipliers
% Initial parameter values
cvC45 = cvC_data(4:5); % empirical coefficient of variations (cv) for C4 and C5
la0 = (1-1/mean(cvC45).^2)/2; % la < 1/2, la < 0 for MH sampling (match cv with averaged empirical cv for C4 and C5)
la01 = -(mC_data(2)*la0)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la1 > 0 (match mode with empirical mean value for C1)
la02 = -(mC_data(1)*la0)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la2 > 0 (match mode with empirical mean value for C2)
la03 = (2*mC_data(3)*la0)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la3 in R such that 2*sqrt(la1*la2)-la3 > 0 (match mode with empirical mean value for C3)
a0 = 1-2*la0; % a > 0
la04 = a0/mC_data(4); % la4 > 0 (match mathematical expectation with empirical mean value for C4)
la05 = a0/mC_data(5); % la5 > 0 (match mathematical expectation with empirical mean value for C5)
% la04 = -2*la0/mC_data(4); % la4 > 0 (match mode with empirical mean value for C4)
% la05 = -2*la0/mC_data(5); % la5 > 0 (match mode with empirical mean value for C5)
b04 = 1/la04; % b4 > 0
b05 = 1/la05; % b5 > 0

% Initial parameter vector
if useRedParam
    lambda0 = [la01 la02 la03 la0]; % reduced parameterization
else
    lambda0 = [la01 la02 la03 la04 la05 la0]; % full parameterization
end

%% Display convergence for initial parameter vector
if displayConvergenceInit
    location = 'SouthEast';
    hmC1 = figure('Name','Convergence mean(C1) for initial parameters');
    clf
    hmC2 = figure('Name','Convergence mean(C2) for initial parameters');
    clf
    hmC3 = figure('Name','Convergence mean(C3) for initial parameters');
    clf
    hmCnorm = figure('Name','Convergence norm(mean(C)) for initial parameters');
    clf
    hmphiC = figure('Name','Convergence mean(phi(C))=mean(log(det([C]))) for initial parameters');
    clf
    herrmC1 = figure('Name','Convergence error mean(C1) for initial parameters');
    clf
    herrmC2 = figure('Name','Convergence error mean(C2) for initial parameters');
    clf
    herrmC3 = figure('Name','Convergence error mean(C3) for initial parameters');
    clf
    herrmCnorm = figure('Name','Convergence error norm(mean(C)) for initial parameters');
    clf
    herrmphiC = figure('Name','Convergence error mean(phi(C))=mean(log(det([C]))) for initial parameters');
    clf
    hresnorm = figure('Name','Convergence squared norm of residual for initial parameters');
    clf
    hfval = figure('Name','Convergence objective function value for initial parameters');
    clf
    
    rng(s); % initialize the random number generator using the default settings contained in s
    Ndraws = 3; % number of draws
    N = 1e7; % number of samples
    Nexponent = floor(log10(N));
    
    legmC1    = cell(1,Ndraws+1);
    legmC2    = cell(1,Ndraws+1);
    legmC3    = cell(1,Ndraws+1);
    legmCnorm = cell(1,Ndraws+1);
    legmphiC  = cell(1,Ndraws+1);
    legerrmC1    = cell(1,Ndraws);
    legerrmC2    = cell(1,Ndraws);
    legerrmC3    = cell(1,Ndraws);
    legerrmCnorm = cell(1,Ndraws);
    legerrmphiC  = cell(1,Ndraws);
    legresnorm   = cell(1,Ndraws);
    legfval      = cell(1,Ndraws);
    
    fprintf('\n');
    fprintf('Convergence analysis\n');
    if useRedParam
        fprintf('initial lambda = (%g, %g, %g, %g, %g, %g)\n',lambda0(1:3),la04,la05,lambda0(4));
    else
        fprintf('initial lambda = (%g, %g, %g, %g, %g, %g)\n',lambda0);
    end
    
    mC45 = gamstat(a0,[b04,b05]);
    nuC4 = psi(a0)+log(b04);
    nuC5 = psi(a0)+log(b05);
    for n=1:Ndraws
        tinit = tic;
        [f0,C_sample] = funoptimlseElasIsotTrans(lambda0,C_data,mC_data,nuC_data,N,MCMCalgo);
        
        mC123 = mean(C_sample,1);
        mC = [mC123 mC45];
        mCnorm = norm(mC);
        phiC123 = log(C_sample(:,1).*C_sample(:,2)-C_sample(:,3).^2);
        nuC123 = mean(phiC123,1);
        nuC = nuC123 + 2*(nuC4 + nuC5);
        fC = [mC nuC];
        resnorm   = norm(fC - fC_data)^2;
        errmCnorm = norm(mC - mC_data)/norm(mC_data);
        err       = abs(fC - fC_data)./abs(fC_data);
        fprintf('draw #%d: fval = %g, resnorm = %g, error(norm(mean(C))) = %g, error = (%g, %g, %g, %g, %g, %g), elapsed time = %f s\n',n,f0,resnorm,errmCnorm,err,toc(tinit));
        
        Ns = (1:N)';
        mC123s = cumsum(C_sample,1) ./ Ns;
        mC1s = mC123s(:,1);
        mC2s = mC123s(:,2);
        mC3s = mC123s(:,3);
        
        mCs  = [mC123s repmat(mC45,N,1)];
        mCnorms = sqrt(sum(mCs.^2,2));
        
        nuC123s = cumsum(phiC123,1) ./ Ns;
        nuCs = nuC123s + 2*(nuC4 + nuC5);
        
        fCs = [mCs nuCs];
        resnorms = sum((fCs - fC_data).^2,2); % squared norm of residual
        % resnorms = sum((fCs - fC_data).^2,2)/norm(fC_data)^2; % relative squared norm of residual
        % resnorms = sum((fCs - fC_data).^2,2)/sum(fC_data.^2,2); % relative squared norm of residual
        % fvals = sum((mCs - mC_data).^2,2)/norm(mC_data)^2 + (abs(nuCs - nuC_data).^2)./abs(nuC_data)^2; % objective function value (relative squared norm of residual)
        fvals = sum((mCs - mC_data).^2,2)/sum(mC_data.^2,2) + (abs(nuCs - nuC_data).^2)./abs(nuC_data)^2; % objective function value (relative squared norm of residual)
        % errmCnorms = sqrt(sum((mCs - mC_data).^2,2))/norm(mC_data); % relative error on norm(mC)
        errmCnorms = sqrt(sum((mCs - mC_data).^2,2)/sum(mC_data.^2,2)); % relative error on norm(mC)
        errs = abs(fCs - fC_data)./abs(fC_data); % relative errors
        
        figure(hmC1)
        semilogx(1:N,mC1s,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legmC1{n} = ['$\widehat{\underline{c}}_1(${\boldmath$\lambda$}$^{\mathrm{init}})$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(hmC2)
        semilogx(1:N,mC2s,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legmC2{n} = ['$\widehat{\underline{c}}_2(${\boldmath$\lambda$}$^{\mathrm{init}})$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(hmC3)
        semilogx(1:N,mC3s,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legmC3{n} = ['$\widehat{\underline{c}}_3(${\boldmath$\lambda$}$^{\mathrm{init}})$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(hmCnorm)
        semilogx(1:N,mCnorms,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legmCnorm{n} = ['$||\widehat{\underline{\mathbf{c}}}(${\boldmath$\lambda$}$^{\mathrm{init}})||$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(hmphiC)
        semilogx(1:N,nuCs,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legmphiC{n} = ['$\widehat{\nu}_{\mathbf{C}}(${\boldmath$\lambda$}$^{\mathrm{init}})$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(herrmC1)
        loglog(1:N,errs(:,1),'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legerrmC1{n} = ['$|\widehat{\underline{c}}_1(${\boldmath$\lambda$}$^{\mathrm{init}})-\underline{c}_1|/|\underline{c}_1|$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(herrmC2)
        loglog(1:N,errs(:,2),'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legerrmC2{n} = ['$|\widehat{\underline{c}}_2(${\boldmath$\lambda$}$^{\mathrm{init}})-\underline{c}_2|/|\underline{c}_2|$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(herrmC3)
        loglog(1:N,errs(:,3),'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legerrmC3{n} = ['$|\widehat{\underline{c}}_3(${\boldmath$\lambda$}$^{\mathrm{init}})-\underline{c}_3|/|\underline{c}_3|$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(herrmCnorm)
        loglog(1:N,errmCnorms,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legerrmCnorm{n} = ['$||\widehat{\underline{\mathbf{c}}}(${\boldmath$\lambda$}$^{\mathrm{init}})-\underline{\mathbf{c}}||/||\underline{\mathbf{c}}||$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(herrmphiC)
        loglog(1:N,errs(:,6),'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legerrmphiC{n} = ['$|\widehat{\nu}_{\mathbf{C}}(${\boldmath$\lambda$}$^{\mathrm{init}})-\nu_{\mathbf{C}}|/|\nu_{\mathbf{C}}|$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(hresnorm)
        loglog(1:N,resnorms,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legresnorm{n} = ['$||\widehat{\mathbf{f}}(${\boldmath$\lambda$}$^{\mathrm{init}}) - \mathbf{f}||^2$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(hfval)
        loglog(1:N,fvals,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legfval{n} = ['$\mathcal{J}(${\boldmath$\lambda$}$^{\mathrm{init}})$, draw $\#' num2str(n) '$'];
        hold on
    end
    figure(hmC1)
    plot([1,N],[mC_data(1),mC_data(1)],'--','Color','k','LineWidth',linewidth)
    legmC1{end} = '$\underline{c}_1$';
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Mean value of $C_1$ [GPa]','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Valeur moyenne de $C_1$ [GPa]','Interpreter',interpreter)
    legend(legmC1{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_mean_C1_lambda_init',formats);
    % mymatlab2tikz(pathname,'convergence_mean_C1_lambda_init.tex');
    
    figure(hmC2)
    plot([1,N],[mC_data(2),mC_data(2)],'--','Color','k','LineWidth',linewidth)
    legmC2{end} = '$\underline{c}_2$';
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Mean value of $C_2$ [GPa]','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Valeur moyenne de $C_2$ [GPa]','Interpreter',interpreter)
    legend(legmC2{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_mean_C2_lambda_init',formats);
    % mymatlab2tikz(pathname,'convergence_mean_C2_lambda_init.tex');
    
    figure(hmC3)
    plot([1,N],[mC_data(3),mC_data(3)],'--','Color','k','LineWidth',linewidth)
    legmC3{end} = '$\underline{c}_3$';
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Mean value of $C_3$ [GPa]','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Valeur moyenne de $C_3$ [GPa]','Interpreter',interpreter)
    legend(legmC3{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_mean_C3_lambda_init',formats);
    % mymatlab2tikz(pathname,'convergence_mean_C3_lambda_init.tex');
    
    figure(hmCnorm)
    plot([1,N],[mCnorm_data,mCnorm_data],'--','Color','k','LineWidth',linewidth)
    legmCnorm{end} = '$||\underline{\mathbf{c}}||$';
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Norm of mean value of \boldmath$C$ [GPa]','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Norme de la valeur moyenne de \boldmath$C$ [GPa]','Interpreter',interpreter)
    legend(legmCnorm{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_norm_mean_C_lambda_init',formats);
    % mymatlab2tikz(pathname,'convergence_norm_mean_C_lambda_init.tex');
    
    figure(hmphiC)
    plot([1,N],[nuC_data,nuC_data],'--','Color','k','LineWidth',linewidth)
    legmphiC{end} = '$\nu_{\mathbf{C}}$';
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Mean value of $\varphi(${\boldmath$C$}$)$','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Valeur moyenne de $\varphi(${\boldmath$C$}$)$','Interpreter',interpreter)
    legend(legmphiC{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_mean_phiC_lambda_init',formats);
    % mymatlab2tikz(pathname,'convergence_mean_phiC_lambda_init.tex');
    
    figure(herrmC1)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Relative error','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Erreur relative','Interpreter',interpreter)
    legend(legerrmC1{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_error_mean_C1_lambda_init',formats);
    % mymatlab2tikz(pathname,'convergence_error_mean_C1_lambda_init.tex');
    
    figure(herrmC2)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Relative error','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Erreur relative','Interpreter',interpreter)
    legend(legerrmC2{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_error_mean_C2_lambda_init',formats);
    % mymatlab2tikz(pathname,'convergence_error_mean_C2_lambda_init.tex');
    
    figure(herrmC3)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Relative error','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Erreur relative','Interpreter',interpreter)
    legend(legerrmC3{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_error_mean_C3_lambda_init',formats);
    % mymatlab2tikz(pathname,'convergence_error_mean_C3_lambda_init.tex');
    
    figure(herrmCnorm)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Relative error','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Erreur relative','Interpreter',interpreter)
    legend(legerrmCnorm{:},'Location','NorthEast','Interpreter',interpreter)
    mysaveas(pathname,'convergence_error_norm_mean_C_lambda_init',formats);
    % mymatlab2tikz(pathname,'convergence_error_norm_mean_C_lambda_init.tex');
    
    figure(herrmphiC)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Relative error','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Erreur relative','Interpreter',interpreter)
    legend(legerrmphiC{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_error_mean_phiC_lambda_init',formats);
    % mymatlab2tikz(pathname,'convergence_error_mean_phiC_lambda_init.tex');
    
    figure(hresnorm)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Squared norm of residual','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Norme au carr\''e du r\''esidu','Interpreter',interpreter)
    legend(legresnorm{:},'Location','NorthEast','Interpreter',interpreter)
    mysaveas(pathname,'convergence_resnorm_lambda_init',formats);
    % mymatlab2tikz(pathname,'convergence_resnorm_lambda_init.tex');
    
    figure(hfval)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Objective function value','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Valeur de la fonction objectif','Interpreter',interpreter)
    legend(legfval{:},'Location','NorthEast','Interpreter',interpreter)
    mysaveas(pathname,'convergence_fval_lambda_init',formats);
    % mymatlab2tikz(pathname,'convergence_fval_lambda_init.tex');
end

% Objective function value for initial parameter vector
N = 1e6; % number of samples
switch lower(MCMCalgo)
    case {'imh','rwmh'}
        [fval0,C0_sample,accept0] = funoptimlseElasIsotTrans(lambda0,C_data,mC_data,nuC_data,N,MCMCalgo,s);
    case 'ss'
        [fval0,C0_sample,neval0] = funoptimlseElasIsotTrans(lambda0,C_data,mC_data,nuC_data,N,MCMCalgo,s);
end

% Parameter estimation
myparallel('start')
t = tic;
% lambda = lseStoLinElasIsotTrans(C_data,lambda0,N,MCMCalgo,s,'useParallel',useParallel,'plotFcn',plotFcn,globalSolver);
if isempty(globalSolver)
    [lambda,fval,exitflag,output] = lseStoLinElasIsotTrans(C_data,lambda0,N,MCMCalgo,s,'display','iter-detailed','useParallel',useParallel,'plotFcn',plotFcn);
else
    [lambda,fval,exitflag,output,solutions] = lseStoLinElasIsotTrans(C_data,lambda0,N,MCMCalgo,s,'display','iter-detailed','useParallel',useParallel,'plotFcn',plotFcn,globalSolver);
end
toc(t)
myparallel('stop')

% Optimal parameter values
la1 = lambda(1); % la1 > 0
la2 = lambda(2); % la2 > 0
la3 = lambda(3); % la3 in R such that 2*sqrt(la1*la2)-la3 > 0
if useRedParam
    la  = lambda(4);    % la < 1/2, la < 0 for MH sampling using a trivariate normal proposal pdf with positive-definite covariance matrix Sigma
    a   = 1-2*la;       % a > 0
    la4 = a/mC_data(4); % la4 > 0
    la5 = a/mC_data(5); % la5 > 0
else
    la4 = lambda(4); % la4 > 0
    la5 = lambda(5); % la5 > 0
    la  = lambda(6); % la < 1/2, la < 0 for MH sampling using a trivariate normal proposal pdf with positive-definite covariance matrix Sigma
    a   = 1-2*la;    % a > 0
end
b4 = 1/la4; % b4 > 0
b5 = 1/la5; % b5 > 0

%% Display convergence for optimal parameter vector
if displayConvergenceOpt
    location = 'SouthEast';
    hmC1 = figure('Name','Convergence mean(C1) for optimal parameters');
    clf
    hmC2 = figure('Name','Convergence mean(C2) for optimal parameters');
    clf
    hmC3 = figure('Name','Convergence mean(C3) for optimal parameters');
    clf
    hmCnorm = figure('Name','Convergence norm(mean(C)) for optimal parameters');
    clf
    hmphiC = figure('Name','Convergence mean(phi(C))=mean(log(det([C]))) for optimal parameters');
    clf
    herrmC1 = figure('Name','Convergence error mean(C1) for optimal parameters');
    clf
    herrmC2 = figure('Name','Convergence error mean(C2) for optimal parameters');
    clf
    herrmC3 = figure('Name','Convergence error mean(C3) for optimal parameters');
    clf
    herrmCnorm = figure('Name','Convergence error norm(mean(C)) for optimal parameters');
    clf
    herrmphiC = figure('Name','Convergence error mean(phi(C))=mean(log(det([C]))) for optimal parameters');
    clf
    hresnorm = figure('Name','Convergence squared norm of residual for optimal parameters');
    clf
    hfval = figure('Name','Convergence objective function value for optimal parameters');
    clf
    
    rng(s); % initialize the random number generator using the default settings contained in s
    Ndraws = 3; % number of draws
    N = 1e7; % number of samples
    Nexponent = floor(log10(N));
    
    legmC1    = cell(1,Ndraws+1);
    legmC2    = cell(1,Ndraws+1);
    legmC3    = cell(1,Ndraws+1);
    legmCnorm = cell(1,Ndraws+1);
    legmphiC  = cell(1,Ndraws+1);
    legerrmC1    = cell(1,Ndraws);
    legerrmC2    = cell(1,Ndraws);
    legerrmC3    = cell(1,Ndraws);
    legerrmCnorm = cell(1,Ndraws);
    legerrmphiC  = cell(1,Ndraws);
    legresnorm   = cell(1,Ndraws);
    legfval      = cell(1,Ndraws);
    
    fprintf('\n');
    fprintf('Convergence analysis\n');
    if useRedParam
        fprintf('optimal lambda = (%g, %g, %g, %g, %g, %g)\n',lambda(1:3),la4,la5,lambda(4));
    else
        fprintf('optimal lambda = (%g, %g, %g, %g, %g, %g)\n',lambda);
    end
    
    mC45 = gamstat(a,[b4,b5]);
    nuC4 = psi(a)+log(b4);
    nuC5 = psi(a)+log(b5);
    for n=1:Ndraws
        topt = tic;
        [f,C_sample] = funoptimlseElasIsotTrans(lambda,C_data,mC_data,nuC_data,N,MCMCalgo);
        
        mC123 = mean(C_sample,1);
        mC = [mC123 mC45];
        mCnorm = norm(mC);
        phiC123 = log(C_sample(:,1).*C_sample(:,2)-C_sample(:,3).^2);
        nuC123 = mean(phiC123,1);
        nuC = nuC123 + 2*(nuC4 + nuC5);
        fC = [mC nuC];
        resnorm   = norm(fC - fC_data)^2;
        errmCnorm = norm(mC - mC_data)/norm(mC_data);
        err       = abs(fC - fC_data)./abs(fC_data);
        fprintf('draw #%d: fval = %g, resnorm = %g, error(norm(mean(C))) = %g, error = (%g, %g, %g, %g, %g, %g), elapsed time = %f s\n',n,f,resnorm,errmCnorm,err,toc(topt));
        
        Ns = (1:N)';
        mC123s = cumsum(C_sample,1) ./ Ns;
        mC1s = mC123s(:,1);
        mC2s = mC123s(:,2);
        mC3s = mC123s(:,3);
        
        mCs  = [mC123s repmat(mC45,N,1)];
        mCnorms = sqrt(sum(mCs.^2,2));
        
        nuC123s = cumsum(phiC123,1) ./ Ns;
        nuCs = nuC123s + 2*(nuC4 + nuC5);
        
        fCs = [mCs nuCs];
        resnorms = sum((fCs - fC_data).^2,2); % squared norm of residual
        % resnorms = sum((fCs - fC_data).^2,2)/norm(fC_data)^2; % relative squared norm of residual
        % resnorms = sum((fCs - fC_data).^2,2)/sum(fC_data.^2,2); % relative squared norm of residual
        % fvals = sum((mCs - mC_data).^2,2)/norm(mC_data)^2 + (abs(nuCs - nuC_data).^2)./abs(nuC_data)^2; % objective function value (relative squared norm of residual)
        fvals = sum((mCs - mC_data).^2,2)/sum(mC_data.^2,2) + (abs(nuCs - nuC_data).^2)./abs(nuC_data)^2; % objective function value (relative squared norm of residual)
        % errmCnorms = sqrt(sum((mCs - mC_data).^2,2))/norm(mC_data); % relative error on norm(mC)
        errmCnorms = sqrt(sum((mCs - mC_data).^2,2)/sum(mC_data.^2,2)); % relative error on norm(mC)
        errs = abs(fCs - fC_data)./abs(fC_data); % relative errors
        
        figure(hmC1)
        semilogx(1:N,mC1s,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legmC1{n} = ['$\widehat{\underline{c}}_1(${\boldmath$\lambda$}$^{\mathrm{opt}})$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(hmC2)
        semilogx(1:N,mC2s,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legmC2{n} = ['$\widehat{\underline{c}}_2(${\boldmath$\lambda$}$^{\mathrm{opt}})$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(hmC3)
        semilogx(1:N,mC3s,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legmC3{n} = ['$\widehat{\underline{c}}_3(${\boldmath$\lambda$}$^{\mathrm{opt}})$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(hmCnorm)
        semilogx(1:N,mCnorms,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legmCnorm{n} = ['$||\widehat{\underline{\mathbf{c}}}(${\boldmath$\lambda$}$^{\mathrm{opt}})||$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(hmphiC)
        semilogx(1:N,nuCs,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legmphiC{n} = ['$\widehat{\nu}_{\mathbf{C}}(${\boldmath$\lambda$}$^{\mathrm{opt}})$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(herrmC1)
        loglog(1:N,errs(:,1),'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legerrmC1{n} = ['$|\widehat{\underline{c}}_1(${\boldmath$\lambda$}$^{\mathrm{opt}})-\underline{c}_1|/|\underline{c}_1|$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(herrmC2)
        loglog(1:N,errs(:,2),'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legerrmC2{n} = ['$|\widehat{\underline{c}}_2(${\boldmath$\lambda$}$^{\mathrm{opt}})-\underline{c}_2|/|\underline{c}_2|$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(herrmC3)
        loglog(1:N,errs(:,3),'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legerrmC3{n} = ['$|\widehat{\underline{c}}_3(${\boldmath$\lambda$}$^{\mathrm{opt}})-\underline{c}_3|/|\underline{c}_3|$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(herrmCnorm)
        loglog(1:N,errmCnorms,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legerrmCnorm{n} = ['$||\widehat{\underline{\mathbf{c}}}(${\boldmath$\lambda$}$^{\mathrm{opt}})-\underline{\mathbf{c}}||/||\underline{\mathbf{c}}||$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(herrmphiC)
        loglog(1:N,errs(:,6),'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legerrmphiC{n} = ['$|\widehat{\nu}_{\mathbf{C}}(${\boldmath$\lambda$}$^{\mathrm{opt}})-\nu_{\mathbf{C}}|/|\nu_{\mathbf{C}}|$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(hresnorm)
        loglog(1:N,resnorms,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legresnorm{n} = ['$||\widehat{\mathbf{f}}(${\boldmath$\lambda$}$^{\mathrm{opt}}) - \mathbf{f}||^2$, draw $\#' num2str(n) '$'];
        hold on
        
        figure(hfval)
        loglog(1:N,fvals,'-','Color',getfacecolor(3+n),'LineWidth',linewidth)
        legfval{n} = ['$\mathcal{J}(${\boldmath$\lambda$}$^{\mathrm{opt}})$, draw $\#' num2str(n) '$'];
        hold on
    end
    figure(hmC1)
    plot([1,N],[mC_data(1),mC_data(1)],'--','Color','k','LineWidth',linewidth)
    legmC1{end} = '$\underline{c}_1$';
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Mean value of $C_1$ [GPa]','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Valeur moyenne de $C_1$ [GPa]','Interpreter',interpreter)
    legend(legmC1{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_mean_C1_lambda_opt',formats);
    % mymatlab2tikz(pathname,'convergence_mean_C1_lambda_opt.tex');
    
    figure(hmC2)
    plot([1,N],[mC_data(2),mC_data(2)],'--','Color','k','LineWidth',linewidth)
    legmC2{end} = '$\underline{c}_2$';
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Mean value of $C_2$ [GPa]','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Valeur moyenne de $C_2$ [GPa]','Interpreter',interpreter)
    legend(legmC2{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_mean_C2_lambda_opt',formats);
    % mymatlab2tikz(pathname,'convergence_mean_C2_lambda_opt.tex');
    
    figure(hmC3)
    plot([1,N],[mC_data(3),mC_data(3)],'--','Color','k','LineWidth',linewidth)
    legmC3{end} = '$\underline{c}_3$';
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Mean value of $C_3$ [GPa]','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Valeur moyenne de $C_3$ [GPa]','Interpreter',interpreter)
    legend(legmC3{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_mean_C3_lambda_opt',formats);
    % mymatlab2tikz(pathname,'convergence_mean_C3_lambda_opt.tex');
    
    figure(hmCnorm)
    plot([1,N],[mCnorm_data,mCnorm_data],'--','Color','k','LineWidth',linewidth)
    legmCnorm{end} = '$||\underline{\mathbf{c}}||$';
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Norm of mean value of \boldmath$C$ [GPa]','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Norme de la valeur moyenne de \boldmath$C$ [GPa]','Interpreter',interpreter)
    legend(legmCnorm{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_norm_mean_C_lambda_opt',formats);
    % mymatlab2tikz(pathname,'convergence_norm_mean_C_lambda_opt.tex');
    
    figure(hmphiC)
    plot([1,N],[nuC_data,nuC_data],'--','Color','k','LineWidth',linewidth)
    legmphiC{end} = '$\nu_{\mathbf{C}}$';
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Mean value of $\varphi(${\boldmath$C$}$)$','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Valeur moyenne de $\varphi(${\boldmath$C$}$)$','Interpreter',interpreter)
    legend(legmphiC{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_mean_phiC_lambda_opt',formats);
    % mymatlab2tikz(pathname,'convergence_mean_phiC_lambda_opt.tex');
    
    figure(herrmC1)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Relative error','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Erreur relative','Interpreter',interpreter)
    legend(legerrmC1{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_error_mean_C1_lambda_opt',formats);
    % mymatlab2tikz(pathname,'convergence_error_mean_C1_lambda_opt.tex');
    
    figure(herrmC2)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Relative error','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Erreur relative','Interpreter',interpreter)
    legend(legerrmC2{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_error_mean_C2_lambda_opt',formats);
    % mymatlab2tikz(pathname,'convergence_error_mean_C2_lambda_opt.tex');
    
    figure(herrmC3)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Relative error','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Erreur relative','Interpreter',interpreter)
    legend(legerrmC3{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_error_mean_C3_lambda_opt',formats);
    % mymatlab2tikz(pathname,'convergence_error_mean_C3_lambda_opt.tex');
    
    figure(herrmCnorm)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Relative error','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Erreur relative','Interpreter',interpreter)
    legend(legerrmCnorm{:},'Location','NorthEast','Interpreter',interpreter)
    mysaveas(pathname,'convergence_error_norm_mean_C_lambda_opt',formats);
    % mymatlab2tikz(pathname,'convergence_error_norm_mean_C_lambda_opt.tex');
    
    figure(herrmphiC)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Relative error','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Erreur relative','Interpreter',interpreter)
    legend(legerrmphiC{:},'Location',location,'Interpreter',interpreter)
    mysaveas(pathname,'convergence_error_mean_phiC_lambda_opt',formats);
    % mymatlab2tikz(pathname,'convergence_error_mean_phiC_lambda_opt.tex');
    
    figure(hresnorm)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Squared norm of residual','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Norme au carr\''e du r\''esidu','Interpreter',interpreter)
    legend(legresnorm{:},'Location','NorthEast','Interpreter',interpreter)
    mysaveas(pathname,'convergence_resnorm_lambda_opt',formats);
    % mymatlab2tikz(pathname,'convergence_resnorm_lambda_opt.tex');
    
    figure(hfval)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xticks(10.^(0:Nexponent))
    xlim tight
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Objective function value','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Valeur de la fonction objectif','Interpreter',interpreter)
    legend(legfval{:},'Location','NorthEast','Interpreter',interpreter)
    mysaveas(pathname,'convergence_fval_lambda_opt',formats);
    % mymatlab2tikz(pathname,'convergence_fval_lambda_opt.tex');
end

%% Pdfs and cdfs
pdf_C4 = @(c4) gampdf(c4,a,b4); % Gamma probability density function of C4
pdf_C5 = @(c5) gampdf(c5,a,b5); % Gamma probability density function of C5
cdf_C4 = @(c4) gamcdf(c4,a,b4); % Gamma cumulative density function of C4
cdf_C5 = @(c5) gamcdf(c5,a,b5); % Gamma cumulative density function of C5

%% Sample generation
rng(s); % initialize the random number generator using the default settings contained in s
N = 1e6; % number of samples
switch lower(MCMCalgo)
    case {'imh','rwmh'}
        [C_sample,accept] = mhsampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),N,MCMCalgo);
    case 'ss'
        [C_sample,neval] = slicesampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),N);
    otherwise
        error(['MCMC algorithm ' MCMC ' not implemented'])
end
C123_sample = C_sample;
C1_sample = C_sample(:,1);
C2_sample = C_sample(:,2);
C3_sample = C_sample(:,3);

C4_sample = gamrnd(a,b4,N,1);
C5_sample = gamrnd(a,b5,N,1);
C_sample(:,4) = C4_sample;
C_sample(:,5) = C5_sample;

kT_sample = C2_sample/2; % [GPa]
NUL_sample = (C3_sample./kT_sample)/(2*sqrt(2));
EL_sample = C1_sample - 4*(NUL_sample.^2).*kT_sample; % [GPa]
GT_sample = C4_sample/2; % [GPa]
GL_sample = C5_sample/2; % [GPa]
ET_sample = 4./(1./kT_sample+1./GT_sample+4*(NUL_sample.^2)./EL_sample); % [GPa]
NUT_sample = (ET_sample./GT_sample)/2-1;

mC123 = mean(C123_sample,1);
vC123 = var(C123_sample,0,1);
% vC123 = N/(N-1)*moment(C123_sample,2,1);
stdC123 = std(C123_sample,0,1);
cvC123 = stdC123./mC123;

[mC45,vC45] = gamstat(a,[b4,b5]);
stdC45 = sqrt(vC45);
cvC45 = stdC45./mC45;

mC   = [mC123 mC45];
vC   = [vC123 vC45];
stdC = [stdC123 stdC45];
cvC  = [cvC123 cvC45];
% stdC = sqrt(vC);
% cvC  = stdC./mC;
% sC = sqrt(norm(vC));
mCnorm = norm(mC);
% dC = sC/mCnorm;

nuC4 = psi(a)+log(b4);
nuC5 = psi(a)+log(b5);
phiC123 = log(C1_sample.*C2_sample-C3_sample.^2);
nuC123 = mean(phiC123,1);
nuC = nuC123 + 2*(nuC4 + nuC5);
fC = [mC nuC];

mC_sample = mean(C_sample,1);
vC_sample = var(C_sample,0,1); % vC_sample = N/(N-1)*moment(C_sample,2,1);
stdC_sample = std(C_sample,0,1);
cvC_sample = stdC_sample./mC_sample;
% sC_sample = sqrt(norm(vC_sample));
mCnorm_sample = norm(mC_sample);
% dC_sample = sC_sample/mCnorm_sample;
phiC_sample = log((C1_sample.*C2_sample-C3_sample.^2).*(C4_sample.^2).*(C5_sample.^2));
% phiC_sample = log(C1_sample.*C2_sample-C3_sample.^2) + 2*log(C4_sample) + 2*log(C5_sample);
nuC_sample = mean(phiC_sample,1);
fC_sample = [mC_sample nuC_sample];

%% Outputs
filenameResults = fullfile(pathname,'results.txt');
fid = fopen(filenameResults,'w');
fprintf(fid,'nb data    = %g\n',size(C_data,1));
fprintf(fid,'nb samples = %g\n',N);
fprintf(fid,'MCMC algo  = %s\n',MCMCalgo);
if useRedParam
    fprintf(fid,'parameterization = reduced\n');
else
    fprintf(fid,'parameterization = full\n');
end

fprintf(fid,'\n');
fprintf(fid,'Initial parameter values\n');
fprintf(fid,'lambda_1 = %.4f\n',la01);
fprintf(fid,'lambda_2 = %.4f\n',la02);
fprintf(fid,'lambda_3 = %.4f\n',la03);
fprintf(fid,'lambda_4 = %.4f\n',la04);
fprintf(fid,'lambda_5 = %.4f\n',la05);
fprintf(fid,'lambda   = %.4f\n',la0);
fprintf(fid,'alpha_4  = alpha_5 = %.4f\n',a0);
fprintf(fid,'beta_4   = %.6f\n',b04);
fprintf(fid,'beta_5   = %.6f\n',b05);
fprintf(fid,'fval     = %.4e\n',fval0);
switch lower(MCMCalgo)
    case {'imh','rwmh'}
        fprintf(fid,'acceptance rate = %.2f%%\n',accept0*100);
        fprintf(fid,'rejection rate  = %.2f%%\n',(1-accept0)*100);
    case 'ss'
        fprintf(fid,'nb function evaluations per sample = %d\n',neval0);
end

fprintf(fid,'\n');
fprintf(fid,'Optimal parameter values\n');
fprintf(fid,'lambda_1 = %.4f\n',la1);
fprintf(fid,'lambda_2 = %.4f\n',la2);
fprintf(fid,'lambda_3 = %.4f\n',la3);
fprintf(fid,'lambda_4 = %.4f\n',la4);
fprintf(fid,'lambda_5 = %.4f\n',la5);
fprintf(fid,'lambda   = %.4f\n',la);
fprintf(fid,'alpha_4  = alpha_5 = %.4f\n',a);
fprintf(fid,'beta_4   = %.6f\n',b4);
fprintf(fid,'beta_5   = %.6f\n',b5);
fprintf(fid,'fval     = %.4e\n',fval);
switch lower(MCMCalgo)
    case {'imh','rwmh'}
        fprintf(fid,'acceptance rate = %.2f%%\n',accept*100);
        fprintf(fid,'rejection rate  = %.2f%%\n',(1-accept)*100);
    case 'ss'
        fprintf(fid,'nb function evaluations per sample = %d\n',neval);
end

for i=1:5
    fprintf(fid,'\n');
    fprintf(fid,'mean(C%u_data)   = %.4f GPa\n',i,mC_data(i));
    fprintf(fid,'mean(C%u)        = %.4f GPa\n',i,mC(i));
    fprintf(fid,'mean(C%u_sample) = %.4f GPa\n',i,mC_sample(i));
    fprintf(fid,'var(C%u_data)    = %.4f (GPa)^2\n',i,vC_data(i));
    fprintf(fid,'var(C%u)         = %.4f (GPa)^2\n',i,vC(i));
    fprintf(fid,'var(C%u_sample)  = %.4f (GPa)^2\n',i,vC_sample(i));
    fprintf(fid,'std(C%u_data)    = %.4f GPa\n',i,stdC_data(i));
    fprintf(fid,'std(C%u)         = %.4f GPa\n',i,stdC(i));
    fprintf(fid,'std(C%u_sample)  = %.4f GPa\n',i,stdC_sample(i));
    fprintf(fid,'cv(C%u_data)     = %.4f\n',i,cvC_data(i));
    fprintf(fid,'cv(C%u)          = %.4f\n',i,cvC(i));
    fprintf(fid,'cv(C%u_sample)   = %.4f\n',i,cvC_sample(i));
end
fprintf(fid,'\n');
fprintf(fid,'norm(mC_data)   = %.4f GPa\n',mCnorm_data);
fprintf(fid,'norm(mC)        = %.4f GPa\n',mCnorm);
fprintf(fid,'norm(mC_sample) = %.4f GPa\n',mCnorm_sample);
fprintf(fid,'\n');
fprintf(fid,'nu(C_data)   = mean(log(det([C_data]))   = %.4f\n',nuC_data);
fprintf(fid,'nu(C)        = mean(log(det([C]))        = %.4f\n',nuC);
fprintf(fid,'nu(C_sample) = mean(log(det([C_sample])) = %.4f\n',nuC_sample);

err_meanCi = abs(mC - mC_data)./abs(mC_data);
err_meanC  = norm(mC - mC_data)/norm(mC_data);
err_nuC    = abs(nuC - nuC_data)/abs(nuC_data);
resnorm    = norm(fC - fC_data)^2;
fval       = norm(mC - mC_data)^2/norm(mC_data)^2 + abs(nuC - nuC_data)^2/abs(nuC_data)^2;

err_meanCi_sample = abs(mC_sample - mC_data)./abs(mC_data);
err_meanC_sample  = norm(mC_sample - mC_data)/norm(mC_data);
err_nuC_sample    = abs(nuC_sample - nuC_data)/abs(nuC_data);
resnorm_sample    = norm(fC_sample - fC_data)^2;
fval_sample       = norm(mC_sample - mC_data)^2/norm(mC_data)^2 + abs(nuC_sample - nuC_data)^2/abs(nuC_data)^2;

fprintf(fid,'\n');
for i=1:5
    fprintf(fid,'relative error on mean(C%u)             = %.4e\n',i,err_meanCi(i));
end
fprintf(fid,'relative error on mean(C) in norm      = %.4e\n',err_meanC);
fprintf(fid,'relative error on mean(log(det([C]))   = %.4e\n',err_nuC);
fprintf(fid,'resnorm = squared norm of residual on C       = %.4e\n',resnorm);
fprintf(fid,'fval = relative squared norm of residual on C = %.4e\n',fval);

fprintf(fid,'\n');
for i=1:5
    fprintf(fid,'relative error on mean(C%u_sample)             = %.4e\n',i,err_meanCi_sample(i));
end
fprintf(fid,'relative error on mean(C_sample) in norm      = %.4e\n',err_meanC_sample);
fprintf(fid,'relative error on mean(log(det([C_sample]))   = %.4e\n',err_nuC_sample);
fprintf(fid,'resnorm_sample = squared norm of residual on C_sample       = %.4e\n',resnorm_sample);
fprintf(fid,'fval_sample = relative squared norm of residual on C_sample = %.4e\n',fval_sample);
fclose(fid);
type(filenameResults) % fprintf('%s', fileread(filenameResults))

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

    figure('Name','Data for C1')
    clf
    bar(C1_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Coefficient $C_1$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    mysaveas(pathname,'data_C1','fig');
    mymatlab2tikz(pathname,'data_C1.tex');
    
    figure('Name','Data for C2')
    clf
    bar(C2_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Coefficient $C_2$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    mysaveas(pathname,'data_C2','fig');
    mymatlab2tikz(pathname,'data_C2.tex');
    
    figure('Name','Data for C3')
    clf
    bar(C3_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Coefficient $C_3$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    mysaveas(pathname,'data_C3','fig');
    mymatlab2tikz(pathname,'data_C3.tex');
    
    figure('Name','Data for C4')
    clf
    bar(C4_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Coefficient $C_4$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    mysaveas(pathname,'data_C4','fig');
    mymatlab2tikz(pathname,'data_C4.tex');
    
    figure('Name','Data for C5')
    clf
    bar(C5_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Coefficient $C_5$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    mysaveas(pathname,'data_C5','fig');
    mymatlab2tikz(pathname,'data_C5.tex');
    
    %% Plot pdfs and cdfs
    npts = 1e3;
    [f1,xi1,bw1] = ksdensity(C1_sample,'NumPoints',npts);
    [f2,xi2,bw2] = ksdensity(C2_sample,'NumPoints',npts);
    [f3,xi3,bw3] = ksdensity(C3_sample,'NumPoints',npts);
    [f4,xi4,bw4] = ksdensity(C4_sample,'NumPoints',npts);
    [f5,xi5,bw5] = ksdensity(C5_sample,'NumPoints',npts);
    
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
    
    c4_min = max(0,mean(C4_data)-5*std(C4_data));
    c4_max = mean(C4_data)+5*std(C4_data);
    c5_min = max(0,mean(C5_data)-5*std(C5_data));
    c5_max = mean(C5_data)+5*std(C5_data);
    
    c4 = linspace(c4_min,c4_max,npts);
    c5 = linspace(c5_min,c5_max,npts);
    
    % Plot pdf of C4
    figure('Name','Probability density function of C4')
    clf
    plot(c4,pdf_C4(c4),'-b','LineWidth',linewidth);
    % hold on
    % plot(C4_data,pdf_C4(C4_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(c4),max(c4)])
    xlabel('$c_4$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_4}(c_4)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C4',formats);
    mymatlab2tikz(pathname,'pdf_C4.tex');
    
    % Plot pdf of C5
    figure('Name','Probability density function of C5')
    clf
    plot(c5,pdf_C5(c5),'-b','LineWidth',linewidth);
    % hold on
    % plot(C5_data,pdf_C5(C5_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(c5),max(c5)])
    xlabel('$c_5$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_5}(c_5)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C5',formats);
    mymatlab2tikz(pathname,'pdf_C5.tex');
    
    % Plot cdf of C4
    figure('name','cumulative distribution function of C4')
    clf
    plot(c4,cdf_C4(c4),'-r','LineWidth',linewidth);
    % hold on
    % plot(C4_data,cdf_C4(C4_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(c4),max(c4)])
    xlabel('$c_4$ [GPa]','Interpreter',interpreter)
    ylabel('$F_{C_4}(c_4)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_C4',formats);
    mymatlab2tikz(pathname,'cdf_C4.tex');
    
    % Plot cdf of C5
    figure('name','cumulative distribution function of C5')
    clf
    plot(c5,cdf_C5(c5),'-r','LineWidth',linewidth);
    % hold on
    % plot(C5_data,cdf_C5(C5_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(c5),max(c5)])
    xlabel('$c_5$ [GPa]','Interpreter',interpreter)
    ylabel('$F_{C_5}(c_5)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_C5',formats);
    mymatlab2tikz(pathname,'cdf_C5.tex');
    
    %% Plot samples
    if (strcmpi(MCMCalgo,'imh') || strcmpi(MCMCalgo,'rwmh')) && accept~=1
        deltas = diff(C_sample); % (N-1)xd
        isAccept = [true; any(deltas~=0,2)]; % new sample row is accepted if any coordinate changed
        C_sample = C_sample(isAccept,:); % accepted samples
        Naccept = size(C_sample,1); % number of accepted samples
        Nmax = min(1e4,Naccept);
    else
        Nmax = min(1e4,N);
    end

    % c1min = min(min(C1_sample),min(C1_data));
    c1min = 0;
    c1max = max(max(C1_sample),max(C1_data));
    c1 = linspace(c1min,c1max,1e2);
    
    % c2min = min(min(C2_sample),min(C2_data));
    c2min = 0;
    c2max = max(max(C2_sample),max(C2_data));
    c2 = linspace(c2min,c2max,1e2);
    
    [C1,C2] = meshgrid(c1,c2);
    C3 = sqrt(C1.*C2);
    
    figure('name','Samples of (C_1,C_2,C_3)')
    clf
    hold on
    % surf(C1,C2,C3,'FaceAlpha',0.5)
    % surf(C1,C2,-C3,'FaceAlpha',0.5)
    mesh(C1,C2,C3)
    mesh(C1,C2,-C3)
    scatter3(C_sample(1:Nmax,1),C_sample(1:Nmax,2),C_sample(1:Nmax,3),'b.')
    scatter3(C1_data,C2_data,C3_data,'r+','LineWidth',linewidth)
    hold off
    % view(-37.5,30) % default view
    view(-30,20)
    grid on
    box on
    set(gca,'Xdir','reverse','Ydir','reverse','FontSize',fontsize)
    xlabel('$c_1$ [GPa]','Interpreter',interpreter);
    ylabel('$c_2$ [GPa]','Interpreter',interpreter);
    zlabel('$c_3$ [GPa]','Interpreter',interpreter);
    % legend('realizations','data','support');
    %legend('réalisations','données','support');
    mysaveas(pathname,'samples_C1_C2_C3',formats);
    % mymatlab2tikz(pathname,'samples_C1_C2_C3.tex');
    
    figure('name','Samples of (C_4,C_5)')
    clf
    scatter(C_sample(1:Nmax,4),C_sample(1:Nmax,5),'b.')
    hold on
    scatter(C4_data,C5_data,'r+','LineWidth',linewidth)
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
    scatter3(ET_sample(1:Nmax),GL_sample(1:Nmax)*1e3,EL_sample(1:Nmax)*1e3,'b.')
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
    scatter(NUL_sample(1:Nmax),NUT_sample(1:Nmax),'b.')
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
