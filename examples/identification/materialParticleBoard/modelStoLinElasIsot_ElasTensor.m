%% Stochastic modeling of random linear elasticity tensor %%
%  with isotropic symmetry                                %%
%%--------------------------------------------------------%%

% clc
clearvars
close all
rng('default');

%% Input data
displaySolution = true;

useRedParam = true; % reduced parameterization
method = 'mle'; % parameter estimation method = 'mle' or 'lse'

filename = 'modelStoLinElasIsot_ElasTensor_';
if useRedParam
    filename = [filename 'ReducedParam_'];
else
    filename = [filename 'FullParam_'];
end
filename = [filename method];
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
E_data = mean_ET_data*1e-3; % [GPa]
NU_data = 0.1+0.2*rand(N_data,1); % artificial data for NU uniformly distributed from 0.1 to 0.3
G_data = E_data./(2*(1+NU_data)); % [GPa]
lambda_data = E_data.*NU_data./((1+NU_data).*(1-2*NU_data)); % [GPa]
C1_data = lambda_data + 2/3*G_data; % [GPa]
C2_data = G_data; % [GPa]
C_data = [C1_data(:) C2_data(:)];

% Empirical estimates
mC_data = mean(C_data,1);
vC_data = var(C_data,0,1); % vC_data = size(C_data,1)/(size(C_data,1)-1)*moment(C_data,2,1);
stdC_data = std(C_data,0,1); % stdC_data = sqrt(vC_data);
cvC_data = stdC_data./mC_data;
% sC_data = sqrt(norm(vC_data));
mCnorm_data = norm(mC_data);
% dC_data = sC_data/mCnorm_data;
phiC_data = log(96*C_data(:,1).*C_data(:,2).^5);
nuC_data = mean(phiC_data,1);
fC_data = [mC_data nuC_data];

%% Parameter estimation for computing Lagrange multipliers
% Initial parameter values
la0 = mean([1-1/cvC_data(1)^2 (1-1/cvC_data(2)^2)/5]); % la < 1/5 (match cv with averaged empirical cv for C1 and C2)
a01 = 1-la0;   % a1 > 0
a02 = 1-5*la0; % a2 > 0
la01 = a01/mC_data(1); % la1 > 0 (match mathematical expectation with empirical mean value for C1)
la02 = a02/mC_data(2); % la2 > 0 (match mathematical expectation with empirical mean value for C2)
% la01 = -la0/mC_data(1);   % la1 > 0 (match mode with empirical mean value for C1)
% la02 = -5*la0/mC_data(2); % la2 > 0 (match mode with empirical mean value for C2)
b01 = 1/la01;  % b1 > 0
b02 = 1/la02;  % b2 > 0

% Initial parameter vector
if useRedParam
    lambda0 = la0; % reduced parameterization
else
    lambda0 = [la01 la02 la0]; % full parameterization
end

% Parameter estimation
t = tic;
switch lower(method)
    case 'mle'
        % Maximum likelihood estimation
        nloglf = @nloglfElasIsot;
        lambda = mleStoLinElasIsot(C_data,lambda0,'display','iter');
        % [lambda,loglfval,loglfval0] = mleStoLinElasIsot(C_data,lambda0,'display','iter');
    case 'lse'
        % Least-squares estimation
        lambda = lseStoLinElasIsot(C_data,lambda0,'display','iter-detailed'); 
        % [lambda,resnorm,residual,exitflag,output] = lseStoLinElasIsot(C_data,lambda0,'display','iter-detailed');
end
toc(t)

% Optimal parameter values
if useRedParam
    la  = lambda(1);     % la < 1/5
    a1  = 1-la;          % a1 > 0
    a2  = 1-5*la;        % a2 > 0
    la1 = a1/mC_data(1); % la1 > 0
    la2 = a2/mC_data(2); % la2 > 0
else
    la1 = lambda(1); % la1 > 0
    la2 = lambda(2); % la2 > 0
    la  = lambda(3); % la < 1/5
    a1  = 1-la;      % a1 > 0
    a2  = 1-5*la;    % a2 > 0
end
b1 = 1/la1; % b1 > 0
b2 = 1/la2; % b2 > 0

[mC,vC] = gamstat([a1,a2],[b1,b2]);
% mC = [a1*b1, a2*b2];
% vC = [a1*b1^2, a2*b2^2];
stdC = sqrt(vC);
cvC = stdC./mC;
% sC = sqrt(norm(vC));
mCnorm = norm(mC);
% dC = sC/mCnorm;
nuC = log(96) + psi(a1)+log(b1) + 5*(psi(a2)+log(b2));
fC = [mC nuC];

if strcmpi(method,'mle')
    if ~exist('loglfval0','var')
        loglfval0 = -nloglf(lambda0,C_data);
    end
    if ~exist('loglfval','var')
        loglfval = -nloglf(lambda,C_data);
    end
end

%% Pdfs and cdfs
k1ln = -a1*log(b1)-gammaln(a1);
k2ln = -a2*log(b2)-gammaln(a2);
% k1ln = (1-la)*log(la1)-gammaln(1-la);
% k2ln = (1-5*la)*log(la2)-gammaln(1-5*la);

pdf_C1 = @(c1) gampdf(c1,a1,b1); % Gamma probability density function of C1
pdf_C2 = @(c2) gampdf(c2,a2,b2); % Gamma probability density function of C2
pdf_C = @(c1,c2) pdf_C1(c1).*pdf_C2(c2); % joint probability density function of C=(C1,C2)
% pdf_C = @(c1,c2) exp(k1ln+k2ln -la*log(c1) -5*la*log(c2) -la1*c1 -la2*c2);
% pdf_C = @(c1,c2) exp(k1ln+k2ln)*c1.^(-la)*c2.^(-5*la)*exp(-la1*c1-la2*c2); % it does not work due to the value of la
pdf_EN = @(e,n) exp(k1ln+k2ln -la*(log(e)-log(3*(1-2*n))) -5*la*(log(e)-log(2*(1+n))) -la1*e./(3*(1-2*n)) -la2*e./(2*(1+n)))...
    .*(e./(2*((1+n).^2).*((1-2*n).^2))); % joint probability density function of (E,N)
% pdf_EN = @(e,n) exp(k1ln+k2ln)*(e./(3*(1-2*n))).^(-la).*(e./(2*(1+n))).^(-5*la)...
%     .*(e./(2*((1+n).^2).*((1-2*n).^2)))...
%     .*exp(-la1*e./(3*(1-2*n))-la2*e./(2*(1+n)));

cdf_C1 = @(c1) gamcdf(c1,a1,b1); % Gamma cumulative density function of C1
cdf_C2 = @(c2) gamcdf(c2,a2,b2); % Gamma cumulative density function of C2
cdf_C = @(c1,c2) cdf_C1(c1).*cdf_C2(c2); % joint cumulative density function of C=(C1,C2)
% cdf_C1 = @(c1) arrayfun(@(x) integral(pdf_C1,0,x), c1);
% cdf_C2 = @(c2) arrayfun(@(x) integral(pdf_C2,0,x), c2);
% cdf_C = @(c1,c2) arrayfun(@(x1,x2) integral2(pdf_C,0,x1,0,x2), c1, c2);
cdf_EN = @(e,n) arrayfun(@(xe,xn) integral2(pdf_EN,0,xe,-1,xn), e, n);

%% Sample generation
N = 1e4; % number of samples
C1_sample = gamrnd(a1,b1,N,1); % [GPa]
C2_sample = gamrnd(a2,b2,N,1); % [GPa]
C_sample = [C1_sample(:) C2_sample(:)];

lambda_sample = C1_sample-2/3*C2_sample; % [GPa]
E_sample = (9*C1_sample.*C2_sample)./(3*C1_sample+C2_sample); % [GPa]
NU_sample = (3*C1_sample-2*C2_sample)./(6*C1_sample+2*C2_sample); % [GPa]

mC_sample = mean(C_sample,1);
vC_sample = var(C_sample,0,1); % vC_sample = N/(N-1)*moment(C_sample,2,1);
stdC_sample = std(C_sample,0,1); % stdC_sample = sqrt(vC_sample);
cvC_sample = stdC_sample./mC_sample;
sC_sample = sqrt(norm(vC_sample));
mCnorm_sample = norm(mC_sample);
dC_sample = sC_sample/mCnorm_sample;
phiC_sample = log(96*C1_sample.*C2_sample.^5);
nuC_sample = mean(phiC_sample,1);
fC_sample = [mC_sample nuC_sample];

%% Outputs
filenameResults = fullfile(pathname,'results.txt');
fid = fopen(filenameResults,'w');
fprintf(fid,'nb data    = %g\n',size(C_data,1));
fprintf(fid,'nb samples = %g\n',N);
fprintf(fid,'estimation method = %s\n',method);
if useRedParam
    fprintf(fid,'parameterization  = reduced\n');
else
    fprintf(fid,'parameterization  = full\n');
end

fprintf(fid,'\n');
fprintf(fid,'Initial parameter values\n');
fprintf(fid,'lambda_1 = %.4f\n',la01);
fprintf(fid,'lambda_2 = %.4f\n',la02);
fprintf(fid,'lambda   = %.4f\n',la0);
fprintf(fid,'alpha_1  = %.4f\n',a01);
fprintf(fid,'beta_1   = %.6f\n',b01);
fprintf(fid,'alpha_2  = %.4f\n',a02);
fprintf(fid,'beta_2   = %.6f\n',b02);
if strcmpi(method,'mle')
    fprintf(fid,'loglf    = %.4f\n',loglfval0);
end

fprintf(fid,'\n');
fprintf(fid,'Optimal parameter values\n');
fprintf(fid,'lambda_1 = %.4f\n',la1);
fprintf(fid,'lambda_2 = %.4f\n',la2);
fprintf(fid,'lambda   = %.4f\n',la);
fprintf(fid,'alpha_1  = %.4f\n',a1);
fprintf(fid,'beta_1   = %.6f\n',b1);
fprintf(fid,'alpha_2  = %.4f\n',a2);
fprintf(fid,'beta_2   = %.6f\n',b2);
if strcmpi(method,'mle')
    fprintf(fid,'loglf    = %.4f\n',loglfval);
end

for i=1:2
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
relresnorm = norm((fC - fC_data)./fC_data)^2;
weightresnorm = norm(mC - mC_data)^2/norm(mC_data)^2 + abs(nuC - nuC_data)^2/abs(nuC_data)^2;

err_meanCi_sample = abs(mC_sample - mC_data)./abs(mC_data);
err_meanC_sample  = norm(mC_sample - mC_data)/norm(mC_data);
err_nuC_sample    = abs(nuC_sample - nuC_data)/abs(nuC_data);
resnorm_sample    = norm(fC_sample - fC_data)^2;
relresnorm_sample = norm((fC_sample - fC_data)./fC_data)^2;
weightresnorm_sample = norm(mC_sample - mC_data)^2/norm(mC_data)^2 + abs(nuC_sample - nuC_data)^2/abs(nuC_data)^2;

fprintf(fid,'\n');
for i=1:2
    fprintf(fid,'relative error on mean(C%u)             = %.4e\n',i,err_meanCi(i));
end
fprintf(fid,'relative error on mean(C) in norm      = %.4e\n',err_meanC);
fprintf(fid,'relative error on mean(log(det([C]))   = %.4e\n',err_nuC);
fprintf(fid,'squared norm of residual on C          = %.4e\n',resnorm);
fprintf(fid,'squared norm of relative residual on C = %.4e\n',relresnorm);
fprintf(fid,'weighted squared norm of residual on C = %.4e\n',weightresnorm);

fprintf(fid,'\n');
for i=1:2
    fprintf(fid,'relative error on mean(C%u_sample)             = %.4e\n',i,err_meanCi_sample(i));
end
fprintf(fid,'relative error on mean(C_sample) in norm      = %.4e\n',err_meanC_sample);
fprintf(fid,'relative error on mean(log(det([C_sample]))   = %.4e\n',err_nuC_sample);
fprintf(fid,'squared norm of residual on C_sample          = %.4e\n',resnorm_sample);
fprintf(fid,'squared norm of relative residual on C_sample = %.4e\n',relresnorm_sample);
fprintf(fid,'weighted squared norm of residual on C_sample = %.4e\n',weightresnorm_sample);
fclose(fid);
type(filenameResults) % fprintf('%s', fileread(filenameResults))

%% Display
if displaySolution
    %% Plot data
    figure('Name','Data for E')
    clf
    bar(E_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Young''s modulus $E$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module d''Young $E$ [GPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_E','fig');
    mymatlab2tikz(pathname,'data_E.tex');
    
    figure('Name','Data for NU')
    clf
    bar(NU_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Poisson''s ratio $\nu$','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Coefficient de Poisson $\nu$','Interpreter',interpreter);
    mysaveas(pathname,'data_NU','fig');
    mymatlab2tikz(pathname,'data_NU.tex');
    
    figure('Name','Data for G')
    clf
    bar(G_data)
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Shear modulus $G$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module de cisaillement $G$ [GPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_G','fig');
    mymatlab2tikz(pathname,'data_G.tex');
    
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
    
    %% Plot log-likelihood function
    if strcmpi(method,'mle')
        npts = 1e2;
        if useRedParam
            loglf = @(La,data) arrayfun(@(la) -nloglf(la,data), La1, La2, La);
            la_delta = 0.5*abs(la);
            la_series = linspace(la-la_delta,la+la_delta,npts);
            loglf = @(La,data) arrayfun(@(la) -nloglf(la,data), La);
            % logL = zeros(size(la));
            % for i=1:length(la)
            %     la_i = la_series(i);
            %     logL(i) = -nloglf(la_i,C_data);
            % end
            logL = loglf(la_series,C_data);
            
            % Plot log-likelihood function loglf for C
            plot(la_series,logL,'-b','LineWidth',linewidth);
            hold on
            scatter(la,loglfval,'o','MarkerEdgeColor','k','MarkerFaceColor','r');
            % scatter(la0,loglfval0,'o','MarkerEdgeColor','k','MarkerFaceColor','g');
            % % plot(la,loglfval,'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','r');
            % % plot(la0,loglfval0,'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
            hold off
            grid on
            box on
            set(gca,'FontSize',fontsize)
            xlabel('$\lambda$','Interpreter',interpreter)
            ylabel('$\mathcal{L}(\lambda_1,\lambda_2,\lambda)$','Interpreter',interpreter)
        else
            loglf = @(La1,La2,La,data) arrayfun(@(la1,la2,la) -nloglf([la1,la2,la],data), La1, La2, La);
            la_delta  = 0.5*abs(la);
            la1_delta = 0.5*abs(la1);
            la2_delta = 0.5*abs(la2);
            la_series = linspace(la-la_delta,la+la_delta,npts);
            la1_series = linspace(la1-la1_delta,la1+la1_delta,npts);
            la2_series = linspace(la2-la2_delta,la2+la2_delta,npts);
            n1 = length(la1_series);
            n2 = length(la2_series);
            n3 = length(la_series);
            % logL12 = zeros(n2,n1);
            % for i=1:n1
            %     la1_i = la1_series(i);
            %     for j=1:n2
            %         la2_j = la2_series(j);
            %         param12_ij = [la1_i,la2_j,la];
            %         logL12(j,i) = -nloglf(param12_ij,C_data);
            %     end
            % end
            % logL13 = zeros(n3,n1);
            % for i=1:n1
            %     la1_i = la1_series(i);
            %     for j=1:n3
            %         la_j = la_series(j);
            %         param13_ij = [la1_i,la2,la_j];
            %         logL13(j,i) = -nloglf(param13_ij,C_data);
            %     end
            % end
            % logL23 = zeros(n3,n2);
            % for i=1:n2
            %     la2_i = la2_series(i);
            %     for j=1:n3
            %         la_j = la_series(j);
            %         param23_ij = [la1,la2_i,la_j];
            %         logL23(j,i) = -nloglf(param23_ij,C_data);
            %     end
            % end
            [La1,La2] = meshgrid(la1_series,la2_series);
            logL12 = loglf(La1,La2,repmat(la,n2,n1),C_data);
            [La1,La] = meshgrid(la1_series,la_series);
            logL13 = loglf(La1,repmat(la2,n3,n1),La,C_data);
            [La2,La] = meshgrid(la2_series,la_series);
            logL23 = loglf(repmat(la1,n3,n2),La2,La,C_data);
            
            % Plot log-likelihood function loglf_la1_la2 for C
            figure('Name','Surface plot: Log-likelihood function for C')
            clf
            surfc(la1_series,la2_series,logL12,'EdgeColor','none');
            colorbar
            hold on
            scatter3(la1,la2,loglfval,'MarkerEdgeColor','k','MarkerFaceColor','r');
            % scatter3(la01,la02,loglfval0,'MarkerEdgeColor','k','MarkerFaceColor','g');
            hold off
            set(gca,'FontSize',fontsize)
            xlabel('$\lambda_1$ [GPa$^{-1}$]','Interpreter',interpreter)
            ylabel('$\lambda_2$ [GPa$^{-1}$]','Interpreter',interpreter)
            zlabel('$\mathcal{L}(\lambda_1,\lambda_2,\lambda)$','Interpreter',interpreter)
            mysaveas(pathname,'loglf_12_3D',formats);
            % mymatlab2tikz(pathname,'loglf_12_3D.tex');
            
            figure('Name','Contour plot: Log-likelihood function for C')
            clf
            contourf(la1_series,la2_series,logL12,50);
            colorbar
            hold on
            scatter(la1,la2,'MarkerEdgeColor','k','MarkerFaceColor','r');
            % scatter(la01,la02,'MarkerEdgeColor','k','MarkerFaceColor','g');
            hold off
            set(gca,'FontSize',fontsize)
            xlabel('$\lambda_1$ [GPa$^{-1}$]','Interpreter',interpreter)
            ylabel('$\lambda_2$ [GPa$^{-1}$]','Interpreter',interpreter)
            zlabel('$\mathcal{L}(\lambda_1,\lambda_2,\lambda)$','Interpreter',interpreter)
            mysaveas(pathname,'loglf_12_2D',formats);
            % mymatlab2tikz(pathname,'loglf_12_2D.tex');
            
            % Plot log-likelihood function loglf_la1_la for C
            figure('Name','Surface plot: Log-likelihood function for C')
            clf
            surfc(la1_series,la_series,logL13,'EdgeColor','none');
            colorbar
            hold on
            scatter3(la1,la,loglfval,'MarkerEdgeColor','k','MarkerFaceColor','r');
            % scatter3(la01,la0,loglfval0,'MarkerEdgeColor','k','MarkerFaceColor','g');
            hold off
            set(gca,'FontSize',fontsize)
            xlabel('$\lambda_1$ [GPa$^{-1}$]','Interpreter',interpreter)
            ylabel('$\lambda$','Interpreter',interpreter)
            zlabel('$\mathcal{L}(\lambda_1,\lambda_2,\lambda)$','Interpreter',interpreter)
            mysaveas(pathname,'loglf_13_3D',formats);
            % mymatlab2tikz(pathname,'loglf_13_3D.tex');
            
            figure('Name','Contour plot: Log-likelihood function for C')
            clf
            contourf(la1_series,la_series,logL13,50);
            colorbar
            hold on
            scatter(la1,la,'MarkerEdgeColor','k','MarkerFaceColor','r');
            % scatter(la01,la0,'MarkerEdgeColor','k','MarkerFaceColor','g');
            hold off
            set(gca,'FontSize',fontsize)
            xlabel('$\lambda_1$ [GPa$^{-1}$]','Interpreter',interpreter)
            ylabel('$\lambda$','Interpreter',interpreter)
            zlabel('$\mathcal{L}(\lambda_1,\lambda_2,\lambda)$','Interpreter',interpreter)
            mysaveas(pathname,'loglf_13_2D',formats);
            % mymatlab2tikz(pathname,'loglf_13_2D.tex');
            
            % Plot log-likelihood function loglf_la2_la for C
            figure('Name','Surface plot: Log-likelihood function for C')
            clf
            surfc(la2_series,la_series,logL23,'EdgeColor','none');
            colorbar
            hold on
            scatter3(la2,la,loglfval,'MarkerEdgeColor','k','MarkerFaceColor','r');
            % scatter3(la02,la0,loglfval0,'MarkerEdgeColor','k','MarkerFaceColor','g');
            hold off
            set(gca,'FontSize',fontsize)
            xlabel('$\lambda_2$ [GPa$^{-1}$]','Interpreter',interpreter)
            ylabel('$\lambda$','Interpreter',interpreter)
            zlabel('$\mathcal{L}(\lambda_1,\lambda_2,\lambda)$','Interpreter',interpreter)
            mysaveas(pathname,'loglf_23_3D',formats);
            % mymatlab2tikz(pathname,'loglf_23_3D.tex');
            
            figure('Name','Contour plot: Log-likelihood function for C')
            clf
            contourf(la2_series,la_series,logL23,50);
            colorbar
            hold on
            scatter(la2,la,'MarkerEdgeColor','k','MarkerFaceColor','r');
            % scatter(la02,la0,'MarkerEdgeColor','k','MarkerFaceColor','g');
            hold off
            set(gca,'FontSize',fontsize)
            xlabel('$\lambda_2$ [GPa$^{-1}$]','Interpreter',interpreter)
            ylabel('$\lambda$','Interpreter',interpreter)
            zlabel('$\mathcal{L}(\lambda_1,\lambda_2,\lambda)$','Interpreter',interpreter)
            mysaveas(pathname,'loglf_23_2D',formats);
            % mymatlab2tikz(pathname,'loglf_23_2D.tex');
        end
    end

    %% Plot pdfs and cdfs
    c1_min = max(0,mean(C1_data)-5*std(C1_data));
    c1_max = mean(C1_data)+5*std(C1_data);
    c2_min = max(0,mean(C2_data)-5*std(C2_data));
    c2_max = mean(C2_data)+5*std(C2_data);
    
    c1 = linspace(c1_min,c1_max,1e2);
    c2 = linspace(c2_min,c2_max,1e2);
    [C1,C2] = meshgrid(c1,c2);
    
    e_min = max(0,mean(E_data)-5*std(E_data));
    e_max = mean(E_data)+5*std(E_data);
    nu_min = max(-1,mean(NU_data)-5*std(NU_data));
    nu_max = min(0.5,mean(NU_data)+5*std(NU_data));
    
    e = linspace(e_min,e_max,1e2);
    nu = linspace(nu_min,nu_max,1e2);
    [E,NU] = meshgrid(e,nu);
    
    % Plot pdf of C1
    figure('Name','Probability density function of C1')
    clf
    plot(c1,pdf_C1(c1),'-b','LineWidth',linewidth);
    % hold on
    % plot(C1_data,pdf_C1(C1_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(c1),max(c1)])
    xlabel('$c_1$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_1}(c_1)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C1',formats);
    mymatlab2tikz(pathname,'pdf_C1.tex');
    
    % Plot pdf of C2
    figure('Name','Probability density function of C2')
    clf
    plot(c2,pdf_C2(c2),'-b','LineWidth',linewidth);
    % hold on
    % plot(C2_data,pdf_C2(C2_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(c2),max(c2)])
    xlabel('$c_2$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_2}(c_2)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C2',formats);
    mymatlab2tikz(pathname,'pdf_C2.tex');
    
    % Plot pdf of C=(C1,C2)
    figure('Name','Probability density function of C=(C1,C2)')
    clf
    surf(C1,C2,pdf_C(C1,C2));
    % hold on
    % plot3(C1_data,C2_data,pdf_C(C1_data,C2_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(c1),max(c1)])
    set(gca,'YLim',[min(c2),max(c2)])
    xlabel('$c_1$ [GPa]','Interpreter',interpreter)
    ylabel('$c_2$ [GPa]','Interpreter',interpreter)
    zlabel('$p_{(C_1,C_2)}(c_1,c_2)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C',formats);
    mymatlab2tikz(pathname,'pdf_C.tex');
    
    % Plot pdf of (E,N)
    figure('Name','Probability density function of (E,N)')
    clf
    surf(E,NU,pdf_EN(E,NU));
    % hold on
    % plot3(E_data,NU_data,pdf_EN(E_data,NU_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(e),max(e)])
    set(gca,'YLim',[min(nu),max(nu)])
    xlabel('$e$ [GPa]','Interpreter',interpreter)
    ylabel('$\nu$','Interpreter',interpreter)
    zlabel('$p_{(E,N)}(e,\nu)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_E_N',formats);
    mymatlab2tikz(pathname,'pdf_E_N.tex');
    
    % Plot cdf of C1
    figure('name','cumulative distribution function of C1')
    clf
    plot(c1,cdf_C1(c1),'-r','LineWidth',linewidth);
    % hold on
    % plot(C1_data,cdf_C1(C1_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(c1),max(c1)])
    xlabel('$c_1$ [GPa]','Interpreter',interpreter)
    ylabel('$F_{C_1}(c_1)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_C1',formats);
    mymatlab2tikz(pathname,'cdf_C1.tex');
    
    % Plot cdf of C2
    figure('name','cumulative distribution function of C2')
    clf
    plot(c2,cdf_C2(c2),'-r','LineWidth',linewidth);
    % hold on
    % plot(C2_data,cdf_C2(C2_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(c2),max(c2)])
    xlabel('$c_2$ [GPa]','Interpreter',interpreter)
    ylabel('$F_{C_2}(c_2)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_C2',formats);
    mymatlab2tikz(pathname,'cdf_C2.tex');
    
    % Plot cdf of C=(C1,C2)
    figure('Name','Cumulative density function of C=(C1,C2)')
    clf
    % p12 = zeros(size(C1));
    % for i=1:length(c1)
    %     for j=1:length(c2)
    %         p12(j,i) = cdf_C(c1(i),c2(j));
    %     end
    % end
    % p12 = cdf_C(C1,C2);
    surf(C1,C2,cdf_C(C1,C2));
    % hold on
    % plot3(C1_data,C2_data,cdf_C(C1_data,C2_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(c1),max(c1)])
    set(gca,'YLim',[min(c2),max(c2)])
    xlabel('$c_1$ [GPa]','Interpreter',interpreter)
    ylabel('$c_2$ [GPa]','Interpreter',interpreter)
    zlabel('$F_{(C_1,C_2)}(c_1,c_2)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_C',formats);
    mymatlab2tikz(pathname,'cdf_C.tex');
    
    % Plot cdf of (E,N)
    figure('Name','Cumulative density function of (E,N)')
    clf
    % pEN = zeros(size(Xe));
    % for i=1:length(e)
    %     for j=1:length(n)
    %         pEN(j,i) = cdf_EN(e(i),n(j));
    %     end
    % end
    % pEN = cdf_EN(E,NU);
    surf(E,NU,cdf_EN(E,NU));
    % hold on
    % plot3(E_data,NU_data,cdf_EN(E_data,NU_data),'k+','LineWidth',linewidth);
    % hold off
    grid on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(e),max(e)])
    set(gca,'YLim',[min(nu),max(nu)])
    xlabel('$e$ [GPa]','Interpreter',interpreter)
    ylabel('$\nu$','Interpreter',interpreter)
    zlabel('$F_{(E,N)}(e,\nu)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_E_N',formats);
    mymatlab2tikz(pathname,'cdf_E_N.tex');

    %% Plot samples
    figure('name','Samples of C=(C_1,C_2)')
    clf
    scatter(C1_sample,C2_sample,'b.')
    hold on
    scatter(C1_data,C2_data,'r+','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$c_1$ [GPa]','Interpreter',interpreter)
    ylabel('$c_2$ [GPa]','Interpreter',interpreter)
    legend('realizations','data')
    %legend('réalisations','données')
    mysaveas(pathname,'samples_C1_C2',formats);
    % mymatlab2tikz(pathname,'samples_C1_C2.tex');
    
    figure('name','Samples of (E,N)')
    clf
    scatter(NU_sample,E_sample,'b.')
    hold on
    scatter(NU_data,E_data,'r+','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    %xlabel('Poisson''s ratio $\nu$','Interpreter',interpreter)
    %ylabel('Young''s modulus $E$ [GPa]','Interpreter',interpreter)
    xlabel('$\nu$','Interpreter',interpreter)
    ylabel('$E$ [GPa]','Interpreter',interpreter)
    legend('realizations','data')
    %legend('réalisations','données')
    mysaveas(pathname,'samples_E_N',formats);
    % mymatlab2tikz(pathname,'samples_E_N.tex');
    
end
