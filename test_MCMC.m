% clc
clearvars
% close all
rng('default');
s = rng;

%% Input data
useRedParam = true; % reduced parameterization
MCMCalgo = 'IMH'; % algorithm for Markov Chain Monte Carlo (MCMC) method = 'IMH', 'RWMH' or 'SS'
filename = 'modelStoLinElasIsotTrans_ElasTensor_';
if useRedParam
    filename = [filename 'Reduced_'];
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
formats = {'epsc'};

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
% vC_data = var(C_data,0,1); % vC_data = size(C_data,1)/(size(C_data,1)-1)*moment(C_data,2,1);
stdC_data = std(C_data,0,1); % stdC_data = sqrt(vC_data);
cvC_data = stdC_data./mC_data;
% sC_data = sqrt(norm(vC_data));
mCnorm_data = norm(mC_data);
% dC_data = sC_data/norm(mC_data);
phiC_data = log((C1_data.*C2_data-C3_data.^2).*(C4_data.^2).*(C5_data.^2));
nuC_data = mean(phiC_data,1);

% Initial guess
rng(s);
cvC45 = mean(cvC_data(4:5));
las = (1-1/(cvC45.^2))/2; % la < 1/2, la < 0 for MH sampling using a trivariate normal proposal pdf with positive-definite covariance matrix Sigma
% las = -100; % la < 1/2, la < 0 for MH sampling using a trivariate normal proposal pdf with positive-definite covariance matrix Sigma
% las = -190:10:-10;

residualnorms        = zeros(length(las),1);
residualnorms_sample = zeros(length(las),1);
funvals              = zeros(length(las),1);
funvals_sample       = zeros(length(las),1);
errors_mCnorm        = zeros(length(las),1);
errors_mCnorm_sample = zeros(length(las),1);
errors               = zeros(length(las),6);
errors_sample        = zeros(length(las),6);
for k=1:length(las)
% Initial parameter values
la = las(k);

la1 = -(mC_data(2)*la)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la1 > 0 (match mode with empirical mean value for C1)
la2 = -(mC_data(1)*la)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la2 > 0 (match mode with empirical mean value for C2)
la3 = (2*mC_data(3)*la)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la3 in R such that 2*sqrt(la1*la2)-la3 > 0 (match mode with empirical mean value for C3)

a = 1-2*la; % a > 0
la4 = a/mC_data(4); % la4 > 0 (match mathematical expectation with empirical mean value for C4)
la5 = a/mC_data(5); % la5 > 0 (match mathematical expectation with empirical mean value for C5)
% la4 = -2*la/mC_data(4); % la4 > 0 (match mode witch empirical mean value for C4)
% la5 = -2*la/mC_data(5); % la5 > 0 (match mode with empirical mean value for C5)
b4 = 1/la4; % b4 > 0
b5 = 1/la5; % b5 > 0

% Initial parameter vector
if useRedParam
    lambda = [la1 la2 la3 la]; % reduced parameterization
else
    lambda = [la1 la2 la3 la4 la5 la]; % full parameterization
end

mC45 = gamstat(a,[b4,b5]);
% [mC45,vC45] = gamstat(a,[b4,b5]);

% mC4 = gamstat(a,b4); % mC4 = a*b4;
% mC5 = gamstat(a,b5); % mC5 = a*b5;
% [mC4,vC4] = gamstat(a,b4); % vC4 = a*b4^2;
% [mC5,vC5] = gamstat(a,b5); % vC5 = a*b5^2;
% mC45 = [mC4 mC5]; % mC45 = a*[b4,b5];
% vC45 = [vC4 vC5]; % vC45 = a*[b4,b5].^2;

% stdC45 = sqrt(v45); % stdC45 = sqrt(a)*[b4,b5];
% cvC45 = stdC45./mC45; % cvC45 = 1/sqrt(a)*[1,1];

% For a Gamma-distributed random variable X~Gamma(a,b), phiX = log(X),
% nuX = E{phiX} = E{log(X)} = psi(a)+log(b)
% For Gamma-distributed random variables C4~Gamma(a,b4) and C5~Gamma(a,b5),
% phiC4 = log(C4) and phiC5 = log(C5),
% nuC4 = E{phiC4} = E{log(C4)} = psi(a)+log(b4)
% nuC5 = E{phiC5} = E{log(C5)} = psi(a)+log(b5)
nuC4 = psi(a)+log(b4);
nuC5 = psi(a)+log(b5);

rng(s); % initialize the random number generator using the default settings contained in s
N = 1e7; % number of samples
Nexponent = floor(log10(N));
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

mC123 = mean(C123_sample,1);
% vC123 = var(C123_sample,0,1); % vC123 = N/(N-1)*moment(C123_sample,2,1);
% stdC123 = std(C123_sample,0,1);
% cvC123 = stdC123./mC123;

mC = [mC123 mC45];
% vC = [vC123 vC45];
% stdC = [stdC123 stdC45];
% cvC = [cvC123 cvC45];
% stdC = sqrt(vC);
% cvC  = stdC./mC;
% sC = sqrt(norm(vC));
mCnorm = norm(mC);
% dC = sC/mCnorm;

mC_sample = mean(C_sample,1);
% vC_sample = var(C_sample,0,1); % vC_sample = N/(N-1)*moment(C_sample,2,1);
% stdC_sample = std(C_sample,0,1);
% cvC_sample = stdC_sample./mC_sample;
% sC_sample = sqrt(norm(vC_sample));
mCnorm_sample = norm(mC_sample);
% dC_sample = sC_sample/mCnorm_sample;

% phiC123 = log(det([C123]) = log(C1*C2-C3^2)
% phiC = log(det([C]) = log((C1*C2-C3^2) * C4^2 * C5^2)
%                     = log(C1*C2-C3^2) + 2*log(C4) + 2*log(C5)
%                     = phiC123         + 2*phiC4   + 2*phiC5
% nuC = E{phiC} = E{log(C1*C2-C3^2)} + 2*E{log(C4)} + 2*E{log(C5)}
%               = E{phiC123}         + 2*E{phiC4}   + 2*E{phiC5}
%               = nuC123             + 2*nuC4       + 2*nuC5
phiC123 = log(C1_sample.*C2_sample-C3_sample.^2);
nuC123 = mean(phiC123,1);
nuC = nuC123 + 2*(nuC4 + nuC5);

phiC_sample = log((C1_sample.*C2_sample-C3_sample.^2).*(C4_sample.^2).*(C5_sample.^2));
% phiC_sample = phiC123 + 2*log(C4_sample) + 2*log(C5_sample);
nuC_sample = mean(phiC_sample,1);

fC = [mC nuC];
fC_sample = [mC_sample nuC_sample];
fC_data = [mC_data nuC_data];

resnorm = norm(fC - fC_data)^2; % squared norm of residual
% resnorm = norm(fC - fC_data)^2/norm(fC_data)^2; % relative squared norm of residual
fval = norm(mC - mC_data)^2/norm(mC_data)^2 + abs((nuC - nuC_data)./nuC_data).^2; % objective function value (relative squared norm of residual)

resnorm_sample = norm(fC_sample - fC_data)^2; % squared norm of residual
% resnorm_sample = norm(fC_sample - fC_data)^2/norm(fC_data)^2; % relative squared norm of residual
fval_sample = norm(mC_sample - mC_data)^2/norm(mC_data)^2 + abs((nuC_sample - nuC_data)./nuC_data).^2; % objective function value (relative squared norm of residual)

errmCnorm = norm(mC - mC_data)/norm(mC_data); % relative error on norm(mC)
errmCnorm_sample = norm(mC_sample - mC_data)/norm(mC_data); % relative error on norm(mC)

err = abs(fC - fC_data)./abs(fC_data); % relative error on mC and nuC
err_sample = abs(fC_sample - fC_data)./abs(fC_data); % relative error on mC and nuC

fprintf('\n');
fprintf('lambda = %g\n',la);
fprintf('nb samples = %d\n',N);
fprintf('MCMC algo  = %s\n',MCMCalgo);
switch lower(MCMCalgo)
    case {'imh','rwmh'}
        fprintf('acceptance rate = %.2f%%\n',accept*100);
        fprintf('rejection rate  = %.2f%%\n',(1-accept)*100);
    case 'ss'
        fprintf('nb function evaluations per sample = %g\n',neval);
end
fprintf('resnorm        = %g\n',resnorm);
fprintf('resnorm_sample = %g\n',resnorm_sample);
fprintf('fval           = %g\n',fval);
fprintf('fval_sample    = %g\n',fval_sample);
fprintf('error norm(mean(C))        = %g\n',errmCnorm);
fprintf('error norm(mean(C_sample)) = %g\n',errmCnorm_sample);
fprintf('error          = (%g, %g, %g, %g, %g, %g)\n',err);
fprintf('error_sample   = (%g, %g, %g, %g, %g, %g)\n',err_sample);

residualnorms(k)        = resnorm;
residualnorms_sample(k) = resnorm_sample;
funvals(k)              = fval;
funvals_sample(k)       = fval_sample;
errors_mCnorm(k)        = errmCnorm;
errors_mCnorm_sample(k) = errmCnorm_sample;
errors(k,:)             = err;
errors_sample(k,:)      = err_sample;

Ns = (1:N)';
mC123s = cumsum(C123_sample,1) ./ Ns;
mC1s = mC123s(:,1);
mC2s = mC123s(:,2);
mC3s = mC123s(:,3);

mCs  = [mC123s repmat(mC45,N,1)];
mCnorms = sqrt(sum(mCs.^2,2));

mCs_sample = cumsum(C_sample,1) ./ Ns;
mC4s = mCs_sample(:,4);
mC5s = mCs_sample(:,5);

nuC123s = cumsum(phiC123,1) ./ Ns;
nuCs = nuC123s + 2*(nuC4 + nuC5);
nuCs_sample = cumsum(phiC_sample,1) ./ Ns;

fCs = [mCs nuCs];
resnorms = sum((fCs - fC_data).^2,2); % squared norm of residual
% resnorms = sum((fCs - fC_data).^2,2)/norm(fC_data)^2; % relative squared norm of residual
% resnorms = sum((fCs - fC_data).^2,2)/sum(fC_data.^2,2); % relative squared norm of residual
% fvals = sum((mCs - mC_data).^2,2)/norm(mC_data)^2 + abs((nuCs - nuC_data)./nuC_data).^2; % objective function value (relative squared norm of residual)
fvals = sum((mCs - mC_data).^2,2)/sum(mC_data.^2,2) + abs((nuCs - nuC_data)./nuC_data).^2; % objective function value (relative squared norm of residual)
% errmCnorms = sqrt(sum((mCs - mC_data).^2,2))/norm(mC_data); % relative error on norm(mC)
errmCnorms = sqrt(sum((mCs - mC_data).^2,2)/sum(mC_data.^2,2)); % relative error on norm(mC)
errs = abs(fCs - fC_data)./abs(fC_data); % relative error

fCs_sample = [mCs_sample nuCs_sample];
resnorms_sample = sum((fCs_sample - fC_data).^2,2); % squared norm of residual
% resnorms_sample = sum((fCs_sample - fC_data).^2,2)/norm(fC_data)^2; % relative squared norm of residual
% resnorms_sample = sum((fCs_sample - fC_data).^2,2)/sum(fC_data.^2,2); % relative squared norm of residual
% fvals_sample = sum((mCs_sample - mC_data).^2,2)/norm(mC_data)^2 + abs((nuCs_sample - nuC_data)./nuC_data).^2; % objective function value (relative squared norm of residual)
fvals_sample = sum((mCs_sample - mC_data).^2,2)/sum(mC_data.^2,2) + abs((nuCs_sample - nuC_data)./nuC_data).^2; % objective function value (relative squared norm of residual)
% errmCnorms_sample = sqrt(sum((mCs_sample - mC_data).^2,2))/norm(mC_data); % relative error on norm(mC)
errmCnorms_sample = sqrt(sum((mCs_sample - mC_data).^2,2)/sum(mC_data.^2,2)); % relative error on norm(mC)
errs_sample = abs(fCs_sample - fC_data)./abs(fC_data); % relative errors

figure('Name','Convergence mean(C1)')
clf
semilogx(1:N,mC1s,'-b','LineWidth',linewidth)
hold on
% semilogx(1:N,mC1s_sample,'-r','LineWidth',linewidth)
plot([1,N],[mC_data(1),mC_data(1)],'--','Color','k','LineWidth',linewidth)
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
legend('$\widehat{\underline{c}}_1(${\boldmath$\lambda$}$^{\mathrm{init}})$','$\underline{c}_1$','Interpreter',interpreter)
title(['$\lambda^{\mathrm{init}} = ' num2str(la) '$'],'Interpreter',interpreter)
% mysaveas(pathname,['convergence_mean_C1_lambda_init_' num2str(la)],formats);
% mymatlab2tikz(pathname,['convergence_mean_C1_lambda_init_' num2str(la) '.tex']);

figure('Name','Convergence mean(C2)')
clf
semilogx(1:N,mC2s,'-b','LineWidth',linewidth)
hold on
% semilogx(1:N,mC2s_sample,'-r','LineWidth',linewidth)
plot([1,N],[mC_data(2),mC_data(2)],'--','Color','k','LineWidth',linewidth)
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xticks(10.^(0:6))
xlim tight
xlabel('Number of samples','Interpreter',interpreter)
ylabel('Mean value of $C_2$ [GPa]','Interpreter',interpreter)
%xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
%ylabel('Valeur moyenne de $C_2$ [GPa]','Interpreter',interpreter)
legend('$\widehat{\underline{c}}_2(${\boldmath$\lambda$}$^{\mathrm{init}})$','$\underline{c}_2$','Interpreter',interpreter)
title(['$\lambda^{\mathrm{init}} = ' num2str(la) '$'],'Interpreter',interpreter)
% mysaveas(pathname,['convergence_mean_C2_lambda_init_' num2str(la)],formats);
% mymatlab2tikz(pathname,['convergence_mean_C2_lambda_init_' num2str(la) '.tex']);

figure('Name','Convergence mean(C3)')
clf
semilogx(1:N,mC3s,'-b','LineWidth',linewidth)
hold on
% semilogx(1:N,mC3s_sample,'-r','LineWidth',linewidth)
plot([1,N],[mC_data(3),mC_data(3)],'--','Color','k','LineWidth',linewidth)
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
legend('$\widehat{\underline{c}}_3(${\boldmath$\lambda$}$^{\mathrm{init}})$','$\underline{c}_3$','Interpreter',interpreter)
title(['$\lambda^{\mathrm{init}} = ' num2str(la) '$'],'Interpreter',interpreter)
% mysaveas(pathname,['convergence_mean_C3_lambda_init_' num2str(la)],formats);
% mymatlab2tikz(pathname,['convergence_mean_C3_lambda_init_' num2str(la) '.tex']);

figure('Name','Convergence norm(mean(C))');
clf
semilogx(1:N,mCnorms,'-b','LineWidth',linewidth)
hold on
% semilogx(1:N,mCnorms_sample,'-r','LineWidth',linewidth)
plot([1,N],[mCnorm_data,mCnorm_data],'--','Color','k','LineWidth',linewidth)
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
legend('$||\widehat{\underline{\mathbf{c}}}(${\boldmath$\lambda$}$^{\mathrm{init}})||$','$||\underline{\mathbf{c}}||$','Interpreter',interpreter)
title(['$\lambda^{\mathrm{init}} = ' num2str(la) '$'],'Interpreter',interpreter)
% mysaveas(pathname,['convergence_norm_mean_C_lambda_init_' num2str(la)],formats);
% mymatlab2tikz(pathname,['convergence_norm_mean_C_lambda_init_' num2str(la) '.tex']);

figure('Name','Convergence mean(phi(C))=mean(log(det([C])))')
clf
semilogx(1:N,nuCs,'-b','LineWidth',linewidth)
hold on
% semilogx(1:N,nuCs_sample,'-r','LineWidth',linewidth)
plot([1,N],[nuC_data,nuC_data],'--','Color','k','LineWidth',linewidth)
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
legend('$\widehat{\nu}_{\mathbf{C}}(${\boldmath$\lambda$}$^{\mathrm{init}})$','$\nu_{\mathbf{C}}$','Interpreter',interpreter)
title(['$\lambda^{\mathrm{init}} = ' num2str(la) '$'],'Interpreter',interpreter)
% mysaveas(pathname,['convergence_mean_phiC_lambda_init_' num2str(la)],formats);
% mymatlab2tikz(pathname,['convergence_mean_phiC_lambda_init_' num2str(la) '.tex']);

figure('Name','Convergence squared norm of residual')
clf
loglog(1:N,resnorms,'-b','LineWidth',linewidth)
% loglog(1:N,resnorms_sample,'-r','LineWidth',linewidth)
grid on
box on
set(gca,'FontSize',fontsize)
xticks(10.^(0:Nexponent))
xlim tight
xlabel('Number of samples','Interpreter',interpreter)
ylabel('Squared norm of residual','Interpreter',interpreter)
%xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
%ylabel('Norme au carr\''e du r\''esidu','Interpreter',interpreter)
title(['$\lambda^{\mathrm{init}} = ' num2str(la) '$'],'Interpreter',interpreter)
% mysaveas(pathname,['convergence_resnorm_lambda_init_' num2str(la)],formats);
% mymatlab2tikz(pathname,['convergence_resnorm_lambda_init_' num2str(la) '.tex']);

figure('Name','Convergence objective function value')
clf
loglog(1:N,fvals,'-b','LineWidth',linewidth)
% loglog(1:N,fvals_sample,'-r','LineWidth',linewidth)
grid on
box on
set(gca,'FontSize',fontsize)
xticks(10.^(0:Nexponent))
xlim tight
xlabel('Number of samples','Interpreter',interpreter)
ylabel('Objective function value $\mathcal{J}(${\boldmath$\lambda$}$^{\mathrm{init}})$','Interpreter',interpreter)
%xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
%ylabel('Valeur de la fonction objectif','Interpreter',interpreter)
title(['$\lambda^{\mathrm{init}} = ' num2str(la) '$'],'Interpreter',interpreter)
% mysaveas(pathname,['convergence_fval_lambda_init_' num2str(la)],formats);
% mymatlab2tikz(pathname,['convergence_fval_lambda_init_' num2str(la) '.tex']);

figure('Name','Convergence relative errors')
clf
loglog(1:N,errs(:,1),'-b','LineWidth',linewidth)
hold on
loglog(1:N,errs(:,2),'-r','LineWidth',linewidth)
loglog(1:N,errs(:,3),'-g','LineWidth',linewidth)
% loglog(1:N,errs(:,4),'-c','LineWidth',linewidth)
% loglog(1:N,errs(:,5),'-m','LineWidth',linewidth)
loglog(1:N,errmCnorms,'Color',[0.9290 0.6940 0.1250],'LineWidth',linewidth)
loglog(1:N,errs(:,6),'-k','LineWidth',linewidth)
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
legend('$|\widehat{\underline{c}}_1(${\boldmath$\lambda$}$^{\mathrm{init}})-\underline{c}_1|/|\underline{c}_1|$',...
    '$|\widehat{\underline{c}}_2(${\boldmath$\lambda$}$^{\mathrm{init}})-\underline{c}_2|/|\underline{c}_2|$',...
    '$|\widehat{\underline{c}}_3(${\boldmath$\lambda$}$^{\mathrm{init}})-\underline{c}_3|/|\underline{c}_3|$',...
    ...% '$|\widehat{\underline{c}}_4(${\boldmath$\lambda$}$^{\mathrm{init}})-\underline{c}_4|/|\underline{c}_4|$',...
    ...% '$|\widehat{\underline{c}}_5(${\boldmath$\lambda$}$^{\mathrm{init}})-\underline{c}_5|/|\underline{c}_5|$',...
    '$||\widehat{\underline{\mathbf{c}}}(${\boldmath$\lambda$}$^{\mathrm{init}})-\underline{\mathbf{c}}||/||\underline{\mathbf{c}}||$',...
    '$|\widehat{\nu}_{\mathbf{C}}(${\boldmath$\lambda$}$^{\mathrm{init}})-\nu_{\mathbf{C}}|/|\nu_{\mathbf{C}}|$',...
    'Interpreter',interpreter)
% legend('$C_1$','$C_2$','$C_3$','$C_4$','$C_5$','$\varphi(${\boldmath$C$}$)$','Interpreter',interpreter)
title(['$\lambda^{\mathrm{init}} = ' num2str(la) '$'],'Interpreter',interpreter)
% mysaveas(pathname,['convergence_relative_errors_lambda_init_' num2str(la)],formats);
% mymatlab2tikz(pathname,['convergence_relative_errors_lambda_init_' num2str(la) '.tex']);

%% Plot samples
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
scatter3(C1_sample,C2_sample,C3_sample,'b.')
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
title(['$\lambda^{\mathrm{init}}$ = ' num2str(la)],'Interpreter',interpreter)
% mysaveas(pathname,['samples_C1_C2_C3_lambda_init_' num2str(la)],formats);
% mymatlab2tikz(pathname,['samples_C1_C2_C3_lambda_init_' num2str(la) '.tex']);

end

if length(las)>1
    figure('Name','Evolution of squared norm of residual w.r.t. lambda')
    clf
    semilogy(las,residualnorms,'-b','LineWidth',linewidth)
    % semilogy(las,residualnorms_sample,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$\lambda^{\mathrm{init}}$','Interpreter',interpreter)
    ylabel('Squared norm of the residual','Interpreter',interpreter)
    %ylabel('Norme au carr\''e du r\''esidu','Interpreter',interpreter)
    % mysaveas(pathname,'evol_resnorm_lambda_init',formats);
    % mymatlab2tikz(pathname,'evol_resnorm_lambda_init.tex');

    figure('Name','Evolution of objective function value w.r.t. lambda')
    clf
    semilogy(las,funvals,'-b','LineWidth',linewidth)
    % semilogy(las,funvals_sample,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$\lambda^{\mathrm{init}}$','Interpreter',interpreter)
    ylabel('Objective function value','Interpreter',interpreter)
    %ylabel('Valeur de la fonction objectif','Interpreter',interpreter)
    % mysaveas(pathname,'evol_fval_lambda_init',formats);
    % mymatlab2tikz(pathname,'evol_fval_lambda_init.tex');
    
    figure('Name','Evolution of relative errors w.r.t. lambda')
    clf
    semilogy(las,errors(:,1),'-b','LineWidth',linewidth)
    hold on
    semilogy(las,errors(:,2),'-r','LineWidth',linewidth)
    semilogy(las,errors(:,3),'-g','LineWidth',linewidth)
    % semilogy(las,errors(:,4),'-c','LineWidth',linewidth)
    % semilogy(las,errors(:,5),'--m','LineWidth',linewidth)
    semilogy(las,errors_mCnorm,'-m','LineWidth',linewidth)
    semilogy(las,errors(:,6),'-k','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$\lambda^{\mathrm{init}}$','Interpreter',interpreter)
    ylabel('Relative error','Interpreter',interpreter)
    legend('$\widehat{\underline{c}}_1$','$\widehat{\underline{c}}_2$','$\widehat{\underline{c}}_3$','$||\widehat{\underline{\mathbf{c}}}||$','$\widehat{\nu}_{\mathbf{C}}$','Interpreter',interpreter)
    % legend('$\widehat{\underline{c}}_1$','$\widehat{\underline{c}}_2$','$\widehat{\underline{c}}_3$','$\widehat{\underline{c}}_4$','$\widehat{\underline{c}}_5$','$||\widehat{\underline{\mathbf{c}}}||$','$\widehat{\nu}_{\mathbf{C}}$','Interpreter',interpreter)
    %ylabel('Erreur relative','Interpreter',interpreter)
    % mysaveas(pathname,'evol_error_lambda_init',formats);
    % mymatlab2tikz(pathname,'evol_error_lambda_init.tex');
end
