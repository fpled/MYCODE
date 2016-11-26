%% Probabilistic modeling of cortical bone microstructure %%
%%--------------------------------------------------------%%

% clc
clear all
close all

pathname = '../Presentation/figures/';
fontsize = 16;
linewidth = 1;
markersize = 36;
interpreter = 'latex';

volume_fraction = 'pore_cort';
% volume_fraction = 'col_ultra';
% volume_fraction = 'HA_foam';

%% Input data
switch volume_fraction
    case 'pore_cort'
        f = [6.47 6.49 4.15 5.26 6.21 8.46 8.80 5.91 20.7 15.5 5.89 7.39]'/100;
        mean_f = 8.44/100;
        disp_f = 57.1/100;
    case 'col_ultra'
        f = [35.4 35.7 35.9 35.2 35.4 35.7 35.2 35.5 35.7 35.4 35.8 35.7]'/100;
        mean_f = 35.5/100;
        disp_f = 1.54/100;
    case 'HA_foam'
        f = [50.4 52.8 54.4 49.2 50.8 53.0 49.3 51.6 53.0 50.5 53.4 53.0]'/100;
        mean_f = 51.5/100;
        disp_f = 7.94/100;
    otherwise
        error('Wrong volume fraction')
end
cens = [0 0 0 0 0 0 0 0 1 1 0 0]';
std_f = mean_f*disp_f;

%% Maximization of Log-likelihood function
fprintf('\nMaximization of Log-likelihood function');
fprintf('\n---------------------------------------\n');

opts = statset('TolFun',1e-6,'TolX',1e-6,'FunValCheck','off');

%% Available information = moment of order 1 + moment of order 2
fprintf('\nInformation = %s\n',strjoin({'moment1','moment2'},' + '));

% All data
% Custom distribution
nloglf = @(lambda,data,cens,freq) length(data)*log(integral(@(x) exp(-lambda(1)*x-lambda(2)*x.^2),0,1))...
    +lambda(1)*sum(data)...
    +lambda(2)*sum(data.^2);
lambda = mle(f,'nloglf',nloglf,'start',[0 0],'options',opts);
lambda_mle_m12 = [log(integral(@(x) exp(-lambda(1)*x-lambda(2)*x.^2),0,1)),lambda(1:end)];
pdf_mle_m12 = @(x) exp(-lambda_mle_m12(1)-lambda_mle_m12(2)*x-lambda_mle_m12(3)*x.^2);

fprintf('\nnb samples = %d',length(f));
for i=1:length(lambda_mle_m12)
    fprintf(['\nlambda_' num2str(i-1) ' = %.3f'],lambda_mle_m12(i));
end
fprintf('\n');

% All data except censored ones
% Custom distribution
nloglf = @(lambda,data,cens,freq) length(data(~cens))*log(integral(@(x) exp(-lambda(1)*x-lambda(2)*x.^2),0,1))...
    +lambda(1)*sum(data(~cens))...
    +lambda(2)*sum(data(~cens).^2);
lambda_cens = mle(f,'nloglf',nloglf,'start',[0 0],'censoring',cens,'options',opts);
lambda_mle_m12_cens = [log(integral(@(x) exp(-lambda_cens(1)*x-lambda_cens(2)*x.^2),0,1)),lambda_cens(1:end)];
pdf_mle_m12_cens = @(x) exp(-lambda_mle_m12_cens(1)-lambda_mle_m12_cens(2)*x-lambda_mle_m12_cens(3)*x.^2);

fprintf('\nnb samples = %d',length(f(~cens)));
for i=1:length(lambda_mle_m12_cens)
    fprintf(['\nlambda_' num2str(i-1) ' = %.3f'],lambda_mle_m12_cens(i));
end
fprintf('\n');

%% Available information = repulsion in 0 + repulsion in 1
fprintf('\nInformation = %s\n',strjoin({'repulsion0','repulsion1'},' + '));

% All data
% Beta distribution
p = mle(f,'distribution','beta');
pdf_mle_r01 = @(x) pdf('beta',x,p(1),p(2));
lambda_mle_r01 = [betaln(p(1),p(2)),1-p(1:end)];

% Custom distribution
% nloglf = @(lambda,data,cens,freq) length(data)*betaln(1-lambda(1),1-lambda(2))...
%     +lambda(1)*sum(log(data))...
%     +lambda(2)*sum(log1p(-data));
% lambda = mle(f,'nloglf',nloglf,'start',[0 0]);
% lambda_mle_r01 = [betaln(1-lambda(1),1-lambda(2)),lambda(1:end)];
% pdf_mle_r01 = @(x) exp(-lambda_mle_r01(1)-lambda_mle_r01(2)*log(x)-lambda_mle_r01(3)*log1p(-x));

fprintf('\nnb samples = %d',length(f));
for i=1:length(lambda_mle_r01)
    fprintf(['\nlambda_' num2str(i-1) ' = %.3f'],lambda_mle_r01(i));
end
fprintf('\n');

% All data except censored ones
% Custom distribution
nloglf = @(lambda,data,cens,freq) length(data(~cens))*betaln(1-lambda(1),1-lambda(2))...
    +lambda(1)*sum(log(data(~cens)))...
    +lambda(2)*sum(log1p(-data(~cens)));
lambda_cens = mle(f,'nloglf',nloglf,'start',[0 0],'censoring',cens,'options',opts);
lambda_mle_r01_cens = [betaln(1-lambda_cens(1),1-lambda_cens(2)),lambda_cens(1:end)];
pdf_mle_r01_cens = @(x) exp(-lambda_mle_r01_cens(1)-lambda_mle_r01_cens(2)*log(x)-lambda_mle_r01_cens(3)*log1p(-x));

fprintf('\nnb samples = %d',length(f(~cens)));
for i=1:length(lambda_mle_r01_cens)
    fprintf(['\nlambda_' num2str(i-1) ' = %.3f'],lambda_mle_r01_cens(i));
end
fprintf('\n');

%% Available information = moment of order 1 + repulsion in 0 + repulsion in 1
fprintf('\nInformation = %s\n',strjoin({'moment1','repulsion0','repulsion1'},' + '));

% All data
% Custom distribution
nloglf = @(lambda,data,cens,freq) length(data)*log(integral(@(x) exp(-lambda(1)*x-lambda(2)*log(x)-lambda(3)*log1p(-x)),0,1))...
    +lambda(1)*sum(data)...
    +lambda(2)*sum(log(data))...
    +lambda(3)*sum(log1p(-data));
lambda = mle(f,'nloglf',nloglf,'start',[0 0 0],'options',opts);
lambda_mle_m1r01 = [log(integral(@(x) exp(-lambda(1)*x-lambda(2)*log(x)-lambda(3)*log1p(-x)),0,1)),lambda(1:end)];
pdf_mle_m1r01 = @(x) exp(-lambda_mle_m1r01(1)-lambda_mle_m1r01(2)*x-lambda_mle_m1r01(3)*log(x)-lambda_mle_m1r01(4)*log1p(-x));

fprintf('\nnb samples = %d',length(f));
for i=1:length(lambda_mle_m1r01)
    fprintf(['\nlambda_' num2str(i-1) ' = %.3f'],lambda_mle_m1r01(i));
end
fprintf('\n');

% All data except censored ones
% Custom distribution
nloglf = @(lambda,data,cens,freq) length(data(~cens))*log(integral(@(x) exp(-lambda(1)*x-lambda(2)*log(x)-lambda(3)*log1p(-x)),0,1))...
    +lambda(1)*sum(data(~cens))...
    +lambda(2)*sum(log(data(~cens)))...
    +lambda(3)*sum(log1p(-data(~cens)));
lambda_cens = mle(f,'nloglf',nloglf,'start',[0 0 0],'censoring',cens,'options',opts);
lambda_mle_m1r01_cens = [log(integral(@(x) exp(-lambda_cens(1)*x-lambda_cens(2)*log(x)-lambda_cens(3)*log1p(-x)),0,1)),lambda_cens(1:end)];
pdf_mle_m1r01_cens = @(x) exp(-lambda_mle_m1r01_cens(1)-lambda_mle_m1r01_cens(2)*x-lambda_mle_m1r01_cens(3)*log(x)-lambda_mle_m1r01_cens(4)*log1p(-x));

fprintf('\nnb samples = %d',length(f(~cens)));
for i=1:length(lambda_mle_m1r01_cens)
    fprintf(['\nlambda_' num2str(i-1) ' = %.3f'],lambda_mle_m1r01(i));
end
fprintf('\n');

%% Minimization of strictly convex function H
fprintf('\nMinimization of strictly convex function');
fprintf('\n----------------------------------------\n');

H = @(lambda) integral(@(x) exp(-lambda(1)-lambda(2)*x-lambda(3)*x.^2),0,1)...
    +lambda(1)+lambda(2)*mean_f+lambda(3)*(std_f^2+mean_f^2);
lambda_mee = fminsearch(@(lambda) H(lambda),[0,0,0]);
pdf_mee = @(x) exp(-lambda_mee(1)-lambda_mee(2)*x-lambda_mee(3)*x.^2);

for i=1:length(lambda_mee)
    fprintf(['\nlambda_' num2str(i-1) ' = %.3f'],lambda_mee(i));
end
fprintf('\n');

%% Plot pdfs
figure('Name','Probability density function')
clf
x = linspace(0,1,1000);
hold on
plot(x,pdf_mle_m12(x),'-b','LineWidth',linewidth);
plot(x,pdf_mle_m12_cens(x),'--b','LineWidth',linewidth);
plot(x,pdf_mle_r01(x),'-g','LineWidth',linewidth);
plot(x,pdf_mle_r01_cens(x),'--g','LineWidth',linewidth);
plot(x,pdf_mle_m1r01(x),'-m','LineWidth',linewidth);
plot(x,pdf_mle_m1r01_cens(x),'--m','LineWidth',linewidth);
plot(x,pdf_mee(x),'-r','LineWidth',linewidth);
plot(f,pdf_mle_m12(f),'b+');
plot(f(~cens),pdf_mle_m12_cens(f(~cens)),'b+');
plot(f,pdf_mle_r01(f),'g+');
plot(f(~cens),pdf_mle_r01_cens(f(~cens)),'g+');
plot(f,pdf_mle_m1r01(f),'m+');
plot(f(~cens),pdf_mle_m1r01_cens(f(~cens)),'m+');
hold off
grid on
box on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[0,1])
xlabel('$x$','Interpreter',interpreter)
ylabel('$p_X(x)$','Interpreter',interpreter)
l = legend(['$c_0e^{-\lambda_1 x-\lambda_2 x^2}$, $\max_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mle_m12)-1) '} \ell(\boldsymbol{\lambda})$'],...
    ['$c_0e^{-\lambda_1 x-\lambda_2 x^2}$, $\max_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mle_m12_cens)-1) '} \ell(\boldsymbol{\lambda})$'],...
    ['$c_0e^{-\lambda_1 \ln(x)-\lambda_2 \ln(1-x)}$, $\max_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mle_r01)-1) '} \ell(\boldsymbol{\lambda})$'],...
    ['$c_0e^{-\lambda_1 \ln(x)-\lambda_2 \ln(1-x)}$, $\max_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mle_r01_cens)-1) '} \ell(\boldsymbol{\lambda})$'],...
    ['$c_0e^{-\lambda_1 x-\lambda_2 \ln(x)-\lambda_3 \ln(1-x)}$, $\max_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mle_m1r01)-1) '} \ell(\boldsymbol{\lambda})$'],...
    ['$c_0e^{-\lambda_1 x-\lambda_2 \ln(x)-\lambda_3 \ln(1-x)}$, $\max_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mle_m1r01_cens)-1) '} \ell(\boldsymbol{\lambda})$'],...
    ['$e^{-\lambda_0-\lambda_1 x-\lambda_2 x^2}$, $\min_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mee)) '} \Hc(\boldsymbol{\lambda})$']);
set(l,'Interpreter','none');
mysaveas(pathname,['pdf_' volume_fraction],'fig');
mymatlab2tikz(pathname,['pdf_' volume_fraction '.tex']);

figure('Name','Probability density function')
clf
switch volume_fraction
    case 'pore_cort'
        x_zoom = linspace(0,0.3,1000);
    case 'col_ultra'
        x_zoom = linspace(0.3,0.4,1000);
    case 'HA_foam'
        x_zoom = linspace(0.35,0.65,1000);
    otherwise
        error('Wrong volume fraction')
end
hold on
plot(x_zoom,pdf_mle_m12(x_zoom),'-b','LineWidth',linewidth);
plot(x_zoom,pdf_mle_m12_cens(x_zoom),'--b','LineWidth',linewidth);
plot(x_zoom,pdf_mle_r01(x_zoom),'-g','LineWidth',linewidth);
plot(x_zoom,pdf_mle_r01_cens(x_zoom),'--g','LineWidth',linewidth);
plot(x_zoom,pdf_mle_m1r01(x_zoom),'-m','LineWidth',linewidth);
plot(x_zoom,pdf_mle_m1r01_cens(x_zoom),'--m','LineWidth',linewidth);
plot(x_zoom,pdf_mee(x_zoom),'-r','LineWidth',linewidth);
plot(f,pdf_mle_m12(f),'b+');
plot(f(~cens),pdf_mle_m12_cens(f(~cens)),'b+');
plot(f,pdf_mle_r01(f),'g+');
plot(f(~cens),pdf_mle_r01_cens(f(~cens)),'g+');
plot(f,pdf_mle_m1r01(f),'m+');
plot(f(~cens),pdf_mle_m1r01_cens(f(~cens)),'m+');
hold off
grid on
box on
set(gca,'FontSize',fontsize)
switch volume_fraction
    case 'pore_cort'
        set(gca,'XLim',[0,0.25])
    case 'col_ultra'
        set(gca,'XLim',[0.32,0.4])
    case 'HA_foam'
        set(gca,'XLim',[0.35,0.65])
    otherwise
        error('Wrong volume fraction')
end
xlabel('$x$','Interpreter',interpreter)
ylabel('$p_X(x)$','Interpreter',interpreter)
l = legend(['$c_0e^{-\lambda_1 x-\lambda_2 x^2}$, $\max_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mle_m12)-1) '} \ell(\boldsymbol{\lambda})$'],...
    ['$c_0e^{-\lambda_1 x-\lambda_2 x^2}$, $\max_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mle_m12_cens)-1) '} \ell(\boldsymbol{\lambda})$'],...
    ['$c_0e^{-\lambda_1 \ln(x)-\lambda_2 \ln(1-x)}$, $\max_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mle_r01)-1) '} \ell(\boldsymbol{\lambda})$'],...
    ['$c_0e^{-\lambda_1 \ln(x)-\lambda_2 \ln(1-x)}$, $\max_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mle_r01_cens)-1) '} \ell(\boldsymbol{\lambda})$'],...
    ['$c_0e^{-\lambda_1 x-\lambda_2 \ln(x)-\lambda_3 \ln(1-x)}$, $\max_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mle_m1r01)-1) '} \ell(\boldsymbol{\lambda})$'],...
    ['$c_0e^{-\lambda_1 x-\lambda_2 \ln(x)-\lambda_3 \ln(1-x)}$, $\max_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mle_m1r01_cens)-1) '} \ell(\boldsymbol{\lambda})$'],...
    ['$e^{-\lambda_0-\lambda_1 x-\lambda_2 x^2}$, $\min_{\boldsymbol{\lambda} \in \Rbb^' num2str(length(lambda_mee)) '} \Hc(\boldsymbol{\lambda})$']);
set(l,'Interpreter','none');
mysaveas(pathname,['pdf_zoom_' volume_fraction],'fig');
mymatlab2tikz(pathname,['pdf_zoom_' volume_fraction '.tex']);

%% Plot log-likelihood function logL
% figure('Name','Maximization of Log-likelihood')
% clf
% N = 100;
% delta = 0.5*abs(lambda_mle_m12(:));
% la1 = linspace(lambda_mle_m12(1)-delta(1),lambda_mle_m12(1)+delta(1),N);
% la2 = linspace(lambda_mle_m12(2)-delta(2),lambda_mle_m12(2)+delta(2),N);
% logL = zeros(N); % Preallocate memory
% for i = 1:N
%     for j = 1:N
%         logL(i,j) = log_likelihood(f,[la1(i),la2(j)],{'moment1','moment2'});
%     end
% end
% [La1,La2] = meshgrid(la1,la2);
% surfc(La1,La2,logL);
% hold on
% plot3(lambda_mle_m12(1),lambda_mle_m12(2),log_likelihood(f,lambda_mle_m12,{'moment1','moment2'}),...
%     'ro','MarkerSize',5,...
%     'MarkerFaceColor','r');
% set(gca,'FontSize',fontsize)
% xlabel('$\lambda_1$','Interpreter',interpreter)
% ylabel('$\lambda_2$','Interpreter',interpreter)
% zlabel('$\ell(\boldsymbol{\lambda})$','Interpreter','none')
% mysaveas(PathName,['max_logL_' volume_fraction],'fig');
% mymatlab2tikz(PathName,['max_logL_' volume_fraction '.tex']);

%% Plot strictly convex function H
% figure('Name','Minimization of strictly convex function')
% clf
% delta = 0.5*abs(lambda_mee(:));
% la1 = linspace(lambda_mee(1)-delta(1),lambda_mee(1)+delta(1),N);
% la2 = linspace(lambda_mee(2)-delta(2),lambda_mee(2)+delta(2),N);
% la3 = linspace(lambda_mee(3)-delta(3),lambda_mee(3)+delta(3),N);
% H = zeros(N); % Preallocate memory
% for i = 1:N
%     for j = 1:N
%         H(i,j) = convex_fun([lambda_mee(1),la2(i),la3(j)],'mean',mean_f,'std',std_f);
%     end
% end
% [La2,La3] = meshgrid(la2,la3);
% surfc(La2,La3,H);
% hold on
% plot3(lambda_mee(2),lambda_mee(3),convex_fun(lambda_mee,'mean',mean_f,'std',std_f),...
%     'ro','MarkerSize',5,...
%     'MarkerFaceColor','r');
% set(gca,'FontSize',fontsize)
% xlabel('$\lambda_1$','Interpreter',interpreter)
% ylabel('$\lambda_2$','Interpreter',interpreter)
% zlabel('$H(\boldsymbol{\lambda})$','Interpreter','none')
% mysaveas(PathName,['min_H_' volume_fraction],'fig');
% mymatlab2tikz(PathName,['min_H_' volume_fraction '.tex']);

%% Statistics
fprintf('\nStatistics');
fprintf('\n----------\n');

norm_mle_m12 = integral(@(x) pdf_mle_m12(x),0,1);
norm_mle_r01 = integral(@(x) pdf_mle_r01(x),0,1);
norm_mle_m1r01 = integral(@(x) pdf_mle_m1r01(x),0,1);
norm_mee = integral(@(x) pdf_mee(x),0,1);
mean_mle_m12 = integral(@(x) x.*pdf_mle_m12(x),0,1);
mean_mle_r01 = integral(@(x) x.*pdf_mle_r01(x),0,1);
mean_mle_m1r01 = integral(@(x) x.*pdf_mle_m1r01(x),0,1);
mean_mee = integral(@(x) x.*pdf_mee(x),0,1);
var_mle_m12 = integral(@(x) (x-mean_mle_m12).^2.*pdf_mle_m12(x),0,1);
var_mle_r01 = integral(@(x) (x-mean_mle_r01).^2.*pdf_mle_r01(x),0,1);
var_mle_m1r01 = integral(@(x) (x-mean_mle_m1r01).^2.*pdf_mle_m1r01(x),0,1);
var_mee = integral(@(x) (x-mean_mee).^2.*pdf_mee(x),0,1);
std_mle_m12 = sqrt(var_mle_m12);
std_mle_r01 = sqrt(var_mle_r01);
std_mle_m1r01 = sqrt(var_mle_m1r01);
std_mee = sqrt(var_mee);
disp_mle_m12 = std_mle_m12/mean_mle_m12;
disp_mle_r01 = std_mle_r01/mean_mle_r01;
disp_mle_m1r01 = std_mle_m1r01/mean_mle_m1r01;
disp_mee = std_mee/mean_mee;

formatVal1 = '%12.6f';
formatVal2 = '%15.6f';
disp('                     +--------------+--------------+-----------------+--------------+')
disp('                     |     Max L    |     Max L    |      Max L      |     Max H    |')
disp('                     | info = m1+m2 | info = r0+r1 | info = m1+r0+r1 |              |')
disp('+--------------------+--------------+--------------+-----------------+--------------+')
fprintf(['| Mean               | ' formatVal1 ' | ' formatVal1 ' | ' formatVal2 ' | ' formatVal1 ' |\n'],mean_mle_m12,mean_mle_r01,mean_mle_m1r01,mean_mee)
fprintf(['| Variance           | ' formatVal1 ' | ' formatVal1 ' | ' formatVal2 ' | ' formatVal1 ' |\n'],var_mle_m12,var_mle_r01,var_mle_m1r01,var_mee)
fprintf(['| Standard deviation | ' formatVal1 ' | ' formatVal1 ' | ' formatVal2 ' | ' formatVal1 ' |\n'],std_mle_m12,std_mle_r01,std_mle_m1r01,std_mee)
fprintf(['| Dispersion         | ' formatVal1 ' | ' formatVal1 ' | ' formatVal2 ' | ' formatVal1 ' |\n'],disp_mle_m12,disp_mle_r01,disp_mle_m1r01,disp_mee)
fprintf(['| Normalization      | ' formatVal1 ' | ' formatVal1 ' | ' formatVal2 ' | ' formatVal1 ' |\n'],norm_mle_m12,norm_mle_r01,norm_mle_m1r01,norm_mee)
disp('+--------------------+--------------+--------------+-----------------+--------------+')
