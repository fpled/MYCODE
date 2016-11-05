%% Stochastic modeling of Young modulus %%
%%--------------------------------------%%

% clc
clear all
close all

filename = 'model_sto_Young';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'results',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc2'};
renderer = 'OpenGL';

fontsize = 16;
linewidth = 1;
markersize = 36;
interpreter = 'latex';

%% Input data
% Data on sample D
E = [4.211 4.057 3.685 3.921 3.839 3.845 3.795...
    3.406 3.389 3.299 3.485 3.319 3.267 3.349 3.307...
    4.684 4.245 4.076 4.407 4.283 4.054 4.226 4.041...
    4.104 4.075 3.556 3.319 3.848 3.707 3.664 3.493 3.550];
mean_E = mean(E);
std_E = std(E);

%% Plot samples
figure('Name','Experimental data')
clf
bar(1:length(E),E)
set(gca,'FontSize',fontsize)
set(gca,'XLim',[0,length(E)+1])
xlabel('\''Echantillon','Interpreter',interpreter);
ylabel('Module d''Young (GPa)','Interpreter',interpreter); 
mysaveas(pathname,'data_E','fig');
mymatlab2tikz(pathname,'data_E.tex');

%% Maximization of Log-likelihood function
fprintf('\nMaximization of Log-likelihood function');
fprintf('\n---------------------------------------\n');

opts = statset('TolFun',1e-6,'TolX',1e-6,'FunValCheck','off');

% Available information = moment of order 1 + repulsion in 0
fprintf('\nInformation = %s\n',strjoin({'moment1','repulsion0'},' + '));

% Gamma distribution
phat = gamfit(E);
% phat = mle(E,'distribution','gam');
% nloglf = @(phat,data,cens,freq) length(data)*log(gamma(phat(1)))...
%     +length(data)*phat(1)*log(phat(2))...
%     +(1-phat(1))*sum(log(data))...
%     +1/phat(2)*sum(data);
% phat = mle(E,'nloglf',nloglf,'start',[2 0],'lowerbound',[2 0],'options',opts);

fprintf('\nnb_samples = %d',length(E));
fprintf('\na = %.4f',phat(1));
fprintf('\nb = %.4f',phat(2));

mE = phat(1)*phat(2);
vE = phat(1)*phat(2)^2;
sE = sqrt(vE);
dE = sE/mE; % 1/sqrt(phat(1))
fprintf('\nmean(E) = %.4f',mE);
fprintf('\nvar(E)  = %.4f',vE);
fprintf('\nstd(E)  = %.4f',sE);
fprintf('\ndisp(E) = %.4f',dE);
fprintf('\n');

%% Plot pdf and cdf
pdf_gam = @(x) gampdf(x,phat(1),phat(2));
% pdf_gam = @(x) pdf('gam',x,phat(1),phat(2));
% pdf_gam = @(x) 1/(phat(2)^phat(1)*gamma(phat(1)))*(x.^(phat(1)-1)).*exp(-x./phat(2));
cdf_gam = @(x) gamcdf(x,phat(1),phat(2));

xmin = max(0,mE-5*sE);
xmax = mE+5*sE;
x = linspace(xmin,xmax,1e3);

figure('Name','Probability density function')
clf
hold on
plot(x,pdf_gam(x),'-b','LineWidth',linewidth);
plot(E,pdf_gam(E),'r+');
hold off
grid on
box on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[xmin,xmax])
xlabel('$e$ (GPa)','Interpreter',interpreter)
ylabel('$p_E(e)$','Interpreter',interpreter)
% l = legend('$p_E(e)$','$(e_i,p_E(e_i))_{i=1}^n$');
% set(l,'Interpreter',interpreter,'Location','northwest');
mysaveas(pathname,'pdf_E','fig');
mymatlab2tikz(pathname,'pdf_E.tex',...
    'extraAxisOptions',{'ylabel style={overlay}'});

figure('Name','Cumulative distribution function')
clf
hold on
plot(x,cdf_gam(x),'-b','LineWidth',linewidth);
plot(E,cdf_gam(E),'r+');
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
mysaveas(pathname,'cdf_E','fig');
mymatlab2tikz(pathname,'cdf_E.tex',...
    'extraAxisOptions',{'ylabel style={overlay}'});

%% Sample generation
N = 1e4; % number of samples
e = gamrnd(phat(1),phat(2),1,N);
% u = randn(1,N);
% e = gaminv(normcdf(u),phat(1),phat(2));

figure('Name','Monte Carlo simulation')
clf
hold on
% plot(1:N,e,'b+')
scatter(1:N,e,markersize,'b+','LineWidth',linewidth);
plot([1 N],[mE mE],'-r','LineWidth',linewidth)
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
ylabel('Module d''Young (GPa)','Interpreter',interpreter)
% l = legend('r\''ealisations $(e_1,\dots,e_N)$','moyenne $m_E$');
% set(l,'Interpreter',interpreter);
mysaveas(pathname,'samples_E','fig');
mymatlab2tikz(pathname,'samples_E.tex');
