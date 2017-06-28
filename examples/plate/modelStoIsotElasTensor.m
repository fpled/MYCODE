%% Stochastic modeling of isotropic elasticity tensor %%
%%----------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');

%% Input data
filename = 'modelStoIsotElasTensor';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','plate',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
fontsize = 16;
linewidth = 1;
markersize = 36;
interpreter = 'latex';
formats = {'fig','epsc2'};
renderer = 'OpenGL';

Num_sample = 20;
load('resultats_ET_GL_EL_nuL_all_samples')

E_data = zeros(1,Num_sample);  
G_data = zeros(1,Num_sample);  
for i = 1:Num_sample
    num_ech = ['C' num2str(i)];
   eval( ['E_data(i) = mean_ET_' num_ech '_rest/2.1;'] ); % Gpa
   eval( ['G_data(i) = mean_GL_' num_ech '_rest*7/1000;'] ); % Gpa 
end
nu_data = E_data./(2*G_data)-1;
lambda_data = E_data.*nu_data./((1+nu_data).*(1-2*nu_data)); 
C1_data = lambda_data + 2/3*G_data; 
C2_data = G_data; 

%% Plot samples
figure('Name','Experimental data E')
clf
bar(1:length(E_data),E_data)
set(gca,'FontSize',fontsize)
set(gca,'XLim',[0,length(E_data)+1])
% xlabel('\''Echantillon','Interpreter',interpreter);
% ylabel('Module d''Young $E$ (GPa)','Interpreter',interpreter); 
xlabel('Sample','Interpreter',interpreter);
ylabel('Young modulus $E$ (GPa)','Interpreter',interpreter);
mysaveas(pathname,'data_E','fig');
mymatlab2tikz(pathname,'data_E.tex');

figure('Name','Experimental data G')
clf
bar(1:length(G_data),G_data)
set(gca,'FontSize',fontsize)
set(gca,'XLim',[0,length(G_data)+1])
% xlabel('\''Echantillon','Interpreter',interpreter);
% ylabel('Module de cisaillement $G$ (GPa)','Interpreter',interpreter); 
xlabel('Sample','Interpreter',interpreter);
ylabel('Shear modulus $G$ (GPa)','Interpreter',interpreter);
mysaveas(pathname,'data_G','fig');
mymatlab2tikz(pathname,'data_G.tex');

%% Maximization of Log-likelihood function
n_data = length(C1_data);
m_data = length(C2_data);
data = [C1_data C2_data];

nloglf = @(lambda,data,cens,freq) -n_data*( (1-lambda(3))*log(lambda(1)) - gammaln(1-lambda(3)) )...
    -m_data*( (1-5*lambda(3))*log(lambda(2)) - gammaln(1-5*lambda(3)) )...
    +lambda(3)*( sum(log(data(1:n_data))) + 5*sum(log(data(n_data+1:end))) )...
    +lambda(1)*sum(data(1:n_data)) + lambda(2)*sum(data(n_data+1:end));
lambda = mle(data,'nloglf',nloglf,'start',[1 1 0],'lowerbound',[0 0 -Inf],'upperbound',[Inf Inf 1/5]);
a1 = 1-lambda(3);
b1 = 1/lambda(1);
a2 = 1-5*lambda(3);
b2 = 1/lambda(2);
% k1ln = (1-lambda(3))*log(lambda(1))-gammaln(1-lambda(3));
% k2ln = (1-5*lambda(3))*log(lambda(2))-gammaln(1-5*lambda(3));
k1ln = -gammaln(a1)-a1*log(b1);
k2ln = -gammaln(a2)-a2*log(b2);

pdf_C1 = @(c1) gampdf(c1,a1,b1); % Gamma probability density function of C1
pdf_C2 = @(c2) gampdf(c2,a2,b2); % Gamma probability density function of C2
% pdf_C = @(c1,c2) exp(k1ln+k2ln)*c1.^(-lambda(3))*c2.^(-5*lambda(3))*exp(-lambda(1)*c1-lambda(2)*c2);
pdf_C = @(c1,c2) pdf_C1(c1).*pdf_C2(c2); % joint probability density function of C=(C1,C2)
pdf_EN = @(e,n) exp(k1ln+k2ln)*(e./(3*(1-2*n))).^(-lambda(3)).*(e./(2*(1+n))).^(-5*lambda(3))...
    .*(e./(2*((1+n).^2).*((1-2*n).^2)))...
    .*exp(-lambda(1)*e./(3*(1-2*n))-lambda(2)*e./(2*(1+n))); % joint probability density function of (E,N)

cdf_C1 = @(c1) gamcdf(c1,a1,b1); % Gamma cumulative density function of C1
cdf_C2 = @(c2) gamcdf(c2,a2,b2); % Gamma cumulative density function of C2
% cdf_C1 = @(c1) integral(@(x) pdf_C1(x),-Inf,c1);
% cdf_C2 = @(c2) integral(@(x) pdf_C2(x),-Inf,c2);
cdf_C = @(c1,c2) cdf_C1(c1).*cdf_C2(c2); % joint cumulative density function of C=(C1,C2)
% cdf_C = @(c1,c2) integral2(@(x1,x2) pdf_C(x1,x2),-Inf,c1,-Inf,c2);
cdf_EN = @(e,n) integral2(@(xe,xn) pdf_EN(xe,xn),0,e,-1,n);

%% Plot pdfs
N = 100;

x1_min = max(0,mean(C1_data)-4*std(C1_data));
x1_max = mean(C1_data)+4*std(C1_data);
x2_min = max(0,mean(C2_data)-8*std(C2_data));
x2_max = mean(C2_data)+8*std(C2_data);

xe_min = max(0,mean(E_data)-10*std(E_data));
xe_max = mean(E_data)+10*std(E_data);
xn_min = max(-1,mean(nu_data)-5*std(nu_data));
xn_max = min(0.5,mean(nu_data)+5*std(nu_data));

% Plot pdf of C1
figure('Name','Probability density function of C1')
clf
x1 = linspace(x1_min,x1_max,N);
y1 = pdf_C1(x1);
plot(x1,y1,'-b');
hold on
plot(C1_data,pdf_C1(C1_data),'r+');
hold off
grid on
box on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[min(x1),max(x1)])
xlabel('$c_1$ (GPa)','Interpreter',interpreter)
ylabel('$p_{C_1}(c_1)$','Interpreter',interpreter)
mysaveas(pathname,'pdf_C1',formats);
mymatlab2tikz(pathname,'pdf_C1.tex');

% Plot pdf of C2
figure('Name','Probability density function of C2')
clf
x2 = linspace(x2_min,x2_max,N);
y2 = pdf_C2(x2);
plot(x2,y2,'-b');
hold on
plot(C2_data,pdf_C2(C2_data),'r+');
hold off
grid on
box on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[min(x2),max(x2)])
xlabel('$c_2$ (GPa)','Interpreter',interpreter)
ylabel('$p_{C_2}(c_2)$','Interpreter',interpreter)
mysaveas(pathname,'pdf_C2',formats);
mymatlab2tikz(pathname,'pdf_C2.tex');

% Plot pdf of C
figure('Name','Probability density function of C=(C1,C2)')
clf
[X1,X2] = meshgrid(x1,x2);
Z = pdf_C(X1,X2);
surfc(X1,X2,Z);
hold on
plot3(C1_data,C2_data,pdf_C(C1_data,C2_data),'r+');
hold off
grid on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[min(x1),max(x1)])
set(gca,'YLim',[min(x2),max(x2)])
xlabel('$c_1$ (GPa)','Interpreter',interpreter)
ylabel('$c_2$ (GPa)','Interpreter',interpreter)
zlabel('$p_{(C_1,C_2)}(c_1,c_2)$','Interpreter',interpreter)
mysaveas(pathname,'pdf_C',formats);
mymatlab2tikz(pathname,'pdf_C.tex');

% Plot pdf of EN
figure('Name','Probability density function of (E,N)')
clf
xe = linspace(xe_min,xe_max,N);
xn = linspace(xn_min,xn_max,N);
[Xe,Xn] = meshgrid(xe,xn);
Z = pdf_EN(Xe,Xn);
surfc(Xe,Xn,Z);
hold on
plot3(E_data,nu_data,pdf_EN(E_data,nu_data),'r+');
hold off
grid on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[min(xe),max(xe)])
set(gca,'YLim',[min(xn),max(xn)])
xlabel('$e$ (GPa)','Interpreter',interpreter)
ylabel('$n$','Interpreter',interpreter)
zlabel('$p_{(E,N)}(e,n)$','Interpreter',interpreter)
mysaveas(pathname,'pdf_EN',formats);
mymatlab2tikz(pathname,'pdf_EN.tex');

%% Plot cdfs
% Plot cdf of C1
figure('name','cumulative distribution function of C1')
clf
z1 = cdf_C1(x1);
plot(x1,z1,'-b');
hold on
plot(C1_data,gamcdf(C1_data,a1,b1),'r+');
hold off
grid on
box on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[min(x1),max(x1)])
xlabel('$c_1$ (GPa)','Interpreter',interpreter)
ylabel('$F_{C_1}(c_1)$','Interpreter',interpreter)
mysaveas(pathname,'cdf_C1',formats);
mymatlab2tikz(pathname,'cdf_C1.tex');

% Plot cdf of C2
figure('name','cumulative distribution function of C2')
clf
z2 = cdf_C2(x2);
plot(x2,z2,'-b');
hold on
plot(C2_data,gamcdf(C2_data,a2,b2),'r+');
hold off
grid on
box on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[min(x2),max(x2)])
xlabel('$c_2$ (GPa)','Interpreter',interpreter)
ylabel('$F_{C_2}(c_2)$','Interpreter',interpreter)
mysaveas(pathname,'cdf_C2',formats);
mymatlab2tikz(pathname,'cdf_C2.tex');

% Plot cdf of C
figure('Name','Cumulative density function of C=(C1,C2)')
clf
Z = cdf_C(X1,X2);
surf(X1,X2,Z);
hold on
plot3(C1_data,C2_data,cdf_C(C1_data,C2_data),'r+');
hold off
grid on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[min(x1),max(x1)])
set(gca,'YLim',[min(x2),max(x2)])
xlabel('$c_1$ (GPa)','Interpreter',interpreter)
ylabel('$c_2$ (GPa)','Interpreter',interpreter)
zlabel('$F_{(C_1,C_2)}(c_1,c_2)$','Interpreter',interpreter)
mysaveas(pathname,'cdf_C',formats);
mymatlab2tikz(pathname,'cdf_C.tex');

% Plot cdf of (E,N)
figure('Name','Cumulative density function of (E,N)')
clf
xxn = xn(2:end-1);
[XE,XN] = meshgrid(xe,xxn);
Z = zeros(size(XE));
for i=1:length(xe)
    for j=1:length(xxn)
        Z(j,i) = cdf_EN(xe(i),xn(j));
    end
end
surf(XE,XN,Z);
hold on
Z_data = zeros(1,length(E_data));
for i=1:length(E_data)
    Z_data(i) = cdf_EN(E_data(i),nu_data(i));
end
plot3(E_data,nu_data,Z_data,'r+');
hold off
grid on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[min(xe),max(xe)])
set(gca,'YLim',[min(xn),max(xn)])
xlabel('$e$ (GPa)','Interpreter',interpreter)
ylabel('$n$','Interpreter',interpreter)
zlabel('$F_{(E,N)}(e,n)$','Interpreter',interpreter)
mysaveas(pathname,'cdf_EN',formats);
mymatlab2tikz(pathname,'cdf_EG.tex');

%%
N = 3000; % nombre de realisations  
C1_Model = gamrnd(a1,b1,N,1);
C2_Model = gamrnd(a2,b2,N,1);
lambda_Model = C1_Model-2/3*C2_Model;
nu_Model = (3*C1_Model-2*C2_Model)./(6*C1_Model+2*C2_Model);
E_Model = (9*C1_Model.*C2_Model)./(3*C1_Model+C2_Model);

figure('name','random samples of (E,nu)')
clf
plot(nu_Model,E_Model,'b+')
grid on
box on
xlabel('Poisson ratio','Fontsize',fontsize,'Interpreter','latex');
ylabel('Young modulus $E$ (GPa)','Fontsize',fontsize,'Interpreter','latex');
mysaveas(pathname,'samples_EN',formats);
% mymatlab2tikz(pathname,'samples_EN.tex');
