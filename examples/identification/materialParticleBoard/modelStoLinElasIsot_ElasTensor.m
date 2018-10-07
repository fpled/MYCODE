%% Stochastic modeling of random linear elasticity tensor %%
%  with isotropic symmetry                                %%
%%--------------------------------------------------------%%

% clc
clearvars
close all
% rng('default');

%% Input data
displaySolution = true;

filename = 'modelStoLinElasIsot_ElasTensor';
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
filenameAna = 'data_ET_GL.mat';
filenameNum = 'data_EL_NUL.mat';
pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathnameIdentification,filenameAna));
load(fullfile(pathnameIdentification,filenameNum));

E_data = mean_ET_data*1e-3; % GPa
G_data = mean_GL_data*1e-3*13; % GPa
NU_data = E_data./(2*G_data)-1;
lambda_data = E_data.*NU_data./((1+NU_data).*(1-2*NU_data)); % GPa
C1_data = lambda_data + 2/3*G_data; % GPa
C2_data = G_data; % GPa
fprintf('\nnb data = %d',length(E_data));

%% Maximum likelihood estimation
lambda = mleStoLinElasTensorIsot(C1_data,C2_data);
fprintf('\nlambda_1 = %.4f',lambda(1));
fprintf('\nlambda_2 = %.4f',lambda(2));
fprintf('\nlambda   = %.4f',lambda(3));
fprintf('\n');

a1 = 1-lambda(3);
b1 = 1/lambda(1);
a2 = 1-5*lambda(3);
b2 = 1/lambda(2);
fprintf('\nalpha_1 = %.4f',a1);
fprintf('\nbeta_1  = %.4f',b1);
fprintf('\n');

fprintf('\nalpha_2 = %.4f',a2);
fprintf('\nbeta_2  = %.4f',b2);
fprintf('\n');

mC1 = a1*b1;
vC1 = a1*b1^2;
sC1 = sqrt(vC1);
dC1 = sC1/mC1; % dC1 = 1/sqrt(a1);
fprintf('\nmean(C1) = %.4f GPa',mC1);
fprintf('\nvar(C1)  = %.4f (GPa)^2',vC1);
fprintf('\nstd(C1)  = %.4f GPa',sC1);
fprintf('\ndisp(C1) = %.4f',dC1);
fprintf('\n');

mC2 = a2*b2;
vC2 = a2*b2^2;
sC2 = sqrt(vC2);
dC2 = sC2/mC2; % dC2 = 1/sqrt(a2);
fprintf('\nmean(C2) = %.4f GPa',mC2);
fprintf('\nvar(C2)  = %.4f (GPa)^2',vC2);
fprintf('\nstd(C2)  = %.4f GPa',sC2);
fprintf('\ndisp(C2) = %.4f',dC2);
fprintf('\n');

k1ln = -gammaln(a1)-a1*log(b1);
k2ln = -gammaln(a2)-a2*log(b2);
% k1ln = (1-lambda(3))*log(lambda(1))-gammaln(1-lambda(3));
% k2ln = (1-5*lambda(3))*log(lambda(2))-gammaln(1-5*lambda(3));

%% Pdfs and cdfs
pdf_C1 = @(c1) gampdf(c1,a1,b1); % Gamma probability density function of C1
pdf_C2 = @(c2) gampdf(c2,a2,b2); % Gamma probability density function of C2
pdf_C = @(c1,c2) pdf_C1(c1).*pdf_C2(c2); % joint probability density function of C=(C1,C2)
% pdf_C = @(c1,c2) exp(k1ln+k2ln)*c1.^(-lambda(3))*c2.^(-5*lambda(3))*exp(-lambda(1)*c1-lambda(2)*c2);
pdf_EN = @(e,n) exp(k1ln+k2ln)*(e./(3*(1-2*n))).^(-lambda(3)).*(e./(2*(1+n))).^(-5*lambda(3))...
    .*(e./(2*((1+n).^2).*((1-2*n).^2)))...
    .*exp(-lambda(1)*e./(3*(1-2*n))-lambda(2)*e./(2*(1+n))); % joint probability density function of (E,N)

cdf_C1 = @(c1) gamcdf(c1,a1,b1); % Gamma cumulative density function of C1
cdf_C2 = @(c2) gamcdf(c2,a2,b2); % Gamma cumulative density function of C2
cdf_C = @(c1,c2) cdf_C1(c1).*cdf_C2(c2); % joint cumulative density function of C=(C1,C2)
% cdf_C1 = @(c1) integral(@(x) pdf_C1(x),-Inf,c1);
% cdf_C2 = @(c2) integral(@(x) pdf_C2(x),-Inf,c2);
% cdf_C = @(c1,c2) integral2(@(x1,x2) pdf_C(x1,x2),-Inf,c1,-Inf,c2);
cdf_EN = @(e,n) integral2(@(xe,xn) pdf_EN(xe,xn),0,e,-1,n);

%% Sample generation
N = 1e3; % number of samples
C1_sample = gamrnd(a1,b1,N,1);
C2_sample = gamrnd(a2,b2,N,1);
lambda_sample = C1_sample-2/3*C2_sample;
E_sample = (9*C1_sample.*C2_sample)./(3*C1_sample+C2_sample);
NU_sample = (3*C1_sample-2*C2_sample)./(6*C1_sample+2*C2_sample);

%% Display
if displaySolution
    %% Plot data
    figure('Name','Data of E')
    clf
    bar(1:length(E_data),E_data)
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[0,length(E_data)+1])
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Young''s modulus $E$ [GPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_E','fig');
    mymatlab2tikz(pathname,'data_E.tex');
    
    figure('Name','Data of G')
    clf
    bar(1:length(G_data),G_data)
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[0,length(G_data)+1])
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Shear modulus $G$ [GPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_G','fig');
    mymatlab2tikz(pathname,'data_G.tex');
    
    %% Plot pdfs and cdfs
    x1_min = max(0,mean(C1_data)-4*std(C1_data));
    x1_max = mean(C1_data)+4*std(C1_data);
    x2_min = max(0,mean(C2_data)-8*std(C2_data));
    x2_max = mean(C2_data)+8*std(C2_data);
    
    xe_min = max(0,mean(E_data)-10*std(E_data));
    xe_max = mean(E_data)+10*std(E_data);
    xn_min = max(-1,mean(NU_data)-5*std(NU_data));
    xn_max = min(0.5,mean(NU_data)+5*std(NU_data));
    
    % Plot pdf of C1
    figure('Name','Probability density function of C1')
    clf
    x1 = linspace(x1_min,x1_max,1e2);
    y1 = pdf_C1(x1);
    plot(x1,y1,'-b');
    hold on
    plot(C1_data,pdf_C1(C1_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(x1),max(x1)])
    xlabel('$c_1$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_1}(c_1)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C1',formats);
    mymatlab2tikz(pathname,'pdf_C1.tex');
    
    % Plot pdf of C2
    figure('Name','Probability density function of C2')
    clf
    x2 = linspace(x2_min,x2_max,1e2);
    y2 = pdf_C2(x2);
    plot(x2,y2,'-b');
    hold on
    plot(C2_data,pdf_C2(C2_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(x2),max(x2)])
    xlabel('$c_2$ [GPa]','Interpreter',interpreter)
    ylabel('$p_{C_2}(c_2)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C2',formats);
    mymatlab2tikz(pathname,'pdf_C2.tex');
    
    % Plot pdf of C=(C1,C2)
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
    xlabel('$c_1$ [GPa]','Interpreter',interpreter)
    ylabel('$c_2$ [GPa]','Interpreter',interpreter)
    zlabel('$p_{(C_1,C_2)}(c_1,c_2)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C',formats);
    mymatlab2tikz(pathname,'pdf_C.tex');
    
    % Plot pdf of (E,N)
    figure('Name','Probability density function of (E,N)')
    clf
    xe = linspace(xe_min,xe_max,1e2);
    xn = linspace(xn_min,xn_max,1e2);
    [Xe,Xn] = meshgrid(xe,xn);
    Z = pdf_EN(Xe,Xn);
    surfc(Xe,Xn,Z);
    hold on
    plot3(E_data,NU_data,pdf_EN(E_data,NU_data),'r+');
    hold off
    grid on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(xe),max(xe)])
    set(gca,'YLim',[min(xn),max(xn)])
    xlabel('$e$ [GPa]','Interpreter',interpreter)
    ylabel('$n$','Interpreter',interpreter)
    zlabel('$p_{(E,N)}(e,n)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_EN',formats);
    mymatlab2tikz(pathname,'pdf_EN.tex');
    
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
    xlabel('$c_1$ [GPa]','Interpreter',interpreter)
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
    xlabel('$c_2$ [GPa]','Interpreter',interpreter)
    ylabel('$F_{C_2}(c_2)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_C2',formats);
    mymatlab2tikz(pathname,'cdf_C2.tex');
    
    % Plot cdf of C=(C1,C2)
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
    xlabel('$c_1$ [GPa]','Interpreter',interpreter)
    ylabel('$c_2$ [GPa]','Interpreter',interpreter)
    zlabel('$F_{(C_1,C_2)}(c_1,c_2)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_C',formats);
    mymatlab2tikz(pathname,'cdf_C.tex');
    
    % Plot cdf of (E,N)
    figure('Name','Cumulative density function of (E,N)')
    clf
    [XE,XN] = meshgrid(xe,xn);
    Z = zeros(size(XE));
    for i=1:length(xe)
        for j=1:length(xn)
            Z(j,i) = cdf_EN(xe(i),xn(j));
        end
    end
    surf(XE,XN,Z);
    hold on
    Z_data = zeros(1,length(E_data));
    for i=1:length(E_data)
        Z_data(i) = cdf_EN(E_data(i),NU_data(i));
    end
    plot3(E_data,NU_data,Z_data,'r+');
    hold off
    grid on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(xe),max(xe)])
    set(gca,'YLim',[min(xn),max(xn)])
    xlabel('$e$ [GPa]','Interpreter',interpreter)
    ylabel('$n$','Interpreter',interpreter)
    zlabel('$F_{(E,N)}(e,n)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_EN',formats);
    mymatlab2tikz(pathname,'cdf_EN.tex');

    %% Plot samples
    figure('name','Samples of (E,N)')
    clf
    scatter(NU_sample,E_sample,'b.')
    hold on
    scatter(NU_data,E_data,'r+')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Poisson''s ratio $\nu$','Interpreter',interpreter);
    ylabel('Young''s modulus $E$ [GPa]','Interpreter',interpreter);
    legend('sample','data')
    mysaveas(pathname,'samples_EN',formats);
    mymatlab2tikz(pathname,'samples_EN.tex');
    
end
