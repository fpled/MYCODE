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
C_data = [C1_data(:) C2_data(:)];

mC_data = mean(C_data,1);
vC_data = var(C_data,0,1);
% vC_data = size(C_data,1)/(size(C_data,1)-1)*moment(C_data,2,1);
sC_data = sqrt(norm(vC_data));
dC_data = sC_data/norm(mC_data);
phiC_data = log(96*C_data(:,1).*C_data(:,2).^5);
nuC_data = mean(phiC_data,1);

%% Parameter estimation for computing Lagrange multipliers
lambda = mleStoLinElasTensorIsot(C_data); % Maximum likelihood estimation
% lambda = lseStoLinElasTensorIsot(C_data); % Least-squares estimation

la1 = lambda(1);
la2 = lambda(2);
la  = lambda(3);

a1 = 1-la;
b1 = 1/la1;
a2 = 1-5*la;
b2 = 1/la2;

mC = [a1*b1 a2*b2];
vC = [a1*b1^2 a2*b2^2];
sC = sqrt(norm(vC));
dC = sC/norm(mC);
nuC = log(96) + psi(a1)+log(b1) + 5*(psi(a2)+log(b2));

k1ln = -gammaln(a1)-a1*log(b1);
k2ln = -gammaln(a2)-a2*log(b2);
% k1ln = (1-la)*log(la1)-gammaln(1-la);
% k2ln = (1-5*la)*log(la2)-gammaln(1-5*la);

%% Pdfs and cdfs
pdf_C1 = @(c1) gampdf(c1,a1,b1); % Gamma probability density function of C1
pdf_C2 = @(c2) gampdf(c2,a2,b2); % Gamma probability density function of C2
pdf_C = @(c1,c2) pdf_C1(c1).*pdf_C2(c2); % joint probability density function of C=(C1,C2)
% pdf_C = @(c1,c2) exp(k1ln+k2ln)*c1.^(-la)*c2.^(-5*la)*exp(-la1*c1-la2*c2);
pdf_EN = @(e,n) exp(k1ln+k2ln)*(e./(3*(1-2*n))).^(-la).*(e./(2*(1+n))).^(-5*la)...
    .*(e./(2*((1+n).^2).*((1-2*n).^2)))...
    .*exp(-la1*e./(3*(1-2*n))-la2*e./(2*(1+n))); % joint probability density function of (E,N)

cdf_C1 = @(c1) gamcdf(c1,a1,b1); % Gamma cumulative density function of C1
cdf_C2 = @(c2) gamcdf(c2,a2,b2); % Gamma cumulative density function of C2
cdf_C = @(c1,c2) cdf_C1(c1).*cdf_C2(c2); % joint cumulative density function of C=(C1,C2)
% cdf_C1 = @(c1) integral(@(x) pdf_C1(x),-Inf,c1);
% cdf_C2 = @(c2) integral(@(x) pdf_C2(x),-Inf,c2);
% cdf_C = @(c1,c2) integral2(@(x1,x2) pdf_C(x1,x2),-Inf,c1,-Inf,c2);
cdf_EN = @(e,n) integral2(@(xe,xn) pdf_EN(xe,xn),0,e,-1,n);

%% Sample generation
N = 1e4; % number of samples
C1_sample = gamrnd(a1,b1,N,1);
C2_sample = gamrnd(a2,b2,N,1);
C_sample = [C1_sample(:) C2_sample(:)];

lambda_sample = C1_sample-2/3*C2_sample;
E_sample = (9*C1_sample.*C2_sample)./(3*C1_sample+C2_sample);
NU_sample = (3*C1_sample-2*C2_sample)./(6*C1_sample+2*C2_sample);

mC_sample = mean(C_sample,1);
vC_sample = var(C_sample,0,1);
% vC_sample = size(C_sample,1)/(size(C_sample,1)-1)*moment(C_sample,2,1);
sC_sample = sqrt(norm(vC_sample));
dC_sample = sC_sample/norm(mC_sample);
phiC_sample = log(96*C1_sample.*C2_sample.^5);
nuC_sample = mean(phiC_sample,1);

%% Outputs
fprintf('\nnb data   = %g',size(C_data,1));
fprintf('\nnb sample = %g',N);
fprintf('\n');

fprintf('\nlambda_1 = %.4f',la1);
fprintf('\nlambda_2 = %.4f',la2);
fprintf('\nlambda   = %.4f',la);
fprintf('\n');
fprintf('\nalpha_1 = %.4f',a1);
fprintf('\nbeta_1  = %.4f',b1);
fprintf('\n');
fprintf('\nalpha_2 = %.4f',a2);
fprintf('\nbeta_2  = %.4f',b2);
fprintf('\n');

for i=1:2
    fprintf('\nmean(C%u)        = %.4f GPa',i,mC(i));
    fprintf('\nmean(C%u_sample) = %.4f GPa',i,mC_sample(i));
    fprintf('\nmean(C%u_data)   = %.4f GPa',i,mC_data(i));
    fprintf('\nvar(C%u)         = %.4f (GPa)^2',i,vC(i));
    fprintf('\nvar(C%u_sample)  = %.4f (GPa)^2',i,vC_sample(i));
    fprintf('\nvar(C%u_data)    = %.4f (GPa)^2',i,vC_data(i));
    fprintf('\nstd(C%u)         = %.4f GPa',i,sqrt(vC(i)));
    fprintf('\nstd(C%u_sample)  = %.4f GPa',i,sqrt(vC_sample(i)));
    fprintf('\nstd(C%u_data)    = %.4f GPa',i,sqrt(vC_data(i)));
    fprintf('\ndisp(C%u)        = %.4f',i,sqrt(vC(i))/mC(i));
    fprintf('\ndisp(C%u_sample) = %.4f',i,sqrt(vC_sample(i))/mC_sample(i));
    fprintf('\ndisp(C%u_data)   = %.4f',i,sqrt(vC_data(i))/mC_data(i));
    fprintf('\n');
end
fprintf('\nnu(C)        = mean(log(det([C]))        = %.4f',i,nuC);
fprintf('\nnu(C_sample) = mean(log(det([C_sample])) = %.4f',i,nuC_sample);
fprintf('\nnu(C_data)   = mean(log(det([C_data]))   = %.4f',i,nuC_data);

err_mean = norm(mC - mC_data)^2/norm(mC_data)^2;
err_nu = (nuC - nuC_data)^2/(nuC_data)^2;
err_mean_sample = norm(mC_sample - mC_data)^2/norm(mC_data)^2;
err_nu_sample = (nuC_sample - nuC_data)^2/(nuC_data)^2;
fprintf('\nerror on mean(C) = %.4e',err_mean);
fprintf('\nerror on mean(C_sample) = %.4e',err_mean_sample);
fprintf('\nerror on mean(log(det([C])) = %.4e',err_nu);
fprintf('\nerror on mean(log(det([C_sample])) = %.4e',err_nu_sample);
fprintf('\n');

%% Display
if displaySolution
    %% Plot data
    figure('Name','Data for E')
    clf
    bar(1:length(E_data),E_data)
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[0,length(E_data)+1])
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Young''s modulus $E$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module d''Young $E$ [GPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_E','fig');
    mymatlab2tikz(pathname,'data_E.tex');
    
    figure('Name','Data for G')
    clf
    bar(1:length(G_data),G_data)
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[0,length(G_data)+1])
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Shear modulus $G$ [GPa]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Module de cisaillement $G$ [GPa]','Interpreter',interpreter);
    mysaveas(pathname,'data_G','fig');
    mymatlab2tikz(pathname,'data_G.tex');
    
    %% Plot pdfs and cdfs
    x1_min = max(0,mean(C1_data)-5*std(C1_data));
    x1_max = mean(C1_data)+5*std(C1_data);
    x2_min = max(0,mean(C2_data)-10*std(C2_data));
    x2_max = mean(C2_data)+10*std(C2_data);
    
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
    surf(X1,X2,Z);
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
    surf(Xe,Xn,Z);
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
    mysaveas(pathname,'pdf_E_N',formats);
    mymatlab2tikz(pathname,'pdf_E_N.tex');
    
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
    mysaveas(pathname,'cdf_E_N',formats);
    mymatlab2tikz(pathname,'cdf_E_N.tex');

    %% Plot samples
    figure('name','Samples of C=(C_1,C_2)')
    clf
    scatter(C1_sample,C2_sample,'b.')
    hold on
    scatter(C1_data,C2_data,'r+')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$C_1$ [GPa]','Interpreter',interpreter);
    ylabel('$C_2$ [GPa]','Interpreter',interpreter);
    %legend('samples','data')
    legend('réalisations','données')
    mysaveas(pathname,'samples_C1_C2',formats);
    mymatlab2tikz(pathname,'samples_C1_C2.tex');
    
    figure('name','Samples of (E,N)')
    clf
    scatter(NU_sample,E_sample,'b.')
    hold on
    scatter(NU_data,E_data,'r+')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$\nu$','Interpreter',interpreter);
    ylabel('$E$ [GPa]','Interpreter',interpreter);
    %legend('samples','data')
    legend('réalisations','données')
    mysaveas(pathname,'samples_E_N',formats);
    mymatlab2tikz(pathname,'samples_E_N.tex');
    
end
