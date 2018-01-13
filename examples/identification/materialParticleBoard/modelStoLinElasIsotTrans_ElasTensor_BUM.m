%% Blockwise updating method %%
%%---------------------------%%
% This method sometimes doesn't work, because the positive-definite of
% covariance matrix SIGMA isn't always guaranteed.
% If there is a error, rerun the code.
% See page 26 in 'Computational Statistics with Matlab'
% See also page 4 in 'Metropolis-Hastings sampling using multivariate gaussian tangents'

% clc
clearvars
close all

%% Input data
displaySolution = true;

filename = 'modelStoElasIsotTransTensor_BUM';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
linewidth = 1;
markersize = 36;
interpreter = 'latex';
formats = {'fig','epsc2'};

% syms c1 c2 c3 lambda lambda1 lambda2 lambda3
% PHI = -lambda1*c1 - lambda2*c2 - lambda3*c3 - lambda*log(c1*c2-c3^2);
% SIGMA = -hessian(PHI,[c1 c2 c3])^(-1)
% gx = gradient(PHI,c3)

%% Data
filenameAna = 'data_ET_GL.mat';
filenameNum = 'data_EL_nuL.mat';
pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathnameIdentification,filenameAna));
load(fullfile(pathnameIdentification,filenameNum));

sample = 'B';
ET_data = zeros(1,27);
GL_data = zeros(1,27);
EL_data = zeros(1,27);
NUL_data = zeros(1,27);
for j=1:27
    sampleNum = [sample num2str(j)];
    ET_data(j) = eval(['mean_ET_' sampleNum '_data;']); % GPa
    GL_data(j) = eval(['mean_GL_' sampleNum '_data;'])*1e-3; % GPa
    EL_data(j) = eval(['mean_EL_' sampleNum '_data;'])*1e-3; % GPa
    NUL_data(j) = eval(['mean_NUL_' sampleNum '_data;']);
end
NUT_data = 0.1+0.2*rand(1,27); % 27 artificial datas for NUT varying from 0.1 to 0.3
GT_data = ET_data./(2*(1+NUT_data));
kT_data = (EL_data.*ET_data)./(2*(1-NUT_data).*EL_data-4*ET_data.*(NUL_data).^2);
C1_data = EL_data + 4*(NUL_data.^2).*kT_data;
C2_data = 2*kT_data;
C3_data = 2*sqrt(2)*kT_data.*NUL_data;
C4_data = 2*GT_data;
C5_data = 2*GL_data;

mc1 = mean(C1_data);
mc2 = mean(C2_data);
mc3 = mean(C3_data);
mc4 = mean(C4_data);
mc5 = mean(C5_data);

%% Sample generation
% Method to calculate lambda which controls the level of fluctuations
lambda = -40; % negative number
lambda1 = -(mc2*lambda)/(-mc3^2+mc1*mc2);
lambda2 = -(mc1*lambda)/(-mc3^2+mc1*mc2);
lambda3 = (2*mc3*lambda)/(-mc3^2+mc1*mc2);
lambda4 = (1-2*lambda)/mc4;
lambda5 = (1-2*lambda)/mc5;

% Parameters of the trivariate normal distribution
Mu_mean = [mc1 mc2 mc3];
Sigma_mean = [-mc1^2/lambda -mc3^2/lambda -(mc1*mc3)/lambda
              -mc3^2/lambda -mc2^2/lambda -(mc2*mc3)/lambda
              -(mc1*mc3)/lambda -(mc2*mc3)/lambda -(mc3^2+mc1*mc2)/(2*lambda)];
f = @(c) -lambda1*c(1)-lambda2*c(2)-lambda3*c(3)-lambda*log(c(1)*c(2)-c(3)^2);

T = 5000;
C_sample = zeros(5,T);
t = 0;
while t<T
    t = t+1;
    c_old = mvnrnd(Mu_mean,Sigma_mean);
    c1_old = c_old(1);
    c2_old = c_old(2);
    c3_old = c_old(3);
    f_old = f(c_old);
    g_old = [-lambda1-(c2_old*lambda)/(-c3_old^2+c1_old*c2_old)
             -lambda2-(c1_old*lambda)/(-c3_old^2+c1_old*c2_old)
             -lambda3+(2*c3_old*lambda)/(-c3_old^2+c1_old*c2_old)];
    sigma_old = [-c1_old^2/lambda -c3_old^2/lambda -(c1_old*c3_old)/lambda
                 -c3_old^2/lambda -c2_old^2/lambda -(c2_old*c3_old)/lambda
                 -(c1_old*c3_old)/lambda -(c2_old*c3_old)/lambda -(c3_old^2+c1_old*c2_old)/(2*lambda)];
    mu_old = (c_old' + sigma_old*g_old)';
    q_old = mvnpdf(c_old,mu_old,sigma_old);
    %
    c_prop = mvnrnd(mu_old,sigma_old);
    c1_prop = c_prop(1);
    c2_prop = c_prop(2);
    c3_prop = c_prop(3);
    log_q_prop = log(q_old);
    f_prop = f(c_prop);
    g_prop = [-lambda1-(c2_prop*lambda)/(-c3_prop^2+c1_prop*c2_prop)
              -lambda2-(c1_prop*lambda)/(-c3_prop^2+c1_prop*c2_prop)
              -lambda3+(2*c3_prop*lambda)/(-c3_prop^2+c1_prop*c2_prop)];
    sigma_prop = [-c1_prop^2/lambda -c3_prop^2/lambda -(c1_prop*c3_prop)/lambda
                  -c3_prop^2/lambda -c2_prop^2/lambda -(c2_prop*c3_prop)/lambda
                  -(c1_prop*c3_prop)/lambda -(c2_prop*c3_prop)/lambda -(c3_prop^2+c1_prop*c2_prop)/(2*lambda)];
    mu_prop = (c_prop' + sigma_prop*g_prop)';
    q_prop = mvnpdf(c_prop,mu_prop,sigma_prop);
    log_q_old = log(q_prop);
    r = exp( (f_prop-f_old)+(log_q_old-log_q_prop) );
    if r>1
        C_sample(1:3,t) = c_prop;
    else
        s=rand;
        if s<r
            C_sample(1:3,t) = c_prop;
        else
            C_sample(1:3,t) = c_old;
        end
    end
end
C_sample(4,:) = gamrnd(1-2*lambda,1/lambda4,T,1);
C_sample(5,:) = gamrnd(1-2*lambda,1/lambda5,T,1);

kT_sample = C_sample(2,:)/2;
NUL_sample = (C_sample(3,:)./kT_sample)/(2*sqrt(2));
EL_sample = C_sample(1,:) - 4*(NUL_sample.^2).*kT_sample;
GT_sample = C_sample(4,:)/2;
GL_sample = C_sample(5,:)/2;
k1 = 1./kT_sample+4*(NUL_sample.^2)./EL_sample;
k2 = 1./GT_sample;
ET_sample = 4./(k1+k2);
NUT_sample = (ET_sample./GT_sample)/2-1;

%% Plot samples
if displaySolution
    figure('name','Samples of (ET,GL,NUT)')
    clf
    scatter3(ET_sample,GL_sample,NUT_sample,'b.')
    hold on
    scatter3(ET_data,GL_data,NUT_data,'r+')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Young modulus $E^T$ (GPa)','Interpreter',interpreter);
    ylabel('Shear modulus $G^L$ (GPa)','Interpreter',interpreter);
    zlabel('Poisson ratio $\nu^T$','Interpreter',interpreter);
    mysaveas(pathname,'samples_ET_GL_NUT',formats);
    mymatlab2tikz(pathname,'samples_ET_GL_NUT.tex');
end

%% Plot pdfs and cdfs
if displaySolution
    pdf_C4 = @(c4) gampdf(c4,1-2*lambda,1/lambda4); % Gamma probability density function of C4
    pdf_C5 = @(c5) gampdf(c5,1-2*lambda,1/lambda5); % Gamma probability density function of C5
    cdf_C4 = @(c4) gamcdf(c4,1-2*lambda,1/lambda4); % Gamma cumulative density function of C4
    cdf_C5 = @(c5) gamcdf(c5,1-2*lambda,1/lambda5); % Gamma cumulative density function of C5
    
    N = 100;
    x4_min = max(0,mean(C4_data)-8*std(C4_data));
    x4_max = mean(C4_data)+8*std(C4_data);
    x5_min = max(0,mean(C5_data)-8*std(C5_data));
    x5_max = mean(C5_data)+8*std(C5_data);
    
    % Plot pdf of C4
    figure('Name','Probability density function of C4')
    clf
    x4 = linspace(x4_min,x4_max,N);
    y4 = pdf_C4(x4);
    plot(x4,y4,'-b');
    hold on
    plot(C4_data,pdf_C4(C4_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(x4),max(x4)])
    xlabel('$c_4$ (GPa)','Interpreter',interpreter)
    ylabel('$p_{C_4}(c_4)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C4',formats);
    mymatlab2tikz(pathname,'pdf_C4.tex');
    
    % Plot pdf of C5
    figure('Name','Probability density function of C5')
    clf
    x5 = linspace(x5_min,x5_max,N);
    y5 = pdf_C5(x5);
    plot(x5,y5,'-b');
    hold on
    plot(C5_data,pdf_C5(C5_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(x5),max(x5)])
    xlabel('$c_5$ (GPa)','Interpreter',interpreter)
    ylabel('$p_{C_5}(c_5)$','Interpreter',interpreter)
    mysaveas(pathname,'pdf_C5',formats);
    mymatlab2tikz(pathname,'pdf_C5.tex');
    
    % Plot cdf of C4
    figure('name','cumulative distribution function of C4')
    clf
    z4 = cdf_C4(x4);
    plot(x4,z4,'-b');
    hold on
    plot(C4_data,cdf_C4(C4_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(x4),max(x4)])
    xlabel('$c_4$ (GPa)','Interpreter',interpreter)
    ylabel('$F_{C_4}(c_4)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_C4',formats);
    mymatlab2tikz(pathname,'cdf_C4.tex');
    
    % Plot cdf of C5
    figure('name','cumulative distribution function of C5')
    clf
    z5 = cdf_C5(x5);
    plot(x5,z5,'-b');
    hold on
    plot(C5_data,cdf_C5(C5_data),'r+');
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[min(x5),max(x5)])
    xlabel('$c_5$ (GPa)','Interpreter',interpreter)
    ylabel('$F_{C_5}(c_5)$','Interpreter',interpreter)
    mysaveas(pathname,'cdf_C5',formats);
    mymatlab2tikz(pathname,'cdf_C5.tex');
end
