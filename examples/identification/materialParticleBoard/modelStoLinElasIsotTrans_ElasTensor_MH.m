%% Markov Chain Monte-Carlo (MCMC) method - Metropolis-Hasting algorithm %%
%%-----------------------------------------------------------------------%%

% clc
clearvars
close all

%% Input data
displaySolution = true;

filename = 'modelStoElasIsotTransTensor_MH';
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
filenameNum = 'data_EL_nuL.mat';
pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathnameIdentification,filenameAna));
load(fullfile(pathnameIdentification,filenameNum));

numSamples = 27;
ET_data = mean_ET_data*1e-3; % GPa
GL_data = mean_GL_data*1e-3; % GPa
EL_data = mean_EL_data*1e-3; % GPa
NUL_data = mean_NUL_data;
NUT_data = 0.1+0.2*rand(numSamples,1); % 27 artificial datas for NUT varying from 0.1 to 0.3
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
proppdf = @(c) mvnpdf([c(1) c(2) c(3)],Mu_mean,Sigma_mean);
proprnd = @(c) mvnrnd(Mu_mean,Sigma_mean);
pdf = @(c) (c(1)*c(2)-c(3)^2)^(-lambda)*exp(-lambda1*c(1)-lambda2*c(2)-lambda3*c(3));
nsamples = 5000;
start = [0 0 0];
smpl = mhsample(start,nsamples,'pdf',pdf,'proppdf',proppdf, 'proprnd',proprnd);
