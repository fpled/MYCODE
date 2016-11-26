%% Nonparametric model for symmetric positive-matrix real-valued random matrices %%
%%-------------------------------------------------------------------------------%%
% [Soize 2000]

% clc
% clear all
% close all
% set(0,'DefaultFigureVisible','off');
rng('default');
% myparallel('start');

%% Input data

filename = 'nonparamtric_model';
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

%% Nonparametric model

n = 5;
B = randn(n); % generate a n-by-n normal random matrix
A = B*B'; % construct a symmetric positive matrix
R = chol(A); % upper Cholesky factor of A so that A = R'*R

delta = 0.1;
% lambda = 1/(2*delta^2)*(1-delta^2*(n-1)+trace(A)^2/trace(A^2));
lambda = 3;
m = n-1+2*lambda;
sigma = sqrt(2/m); % sigma = sqrt(2/(n-1+2*lambda));

N = 1e2;

% Case lambda is a positive integer

C = A/m;
Ar = zeros(n,n,N);
parfor i=1:N
    Ar(:,:,i) = wishrnd(C,m); % generate a n-by-n Wishart random matrix Ar with n-by-n covariance matrix C and with m degrees of freedom
end
Ar = sum(Ar,end)/N;

X = mvnrnd(0,C,m); % generate a m-by-n matrix X of normal random vectors X_j, centered (with zero mean vector) and a common n-by-n covariance matrix C
Ar = X'*X;
% Case lambda is not an integer



