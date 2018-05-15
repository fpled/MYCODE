%% Nonparametric model for symmetric positive-definite real-valued random matrices %%
%%---------------------------------------------------------------------------------%%
% [Soize, 2000]

% clc
clearvars
close all
% rng('default');
myparallel('start');

%% Input data
displayCv = true;

filename = 'nonparamtricModel';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
linewidth = 1;
markersize = 36;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Nonparametric model
n = 3;
R = randn(n); % generate a n-by-n normal random matrix
A = R*R'; % construct a symmetric positive (semi-)definite matrix
L = chol(A); % upper Cholesky factor of A so that A = L'*L

N = 1e3; % nb samples

gam = 2; % existence of second-order moments of the inverse of random matrix
delta = 0.1; % dispersion parameter
% lambda = 1/(2*delta^2)*(1-delta^2*(n-1)+trace(A)^2/trace(A^2));
lambda = round(1/(2*delta^2)*(1-delta^2*(n-1)+trace(A)^2/trace(A^2)));
lambda_inf = max(0,(gam-1)/n+(3-n)/2);
if lambda<=lambda_inf
    error('Parameter lambda = %g should be > %g',lambda,lambda_inf)
end
m = n-1+2*lambda;
sigma = sqrt(2/m);

% Univariate normal and Gamma distributions
a = (n-(1:n)'+2*lambda)/2;
a_ngam = zeros(n,n,N);
t = tic;
parfor i=1:N
    U = randn(n*(n-1)/2,1); % generate a (n*(n-1)/2)-by-1 random vector U of n*(n-1)/2 standard normal random variables U_l, centered (with zero mean) and reduced (with unit variance)
    X = sigma*U; % (n*(n-1)/2)-by-1 random vector X of n*(n-1)/2 normal random variables X_j, centered (with zero mean) and a common standard deviation sigma
    LL = triu(ones(n),1);
    LL(LL==1) = 1/sqrt(2)*X; % non diagonal part of random matrix LL
    Y = gamrnd(a,1); % generate a n-by-1 random vector Y of n Gamma random variables Y_l with shape parameter a=(n-l+2*lambda)/2 and scale parameter b=1
    LL = LL + diag(sigma*sqrt(Y)); % add diagonal part of random matrix LL
    G = LL'*LL;
    a_ngam(:,:,i) = L'*G*L;
end
err_mean_ngam = zeros(N,1);
parfor i=1:N
    err_mean_ngam(i) = norm(mean(a_ngam(:,:,1:i),3)-A,'fro')/norm(A,'fro');
end
time_ngam = toc(t);

% Check if lambda is an integer
if mod(lambda,1)==0
    C = A/m;  
    % Wishart distribution
    a_wish = zeros(n,n,N);
    t = tic;
    parfor i=1:N
        a_wish(:,:,i) = wishrnd(C,m); % generate a n-by-n Wishart random matrix a_wish with n-by-n covariance matrix C and with m degrees of freedom
    end
    err_mean_wish = zeros(N,1);
    parfor i=1:N
        err_mean_wish(i) = norm(mean(a_wish(:,:,1:i),3)-A,'fro')/norm(A,'fro');
    end
    time_wish = toc(t);
    
    % Multivariate normal distribution
    a_mvn = zeros(n,n,N);
    t = tic;
    parfor i=1:N
        X = mvnrnd(zeros(1,n),C,m)'; % generate a n-by-m matrix X of m normal random vectors X_j, centered (with zero mean vector) and a common n-by-n covariance matrix C
        a_mvn(:,:,i) = X*X';
    end
    err_mean_mvn = zeros(N,1);
    parfor i=1:N
        err_mean_mvn(i) = norm(mean(a_mvn(:,:,1:i),3)-A,'fro')/norm(A,'fro');
    end
    time_mvn = toc(t);
    
    % Multivariate standard normal distribution
    a_stdmvn = zeros(n,n,N);
    t = tic;
    parfor i=1:N
        U = mvnrnd(zeros(1,n),eye(n),m)'; % generate a n-by-m matrix U of m standard normal random vectors U_j, centered (with zero mean vector) and reduced (with identity covariance matrix)
        X = 1/sqrt(m)*L'*U;
        a_stdmvn(:,:,i) = X*X';
    end
    err_mean_stdmvn = zeros(N,1);
    parfor i=1:N
        err_mean_stdmvn(i) = norm(mean(a_stdmvn(:,:,1:i),3)-A,'fro')/norm(A,'fro');
    end
    time_stdmvn = toc(t);
end

%% Statistical outputs
fprintf('\nNb samples = %g\n',N);

fprintf('expect(A) =\n');
disp(A);

fprintf('Univariate normal and Gamma distributions\n');
fprintf('mean(A) =\n');
disp(mean(a_ngam,3));
fprintf('error = %e\n',err_mean_ngam(end));
fprintf('elapsed time = %f s\n',time_ngam);

if mod(lambda,1)==0 % check if lambda is an integer
    fprintf('\nWishart distribution\n');
    fprintf('mean(A) =\n');
    disp(mean(a_wish,3));
    fprintf('error = %e\n',err_mean_wish(end));
    fprintf('elapsed time = %f s\n',time_wish);
    
    fprintf('\nMultivariate normal distribution\n');
    fprintf('mean(A) =\n');
    disp(mean(a_mvn,3));
    fprintf('error = %e\n',err_mean_mvn(end));
    fprintf('elapsed time = %f s\n',time_mvn);
    
    fprintf('\nMultivariate standard normal distribution\n');
    fprintf('mean(A) =\n');
    disp(mean(a_stdmvn,3));
    fprintf('error = %e\n',err_mean_stdmvn(end));
    fprintf('elapsed time = %f s\n',time_stdmvn);
end

%% Display convergence
if displayCv
    figure('Name','Convergence empirical mean')
    clf
    loglog(1:N,err_mean_ngam,'-b','LineWidth',linewidth)
    hold on
    if mod(lambda,1)==0 % check if lambda is an integer
        loglog(1:N,err_mean_wish,'-r','LineWidth',linewidth)
        loglog(1:N,err_mean_mvn,'-g','LineWidth',linewidth)
        loglog(1:N,err_mean_stdmvn,'-m','LineWidth',linewidth)
    end
    loglog(1:N,1./sqrt(1:N),'-k','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    % xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    % ylabel('Erreur relative','Interpreter',interpreter)
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Relative error','Interpreter',interpreter)
    if mod(lambda,1)==0
        l = legend('Univariate normal and Gamma distributions','Wishart distribution','Multivariate normal distribution','Multivariate standard normal distribution','$1/\sqrt{N}$');
    else
        l = legend('Univariate normal and Gamma distributions','$1/\sqrt{N}$');
    end
    set(l,'Interpreter',interpreter);
    mysaveas(pathname,'convergence_empirical_mean','fig');
    mymatlab2tikz(pathname,'convergence_empirical_mean.tex');
end

myparallel('stop');
