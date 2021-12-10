% clc
clearvars
close all
myparallel('start');

%% Inputs
displayGaussianFields = false;

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','Shinozuka');
if ~exist(pathname,'dir')
    mkdir(pathname);
end
fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

Dim = 2; % space dimension Dim = 2, 3

% n3D = 6; % size of 3D elasticity matrix
% n = Dim*(Dim+1)/2; % size of elasticity matrix
% nU = n*(n+1)/2; % number of Gaussian random fields
nU = 1; % number of Gaussian random fields
N = 1; % number of independent realizations for each Gaussian random field
nV = nU*N; % number of independent realizations for all Gaussian random fields

nu = 2^3; % one-dimensional order (number of terms in each spatial dimension) of the spectral representation
order = nu^Dim; % Dim-dimensional order (number of terms) of the spectral representation

%% Domains and meshes
L = 1e-3; % [m]
a = L/2;
if Dim==1
    D = DOMAIN(1,0.0,L);
    elemtype = 'SEG2';
elseif Dim==2
    e = 1;
    D = DOMAIN(2,[0.0,0.0],[L,L]);
    c = LIGNE([0.0,L/2],[a,L/2]);
    % elemtype = 'TRI3';
    elemtype = 'QUA4';
    % elemtype = 'TRI6';
elseif Dim==3
    e = 0.1e-3;
    D = DOMAIN(3,[0.0,0.0,0.0],[L,L,e]);
    c = QUADRANGLE([0.0,L/2,0.0],[a,L/2,0.0],[a,L/2,e],[0.0,L/2,e]);
    % elemtype = 'TET4';
    elemtype = 'CUB8';
    % elemtype = 'TET10';
end
% option = 'DEFO'; % plane strain
option = 'CONT'; % plane stress
if Dim==1
    cl = 2e-6;
elseif Dim==2
    % cl = 1e-5;
    cl = 2e-6;
elseif Dim==3
    cl = 2e-5;
    % cl = 7.5e-6;
end
nbelem = repmat(round(L/cl),1,Dim);
S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'option',option);
% S = build_model(D,'cl',cl,'elemtype',elemtype,'option',option,'filename',fullfile(pathname,'gmsh_domain'));
% S = gmshdomainwithedgecrack(D,C,cl,cl,fullfile(pathname,'gmsh_domain_single_edge_crack'));
% S = gmsh(D,C,cl,cl,fullfile(pathname,'gmsh_domain_single_edge_crack'));

S = final(S);

% Gauss point coordinates
% x = calc_gausscoord(S,'rigi');
% Node coordinates
node = getnode(S);
x = getcoord(node);

nx = size(x,1); % number of points

lcorr = repmat(L/50,Dim,1); % spatial correlation lengths

fprintf('\nNumber of points  = %d',nx);
fprintf('\nNumber of fields  = %d',nU);
fprintf('\nNumber of samples = %d for each Gaussian random field',N);
fprintf('\nNumber of samples = %d for all Gaussian random fields',nV);
fprintf('\nNumber of terms   = %d in the spectral representation',order);

%% Standard Shinozuka method
fprintf('\nStandard Shinozuka method\n');

rng('default');
srng = rng; % get current random number generator settings
X = rand(2,order,nV);
Phi = X(1,:,:)*2*pi; % random phase shifts Phi uniformly distributed on [0,2*pi]
Psi = X(2,:,:); % random variables Psi uniformly distributed on [0,1]
Z = sqrt(-log(Psi)); % random amplitudes Z with values in [0,+Inf[
clear X Psi

supp = 2*pi./lcorr; % support of power spectral density (PSD) functions

beta = 1:nu;
tau = -1+(beta-1/2)*2/nu; % discretization of normalized wave number tau in [-1,1] for each spatial dimension
q = 1-abs(tau); % triangular function (tri)
if verLessThan('matlab','9.1') % compatibility (<R2016b)
    k = pi*bsxfun(@rdivide,tau,lcorr);
    s = bsxfun(@times,lcorr,q)/pi; % power spectral density (PSD) functions
else
    k = pi*tau./lcorr; % discretization of wave number k in [-pi/lcorr(j),pi/lcorr(j)] for each spatial dimension j
    s = lcorr.*q/pi; % power spectral density (PSD) functions
end
c = sqrt(supp.*s/nu);

tVscalar = tic;
Vscalar = zeros(nx,nV);
parfor i=1:nV
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    Vscalar(:,i) = shinozuka_scalar(c,z,phi,k,x);
end
timeVscalar = toc(tVscalar);
fprintf('\nScalar implementation     : elapsed time = %f s',timeVscalar);

tVcpp = tic;
Vcpp = zeros(nx,nV);
parfor i=1:nV
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    Vcpp(:,i) = shinozuka_cpp(c,z,phi,k,x);
end
timeVcpp = toc(tVcpp);
errVcpp = norm(Vcpp-Vscalar)/norm(Vscalar);
fprintf('\nC++ implementation        : elapsed time = %f s, relative error = %e',timeVcpp,errVcpp);

tVsum = tic;
Vsum = zeros(nx,nV);
parfor i=1:nV
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    Vsum(:,i) = shinozuka_sum(c,z,phi,k,x);
end
timeVsum = toc(tVsum);
errVsum = norm(Vsum-Vscalar)/norm(Vscalar);
fprintf('\nSummation implementation  : elapsed time = %f s, relative error = %e',timeVsum,errVsum);

tVvec = tic;
Vvec = zeros(nx,nV);
parfor i=1:nV
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    Vvec(:,i) = shinozuka_vec(c,z,phi,k,x);
end
timeVvec = toc(tVvec);
errVvec = norm(Vvec-Vscalar)/norm(Vscalar);
fprintf('\nVectorized implementation : elapsed time = %f s, relative error = %e',timeVvec,errVvec);
fprintf('\n');

%% Randomized Shinozuka method
fprintf('\nRandomized Shinozuka method\n');

rng(srng) % set same random number generator settings as for standard Shinozuka method
X = rand(2+Dim,order,nV);
Phi = X(1,:,:)*2*pi; % random phase shifts Phi uniformly distributed on [0,2*pi]
Psi = X(2,:,:); % random variables Psi uniformly distributed on [0,1]
Z = sqrt(-log(Psi)); % random amplitudes Z with values in [0,+Inf[
Tau = (-1+2*X(2+(1:Dim),:,:)); % random wave numbers Tau uniformly distributed on [-1,1]
clear X Psi

tWscalar = tic;
Wscalar = zeros(nx,nV);
parfor i=1:nV
    tau = Tau(:,:,i);
    q = 1-abs(tau); % triangular function (tri)
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        k = pi*bsxfun(@rdivide,tau,lcorr);
        s = bsxfun(@times,lcorr,q)/pi; % power spectral density (PSD) functions
    else
        k = pi*tau./lcorr;
        s = lcorr.*q/pi; % power spectral density (PSD) functions
    end
    c = sqrt(prod(supp.*s,1)/order);
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    Wscalar(:,i) = shinozukaRand_scalar(c,z,phi,k,x);
end
timeWscalar = toc(tWscalar);
fprintf('\nScalar implementation     : elapsed time = %f s',timeWscalar);

tWcpp = tic;
Wcpp = zeros(nx,nV);
parfor i=1:nV
    tau = Tau(:,:,i);
    q = 1-abs(tau); % triangular function (tri)
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        k = pi*bsxfun(@rdivide,tau,lcorr);
        s = bsxfun(@times,lcorr,q)/pi; % power spectral density (PSD) functions
    else
        k = pi*tau./lcorr;
        s = lcorr.*q/pi; % power spectral density (PSD) functions
    end
    c = sqrt(prod(supp.*s,1)/order);
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    Wcpp(:,i) = shinozukaRand_cpp(c,z,phi,k,x);
end
timeWcpp = toc(tWcpp);
errWcpp = norm(Wcpp-Wscalar)/norm(Wscalar);
fprintf('\nC++ implementation        : elapsed time = %f s, relative error = %e',timeWcpp,errWcpp);

tWsum = tic;
Wsum = zeros(nx,nV);
parfor i=1:nV
    tau = Tau(:,:,i);
    q = 1-abs(tau); % triangular function (tri)
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        k = pi*bsxfun(@rdivide,tau,lcorr);
        s = bsxfun(@times,lcorr,q)/pi; % power spectral density (PSD) functions
    else
        k = pi*tau./lcorr;
        s = lcorr.*q/pi; % power spectral density (PSD) functions
    end
    c = sqrt(prod(supp.*s,1)/order);
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    Wsum(:,i) = shinozukaRand_sum(c,z,phi,k,x);
end
timeWsum = toc(tWsum);
errWsum = norm(Wsum-Wscalar)/norm(Wscalar);
fprintf('\nSummation implementation  : elapsed time = %f s, relative error = %e',timeWsum,errWsum);

tWvec = tic;
Wvec = zeros(nx,nV);
parfor i=1:nV
    tau = Tau(:,:,i);
    q = 1-abs(tau); % triangular function (tri)
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        k = pi*bsxfun(@rdivide,tau,lcorr);
        s = bsxfun(@times,lcorr,q)/pi; % power spectral density (PSD) functions
    else
        k = pi*tau./lcorr;
        s = lcorr.*q/pi; % power spectral density (PSD) functions
    end
    c = sqrt(prod(supp.*s,1)/order);
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    Wvec(:,i) =  sqrt(2) * cos(phi + x*k) * (c' .* z');
end
timeWvec = toc(tWvec);
errWvec = norm(Wvec-Wscalar)/norm(Wscalar);
fprintf('\nVectorized implementation : elapsed time = %f s, relative error = %e',timeWvec,errWvec);

%% Display one realization of a Gaussian random field
V = reshape(Vcpp,nx,N,nU);
W = reshape(Wcpp,nx,N,nU);
if displayGaussianFields
    if nx==getnbnode(S)
        figure('Name',['Gaussian field - standard Shinozuka with order = ' num2str(order)])
        clf
        plot_sol(S,V(:,1,1));
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['gaussian_field_shinozuka_std_order_' num2str(order)],formats,renderer);

        figure('Name','Gaussian field - randomized Shinozuka with order = ' num2str(order)])
        clf
        plot_sol(S,W(:,1,1));
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['gaussian_field_shinozuka_rand_order_' num2str(order)],formats,renderer);
    else
        Ve = cell(getnbgroupelem(S),1);
        We = cell(getnbgroupelem(S),1);
        nbgauss = 0;
        for i=1:getnbgroupelem(S)
            elem = getgroupelem(S,i);
            nbelem = getnbelem(elem);
            gauss = calc_gauss(elem,'rigi');
            rep = nbgauss + (1:nbelem*gauss.nbgauss);
            Vi = reshape(V(rep,:,:),N,nU,nbelem,gauss.nbgauss);
            Wi = reshape(W(rep,:,:),N,nU,nbelem,gauss.nbgauss);
            Ve{i} = MYDOUBLEND(Vi);
            We{i} = MYDOUBLEND(Wi);
            nbgauss = rep(end);
        end

        Ve = FEELEMFIELD(Ve,'storage','gauss','type','scalar','ddl',DDL('V'));
        We = FEELEMFIELD(We,'storage','gauss','type','scalar','ddl',DDL('W'));

        figure('Name',['Gaussian field - standard Shinozuka with order = ' num2str(order)])
        clf
        plot(Ve(1),S);
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['gaussian_field_shinozuka_std_order_' num2str(order)],formats,renderer);

        figure('Name',['Gaussian field - randomized Shinozuka with order = ' num2str(order)])
        clf
        plot(We(1),S);
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['gaussian_field_shinozuka_rand_order_' num2str(order)],formats,renderer);
    end
end

myparallel('stop');
