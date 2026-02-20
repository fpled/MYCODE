% clc
clearvars
close all
myparallel('start');

%% Inputs
displayGaussianField = false;

Dim = 2; % space dimension Dim = 2, 3

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results',['shinozuka_' num2str(Dim) 'D']);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

storage = 'node'; % storage at nodal points
% storage = 'gauss'; % storage at gauss points

% n = Dim*(Dim+1)/2; % size of elasticity matrix
% nU = n*(n+1)/2; % number of Gaussian random fields
nU = 1; % number of Gaussian random fields
N = 1; % number of independent realizations for each Gaussian random field
nV = nU*N; % number of independent realizations for all Gaussian random fields

% nu = 4; % nu = 2^2; % one-dimensional order (number of terms in each spatial dimension) of the spectral representation
nu = 8; % nu = 2^3;
% nu = 16; % nu = 2^4;
% nu = 20; % nu = L/lcorr = 1e-3/5e-5; % minimal order to avoid spatial periodicity of gaussian random fields computed using standard Shinozuka, 
                                       % such that the domain extent <= one half period in each spatial direction:
                                       % L_j <= nu * lcorr_j for all spatial dimensions j in 1,...,Dim
% nu = 32; % nu = 2^5;
% nu = 48;
% nu = 64; % nu = 2^6;
order = nu^Dim; % Dim-dimensional order (number of terms) of the spectral representation

%% Domains and meshes
L = 1e-3; % domain size [m]
lcorr = repmat(L/20,Dim,1); % spatial correlation lengths

a = L/2;
b = L/2;
if Dim==1
    D = DOMAIN(1,0.0,L);
    elemtype = 'SEG2';
elseif Dim==2
    e = 1;
    D = DOMAIN(2,[0.0,0.0],[L,L]);
    C = LINE([0.0,b],[a,b]);
    % elemtype = 'TRI3';
    elemtype = 'QUA4';
    % elemtype = 'TRI6';
elseif Dim==3
    % e = 0.1e-3;
    e = 1e-3;
    D = DOMAIN(3,[0.0,0.0,0.0],[L,L,e]);
    C = QUADRANGLE([0.0,b,0.0],[a,b,0.0],[a,b,e],[0.0,b,e]);
    % elemtype = 'TET4';
    elemtype = 'CUB8';
    % elemtype = 'TET10';
end
% option = 'DEFO'; % plane strain
option = 'CONT'; % plane stress
if Dim==1
    cl = 5e-6; % cl = L/200;
elseif Dim==2
    % cl = 1e-5; % cl = L/100;
    cl = 5e-6; % cl = L/200;
elseif Dim==3
    % cl = 2e-5; % cl = L/50;
    cl = 1e-5; % cl = L/100;
    % cl = 5e-6; % cl = L/200;
end
nbelem = repmat(round(L/cl),1,Dim);
S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'option',option);
% S = build_model(D,'cl',cl,'elemtype',elemtype,'option',option,'filename',fullfile(pathname,'gmsh_domain'));
% S = gmshDomainWithSingleEdgeCrack(D,C,cl,cl,fullfile(pathname,'gmsh_domain_single_edge_crack'));
% S = gmsh(D,C,cl,cl,fullfile(pathname,'gmsh_domain_single_edge_crack'));

S = final(S);
% S = final(S,'duplicate');

switch storage
    case 'node'
        % Node coordinates
        node = getnode(S);
        x = getcoord(node);
    case 'gauss'
        % Gauss point coordinates
        x = calc_gausscoord(S,'rigi');
    otherwise
        error('Wrong storage');
end
nx = size(x,1); % number of points

fprintf('\n');
fprintf('Number of points  = %d\n',nx);
fprintf('Number of fields  = %d\n',nU);
fprintf('Number of samples = %d for each Gaussian random field\n',N);
fprintf('Number of samples = %d for all Gaussian random fields\n',nV);
fprintf('Number of terms   = %d in the spectral representation\n',order);

rng('default');
srng = rng; % get current random number generator settings

%% Standard Shinozuka
fprintf('\n');
fprintf('Standard Shinozuka\n');
fprintf('\n');

X = rand(2,order,nV);
Phi = X(1,:,:)*(2*pi); % random phase shifts Phi uniformly distributed on [0,2*pi]
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
fprintf('Scalar implementation     : elapsed time = %f s\n',timeVscalar);

tVcpp = tic;
Vcpp = zeros(nx,nV);
parfor i=1:nV
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    Vcpp(:,i) = shinozuka_cpp(c,z,phi,k,x);
end
timeVcpp = toc(tVcpp);
errVcpp = norm(Vcpp-Vscalar)/norm(Vscalar);
fprintf('C++ implementation        : elapsed time = %f s, relative error = %e\n',timeVcpp,errVcpp);

tVsum = tic;
Vsum = zeros(nx,nV);
parfor i=1:nV
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    Vsum(:,i) = shinozuka_sum(c,z,phi,k,x);
end
timeVsum = toc(tVsum);
errVsum = norm(Vsum-Vscalar)/norm(Vscalar);
fprintf('Summation implementation  : elapsed time = %f s, relative error = %e\n',timeVsum,errVsum);

tVvec = tic;
Vvec = zeros(nx,nV);
parfor i=1:nV
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    Vvec(:,i) = shinozuka_vec(c,z,phi,k,x);
end
timeVvec = toc(tVvec);
errVvec = norm(Vvec-Vscalar)/norm(Vscalar);
fprintf('Vectorized implementation : elapsed time = %f s, relative error = %e\n',timeVvec,errVvec);

%% Randomized Shinozuka
fprintf('\n');
fprintf('Randomized Shinozuka\n');
fprintf('\n');

rng(srng) % set same random number generator settings as for standard Shinozuka
X = rand(2+Dim,order,nV);
Phi = X(1,:,:)*(2*pi); % random phase shifts Phi uniformly distributed on [0,2*pi]
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
fprintf('Scalar implementation     : elapsed time = %f s\n',timeWscalar);

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
fprintf('C++ implementation        : elapsed time = %f s, relative error = %e\n',timeWcpp,errWcpp);

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
fprintf('Summation implementation  : elapsed time = %f s, relative error = %e\n',timeWsum,errWsum);

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
fprintf('Vectorized implementation : elapsed time = %f s, relative error = %e\n',timeWvec,errWvec);

V = reshape(Vcpp,nx,nU,N);
W = reshape(Wcpp,nx,nU,N);

%% Display one realization of a Gaussian random field
fpos = get(groot,'DefaultFigurePosition');
fposStd  = [fpos(1)-fpos(3)/2 fpos(2:4)];
fposRand = [fpos(1)+fpos(3)/2 fpos(2:4)];
if displayGaussianField
    switch storage
        case 'node'
            % Standard Shinozuka
            figure('Name',['Gaussian field - standard Shinozuka (order ' num2str(order) ')'],...
                'Position',fposStd)
            clf
            plot_sol(S,V(:,1,1));
            colorbar
            set(gca,'FontSize',fontsize)
            mysaveas(pathname,['gaussian_field_std_order_' num2str(order)],formats,renderer);
            
            % Randomized Shinozuka
            figure('Name',['Gaussian field - randomized Shinozuka (order ' num2str(order) ')'],...
                'Position',fposRand)
            clf
            plot_sol(S,W(:,1,1));
            colorbar
            set(gca,'FontSize',fontsize)
            mysaveas(pathname,['gaussian_field_rand_order_' num2str(order)],formats,renderer);
            
        case 'gauss'
            Ve = cell(getnbgroupelem(S),1);
            We = cell(getnbgroupelem(S),1);
            nbgauss = 0;
            for i=1:getnbgroupelem(S)
                elem = getgroupelem(S,i);
                nbelem = getnbelem(elem);
                gauss = calc_gauss(elem,'rigi');
                rep = nbgauss + (1:nbelem*gauss.nbgauss);
                Vi = reshape(V(rep,:,1)',nU,1,nbelem,gauss.nbgauss);
                Wi = reshape(W(rep,:,1)',nU,1,nbelem,gauss.nbgauss);
                Ve{i} = MYDOUBLEND(Vi);
                We{i} = MYDOUBLEND(Wi);
                nbgauss = rep(end);
            end
            
            Ve = FEELEMFIELD(Ve,'storage','gauss','type','scalar','ddl',DDL('V'));
            We = FEELEMFIELD(We,'storage','gauss','type','scalar','ddl',DDL('W'));
            
            % Standard Shinozuka
            figure('Name',['Gaussian field - standard Shinozuka (order ' num2str(order) ')'],...
                'Position',fposStd)
            clf
            plot(Ve(1),S);
            colorbar
            set(gca,'FontSize',fontsize)
            mysaveas(pathname,['gaussian_field_std_order_' num2str(order)],formats,renderer);
            
            % Randomized Shinozuka
            figure('Name',['Gaussian field - randomized Shinozuka (order ' num2str(order) ')'],...
                'Position',fposRand)
            clf
            plot(We(1),S);
            colorbar
            set(gca,'FontSize',fontsize)
            mysaveas(pathname,['gaussian_field_rand_order_' num2str(order)],formats,renderer);
        otherwise
            error('Wrong storage');
    end
end

myparallel('stop');
