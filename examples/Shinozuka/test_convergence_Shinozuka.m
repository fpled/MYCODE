% clc
clearvars
close all
% rng('default');
myparallel('start');

%% Inputs
displayGaussianGerms = false;

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
fprintf('\n');

s = rng('default');

%% Standard Shinozuka method
fprintf('\nStandard Shinozuka method\n');
tShinozuka = tic;

V = shinozuka(x,lcorr,nU,N,'order',nu,'state',s);

timeShinozuka = toc(tShinozuka);
fprintf('\nelapsed time = %f s\n',timeShinozuka);

%% Randomized Shinozuka method
fprintf('\nRandomized Shinozuka method\n');
tShinozukaRand = tic;

W = shinozukaRand(x,lcorr,nU,N,'order',order,'state',s);

timeShinozukaRand = toc(tShinozukaRand);
fprintf('\nelapsed time = %f s\n',timeShinozukaRand);

%% Display Gaussian germs
if displayGaussianGerms
    if nx==getnbnode(S)
        figure('Name','Gaussian germ - Standard Shinozuka method')
        clf
        plot_sol(S,V(:,1,1));
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,'gaussian_germ_Shinozuka_std',formats,renderer);

        figure('Name','Gaussian germ - Randomized Shinozuka method')
        clf
        plot_sol(S,W(:,1,1));
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,'gaussian_germ_Shinozuka_rand',formats,renderer);
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

        figure('Name','Gaussian germ - Standard Shinozuka method')
        clf
        plot(Ve(1),S);
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,'gaussian_germ_Shinozuka_std',formats,renderer);

        figure('Name','Gaussian germ - Randomized Shinozuka method')
        clf
        plot(We(1),S);
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,'gaussian_germ_Shinozuka_rand',formats,renderer);
    end
end

myparallel('stop');
