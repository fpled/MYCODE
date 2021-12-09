% clc
clearvars
close all
% rng('default');
myparallel('start');

%% Inputs
displayGaussianGerms = true;
displayCorrelationStructure = true;
displayPdf = true; % display probability density function

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','Shinozuka');
if ~exist(pathname,'dir')
    mkdir(pathname);
end
% sprintf('_%g',lcorr)
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
% N = 10.^(1:5); % number of independent realizations for each Gaussian random field
N = 1e3;
nV = nU*N; % number of independent realizations for all Gaussian random fields

nu = 2.^(2:6); % one-dimensional order (number of terms in each spatial dimension) of the spectral representation
% nu = 2^3;
order = nu.^Dim; % Dim-dimensional order (number of terms) of the spectral representation

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

%% Analytical computation of normalized autocorrelation function indexed by the center point as reference point
idxm = find(x(:,2)==L/2);
idym = find(x(:,1)==L/2);
idxc = intersect(idxm,idym);
xc = x(idxc,:);
xm = x(idxm,1);
ym = x(idym,2);
[xm, repxm] = sort(xm);
[ym, repym] = sort(ym);

% Normalized autocorrelation function at center point
if verLessThan('matlab','9.1') % compatibility (<R2016b)
    corrAna = prod(sinc(bsxfun(@rdivide,x-xc,2*lcorr')).^2,2);
else
    corrAna = prod(sinc((x-xc)./(2*lcorr')).^2,2);
end
corrAnaX = corrAna(idxm);
corrAnaY = corrAna(idym);
corrAnaX = corrAnaX(repxm);
corrAnaY = corrAnaY(repym);

%% Numerical computation of normalized autocorrelation function with the standard and randomized Shinozuka methods
s = rng; % get current random number generator settings
for k=1:length(nu)
    nuk = nu(k); % one-dimensional order (number of terms in each spatial dimension) of the spectral representation
    orderk = order(k); % Dim-dimensional order (number of terms) of the spectral representation

    fprintf('\nNumber of terms   = %d in the spectral representation',orderk);
    fprintf('\n');

    %% Standard Shinozuka method
    fprintf('\nStandard Shinozuka method\n');

    %% Gaussian germs
    fprintf('\nComputation of Gaussian germs\n');
    tShinozuka = tic;

    V = shinozuka(x,lcorr,nU,max(N),'order',nuk,'state',s);

    timeShinozuka = toc(tShinozuka);
    fprintf('\nelapsed time = %f s\n',timeShinozuka);

    %% Normalized autocorrelation function at center point
    fprintf('\nComputation of autocorrelation function\n');
    tCorrV = tic;

    corrV = corr(V',V(idxc,:)'); % Dim-dimensional autocorrelation function at center point
    corrVX = corrV(idxm); % one-dimensional autocorrelation function at center point along x axis
    corrVY = corrV(idym); % one-dimensional autocorrelation function at center point along y axis
    corrVX = corrVX(repxm);
    corrVY = corrVY(repym);

    errCorrV = norm(corrV - corrAna)/norm(corrAna);
    errCorrVX = norm(corrVX - corrAnaX)/norm(corrAnaX);
    errCorrVY = norm(corrVY - corrAnaY)/norm(corrAnaY);

    timeCorrV = toc(tCorrV);
    fprintf('\nelapsed time = %f s\n',timeCorrV);

    % Save variables
    save(fullfile(pathname,['gaussian_germ_Shinozuka_std_order_' num2str(orderk) '.mat']),...
        'V','corrV','corrVX','corrVY','timeShinozuka','timeCorrV');

    %% Randomized Shinozuka method
    fprintf('\nRandomized Shinozuka method\n');

    %% Gaussian germs
    fprintf('\nComputation of Gaussian germs\n');
    tShinozukaRand = tic;

    W = shinozukaRand(x,lcorr,nU,max(N),'order',orderk,'state',s);

    timeShinozukaRand = toc(tShinozukaRand);
    fprintf('\nelapsed time = %f s\n',timeShinozukaRand);

    %% Normalized autocorrelation function at center point
    fprintf('\nComputation of autocorrelation function\n');
    tCorrW = tic;
    
    corrW = corr(W',W(idxc,:)'); % Dim-dimensional autocorrelation function at center point
    corrWX = corrW(idxm); % one-dimensional autocorrelation function at center point along x axis
    corrWY = corrW(idym); % one-dimensional autocorrelation function at center point along y axis
    corrWX = corrWX(repxm);
    corrWY = corrWY(repym);

    errCorrW = norm(corrW - corrAna)/norm(corrAna);
    errCorrWX = norm(corrWX - corrAnaX)/norm(corrAnaX);
    errCorrWY = norm(corrWY - corrAnaY)/norm(corrAnaY);

    timeCorrW = toc(tCorrW);
    fprintf('\nelapsed time = %f s\n',timeCorrW);

    % Save variables
    save(fullfile(pathname,['gaussian_germ_Shinozuka_rand_order_' num2str(orderk) '.mat']),...
        'W','corrW','corrWX','corrWY','timeShinozukaRand','timeCorrW');
end

%% Display one realization of a Gaussian germ
if displayGaussianGerms
    for k=1:length(nu)
    orderk = order(k);
    load(fullfile(pathname,['gaussian_germ_Shinozuka_std_order_' num2str(orderk) '.mat']),'V');
    load(fullfile(pathname,['gaussian_germ_Shinozuka_rand_order_' num2str(orderk) '.mat']),'W');
    
    if nx==getnbnode(S)
        figure('Name',['Gaussian germ - Standard Shinozuka method with order ' num2str(orderk)])
        clf
        plot_sol(S,V(:,1,1));
        xlabel('$x$ [m]','Interpreter',interpreter)
        if Dim>=2
            ylabel('$y$ [m]','Interpreter',interpreter)
        end
        if Dim==3
            zlabel('$z$ [m]','Interpreter',interpreter)
        end
        axis on
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['gaussian_germ_Shinozuka_std_order_' num2str(orderk)],formats,renderer);

        figure('Name',['Gaussian germ - Randomized Shinozuka method with order ' num2str(orderk)])
        clf
        plot_sol(S,W(:,1,1));
        xlabel('$x$ [m]','Interpreter',interpreter)
        if Dim>=2
            ylabel('$y$ [m]','Interpreter',interpreter)
        end
        if Dim==3
            zlabel('$z$ [m]','Interpreter',interpreter)
        end
        axis on
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['gaussian_germ_Shinozuka_rand_order_' num2str(orderk)],formats,renderer);
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

        figure('Name',['Gaussian germ - Standard Shinozuka method with order ' num2str(orderk)])
        clf
        plot(Ve(1),S);
        xlabel('$x$ [m]','Interpreter',interpreter)
        if Dim>=2
            ylabel('$y$ [m]','Interpreter',interpreter)
        end
        if Dim==3
            zlabel('$z$ [m]','Interpreter',interpreter)
        end
        axis on
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['gaussian_germ_Shinozuka_std_order_' num2str(orderk)],formats,renderer);

        figure('Name',['Gaussian germ - Randomized Shinozuka method with order ' num2str(orderk)])
        clf
        plot(We(1),S);
        xlabel('$x$ [m]','Interpreter',interpreter)
        if Dim>=2
            ylabel('$y$ [m]','Interpreter',interpreter)
        end
        if Dim==3
            zlabel('$z$ [m]','Interpreter',interpreter)
        end
        axis on
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['gaussian_germ_Shinozuka_rand_order_' num2str(orderk)],formats,renderer);
    end
    end
end

%% Display correlation structure
if displayCorrelationStructure
    figure('Name','Autocorrelation function')
    clf
    plot(corrAna,S);
    xlabel('$x$ [m]','Interpreter',interpreter)
    ylabel('$y$ [m]','Interpreter',interpreter)
    axis on
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'autocorrelation_function',formats,renderer);

    for k=1:length(nu)
        orderk = order(k);
        load(fullfile(pathname,['gaussian_germ_Shinozuka_std_order_' num2str(orderk) '.mat']),'corrV');
        load(fullfile(pathname,['gaussian_germ_Shinozuka_rand_order_' num2str(orderk) '.mat']),'corrW');

        figure('Name',['Autocorrelation function - Standard Shinozuka method with order ' num2str(orderk)])
        clf
        plot(corrV,S);
        xlabel('$x$ [m]','Interpreter',interpreter)
        ylabel('$y$ [m]','Interpreter',interpreter)
        axis on
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['autocorrelation_function_Shinozuka_std_order_' num2str(orderk)],formats,renderer);

        figure('Name',['Autocorrelation function - Randomized Shinozuka method with order ' num2str(orderk)])
        clf
        plot(corrW,S);
        xlabel('$x$ [m]','Interpreter',interpreter)
        ylabel('$y$ [m]','Interpreter',interpreter)
        axis on
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['autocorrelation_function_Shinozuka_rand_order_' num2str(orderk)],formats,renderer);
    end

    figure('Name','Autocorrelation function along x axis')
    plot(xm,corrAnaX,'r-');
    hold on
    leg = cell(1+2*length(nu),1);
    leg{1} = 'Reference';
    for k=1:length(nu)
        orderk = order(k);
        load(fullfile(pathname,['gaussian_germ_Shinozuka_std_order_' num2str(orderk) '.mat']),'corrVX');
        load(fullfile(pathname,['gaussian_germ_Shinozuka_rand_order_' num2str(orderk) '.mat']),'corrWX');
        plot(xm,corrVX,'-','Color',getfacecolor(k+1));
        plot(xm,corrWX,'--','Color',getfacecolor(k+1));
        leg{2*k} = ['Standard Shinozuka - order = ' num2str(orderk)];
        leg{2*k+1} = ['Randomized Shinozuka - order = ' num2str(orderk)];
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$x$ [m]','Interpreter',interpreter)
    ylabel('$\rho(x - x_c)$','Interpreter',interpreter)
    legend(leg{:})
    mysaveas(pathname,'autocorrelation_function_x_axis',formats,renderer);

    figure('Name','Autocorrelation function along y axis')
    plot(ym,corrAnaY,'r--');
    hold on
    leg = cell(1+2*length(nu),1);
    leg{1} = 'Reference';
    for k=1:length(nu)
        orderk = order(k);
        load(fullfile(pathname,['gaussian_germ_Shinozuka_std_order_' num2str(orderk) '.mat']),'corrVY');
        load(fullfile(pathname,['gaussian_germ_Shinozuka_rand_order_' num2str(orderk) '.mat']),'corrWY');
        plot(ym,corrVY,'-','Color',getfacecolor(i+1));
        plot(ym,corrWY,'--','Color',getfacecolor(i+1));
        leg{2*k} = ['Standard Shinozuka - order = ' num2str(orderk)];
        leg{2*k+1} = ['Randomized Shinozuka - order = ' num2str(orderk)];
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$y$ [m]','Interpreter',interpreter)
    ylabel('$\rho(y - y_c)$','Interpreter',interpreter)
    legend(leg{:})
    mysaveas(pathnamei,'autocorrelation_function_y_axis',formats,renderer);
end

myparallel('stop');
