% clc
clearvars
close all
myparallel('start');

%% Inputs
computeGaussianField = true;
computeAutocorrelation = true;
displayGaussianFields = false;
displayAutocorrelation = false;

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','shinozuka');
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
storage = 'node'; % storage at nodal points
% storage = 'gauss'; % storage at gauss points

% n3D = 6; % size of 3D elasticity matrix
% n = Dim*(Dim+1)/2; % size of elasticity matrix
% nU = n*(n+1)/2; % number of Gaussian random fields
nU = 1; % number of Gaussian random fields

% N = [1e1 5e1 1e2 2.5e2 5e2 7.5e2 1e3 2.5e3 5e3 7.5e3 1e4]; % number of independent realizations for each Gaussian random field
N = [1e1 5e1 1e2];

nV = nU*N; % number of independent realizations for all Gaussian random fields

% nu = [4 8 16 20 24 28 32 36 40 44 48 56 64]; % one-dimensional order (number of terms in each spatial dimension) of the spectral representation
% nu = [4 8 16 20 24];
% nu = [64 68];
% nu = [124 128];
order = nu.^Dim; % Dim-dimensional order (number of terms) of the spectral representation

%% Domains and meshes
L = 1e-3; % [m]
a = L/2;
b = L/2;
if Dim==1
    D = DOMAIN(1,0.0,L);
    elemtype = 'SEG2';
elseif Dim==2
    e = 1;
    D = DOMAIN(2,[0.0,0.0],[L,L]);
    C = LIGNE([0.0,b],[a,b]);
    % elemtype = 'TRI3';
    elemtype = 'QUA4';
    % elemtype = 'TRI6';
elseif Dim==3
    e = 0.1e-3;
    D = DOMAIN(3,[0.0,0.0,0.0],[L,L,e]);
    C = QUADRANGLE([0.0,b,0.0],[a,b,0.0],[a,b,e],[0.0,b,e]);
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

lcorr = repmat(L/50,Dim,1); % spatial correlation lengths
% lcorr = repmat(7.5e-6,Dim,1);
% lcorr = repmat(4e-6,Dim,1);

fprintf('\nNumber of points  = %d',nx);
fprintf('\nNumber of fields  = %d',nU);
fprintf('\nNumber of samples =%s for each Gaussian random field',sprintf(' %g',N));
fprintf('\nNumber of samples =%s for all Gaussian random fields',sprintf(' %g',nV));
fprintf('\nNumber of terms   =%s in the spectral representation',sprintf(' %g',order));
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
rng('default');
s = rng; % get current random number generator settings
if computeAutocorrelation
    errCorrV = zeros(length(order),length(N));
    errCorrVX = zeros(length(order),length(N));
    errCorrVY = zeros(length(order),length(N));
    errCorrW = zeros(length(order),length(N));
    errCorrWX = zeros(length(order),length(N));
    errCorrWY = zeros(length(order),length(N));
end
for k=1:length(order)
    nuk = nu(k); % one-dimensional order (number of terms in each spatial dimension) of the spectral representation
    orderk = order(k); % Dim-dimensional order (number of terms) of the spectral representation

    fprintf('\nOrder = %d (number of terms in the spectral representation)',orderk);
    fprintf('\n');

    %% Gaussian random fields
    if computeGaussianField
        %% Standard Shinozuka method
        fprintf('\nStandard Shinozuka method\n');
        tGaussShinozukaStd = tic;

        V = shinozuka(x,lcorr,nU,max(N),'order',nuk,'state',s);

        timeGaussShinozukaStd = toc(tGaussShinozukaStd);
        fprintf('elapsed time = %f s\n',timeGaussShinozukaStd);

        %% Randomized Shinozuka method
        fprintf('\nRandomized Shinozuka method\n');
        tGaussShinozukaRand = tic;

        W = shinozukaRand(x,lcorr,nU,max(N),'order',orderk,'state',s);

        timeGaussShinozukaRand = toc(tGaussShinozukaRand);
        fprintf('elapsed time = %f s\n',timeGaussShinozukaRand);

        save(fullfile(pathname,['gauss_shinozuka_std_order_' num2str(orderk) '.mat']),'V','timeGaussShinozukaStd');
        save(fullfile(pathname,['gauss_shinozuka_rand_order_' num2str(orderk) '.mat']),'W','timeGaussShinozukaRand');
    elseif (displayGaussianFields || computeAutocorrelation)
        load(fullfile(pathname,['gauss_shinozuka_std_order_' num2str(orderk) '.mat']),'V','timeGaussShinozukaStd');
        load(fullfile(pathname,['gauss_shinozuka_rand_order_' num2str(orderk) '.mat']),'W','timeGaussShinozukaRand');
    end

    %% Normalized autocorrelation function at center point
    if computeAutocorrelation
        %% Standard Shinozuka method
        fprintf('\nStandard Shinozuka method\n');
        tCorrShinozukaStd = tic;

        fprintf('Computing autocorrelation function\n');
        corrV = zeros(nx,length(N));
        corrVX = zeros(length(idxm),length(N));
        corrVY = zeros(length(idym),length(N));
        for j=1:length(N)
            idN = 1:N(j);
            corrV(:,j) = corr(V(:,idN)',V(idxc,idN)'); % Dim-dimensional autocorrelation function at center point
            corrVX(:,j) = corrV(idxm,j); % one-dimensional autocorrelation function at center point along x axis
            corrVY(:,j) = corrV(idym,j); % one-dimensional autocorrelation function at center point along y axis
            corrVX(:,j) = corrVX(repxm,j);
            corrVY(:,j) = corrVY(repym,j);

            errCorrV(k,j) = norm(corrV(:,j) - corrAna)/norm(corrAna);
            errCorrVX(k,j) = norm(corrVX(:,j) - corrAnaX)/norm(corrAnaX);
            errCorrVY(k,j) = norm(corrVY(:,j) - corrAnaY)/norm(corrAnaY);
        end
        timeCorrShinozukaStd = toc(tCorrShinozukaStd);
        fprintf('elapsed time = %f s\n',timeCorrShinozukaStd);

        %% Randomized Shinozuka method
        fprintf('\nRandomized Shinozuka method\n');
        tCorrShinozukaRand = tic;

        fprintf('Computing autocorrelation function\n');
        corrW = zeros(nx,length(N));
        corrWX = zeros(length(idxm),length(N));
        corrWY = zeros(length(idym),length(N));
        for j=1:length(N)
            idN = 1:N(j);
            corrW(:,j) = corr(W(:,idN)',W(idxc,idN)'); % Dim-dimensional autocorrelation function at center point
            corrWX(:,j) = corrW(idxm,j); % one-dimensional autocorrelation function at center point along x axis
            corrWY(:,j) = corrW(idym,j); % one-dimensional autocorrelation function at center point along y axis
            corrWX(:,j) = corrWX(repxm,j);
            corrWY(:,j) = corrWY(repym,j);

            errCorrW(k,j) = norm(corrW(:,j) - corrAna)/norm(corrAna);
            errCorrWX(k,j) = norm(corrWX(:,j) - corrAnaX)/norm(corrAnaX);
            errCorrWY(k,j) = norm(corrWY(:,j) - corrAnaY)/norm(corrAnaY);
        end
        timeCorrShinozukaRand = toc(tCorrShinozukaRand);
        fprintf('elapsed time = %f s\n',timeCorrShinozukaRand);

        save(fullfile(pathname,['autocorr_shinozuka_std_order_' num2str(orderk) '.mat']),...
            'corrV','corrVX','corrVY','timeCorrShinozukaStd');
        save(fullfile(pathname,['autocorr_shinozuka_rand_order_' num2str(orderk) '.mat']),...
            'corrW','corrWX','corrWY','timeCorrShinozukaRand');
    else
        load(fullfile(pathname,['autocorr_shinozuka_std_order_' num2str(orderk) '.mat']),...
            'corrV','corrVX','corrVY','timeCorrShinozukaStd');
        load(fullfile(pathname,['autocorr_shinozuka_rand_order_' num2str(orderk) '.mat']),...
            'corrW','corrWX','corrWY','timeCorrShinozukaRand');
    end
end
if computeAutocorrelation
    save(fullfile(pathname,'error_autocorr_shinozuka_std.mat'),...
        'errCorrV','errCorrVX','errCorrVY');
    save(fullfile(pathname,'error_autocorr_shinozuka_rand.mat'),...
        'errCorrW','errCorrWX','errCorrWY');
else
    load(fullfile(pathname,'error_autocorr_shinozuka_std.mat'),...
        'errCorrV','errCorrVX','errCorrVY');
    load(fullfile(pathname,'error_autocorr_shinozuka_rand.mat'),...
        'errCorrW','errCorrWX','errCorrWY');
end

%% Display one realization of a Gaussian random field
if displayGaussianFields
    for k=1:length(nu)
        orderk = order(k);
        load(fullfile(pathname,['gauss_shinozuka_std_order_' num2str(orderk) '.mat']),'V');
        load(fullfile(pathname,['gauss_shinozuka_rand_order_' num2str(orderk) '.mat']),'W');

        switch storage
            case 'node'
                figure('Name',['Gaussian field - standard Shinozuka with order ' num2str(orderk)])
                clf
                plot_sol(S,V(:,1,1));
                colorbar
                set(gca,'FontSize',fontsize)
                mysaveas(pathname,['gaussian_field_shinozuka_std_order_' num2str(orderk)],formats,renderer);

                figure('Name',['Gaussian field - randomized Shinozuka with order ' num2str(orderk)])
                clf
                plot_sol(S,W(:,1,1));
                colorbar
                set(gca,'FontSize',fontsize)
                mysaveas(pathname,['gaussian_field_shinozuka_rand_order_' num2str(orderk)],formats,renderer);
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

                figure('Name',['Gaussian field - standard Shinozuka (order ' num2str(orderk) ')'])
                clf
                plot(Ve(1),S);
                colorbar
                set(gca,'FontSize',fontsize)
                mysaveas(pathname,['gaussian_field_shinozuka_std_order_' num2str(orderk)],formats,renderer);

                figure('Name',['Gaussian field - randomized Shinozuka (order ' num2str(orderk) ')'])
                clf
                plot(We(1),S);
                colorbar
                set(gca,'FontSize',fontsize)
                mysaveas(pathname,['gaussian_field_shinozuka_rand_order_' num2str(orderk)],formats,renderer);
            otherwise
                error('Wrong storage');
        end
    end
end

%% Display normalized autocorrelation functions
if displayAutocorrelation
    %% Autocorrelation function
    figure('Name','Autocorrelation function')
    clf
    plot(corrAna,S);
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
    mysaveas(pathname,'autocorrelation_function',formats,renderer);

    %% Autocorrelation functions computed with standard and randomized Shinozuka methods
    for k=1:length(nu)
        orderk = order(k);
        load(fullfile(pathname,['autocorr_shinozuka_std_order_' num2str(orderk) '.mat']),'corrV');
        load(fullfile(pathname,['autocorr_shinozuka_rand_order_' num2str(orderk) '.mat']),'corrW');

        pathnameShinozukaStd = fullfile(pathname,['autocorr_shinozuka_std_order_' num2str(orderk)]);
        pathnameShinozukaRand = fullfile(pathname,['autocorr_shinozuka_rand_order_' num2str(orderk)]);
        if ~exist(pathnameShinozukaStd,'dir')
            mkdir(pathnameShinozukaStd);
        end
        if ~exist(pathnameShinozukaRand,'dir')
            mkdir(pathnameShinozukaRand);
        end

        close all
        for j=1:length(N)
            figure('Name',['Autocorrelation function - standard Shinozuka (order ' num2str(orderk) ', ' num2str(N(j)) ' samples)'])
            clf
            plot(corrV(:,j),S);
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
            
            mysaveas(pathnameShinozukaStd,['autocorr_shinozuka_std_order_' num2str(orderk) '_N_' num2str(N(j))],formats,renderer);
        end
        for j=1:length(N)
            figure('Name',['Autocorrelation function - randomized Shinozuka (order ' num2str(orderk) ', ' num2str(N(j)) ' samples)'])
            clf
            plot(corrW(:,j),S);
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
            pathnamek = fullfile(pathname,['autocorr_shinozuka_rand_order_' num2str(orderk)]);
            mysaveas(pathnameShinozukaRand,['autocorr_shinozuka_rand_order_' num2str(orderk) '_N_' num2str(N(j))],formats,renderer);
        end
    end

    %% Autocorrelation functions along each axis 
    close all
    color = distinguishable_colors(1+length(nu));

    figure('Name','Autocorrelation function along x axis - standard Shinozuka')
    clf
    plot(xm,corrAnaX,'LineStyle','-','Color',color(1,:),'LineWidth',linewidth);
    hold on
    leg = cell(1+length(nu),1);
    leg{1} = 'reference';
    for k=1:length(nu)
        orderk = order(k);
        load(fullfile(pathname,['autocorr_shinozuka_std_order_' num2str(orderk) '.mat']),'corrVX');
        plot(xm,corrVX(:,end),'LineStyle','--','Color',color(k+1,:),'LineWidth',linewidth);
        leg{k+1} = ['order = ' num2str(orderk)];
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$x$ [m]','Interpreter',interpreter)
    ylabel('$\rho(x - x_c)$','Interpreter',interpreter)
    legend(leg{:})
    mysaveas(pathname,'autocorr_x_axis_shinozuka_std',formats,renderer);
    mymatlab2tikz(pathname,'autocorr_x_axis_shinozuka_std.tex');

    figure('Name','Autocorrelation function along y axis - standard Shinozuka')
    clf
    plot(ym,corrAnaY,'LineStyle','-','Color',color(1,:),'LineWidth',linewidth);
    hold on
    leg = cell(1+length(nu),1);
    leg{1} = 'reference';
    for k=1:length(nu)
        orderk = order(k);
        load(fullfile(pathname,['autocorr_shinozuka_std_order_' num2str(orderk) '.mat']),'corrVY');
        plot(ym,corrVY(:,end),'LineStyle','--','Color',color(k+1,:),'LineWidth',linewidth);
        leg{k+1} = ['order = ' num2str(orderk)];
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$y$ [m]','Interpreter',interpreter)
    ylabel('$\rho(y - y_c)$','Interpreter',interpreter)
    legend(leg{:})
    mysaveas(pathname,'autocorr_y_axis_shinozuka_std',formats,renderer);
    mymatlab2tikz(pathname,'autocorr_y_axis_shinozuka_std.tex');

    figure('Name','Autocorrelation function along x axis - randomized Shinozuka')
    clf
    plot(xm,corrAnaX,'LineStyle','-','Color',color(1,:),'LineWidth',linewidth);
    hold on
    leg = cell(1+length(nu),1);
    leg{1} = 'reference';
    for k=1:length(nu)
        orderk = order(k);
        load(fullfile(pathname,['autocorr_shinozuka_rand_order_' num2str(orderk) '.mat']),'corrWX');
        plot(xm,corrWX(:,end),'LineStyle','--','Color',color(k+1,:),'LineWidth',linewidth);
        leg{k+1} = ['order = ' num2str(orderk)];
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$x$ [m]','Interpreter',interpreter)
    ylabel('$\rho(x - x_c)$','Interpreter',interpreter)
    legend(leg{:})
    mysaveas(pathname,'autocorr_x_axis_shinozuka_rand',formats,renderer);
    mymatlab2tikz(pathname,'autocorr_x_axis_shinozuka_rand.tex');

    figure('Name','Autocorrelation function along y axis - randomized Shinozuka')
    clf
    plot(ym,corrAnaY,'LineStyle','-','Color',color(1,:),'LineWidth',linewidth);
    hold on
    leg = cell(1+length(nu),1);
    leg{1} = 'reference';
    for k=1:length(nu)
        orderk = order(k);
        load(fullfile(pathname,['autocorr_shinozuka_rand_order_' num2str(orderk) '.mat']),'corrWY');
        plot(ym,corrWY(:,end),'LineStyle','--','Color',color(k+1,:),'LineWidth',linewidth);
        leg{k+1} = ['order = ' num2str(orderk)];
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$y$ [m]','Interpreter',interpreter)
    ylabel('$\rho(y - y_c)$','Interpreter',interpreter)
    legend(leg{:})
    mysaveas(pathname,'autocorr_y_axis_shinozuka_rand',formats,renderer);
    mymatlab2tikz(pathname,'autocorr_y_axis_shinozuka_rand.tex');

    %% Relative errors on autocorrelation functions
    figure('Name','Relative error on autocorrelation function - standard Shinozuka')
    clf
    leg = cell(length(nu),1);
    for k=1:length(nu)
        loglog(N,errCorrV(k,:),'LineStyle','-','Color',color(k,:),'LineWidth',linewidth);
        hold on
        leg{k} = ['order = ' num2str(order(k))];
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Realizations')
    ylabel('Relative error')
    legend(leg{:})
    mysaveas(pathname,'error_autocorr_shinozuka_std',formats,renderer);
    mymatlab2tikz(pathname,'error_autocorr_shinozuka_std.tex');

    figure('Name','Relative error on autocorrelation function along x axis - standard Shinozuka')
    clf
    leg = cell(length(nu),1);
    for k=1:length(nu)
        loglog(N,errCorrVX(k,:),'LineStyle','-','Color',color(k,:),'LineWidth',linewidth);
        hold on
        leg{k} = ['order = ' num2str(order(k))];
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Realizations')
    ylabel('Relative error')
    legend(leg{:})
    mysaveas(pathname,'error_autocorr_x_axis_shinozuka_std',formats,renderer);
    mymatlab2tikz(pathname,'error_autocorr_x_axis_shinozuka_std.tex');

    figure('Name','Relative error on autocorrelation function along y axis - standard Shinozuka')
    clf
    leg = cell(length(nu),1);
    for k=1:length(nu)
        loglog(N,errCorrVY(k,:),'LineStyle','-','Color',color(k,:),'LineWidth',linewidth);
        hold on
        leg{k} = ['order = ' num2str(order(k))];
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Realizations')
    ylabel('Relative error')
    legend(leg{:})
    mysaveas(pathname,'error_autocorr_y_axis_shinozuka_std',formats,renderer);
    mymatlab2tikz(pathname,'error_autocorr_y_axis_shinozuka_std.tex');

    figure('Name','Relative error on autocorrelation function - randomized Shinozuka')
    clf
    leg = cell(length(nu),1);
    for k=1:length(nu)
        loglog(N,errCorrW(k,:),'LineStyle','-','Color',color(k,:),'LineWidth',linewidth);
        hold on
        leg{k} = ['order = ' num2str(order(k))];
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Realizations')
    ylabel('Relative error')
    legend(leg{:})
    mysaveas(pathname,'error_autocorr_shinozuka_rand',formats,renderer);
    mymatlab2tikz(pathname,'error_autocorr_shinozuka_rand.tex');

    figure('Name','Relative error on autocorrelation function along x axis - randomized Shinozuka')
    clf
    leg = cell(length(nu),1);
    for k=1:length(nu)
        loglog(N,errCorrWX(k,:),'LineStyle','-','Color',color(k,:),'LineWidth',linewidth);
        hold on
        leg{k} = ['order = ' num2str(order(k))];
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Realizations')
    ylabel('Relative error')
    legend(leg{:})
    mysaveas(pathname,'error_autocorr_x_axis_shinozuka_rand',formats,renderer);
    mymatlab2tikz(pathname,'error_autocorr_x_axis_shinozuka_rand.tex');
    
    figure('Name','Relative error on autocorrelation function along y axis - randomized Shinozuka')
    clf
    leg = cell(length(nu),1);
    for k=1:length(nu)
        loglog(N,errCorrWY(k,:),'LineStyle','-','Color',color(k,:),'LineWidth',linewidth);
        hold on
        leg{k} = ['order = ' num2str(order(k))];
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Realizations')
    ylabel('Relative error')
    legend(leg{:})
    mysaveas(pathname,'error_autocorr_y_axis_shinozuka_rand',formats,renderer);
    mymatlab2tikz(pathname,'error_autocorr_y_axis_shinozuka_rand.tex');

    figure('Name','Relative error on autocorrelation function - standard Shinozuka')
    clf
    surf(N,order,errCorrV,'EdgeColor','none')
    colorbar
    set(gca,'FontSize',fontsize)
    set(gca,'Xdir','reverse','Ydir','reverse')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'ZScale','log')
    xlabel('Realizations')
    ylabel('Order')
    zlabel('Relative error')
    mysaveas(pathname,'error_autocorr_shinozuka_std_surf',formats,renderer);
    
    figure('Name','Relative error on autocorrelation function along x axis - standard Shinozuka')
    clf
    surf(N,order,errCorrVX,'EdgeColor','none')
    colorbar
    set(gca,'FontSize',fontsize)
    set(gca,'Xdir','reverse','Ydir','reverse')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'ZScale','log')
    xlabel('Realizations')
    ylabel('Order')
    zlabel('Relative error')
    mysaveas(pathname,'error_autocorr_x_axis_shinozuka_std_surf',formats,renderer);

    figure('Name','Relative error on autocorrelation function along y axis - standard Shinozuka')
    clf
    surf(N,order,errCorrVY,'EdgeColor','none')
    colorbar
    set(gca,'FontSize',fontsize)
    set(gca,'Xdir','reverse','Ydir','reverse')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'ZScale','log')
    xlabel('Realizations')
    ylabel('Order')
    zlabel('Relative error')
    mysaveas(pathname,'error_autocorr_y_axis_shinozuka_std_surf',formats,renderer);

    figure('Name','Relative error on autocorrelation function - randomized Shinozuka')
    clf
    surf(N,order,errCorrW,'EdgeColor','none')
    colorbar
    set(gca,'FontSize',fontsize)
    set(gca,'Xdir','reverse','Ydir','reverse')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'ZScale','log')
    xlabel('Realizations')
    ylabel('Order')
    zlabel('Relative error')
    mysaveas(pathname,'error_autocorr_shinozuka_rand_surf',formats,renderer);

    figure('Name','Relative error on autocorrelation function along x axis - randomized Shinozuka')
    clf
    surf(N,order,errCorrWX,'EdgeColor','none')
    colorbar
    set(gca,'FontSize',fontsize)
    set(gca,'Xdir','reverse','Ydir','reverse')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'ZScale','log')
    xlabel('Realizations')
    ylabel('Order')
    zlabel('Relative error')
    mysaveas(pathname,'error_autocorr_x_axis_shinozuka_rand_surf',formats,renderer);

    figure('Name','Relative error on autocorrelation function along x axis - randomized Shinozuka')
    clf
    surf(N,order,errCorrWY,'EdgeColor','none')
    colorbar
    set(gca,'FontSize',fontsize)
    set(gca,'Xdir','reverse','Ydir','reverse')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'ZScale','log')
    xlabel('Realizations')
    ylabel('Order')
    zlabel('Relative error')
    mysaveas(pathname,'error_autocorr_y_axis_shinozuka_rand_surf',formats,renderer);
end

myparallel('stop');
