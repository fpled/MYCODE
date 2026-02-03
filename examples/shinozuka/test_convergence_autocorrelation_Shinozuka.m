% clc
clearvars
close all
myparallel('start');

%% Inputs
computeGaussianField = true;
computeAutocorrelation = true;
displayGaussianFields = false;
displayAutocorrelation = false;

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

% N = [1e1 5e1 1e2 5e2 1e3 5e3 1e4]; % number of independent realizations for each Gaussian random field
N = [1e1 1e2 1e3 1e4];

nV = nU*N; % number of independent realizations for all Gaussian random fields

% nu = [4 8 16 20 24 28 32 36 40 44 48 56 64]; % one-dimensional order (number of terms in each spatial dimension) of the spectral representation
nu = [4 8 16 20 32 48 64];
order = nu.^Dim; % Dim-dimensional order (number of terms) of the spectral representation

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
fprintf('Number of samples =%s for each Gaussian random field\n',sprintf(' %g',N));
fprintf('Number of samples =%s for all Gaussian random fields\n',sprintf(' %g',nV));
fprintf('Number of terms   =%s in the spectral representation\n',sprintf(' %g',order));

%% Analytical computation of normalized autocorrelation function indexed by the center point as reference point
switch Dim
    case 1
        idxc = find(x(:,1)==L/2);
    case 2
        idxm = find(x(:,2)==L/2);
        idym = find(x(:,1)==L/2);
        idxc = intersect(idxm,idym); % idxc = find(x(:,1)==L/2 & x(:,2)==L/2);
    case 3
        idxm = find(x(:,2)==L/2 & x(:,3)==L/2);
        idym = find(x(:,1)==L/2 & x(:,3)==L/2);
        idzm = find(x(:,1)==L/2 & x(:,2)==L/2);
        idxc = intersect(intersect(idxm,idym),idzm); % idxc = find(x(:,1)==L/2 & x(:,2)==L/2 & x(:,3)==L/2);
end
xc = x(idxc,:);
if Dim>=2
    xm = x(idxm,1);
    ym = x(idym,2);
    [xm, repxm] = sort(xm);
    [ym, repym] = sort(ym);
    if Dim==3
        zm = x(idzm,3);
        [zm, repzm] = sort(zm);
    end
end

% Normalized autocorrelation function at center point
if verLessThan('matlab','9.1') % compatibility (<R2016b)
    corrAna = prod(sinc(bsxfun(@rdivide,x-xc,2*lcorr')).^2,2);
else
    corrAna = prod(sinc((x-xc)./(2*lcorr')).^2,2);
end
if Dim>=2
    corrAnaX = corrAna(idxm);
    corrAnaY = corrAna(idym);
    corrAnaX = corrAnaX(repxm);
    corrAnaY = corrAnaY(repym);
    if Dim==3
        corrAnaZ = corrAna(idzm);
        corrAnaZ = corrAnaZ(repzm);
    end
end

%% Numerical computation of normalized autocorrelation function with the standard and randomized Shinozuka methods
rng('default');
s = rng; % get current random number generator settings
if computeAutocorrelation
    errCorrV = zeros(length(order),length(N));
    errCorrW = zeros(length(order),length(N));
    if Dim>=2
        errCorrVX = zeros(length(order),length(N));
        errCorrVY = zeros(length(order),length(N));
        errCorrWX = zeros(length(order),length(N));
        errCorrWY = zeros(length(order),length(N));
        if Dim==3
            errCorrVZ = zeros(length(order),length(N));
            errCorrWZ = zeros(length(order),length(N));
        end
    end
end
for k=1:length(order)
    nuk = nu(k); % one-dimensional order (number of terms in each spatial dimension) of the spectral representation
    orderk = order(k); % Dim-dimensional order (number of terms) of the spectral representation
    
    if computeGaussianField || computeAutocorrelation
        fprintf('\n');
        fprintf('Order = %d (number of terms in the spectral representation)\n',orderk);
    end
    
    %% Gaussian random fields
    if computeGaussianField
        %% Standard Shinozuka
        fprintf('\n');
        fprintf('Standard Shinozuka\n');
        tGaussStd = tic;
        
        V = shinozuka(x,lcorr,nU,max(N),'order',nuk,'state',s);
        
        timeGaussStd = toc(tGaussStd);
        fprintf('elapsed time = %f s\n',timeGaussStd);
        
        %% Randomized Shinozuka
        fprintf('\n');
        fprintf('Randomized Shinozuka\n');
        tGaussRand = tic;
        
        W = shinozukaRand(x,lcorr,nU,max(N),'order',orderk,'state',s);
        
        timeGaussRand = toc(tGaussRand);
        fprintf('elapsed time = %f s\n',timeGaussRand);
        
        save(fullfile(pathname,['gauss_std_order_' num2str(orderk) '.mat']),'V','timeGaussStd');
        save(fullfile(pathname,['gauss_rand_order_' num2str(orderk) '.mat']),'W','timeGaussRand');
    elseif computeAutocorrelation
        load(fullfile(pathname,['gauss_std_order_' num2str(orderk) '.mat']),'V','timeGaussStd');
        load(fullfile(pathname,['gauss_rand_order_' num2str(orderk) '.mat']),'W','timeGaussRand');
    end
    
    %% Normalized autocorrelation function at center point
    if computeAutocorrelation
        %% Standard Shinozuka
        fprintf('\n');
        fprintf('Standard Shinozuka\n');
        tCorrStd = tic;
        
        fprintf('Computing autocorrelation function\n');
        corrV = zeros(nx,length(N));
        if Dim>=2
            corrVX = zeros(length(idxm),length(N));
            corrVY = zeros(length(idym),length(N));
            if Dim==3
                corrVZ = zeros(length(idzm),length(N));
            end
        end
        for j=1:length(N)
            idN = 1:N(j);
            corrV(:,j) = corr(V(:,idN)',V(idxc,idN)'); % Dim-dimensional autocorrelation function at center point
            errCorrV(k,j) = norm(corrV(:,j) - corrAna)/norm(corrAna);
            if Dim>=2
                corrVX(:,j) = corrV(idxm,j); % one-dimensional autocorrelation function at center point along x axis
                corrVY(:,j) = corrV(idym,j); % one-dimensional autocorrelation function at center point along y axis
                corrVX(:,j) = corrVX(repxm,j);
                corrVY(:,j) = corrVY(repym,j);
                errCorrVX(k,j) = norm(corrVX(:,j) - corrAnaX)/norm(corrAnaX);
                errCorrVY(k,j) = norm(corrVY(:,j) - corrAnaY)/norm(corrAnaY);
                if Dim==3
                    corrVZ(:,j) = corrV(idzm,j); % one-dimensional autocorrelation function at center point along z axis
                    corrVZ(:,j) = corrVZ(repzm,j);
                    errCorrVZ(k,j) = norm(corrVZ(:,j) - corrAnaZ)/norm(corrAnaZ);
                end
            end
        end
        timeCorrStd = toc(tCorrStd);
        fprintf('elapsed time = %f s\n',timeCorrStd);
        
        %% Randomized Shinozuka
        fprintf('\n');
        fprintf('Randomized Shinozuka\n');
        tCorrRand = tic;
        
        fprintf('Computing autocorrelation function\n');
        corrW = zeros(nx,length(N));
        if Dim>=2
            corrWX = zeros(length(idxm),length(N));
            corrWY = zeros(length(idym),length(N));
            if Dim==3
                corrWZ = zeros(length(idzm),length(N));
            end
        end
        for j=1:length(N)
            idN = 1:N(j);
            corrW(:,j) = corr(W(:,idN)',W(idxc,idN)'); % Dim-dimensional autocorrelation function at center point
            errCorrW(k,j) = norm(corrW(:,j) - corrAna)/norm(corrAna);
            if Dim>=2
                corrWX(:,j) = corrW(idxm,j); % one-dimensional autocorrelation function at center point along x axis
                corrWY(:,j) = corrW(idym,j); % one-dimensional autocorrelation function at center point along y axis
                corrWX(:,j) = corrWX(repxm,j);
                corrWY(:,j) = corrWY(repym,j);
                errCorrWX(k,j) = norm(corrWX(:,j) - corrAnaX)/norm(corrAnaX);
                errCorrWY(k,j) = norm(corrWY(:,j) - corrAnaY)/norm(corrAnaY);
                if Dim==3
                    corrWZ(:,j) = corrW(idzm,j); % one-dimensional autocorrelation function at center point along z axis
                    corrWZ(:,j) = corrWZ(repzm,j);
                    errCorrWZ(k,j) = norm(corrWZ(:,j) - corrAnaZ)/norm(corrAnaZ);
                end
            end
        end
        timeCorrRand = toc(tCorrRand);
        fprintf('elapsed time = %f s\n',timeCorrRand);
        
        save(fullfile(pathname,['autocorr_std_order_' num2str(orderk) '.mat']),'corrV','timeCorrStd');
        save(fullfile(pathname,['autocorr_rand_order_' num2str(orderk) '.mat']),'corrW','timeCorrRand');
        if Dim>=2
            save(fullfile(pathname,['autocorr_std_order_' num2str(orderk) '.mat']),'corrVX','corrVY','-append');
            save(fullfile(pathname,['autocorr_rand_order_' num2str(orderk) '.mat']),'corrWX','corrWY','-append');
            if Dim==3
                save(fullfile(pathname,['autocorr_std_order_' num2str(orderk) '.mat']),'corrVZ','-append');
                save(fullfile(pathname,['autocorr_rand_order_' num2str(orderk) '.mat']),'corrWZ','-append');
            end
        end
    end
end
if computeAutocorrelation
    save(fullfile(pathname,'error_autocorr_std.mat'),'errCorrV');
    save(fullfile(pathname,'error_autocorr_rand.mat'),'errCorrW');
    if Dim>=2
        save(fullfile(pathname,'error_autocorr_std.mat'),'errCorrVX','errCorrVY','-append');
        save(fullfile(pathname,'error_autocorr_rand.mat'),'errCorrWX','errCorrWY','-append');
        if Dim==3
            save(fullfile(pathname,'error_autocorr_std.mat'),'errCorrVZ','-append');
            save(fullfile(pathname,'error_autocorr_rand.mat'),'errCorrWZ''-append');
        end
    end
else
    load(fullfile(pathname,'error_autocorr_std.mat'),'errCorrV');
    load(fullfile(pathname,'error_autocorr_rand.mat'),'errCorrW');
    if Dim>=2
        load(fullfile(pathname,'error_autocorr_std.mat'),'errCorrVX','errCorrVY');
        load(fullfile(pathname,'error_autocorr_rand.mat'),'errCorrWX','errCorrWY');
        if Dim==3
            load(fullfile(pathname,'error_autocorr_std.mat'),'errCorrVZ');
            load(fullfile(pathname,'error_autocorr_rand.mat'),'errCorrWZ');
        end
    end
end

fpos = get(groot,'DefaultFigurePosition');
fposStd = [fpos(1)-fpos(3)/2 fpos(2:4)];
fposRand = [fpos(1)+fpos(3)/2 fpos(2:4)];

%% Display one realization of a Gaussian random field
if displayGaussianFields
    for k=1:length(order)
        orderk = order(k);
        load(fullfile(pathname,['gauss_std_order_' num2str(orderk) '.mat']),'V');
        load(fullfile(pathname,['gauss_rand_order_' num2str(orderk) '.mat']),'W');
        
        switch storage
            case 'node'
                % Standard Shinozuka
                figure('Name',['Gaussian field - standard Shinozuka (order ' num2str(orderk) ')'],...
                    'Position',fposStd)
                clf
                plot_sol(S,V(:,1,1));
                colorbar
                set(gca,'FontSize',fontsize)
                mysaveas(pathname,['gaussian_field_std_order_' num2str(orderk)],formats,renderer);
                
                % Randomized Shinozuka
                figure('Name',['Gaussian field - randomized Shinozuka (order ' num2str(orderk) ')'],...
                    'Position',fposRand)
                clf
                plot_sol(S,W(:,1,1));
                colorbar
                set(gca,'FontSize',fontsize)
                mysaveas(pathname,['gaussian_field_rand_order_' num2str(orderk)],formats,renderer);
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
                figure('Name',['Gaussian field - standard Shinozuka (order ' num2str(orderk) ')'],...
                    'Position',fposStd)
                clf
                plot(Ve(1),S);
                colorbar
                set(gca,'FontSize',fontsize)
                mysaveas(pathname,['gaussian_field_std_order_' num2str(orderk)],formats,renderer);
                
                % Randomized Shinozuka
                figure('Name',['Gaussian field - randomized Shinozuka (order ' num2str(orderk) ')'],...
                    'Position',fposRand)
                clf
                plot(We(1),S);
                colorbar
                set(gca,'FontSize',fontsize)
                mysaveas(pathname,['gaussian_field_rand_order_' num2str(orderk)],formats,renderer);
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
    for k=1:length(order)
        orderk = order(k);
        load(fullfile(pathname,['autocorr_std_order_' num2str(orderk) '.mat']),'corrV');
        load(fullfile(pathname,['autocorr_rand_order_' num2str(orderk) '.mat']),'corrW');
        
        pathnameShinozukaStd = fullfile(pathname,['autocorr_std_order_' num2str(orderk)]);
        pathnameShinozukaRand = fullfile(pathname,['autocorr_rand_order_' num2str(orderk)]);
        if ~exist(pathnameShinozukaStd,'dir')
            mkdir(pathnameShinozukaStd);
        end
        if ~exist(pathnameShinozukaRand,'dir')
            mkdir(pathnameShinozukaRand);
        end
        
        close all
        
        % Standard Shinozuka
        for j=1:length(N)
            figure('Name',['Autocorrelation function - standard Shinozuka (order ' num2str(orderk) ', ' num2str(N(j)) ' samples)'],...
                'Position',fposStd)
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
            mysaveas(pathnameShinozukaStd,['autocorr_std_order_' num2str(orderk) '_N_' num2str(N(j))],formats,renderer);
        end
        
        % Randomized Shinozuka
        for j=1:length(N)
            figure('Name',['Autocorrelation function - randomized Shinozuka (order ' num2str(orderk) ', ' num2str(N(j)) ' samples)'],...
                'Position',fposRand)
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
            mysaveas(pathnameShinozukaRand,['autocorr_rand_order_' num2str(orderk) '_N_' num2str(N(j))],formats,renderer);
        end
    end
    
    %% One-dimensional autocorrelation functions along each axis
    if Dim>=2
        close all
        colors = distinguishable_colors(1+length(order));
        
        % Standard Shinozuka - x axis
        figure('Name','Autocorrelation function along x axis - standard Shinozuka',...
            'Position',fposStd)
        clf
        plot(xm,corrAnaX,'LineStyle','-','Color',colors(1,:),'LineWidth',linewidth);
        hold on
        leg = cell(1+length(order),1);
        leg{1} = 'reference';
        for k=1:length(order)
            orderk = order(k);
            load(fullfile(pathname,['autocorr_std_order_' num2str(orderk) '.mat']),'corrVX');
            plot(xm,corrVX(:,end),'LineStyle','--','Color',colors(k+1,:),'LineWidth',linewidth);
            leg{k+1} = ['order = ' num2str(orderk)];
        end
        hold off
        grid on
        box on
        set(gca,'FontSize',fontsize)
        xlabel('$x$ [m]','Interpreter',interpreter)
        ylabel('$\rho(x - x_c)$','Interpreter',interpreter)
        legend(leg{:})
        mysaveas(pathname,'autocorr_x_axis_std',formats,renderer);
        mymatlab2tikz(pathname,'autocorr_x_axis_std.tex');
        
        % Standard Shinozuka - y axis
        figure('Name','Autocorrelation function along y axis - standard Shinozuka',...
            'Position',fposStd)
        clf
        plot(ym,corrAnaY,'LineStyle','-','Color',colors(1,:),'LineWidth',linewidth);
        hold on
        leg = cell(1+length(order),1);
        leg{1} = 'reference';
        for k=1:length(order)
            orderk = order(k);
            load(fullfile(pathname,['autocorr_std_order_' num2str(orderk) '.mat']),'corrVY');
            plot(ym,corrVY(:,end),'LineStyle','--','Color',colors(k+1,:),'LineWidth',linewidth);
            leg{k+1} = ['order = ' num2str(orderk)];
        end
        hold off
        grid on
        box on
        set(gca,'FontSize',fontsize)
        xlabel('$y$ [m]','Interpreter',interpreter)
        ylabel('$\rho(y - y_c)$','Interpreter',interpreter)
        legend(leg{:})
        mysaveas(pathname,'autocorr_y_axis_std',formats,renderer);
        mymatlab2tikz(pathname,'autocorr_y_axis_std.tex');
        
        % Standard Shinozuka - z axis
        if Dim==3
            figure('Name','Autocorrelation function along z axis - standard Shinozuka',...
                'Position',fposStd)
            clf
            plot(zm,corrAnaZ,'LineStyle','-','Color',colors(1,:),'LineWidth',linewidth);
            hold on
            leg = cell(1+length(order),1);
            leg{1} = 'reference';
            for k=1:length(order)
                orderk = order(k);
                load(fullfile(pathname,['autocorr_std_order_' num2str(orderk) '.mat']),'corrVZ');
                plot(zm,corrVZ(:,end),'LineStyle','--','Color',colors(k+1,:),'LineWidth',linewidth);
                leg{k+1} = ['order = ' num2str(orderk)];
            end
            hold off
            grid on
            box on
            set(gca,'FontSize',fontsize)
            xlabel('$z$ [m]','Interpreter',interpreter)
            ylabel('$\rho(z - z_c)$','Interpreter',interpreter)
            legend(leg{:})
            mysaveas(pathname,'autocorr_z_axis_std',formats,renderer);
            mymatlab2tikz(pathname,'autocorr_z_axis_std.tex');
        end
        
        % Randomized Shinozuka - x axis
        figure('Name','Autocorrelation function along x axis - randomized Shinozuka',...
            'Position',fposRand)
        clf
        plot(xm,corrAnaX,'LineStyle','-','Color',colors(1,:),'LineWidth',linewidth);
        hold on
        leg = cell(1+length(order),1);
        leg{1} = 'reference';
        for k=1:length(order)
            orderk = order(k);
            load(fullfile(pathname,['autocorr_rand_order_' num2str(orderk) '.mat']),'corrWX');
            plot(xm,corrWX(:,end),'LineStyle','--','Color',colors(k+1,:),'LineWidth',linewidth);
            leg{k+1} = ['order = ' num2str(orderk)];
        end
        hold off
        grid on
        box on
        set(gca,'FontSize',fontsize)
        xlabel('$x$ [m]','Interpreter',interpreter)
        ylabel('$\rho(x - x_c)$','Interpreter',interpreter)
        legend(leg{:})
        mysaveas(pathname,'autocorr_x_axis_rand',formats,renderer);
        mymatlab2tikz(pathname,'autocorr_x_axis_rand.tex');
        
        % Randomized Shinozuka - y axis
        figure('Name','Autocorrelation function along y axis - randomized Shinozuka',...
            'Position',fposRand)
        clf
        plot(ym,corrAnaY,'LineStyle','-','Color',colors(1,:),'LineWidth',linewidth);
        hold on
        leg = cell(1+length(order),1);
        leg{1} = 'reference';
        for k=1:length(order)
            orderk = order(k);
            load(fullfile(pathname,['autocorr_rand_order_' num2str(orderk) '.mat']),'corrWY');
            plot(ym,corrWY(:,end),'LineStyle','--','Color',colors(k+1,:),'LineWidth',linewidth);
            leg{k+1} = ['order = ' num2str(orderk)];
        end
        hold off
        grid on
        box on
        set(gca,'FontSize',fontsize)
        xlabel('$y$ [m]','Interpreter',interpreter)
        ylabel('$\rho(y - y_c)$','Interpreter',interpreter)
        legend(leg{:})
        mysaveas(pathname,'autocorr_y_axis_rand',formats,renderer);
        mymatlab2tikz(pathname,'autocorr_y_axis_rand.tex');
        
        % Randomized Shinozuka - z axis
        if Dim==3
            figure('Name','Autocorrelation function along z axis - randomized Shinozuka',...
                'Position',fposRand)
            clf
            plot(yz,corrAnaZ,'LineStyle','-','Color',colors(1,:),'LineWidth',linewidth);
            hold on
            leg = cell(1+length(order),1);
            leg{1} = 'reference';
            for k=1:length(order)
                orderk = order(k);
                load(fullfile(pathname,['autocorr_rand_order_' num2str(orderk) '.mat']),'corrWZ');
                plot(zm,corrWZ(:,end),'LineStyle','--','Color',colors(k+1,:),'LineWidth',linewidth);
                leg{k+1} = ['order = ' num2str(orderk)];
            end
            hold off
            grid on
            box on
            set(gca,'FontSize',fontsize)
            xlabel('$z$ [m]','Interpreter',interpreter)
            ylabel('$\rho(z - z_c)$','Interpreter',interpreter)
            legend(leg{:})
            mysaveas(pathname,'autocorr_z_axis_rand',formats,renderer);
            mymatlab2tikz(pathname,'autocorr_z_axis_rand.tex');
        end
    end
    
    %% Relative errors on autocorrelation functions
    % Standard Shinozuka - line plot
    figure('Name','Relative error on autocorrelation function - standard Shinozuka',...
        'Position',fposStd)
    clf
    leg = cell(length(order),1);
    for k=1:length(order)
        loglog(N,errCorrV(k,:),'LineStyle','-','Color',colors(k,:),'LineWidth',linewidth);
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
    mysaveas(pathname,'error_autocorr_std',formats,renderer);
    mymatlab2tikz(pathname,'error_autocorr_std.tex');
    
    % Randomized Shinozuka - line plot
    figure('Name','Relative error on autocorrelation function - randomized Shinozuka',...
        'Position',fposRand)
    clf
    leg = cell(length(order),1);
    for k=1:length(order)
        loglog(N,errCorrW(k,:),'LineStyle','-','Color',colors(k,:),'LineWidth',linewidth);
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
    mysaveas(pathname,'error_autocorr_rand',formats,renderer);
    mymatlab2tikz(pathname,'error_autocorr_rand.tex');
    
    % Standard Shinozuka - surface plot
    figure('Name','Relative error on autocorrelation function - standard Shinozuka',...
        'Position',fposStd)
    clf
    surf(N,order,errCorrV,'EdgeColor','none')
    colorbar
    set(gca,'ColorScale','log')
    set(gca,'FontSize',fontsize)
    set(gca,'Xdir','reverse','Ydir','reverse')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'ZScale','log')
    xlabel('Realizations')
    ylabel('Order')
    zlabel('Relative error')
    mysaveas(pathname,'error_autocorr_std_surf',formats,renderer);
    
    % Randomized Shinozuka - surface plot
    figure('Name','Relative error on autocorrelation function - randomized Shinozuka',...
        'Position',fposRand)
    clf
    surf(N,order,errCorrW,'EdgeColor','none')
    colorbar
    set(gca,'ColorScale','log')
    set(gca,'FontSize',fontsize)
    set(gca,'Xdir','reverse','Ydir','reverse')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'ZScale','log')
    xlabel('Realizations')
    ylabel('Order')
    zlabel('Relative error')
    mysaveas(pathname,'error_autocorr_rand_surf',formats,renderer);
    
    %% Relative errors on one-dimensional autocorrelation functions along each axis
    if Dim>=2
        % Standard Shinozuka - x axis - line plot
        figure('Name','Relative error on autocorrelation function along x axis - standard Shinozuka',...
            'Position',fposStd)
        clf
        leg = cell(length(order),1);
        for k=1:length(order)
            loglog(N,errCorrVX(k,:),'LineStyle','-','Color',colors(k,:),'LineWidth',linewidth);
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
        mysaveas(pathname,'error_autocorr_x_axis_std',formats,renderer);
        mymatlab2tikz(pathname,'error_autocorr_x_axis_std.tex');
        
        % Standard Shinozuka - y axis - line plot
        figure('Name','Relative error on autocorrelation function along y axis - standard Shinozuka',...
            'Position',fposStd)
        clf
        leg = cell(length(order),1);
        for k=1:length(order)
            loglog(N,errCorrVY(k,:),'LineStyle','-','Color',colors(k,:),'LineWidth',linewidth);
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
        mysaveas(pathname,'error_autocorr_y_axis_std',formats,renderer);
        mymatlab2tikz(pathname,'error_autocorr_y_axis_std.tex');
        
        % Standard Shinozuka - z axis - line plot
        if Dim==3
            figure('Name','Relative error on autocorrelation function along z axis - standard Shinozuka',...
                'Position',fposStd)
            clf
            leg = cell(length(order),1);
            for k=1:length(order)
                loglog(N,errCorrVZ(k,:),'LineStyle','-','Color',colors(k,:),'LineWidth',linewidth);
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
            mysaveas(pathname,'error_autocorr_z_axis_std',formats,renderer);
            mymatlab2tikz(pathname,'error_autocorr_z_axis_std.tex');
        end
        
        % Randomized Shinozuka - x axis - line plot
        figure('Name','Relative error on autocorrelation function along x axis - randomized Shinozuka',...
            'Position',fposRand)
        clf
        leg = cell(length(order),1);
        for k=1:length(order)
            loglog(N,errCorrWX(k,:),'LineStyle','-','Color',colors(k,:),'LineWidth',linewidth);
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
        mysaveas(pathname,'error_autocorr_x_axis_rand',formats,renderer);
        mymatlab2tikz(pathname,'error_autocorr_x_axis_rand.tex');
        
        % Randomized Shinozuka - y axis - line plot
        figure('Name','Relative error on autocorrelation function along y axis - randomized Shinozuka',...
            'Position',fposRand)
        clf
        leg = cell(length(order),1);
        for k=1:length(order)
            loglog(N,errCorrWY(k,:),'LineStyle','-','Color',colors(k,:),'LineWidth',linewidth);
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
        mysaveas(pathname,'error_autocorr_y_axis_rand',formats,renderer);
        mymatlab2tikz(pathname,'error_autocorr_y_axis_rand.tex');
        
        % Randomized Shinozuka - z axis - line plot
        if Dim==3
            figure('Name','Relative error on autocorrelation function along z axis - randomized Shinozuka',...
                'Position',fposRand)
            clf
            leg = cell(length(order),1);
            for k=1:length(order)
                loglog(N,errCorrWZ(k,:),'LineStyle','-','Color',colors(k,:),'LineWidth',linewidth);
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
            mysaveas(pathname,'error_autocorr_z_axis_rand',formats,renderer);
            mymatlab2tikz(pathname,'error_autocorr_z_axis_rand.tex');
        end
        
        % Standard Shinozuka - x axis - surface plot
        figure('Name','Relative error on autocorrelation function along x axis - standard Shinozuka',...
            'Position',fposStd)
        clf
        surf(N,order,errCorrVX,'EdgeColor','none')
        colorbar
        set(gca,'ColorScale','log')
        set(gca,'FontSize',fontsize)
        set(gca,'Xdir','reverse','Ydir','reverse')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        set(gca,'ZScale','log')
        xlabel('Realizations')
        ylabel('Order')
        zlabel('Relative error')
        mysaveas(pathname,'error_autocorr_x_axis_std_surf',formats,renderer);
        
        % Standard Shinozuka - y axis - surface plot
        figure('Name','Relative error on autocorrelation function along y axis - standard Shinozuka',...
            'Position',fposStd)
        clf
        surf(N,order,errCorrVY,'EdgeColor','none')
        colorbar
        set(gca,'ColorScale','log')
        set(gca,'FontSize',fontsize)
        set(gca,'Xdir','reverse','Ydir','reverse')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        set(gca,'ZScale','log')
        xlabel('Realizations')
        ylabel('Order')
        zlabel('Relative error')
        mysaveas(pathname,'error_autocorr_y_axis_std_surf',formats,renderer);
        
        % Standard Shinozuka - z axis - surface plot
        if Dim==3
            figure('Name','Relative error on autocorrelation function along z axis - standard Shinozuka',...
                'Position',fposStd)
            clf
            surf(N,order,errCorrVZ,'EdgeColor','none')
            colorbar
            set(gca,'ColorScale','log')
            set(gca,'FontSize',fontsize)
            set(gca,'Xdir','reverse','Ydir','reverse')
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            set(gca,'ZScale','log')
            xlabel('Realizations')
            ylabel('Order')
            zlabel('Relative error')
            mysaveas(pathname,'error_autocorr_z_axis_std_surf',formats,renderer);
        end
        
        % Randomized Shinozuka - x axis - surface plot
        figure('Name','Relative error on autocorrelation function along x axis - randomized Shinozuka',...
            'Position',fposRand)
        clf
        surf(N,order,errCorrWX,'EdgeColor','none')
        colorbar
        set(gca,'ColorScale','log')
        set(gca,'FontSize',fontsize)
        set(gca,'Xdir','reverse','Ydir','reverse')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        set(gca,'ZScale','log')
        xlabel('Realizations')
        ylabel('Order')
        zlabel('Relative error')
        mysaveas(pathname,'error_autocorr_x_axis_rand_surf',formats,renderer);
        
        % Randomized Shinozuka - y axis - surface plot
        figure('Name','Relative error on autocorrelation function along y axis - randomized Shinozuka',...
            'Position',fposRand)
        clf
        surf(N,order,errCorrWY,'EdgeColor','none')
        colorbar
        set(gca,'ColorScale','log')
        set(gca,'FontSize',fontsize)
        set(gca,'Xdir','reverse','Ydir','reverse')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        set(gca,'ZScale','log')
        xlabel('Realizations')
        ylabel('Order')
        zlabel('Relative error')
        mysaveas(pathname,'error_autocorr_y_axis_rand_surf',formats,renderer);
        
        % Randomized Shinozuka - z axis - surface plot
        if Dim==3
            figure('Name','Relative error on autocorrelation function along z axis - randomized Shinozuka',...
                'Position',fposRand)
            clf
            surf(N,order,errCorrWZ,'EdgeColor','none')
            colorbar
            set(gca,'ColorScale','log')
            set(gca,'FontSize',fontsize)
            set(gca,'Xdir','reverse','Ydir','reverse')
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            set(gca,'ZScale','log')
            xlabel('Realizations')
            ylabel('Order')
            zlabel('Relative error')
            mysaveas(pathname,'error_autocorr_z_axis_rand_surf',formats,renderer);
        end
    end
end

myparallel('stop');
