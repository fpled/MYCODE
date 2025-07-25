% clc
clearvars
close all
myparallel('start');

%% Inputs
displayGaussianField = true;
displayElasticityField = true;

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','shinozuka');
if ~exist(pathname,'dir')
    mkdir(pathname);
end
fontsize = 12;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

Dim = 2; % space dimension Dim = 2, 3
% storage = 'node'; % storage at nodal points
storage = 'gauss'; % storage at gauss points

n3D = 6; % size of 3D elasticity matrix
n = Dim*(Dim+1)/2; % size of elasticity matrix
nU = n*(n+1)/2; % number of Gaussian random fields
% nU = 1; % number of Gaussian random fields
N = 4; % number of independent realizations for each Gaussian random field
nV = nU*N; % number of independent realizations for all Gaussian random fields

nu = 2^3; % one-dimensional order (number of terms in each spatial dimension) of the spectral representation
order = nu^Dim; % Dim-dimensional order (number of terms) of the spectral representation

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
    C = LINE([0.0,b],[a,b]);
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
% S = gmshDomainWithSingleEdgeCrack(D,C,cl,cl,fullfile(pathname,'gmsh_domain_single_edge_crack'));
% S = gmsh(D,C,cl,cl,fullfile(pathname,'gmsh_domain_single_edge_crack'));

S_scalar = final(S);
% S_scalar = final(S_scalar,'duplicate');

%% Materials
% Material symmetry
symmetry = 'isotropic';
% Option
option = 'DEFO'; % plane strain [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
% option = 'CONT'; % plane stress [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
switch lower(symmetry)
    case 'isotropic' % isotropic material
        % Lame coefficients
        % lambda = 121.1538e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
        % mu = 80.7692e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
        % lambda = 121.154e9; % [Hesch, Weinberg, 2014, IJNME]
        % mu = 80.769e9; % [Hesch, Weinberg, 2014, IJNME]
        lambda = 121.15e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
        mu = 80.77e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
        % Young modulus and Poisson ratio
        if Dim==2
            switch lower(option)
                case 'defo'
                    E = mu*(3*lambda+2*mu)/(lambda+mu); % E = 210e9;
                    NU = lambda/(lambda+mu)/2; % NU = 0.3;
                case 'cont'
                    E = 4*mu*(lambda+mu)/(lambda+2*mu);
                    NU = lambda/(lambda+2*mu);
            end
            % E = 210e9; NU = 0.2; % [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
            % E = 210e9; NU = 0.3; % [Molnar, Gravouil, 2015, FEAD], [Gerasimov, De Lorenzis, 2016, CMAME], [Zhou, Rabczuk, Zhuang, 2018, AES], [Wu, Nguyen, 2018, JMPS], [Gerasimov, De Lorenzis, 2019, CMAME], [Wu, Nguyen, Zhou, Huang, 2020, CMAME], [Kristensen, Martinez-Paneda, 2020, TAFM]
            % kappa = 121030e6; NU=0.227; lambda=3*kappa*NU/(1+NU); mu = 3*kappa*(1-2*NU)/(2*(1+NU)); E = 3*kappa*(1-2*NU); % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
        elseif Dim==3
            E = mu*(3*lambda+2*mu)/(lambda+mu);
            NU = lambda/(lambda+mu)/2;
        end
        
    case 'anisotropic' % anisotropic material
        if Dim==2
            switch lower(option)
                case 'defo'
                    % [Nguyen, Yvonnet, Waldmann, He, 2020,IJNME]
                    % Elasticity matrix in reference material coordinate system [Pa]
                    Cmat = 1e9*[65 20 0;
                        20 260 0;
                        0 0 30];
                    theta = deg2rad(ang); % clockwise material orientation angle around z-axis [rad]
                    c = cos(theta);
                    s = sin(theta);
                    % Transition matrix for elasticity matrix from material coordinate system to global coordinate system
                    P = [c^2 s^2 -c*s;
                        s^2 c^2 c*s;
                        2*c*s -2*c*s c^2-s^2];
                    % Elasticity matrix in global coordinate system [Pa]
                    Cmat = P'*Cmat*P;
                case 'cont'
                    error('Not implemented yet')
            end

            if isotropicTest
                lambda = 121.15e9;
                mu = 80.77e9;
                if strcmpi(option,'cont')
                    E = mu*(3*lambda+2*mu)/(lambda+mu);
                    NU = lambda/(lambda+mu)/2;
                    lambda = E*NU/(1-NU^2); % first Lamé coefficient
                end
                Cmat = [lambda+2*mu,lambda,0;...
                    lambda,lambda+2*mu,0;...
                    0,0,mu];
            end

        elseif Dim==3
            error('Not implemented yet')
        end
    otherwise
        error('Wrong material symmetry class');
end
% Density
RHO = 1;

% Material
switch lower(symmetry)
    case 'isotropic' % isotropic material model for isotropic material only
        mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e);
    case 'anisotropic' % anisotropic material model for all symmetry classes
        mat = ELAS_ANISOT('C',Cmat,'RHO',RHO,'DIM3',e);
    otherwise
        error('Wrong material symmetry class');
end
mat = setnumber(mat,1);
S = setoption(S,option);
S = setmaterial(S,mat);

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

fprintf('\nNumber of points  = %d',nx);
fprintf('\nNumber of fields  = %d',nU);
fprintf('\nNumber of samples = %d for each Gaussian random field',N);
fprintf('\nNumber of samples = %d for all Gaussian random fields',nV);
fprintf('\nNumber of terms   = %d in the spectral representation',order);

rng('default');
s = rng; % get current random number generator settings

%% Gaussian random fields
% Standard Shinozuka method
fprintf('\nStandard Shinozuka method\n');
tGaussShinozukaStd = tic;

V = shinozuka(x,lcorr,nU,N,'order',nu,'state',s);

timeGaussShinozukaStd = toc(tGaussShinozukaStd);
fprintf('elapsed time = %f s\n',timeGaussShinozukaStd);

% Randomized Shinozuka method
fprintf('\nRandomized Shinozuka method\n');
tGaussShinozukaRand = tic;

W = shinozukaRand(x,lcorr,nU,N,'order',order,'state',s);

timeGaussShinozukaRand = toc(tGaussShinozukaRand);
fprintf('elapsed time = %f s\n',timeGaussShinozukaRand);

%% Non-Gaussian random elasticity fields
delta = 0.1; % dispersion coefficient
mC = double(calc_opmat(S)); % mean elasticity matrix
mL = chol(mC); % upper triangular matrix of the Cholesky factor of mean elasticity matrix

% Standard Shinozuka method
fprintf('\nStandard Shinozuka method\n');
fprintf('Computing non-Gaussian random elasticity fields\n');
tElasShinozukaStd = tic;

CV = randAnisotElasField(delta,mL,shiftdim(V,1));

timeElasShinozukaStd = toc(tElasShinozukaStd);
fprintf('elapsed time = %f s\n',timeElasShinozukaStd);

% Randomized Shinozuka method
fprintf('\nRandomized Shinozuka method\n');
fprintf('Computing non-Gaussian random elasticity fields\n');
tElasShinozukaRand = tic;

CW = randAnisotElasField(delta,mL,shiftdim(W,1));

timeElasShinozukaRand = toc(tElasShinozukaRand);
fprintf('elapsed time = %f s\n',timeElasShinozukaRand);

%% Display one realization of Gaussian random fields
if displayGaussianField
    switch storage
        case 'node'
            figure('Name',['Gaussian fields - standard Shinozuka (order ' num2str(order) ')'])
            clf
            iU = 0;
            for i=1:n
                for j=i:n
                    iU = iU+1;
                    subplot(n,n,(i-1)*n+j);
                    plot_sol(S_scalar,V(:,iU,1));
                    colorbar
                    set(gca,'FontSize',fontsize)
                    title(['$U_{' num2str(iU) '}$'],'FontSize',12,'Interpreter',interpreter)
                end
            end
            mysaveas(pathname,['gaussian_fields_shinozuka_std_order_' num2str(order)],formats,renderer);

            figure('Name',['Gaussian fields - randomized Shinozuka (order ' num2str(order) ')'])
            clf
            iU = 0;
            for i=1:n
                for j=i:n
                    iU = iU+1;
                    subplot(n,n,(i-1)*n+j);
                    plot_sol(S_scalar,W(:,iU,1));
                    colorbar
                    set(gca,'FontSize',fontsize)
                    title(['$U_{' num2str(iU) '}$'],'FontSize',12,'Interpreter',interpreter)
                end
            end
            mysaveas(pathname,['gaussian_fields_shinozuka_rand_order_' num2str(order)],formats,renderer);
        case 'gauss'
            Ve = cell(getnbgroupelem(S_scalar),1);
            We = cell(getnbgroupelem(S_scalar),1);
            nbgauss = 0;
            for i=1:getnbgroupelem(S_scalar)
                elem = getgroupelem(S_scalar,i);
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

            figure('Name',['Gaussian fields - standard Shinozuka with order = ' num2str(order)])
            clf
            iU = 0;
            for i=1:n
                for j=i:n
                    iU = iU+1;
                    subplot(n,n,(i-1)*n+j);
                    plot(Ve(iU,1),S_scalar);
                    colorbar
                    set(gca,'FontSize',fontsize)
                    title(['$U_{' num2str(iU) '}$'],'FontSize',12,'Interpreter',interpreter)
                end
            end
            mysaveas(pathname,['gaussian_fields_shinozuka_std_order_' num2str(order)],formats,renderer);

            figure('Name',['Gaussian fields - randomized Shinozuka with order = ' num2str(order)])
            clf
            iU = 0;
            for i=1:n
                for j=i:n
                    iU = iU+1;
                    subplot(n,n,(i-1)*n+j);
                    plot(We(iU,1),S_scalar);
                    colorbar
                    set(gca,'FontSize',fontsize)
                    title(['$U_{' num2str(iU) '}$'],'FontSize',12,'Interpreter',interpreter)
                end
            end
            mysaveas(pathname,['gaussian_fields_shinozuka_rand_order_' num2str(order)],formats,renderer);
        otherwise
            error('Wrong storage');
    end
end

%% Display one realization of non-Gaussian random elasticity field
if displayElasticityField
    switch storage
        case 'node'
            figure('Name',['Elasticity field - standard Shinozuka (order ' num2str(order) ')'])
            clf
            for i=1:n
                for j=i:n
                    subplot(n,n,(i-1)*n+j);
                    CVij = reshape(CV(i,j,1,:),nx,1);
                    plot_sol(S_scalar,CVij);
                    colorbar
                    set(gca,'FontSize',fontsize)
                    title(['$C_{' num2str(i) num2str(j) '}$'],'FontSize',12,'Interpreter',interpreter)
                end
            end
            mysaveas(pathname,['elasticity_field_shinozuka_std_order_' num2str(order)],formats,renderer);

            figure('Name',['Elasticity field - randomized Shinozuka (order ' num2str(order) ')'])
            clf
            for i=1:n
                for j=i:n
                    subplot(n,n,(i-1)*n+j);
                    CWij = reshape(CW(i,j,1,:),nx,1);
                    plot_sol(S_scalar,CWij);
                    colorbar
                    set(gca,'FontSize',fontsize)
                    title(['$C_{' num2str(i) num2str(j) '}$'],'FontSize',12,'Interpreter',interpreter)
                end
            end
            mysaveas(pathname,['elasticity_field_shinozuka_rand_order_' num2str(order)],formats,renderer);

        case 'gauss'
            CVe = cell(getnbgroupelem(S_scalar),1);
            CWe = cell(getnbgroupelem(S_scalar),1);
            nbgauss = 0;
            for i=1:getnbgroupelem(S_scalar)
                elem = getgroupelem(S_scalar,i);
                nbelem = getnbelem(elem);
                gauss = calc_gauss(elem,'rigi');
                rep = nbgauss + (1:nbelem*gauss.nbgauss);
                CVi = reshape(CV(:,:,1,rep),n,n,nbelem,gauss.nbgauss);
                CWi = reshape(CW(:,:,1,rep),n,n,nbelem,gauss.nbgauss);
                CVe{i} = MYDOUBLEND(CVi);
                CWe{i} = MYDOUBLEND(CWi);
                nbgauss = rep(end);
                syscoordgauss = getsyscoordlocal(elem);
                fieldddl = DDL(DDLTENS4('C',syscoordgauss));
            end

            CVe = FEELEMFIELD(CVe,'storage','gauss','type','scalar','ddl',fieldddl);
            CWe = FEELEMFIELD(CWe,'storage','gauss','type','scalar','ddl',fieldddl);

            figure('Name',['Gaussian field - standard Shinozuka with order = ' num2str(order)])
            clf
            for i=1:n
                for j=i:n
                    subplot(n,n,(i-1)*n+j)
                    plot(CVe(i,j),S_scalar);
                    colorbar
                    set(gca,'FontSize',fontsize)
                    title(['$C_{' num2str(i) num2str(j) '}$'],'FontSize',12,'Interpreter',interpreter)
                end
            end
            mysaveas(pathname,['gaussian_field_shinozuka_std_order_' num2str(order)],formats,renderer);

            figure('Name',['Gaussian field - randomized Shinozuka with order = ' num2str(order)])
            clf
            for i=1:n
                for j=i:n
                    subplot(n,n,(i-1)*n+j)
                    plot(CWe(i,j),S_scalar);
                    colorbar
                    set(gca,'FontSize',fontsize)
                    title(['$C_{' num2str(i) num2str(j) '}$'],'FontSize',12,'Interpreter',interpreter)
                end
            end
            mysaveas(pathname,['gaussian_field_shinozuka_rand_order_' num2str(order)],formats,renderer);
        otherwise
            error('Wrong storage');
    end
end

myparallel('stop');
