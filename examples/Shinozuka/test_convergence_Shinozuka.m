% clc
clearvars
close all
% rng('default');
myparallel('start');

%% Inputs
displayGaussianGerms = true;
displayCorrelationStructure = true;
displayPdf = true; % display probability density function

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
% N = [20 100 500 2500 12500]; % number of independent realizations for each Gaussian random field
N = 100;
nV = nU*N; % number of independent realizations for all Gaussian random fields

% nu = [2^2 2^3 2^4 2^5]; % one-dimensional order (number of terms in each spatial dimension) of the spectral representation
nu = 2^3;
order = nu.^Dim; % Dim-dimensional order (number of terms) of the spectral representation

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','Shinozuka');
if ~exist(pathname,'dir')
    mkdir(pathname);
end

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
    % cl = 1.25e-6;
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

s = rng('default');

%% Analytic autocorrelation with the middle point as a reference
if displayCorrelationStructure
    IDLineX = find(x(:,2)==L/2);
    IDLineY = find(x(:,1)==L/2);
    XMid = x(IDLineX,1);
    YMid = x(IDLineY,2);
    [XMid, reOrderX] = sort(XMid);
    [YMid, reOrderY] = sort(YMid);
    nX = size(IDLineX,1);
    nY = size(IDLineY,1);
    IDCenterX = floor((nX+1)/2);
    IDCenterY = floor((nY+1)/2);
    IDCenter = intersect(IDLineX,IDLineY);
    xCenter = x(IDCenter,:);

    % Expected correlation function
    corrAna = prod(sinc((x-xCenter)/2./lcorr').^2,2); % vector of correlation to the central point
    corrAnaX = corrAna(IDLineX);
    corrAnaX = corrAnaX(reOrderX);
    corrAnaY = corrAna(IDLineY);
    corrAnaY = corrAnaY(reOrderY);

    figure('Name','Expected autocorrelation function')
    clf
    plot(corrAna,S);
    xlabel('$x$ [m]')
    ylabel('$y$ [m]')
    axis on
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Autocorrelation_Analytic_Shinozuka',formats,renderer);
end

%% Computation
for Ni = N
    nVi = nU*Ni; % number of independent realizations for all Gaussian random fields
    for nui = nu
        close all
        orderi = nui^Dim; % Dim-dimensional order (number of terms) of the spectral representation

        pathnamei = fullfile(pathname,['N' num2str(Ni) 'Order' num2str(orderi)]);
        if ~exist(pathnamei,'dir')
            mkdir(pathnamei);
        end
        pathResultsi = fullfile(pathnamei,'results.txt');
        fid = fopen(pathResultsi,'wt');

        fprintf('\nNumber of points  = %d',nx);
        fprintf('\nNumber of fields  = %d',nU);
        fprintf('\nNumber of samples = %d for each Gaussian random field',Ni);
        fprintf('\nTotal number of realizations = %d',nVi);
        fprintf('\nOrder of the approximation = %d',orderi);
        fprintf('\n');

        fprintf(fid,'\nNumber of points  = %d',nx);
        fprintf(fid,'\nNumber of fields  = %d',nU);
        fprintf(fid,'\nNumber of samples = %d for each Gaussian random field',Ni);
        fprintf(fid,'\nTotal number of realizations = %d',nVi);
        fprintf(fid,'\nOrder of the approximation = %d',orderi);
        fprintf(fid,'\n');

        %% Standard Shinozuka method
        fprintf('\nStandard Shinozuka method\n');
        tShinozuka = tic;

        V = shinozuka(x,lcorr,nU,Ni,'order',nui,'state',s);

        timeShinozuka = toc(tShinozuka);
        fprintf('\nelapsed time = %f s\n',timeShinozuka);

        fprintf(fid,'\nStandard Shinozuka method');
        fprintf(fid,'\nelapsed time = %f s\n',timeShinozuka);

        %% Randomized Shinozuka method
        fprintf('\nRandomized Shinozuka method\n');
        tShinozukaRand = tic;

        W = shinozukaRand(x,lcorr,nU,Ni,'order',orderi,'state',s);

        timeShinozukaRand = toc(tShinozukaRand);
        fprintf('\nelapsed time = %f s\n',timeShinozukaRand);

        fprintf(fid,'\nRandomized Shinozuka method');
        fprintf(fid,'\nelapsed time = %f s\n',timeShinozukaRand);

        %% Display Gaussian germs
        if displayGaussianGerms
            if nx==getnbnode(S)
                figure('Name','Gaussian germ - Standard Shinozuka method')
                clf
                plot_sol(S,V(:,1,1));
                xlabel('$x$ [m]')
                ylabel('$y$ [m]')
                axis on
                colorbar
                set(gca,'FontSize',fontsize)
                mysaveas(pathnamei,'gaussian_germ_Shinozuka_std',formats,renderer);

                figure('Name','Gaussian germ - Randomized Shinozuka method')
                clf
                plot_sol(S,W(:,1,1));
                xlabel('$x$ [m]')
                ylabel('$y$ [m]')
                axis on
                colorbar
                set(gca,'FontSize',fontsize)
                mysaveas(pathnamei,'gaussian_germ_Shinozuka_rand',formats,renderer);
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
                xlabel('$x$ [m]')
                ylabel('$y$ [m]')
                axis on
                colorbar
                set(gca,'FontSize',fontsize)
                mysaveas(pathname,'gaussian_germ_Shinozuka_std',formats,renderer);

                figure('Name','Gaussian germ - Randomized Shinozuka method')
                clf
                plot(We(1),S);
                xlabel('$x$ [m]')
                ylabel('$y$ [m]')
                axis on
                colorbar
                set(gca,'FontSize',fontsize)
                mysaveas(pathname,'gaussian_germ_Shinozuka_rand',formats,renderer);
            end
        end

        %% Display correlation structure
        if displayCorrelationStructure
            % Correlation functions on the mesh
            corrV = autocorrVec(V,IDCenter); % vector of correlation to the central point
            errCorrV = norm(corrV - corrAna)/norm(corrAna);
            % corrVold = corr(V'); % correlation matrix
            % corrVold = corrVold(:,IDCenter); % vector of correlation to the central point
            % errCorrVV = norm(corrV - corrVold)/norm(corrVold);

            corrW = autocorrVec(W,IDCenter);
            errCorrW = norm(corrW - corrAna)/norm(corrAna);
            % corrWold = corr(W');
            % corrWold = corrWold(:,IDCenter);
            % errCorrWW = norm(corrW - corrWold)/norm(corrWold);

            % Correlation functions on the central lines
            corrVX = corrV(IDLineX);
            corrVX = corrVX(reOrderX);
            corrVY = corrV(IDLineY);
            corrVY = corrVY(reOrderY);
            errCorrVX = norm(corrVX - corrAnaX)/norm(corrAnaX);
            errCorrVY = norm(corrVY - corrAnaY)/norm(corrAnaY);

            corrWX = corrW(IDLineX);
            corrWX = corrWX(reOrderX);
            corrWY = corrW(IDLineY);
            corrWY = corrWY(reOrderY);
            errCorrWX = norm(corrWX - corrAnaX)/norm(corrAnaX);
            errCorrWY = norm(corrWY - corrAnaY)/norm(corrAnaY);

            fprintf('\nRelative error to the expected correlation function using the standard Shinozuka method');
            fprintf('\nRelative error on the entire mesh = %f',errCorrV);
            fprintf('\nRelative error on the central x axis = %f',errCorrVX);
            fprintf('\nRelative error on the central y axis = %f\n',errCorrVY);

            fprintf('\nRelative error to the expected correlation function using the randomized Shinozuka method');
            fprintf('\nRelative error on the entire mesh = %f',errCorrW);
            fprintf('\nRelative error on the central x axis = %f',errCorrWX);
            fprintf('\nRelative error on the central y axis = %f\n',errCorrWY);

            fprintf(fid,'\nRelative error to the expected correlation function using the standard Shinozuka method');
            fprintf(fid,'\nRelative error on the entire mesh = %f',errCorrV);
            fprintf(fid,'\nRelative error on the central x axis = %f',errCorrVX);
            fprintf(fid,'\nRelative error on the central y axis = %f\n',errCorrVY);

            fprintf(fid,'\nRelative error to the expected correlation function using the randomized Shinozuka method');
            fprintf(fid,'\nRelative error on the entire mesh = %f',errCorrW);
            fprintf(fid,'\nRelative error on the central x axis = %f',errCorrWX);
            fprintf(fid,'\nRelative error on the central y axis = %f\n',errCorrWY);

            % Plotting

            figure('Name','Autocorrelation function with the standard Shinozuka method')
            plot(corrV,S);
            xlabel('$x$ [m]')
            ylabel('$y$ [m]')
            axis on
            colorbar
            set(gca,'FontSize',fontsize)
            mysaveas(pathnamei,'Autocorrelation_Shinozuka_std',formats,renderer);

            figure('Name','Autocorrelation function with the randomized Shinozuka method')
            plot(corrW,S);
            xlabel('$x$ [m]')
            ylabel('$y$ [m]')
            axis on
            colorbar
            set(gca,'FontSize',fontsize)
            mysaveas(pathnamei,'Autocorrelation_Shinozuka_rand',formats,renderer);

            figure('Name','Autocorrelation on central x axis using the standard Shinozuka method')
            hold on
            plot(XMid,corrVX,'b');
            plot(XMid,corrAnaX,'r--');
            xlabel('$x$ [m]')
            ylabel('$R(x - x_c)$')
            grid on
            grid minor
            % legend('Numerical computation','Analytical computation')
            set(gca,'FontSize',fontsize)
            mysaveas(pathnamei,'autocorrelation_X_Shinozuka_std',formats,renderer);

            figure('Name','Autocorrelation on central y axis using the standard Shinozuka method')
            hold on
            plot(YMid,corrVY,'b');
            plot(YMid,corrAnaY,'r--');
            xlabel('$y$ [m]')
            ylabel('$R(y - y_c)$')
            grid on
            grid minor
            % legend('Numerical computation','Analytical computation')
            set(gca,'FontSize',fontsize)
            mysaveas(pathnamei,'autocorrelation_Y_Shinozuka_std',formats,renderer);

            figure('Name','Autocorrelation on central x axis using the randomized Shinozuka method')
            hold on
            plot(XMid,corrWX,'b');
            plot(XMid,corrAnaX,'r--');
            xlabel('$x$ [m]')
            ylabel('$R(x - x_c)$')
            grid on
            grid minor
            % legend('Numerical computation','Analytical computation')
            set(gca,'FontSize',fontsize)
            mysaveas(pathnamei,'autocorrelation_X_Shinozuka_rand',formats,renderer);

            figure('Name','Autocorrelation on central y axis using the randomized Shinozuka method')
            hold on
            plot(YMid,corrWY,'b');
            plot(YMid,corrAnaY,'r--');
            xlabel('$y$ [m]')
            ylabel('$R(y - y_c)$')
            grid on
            grid minor
            % legend('Numerical computation','Analytical computation')
            set(gca,'FontSize',fontsize)
            mysaveas(pathnamei,'autocorrelation_Y_Shinozuka_rand',formats,renderer);
        end

        %% Probability density functions
        if displayPdf
            % Pdf of the middle point
            IDLineX = find(x(:,2)==L/2);
            IDLineY = find(x(:,1)==L/2);
            IDCenter = intersect(IDLineX,IDLineY);

            nPts = 201;
            [FV,XIV] = ksdensity(V(IDCenter,:));
            [FW,XIW] = ksdensity(W(IDCenter,:));

            gauss = @(x,mu,sig) exp(-((x-mu).^2)./(2*sig.^2))/sig/sqrt(2*pi);
            gaussV = gauss(XIV,0,1);
            gaussW = gauss(XIW,0,1);

            errGaussV = norm(FV - gaussV)/norm(gaussV);
            errGaussW = norm(FW - gaussW)/norm(gaussW);

            fprintf('\nRelative error to the Gaussian of the distribution at the center point');
            fprintf('\nRelative error with the standard method = %f',errGaussV);
            fprintf('\nRelative error with the randomized method = %f\n',errGaussW);

            fprintf(fid,'\nRelative error to the Gaussian of the distribution at the center point');
            fprintf(fid,'\nRelative error with the standard method = %f',errGaussV);
            fprintf(fid,'\nRelative error with the randomized method = %f\n',errGaussW);

            % Plotting
            figure('Name','ksdensity of the center node using the standard Shinozuka method')
            hold on
            plot(XIV,FV,'b')
            plot(XIV,gaussV,'r--')
            xlabel('$V(x_c, y_c)$')
            ylabel('$f(V(x_c, y_c))$')
            grid on
            grid minor
            % legend('Numerical estimation','Analytical computation')
            set(gca,'FontSize',fontsize)
            mysaveas(pathnamei,'pdf_Shinozuka_std',formats,renderer);

            figure('Name','ksdensity of the center node using the randomized Shinozuka method')
            hold on
            plot(XIW,FW,'b')
            plot(XIW,gaussW,'r--')
            xlabel('$W(x_c, y_c)$')
            ylabel('$f(W(x_c, y_c))$')
            grid on
            grid minor
            % legend('Numerical estimation','Analytical computation')
            set(gca,'FontSize',fontsize)
            mysaveas(pathnamei,'pdf_Shinozuka_rand',formats,renderer);
        end
        fclose(fid);
        fprintf(['\n' repmat('-',1,76) '\n'])
    end
end
myparallel('stop');