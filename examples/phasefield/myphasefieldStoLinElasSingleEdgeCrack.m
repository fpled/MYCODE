%
% clc
clearvars
close all

%% Inputs
setProblem = true;
solveProblem = true;
printInfos = true;
displayModel = false;
displaySolution = false;
makeMovie = false;
saveParaview = false;
displayExamples = false;

test = false; % fine meshes
% test = true; % coarse meshes

% Parameters for the Monte Carlo estimation and the parallelization
% NMC = 500; % number of independent realizations for Gaussian random fields, number of Monte-carlo simulations
NMC = 50;
% NMC = 3;
nbEx = 10; % maximum number of realizations to store as examples to display

maxNumCompThreads('automatic'); % set the number of maximum computational threads to the number of computational cores
nbCoresMax = maxNumCompThreads; % now, return the number of maximum computational threads, ie the number of computational cores
% For safety reasons, on large-scale computing stations some cores are not used
nbCoresMax = nbCoresMax - 2*(nbCoresMax>16) - 2*(nbCoresMax>32);
nbIterations = ceil(NMC/nbCoresMax);
nbWorkers = ceil(NMC/nbIterations); % number of workers in the parallel pool
% for 500 realizations :
% | available cores | 96 | 72 | 64 | 48 | 32 |
% |    nbWorkers    | 84 | 63 | 56 | 42 | 30 |

% For computation time estimation
% NMC = 1;
% maxNumCompThreads(1); % Mono-thread computation
% nbWorkers = 1;

% Domain parameters
Dim = 2; % space dimension (2 only for now)
L = 1e-3; % [m] square domain size
a = L/2; % [m] initial crack size
e = 1*(Dim==2) + 0.1e-3*(Dim==3); % [m] domain thickness

% Mesh parameters
switch Dim
    case 2
        elemtype = 'TRI3';
        %         elemtype = 'QUA4';
        %         elemtype = 'TRI6';
    case 3

end

% Characteristic lengths of the elements
% clC = 5e-6*(Dim==2) + 7.5e-6*(Dim==3); % close to the crack
% clD = 5e-6*(Dim==2) + 7.5e-6*(Dim==3); % in the rest of the domain
clC = 2.5e-6*(Dim==2) + 7.5e-6*(Dim==3);
clD = 2.5e-6*(Dim==2) + 7.5e-6*(Dim==3);
if test
    clC = 1e-5*(Dim==2) + 2e-5*(Dim==3);
    clD = 1e-5*(Dim==2) + 2e-5*(Dim==3);
    clC = 1e-4*(Dim==2) + 2e-4*(Dim==3);
    clD = 1e-4*(Dim==2) + 2e-4*(Dim==3);
end

% Mechanical loading
loading = 'Shear';
% loading = 'Tension';

% Mean material properties
% Mean values of Lame coefficients
lambda = 121.15e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
mu = 80.77e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
option = 'DEFO'; % plane strain [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
% option = 'CONT'; % plane stress [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
RHO = 1; % density
symmetry = 'Isotropic'; % symmetry class of random elasticity matrix
% symmetry = 'Anisotropic';

% Phase field parameters
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
PFmodel = 'AnisotropicMiehe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'
g = @(d) (1-d).^2; % energetic degradation function
k = 1e-10; % small artificial residual stiffness
coeff_gc = 1.0;
coeff_l = 1.0;

% Parameters for the Shinozuka method
% nu = 28; % one-dimensional order (number of terms in each spatial dimension) of the spectral representation
% Let the order be computed automatically to avoid periodicity in the
% domain

% lcorr = L/50; % [m] correlation length(s). Set to Inf for homogeneous parameters
% lcorr = 2e-5; % [m] correlation length(s)
lcorr = 1e-5;
% lcorr = 5e-6;
% lcorr = Inf; % no correlation structure

% Statistical fluctuation parameters
deltaC = 0.1; % coefficient of variation of the material (for bulk modulus for an isotropic material)
deltaP1 = 0; % coefficient of variation of fracture toughness
deltaP2 = 0; % coefficient of variation of regularization length
rhoPF = 0; % correlation coefficient between fracture toughness 'gc' and regularization length 'l'


% Display settings
fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Check validity of inputs
% to get rid of 'otherwise' in the 'switch' instructions
fileName = mfilename(); % name of the current file
% validDims = [2 3];
validDims = 2;
mustBeMember(Dim,validDims);

switch Dim
    case 2
        validOptions = {'DEFO','CONT'};
        option = validatestring(option,validOptions,fileName,'option');
        validElemTypes = {'TRI3','TRI6','QUA4'};
    case 3
        validElemTypes = {'TET4','TET10','CUB8'};
end
elemtype = validatestring(elemtype,validElemTypes,fileName,'elemtype');

% validSymmetries = {'Isotropic','Anisotropic'};
validSymmetries = {'Isotropic'};
symmetry = validatestring(symmetry,validSymmetries,fileName,'symmetry');

validLoadings = {'Shear','Tension'};
loading = validatestring(loading,validLoadings,fileName,'loading');

validPFmodels = {'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'};
PFmodel = validatestring(PFmodel,validPFmodels,fileName,'PFmodel');

if rhoPF~=0 % statistical correlation between phase field parameters
    mustBePositive(deltaP1)
    mustBePositive(deltaP2)
end

%% Computation of dependent parameters
% Phase field parameters
[gc,l] = paramMatPhaseSingleEdgeCrack(symmetry,test); % mean phase field parameters
gc = coeff_gc*gc;
l = coeff_l*l;

% Number of independant Gaussian random variables/fields to generate for the random elasticity matrix
nbGermsElas = 0;
if (deltaC~=0)
    switch lower(symmetry)
        case 'isotropic'
            nbGermsElas = 2;
        case 'anisotropic'
            n = Dim*(Dim+1)/2; % size of elasticity matrix (n2D = 3, n3D = 6)
            nbGermsElas = n*(n+1)/2;
    end
end

% Number of independent Gaussian variables/fields to generate for the phase field parameters
nbGermsPF = (deltaP1~=0) + (deltaP2~=0);
nbGermsTot = nbGermsPF + nbGermsElas; % total number of independent Gaussian variables/fields to generate

% Replace the dot by a comma when printing (to avoid file extension issues)
clDString = strrep(num2str(clD),'.',',');
clCString = strrep(num2str(clC),'.',',');
gcString = strrep(num2str(gc),'.',',');
lString = strrep(num2str(l),'.',',');
lcorrString = strrep(num2str(lcorr),'.',',');
CVElasString = strrep(num2str(deltaC),'.',',');
CVPFString = strrep([num2str(deltaP1) '_' num2str(deltaP2)],'.',',');
rhoPFString = strrep(num2str(rhoPF),'.',',');

%% Results directory
pathname = fullfile(getfemobjectoptions('path'),'MYCODE','results');
if test
    pathname = fullfile(pathname,'PFShinozuka_test');
else
    pathname = fullfile(pathname,'PFShinozuka');
end
pathname = fullfile(pathname,'StoLinElasSingleEdgeCrack');

dirName = [num2str(Dim) 'D'];
if Dim==2, dirName = [dirName option]; end
dirName = [dirName symmetry elemtype];
subDirName = [loading PFmodel];

pbPathname = fullfile(pathname,dirName); % pathname for the common problem file, containing the same mesh
pathname = fullfile(pathname,dirName,subDirName); % pathname for the solution files and folders

figPathname = fullfile(pathname,'Figures');
visualPathname = fullfile(pathname,'Visualization');

if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Results files
% File storing the common problem : the time scheme and the initial model
% ie the domain, the initial mesh, material, boundary condition and phase field
pbFileName = ['pb_clD' clDString '_clC' clCString]; % problem file

% Appendix to add to the files' name
fileAppend = ['_clD' clDString '_clC' clCString '_gc' gcString '_l' lString]; % appendix to add to files
if lcorr<Inf, fileAppend = [fileAppend '_lcorr' lcorrString]; end
if deltaC~=0, fileAppend = [fileAppend '_CVElas' CVElasString]; end
if (deltaP1~=0)||(deltaP2~=0), fileAppend = [fileAppend '_CVPF' CVPFString]; end
if rhoPF~=0, fileAppend = [fileAppend, '_rhoPF' rhoPFString]; end

exFileAppend = [fileAppend '_real']; % appendix of examples/realizations files
solFileAppend = fileAppend; % appendix of results/solutions files
if NMC>1, solFileAppend = [solFileAppend '_' num2str(NMC) 'Samples']; end
solFileName = ['sol' solFileAppend]; % solution file

if printInfos
    fprintf(['\npathname           = ' pathname])
    fprintf(['\nproblem file name  = ' pbFileName])
    fprintf(['\nsolution file name = ' solFileName])
    % fprintf(['\nnb iterations      = ' num2str(nbIterations)])
    fprintf('\n')
end

%% Set problem
if setProblem
    %% Time scheme and snapshots indices
    [T,idSnap] = setTimeSchemeSingleEdgeCrack(Dim,symmetry,loading,test);

    %% Construction of domain
    switch Dim
        case 2
            D = DOMAIN(2,[0.0,0.0],[L,L]);
            C = LIGNE([0.0,L/2],[a,L/2]);
        case 3
            D = DOMAIN(3,[0.0,0.0,0.0],[L,L,e]);
            C = QUADRANGLE([0.0,L/2,0.0],[a,L/2,0.0],[a,L/2,e],[0.0,L/2,e]);
    end

    %% Construction of (initial) mesh
    % Phase field model
    switch upper(elemtype)
        case {'TRI3','TRI6'}
            S_phaseInit = gmshdomainwithedgecrack(D,C,clD,clC,fullfile(pbPathname,'gmsh_domain_single_edge_crack'));
        case 'QUA4'
            S_phaseInit = gmshdomainwithedgecrack(D,C,clD,clC,fullfile(pbPathname,'gmsh_domain_single_edge_crack'),2,'recombine');
        case {'TET4','TET10'}
        case 'CUB8'
    end
    SInit = S_phaseInit; % displacement field model

    % Declaration of material class to get the correct number of dof
    mat_phaseInit = FOUR_ISOT;
    mat_phaseInit = setnumber(mat_phaseInit,1);
    S_phaseInit = setmaterial(S_phaseInit,mat_phaseInit);

    matInit = ELAS_ISOT;
    matInit = setnumber(matInit,1);
    SInit = setmaterial(SInit,matInit);

    %% Dirichlet boundary conditions
    % Phase field
    S_phaseInit = final(S_phaseInit,'duplicate');
    S_phaseInit = addcl(S_phaseInit,C,'T',1);

    % Displacement field
    BU = LIGNE([0.0,L],[L,L]);
    BL = LIGNE([0.0,0.0],[L,0.0]);
    BRight = LIGNE([L,0.0],[L,L]);
    BLeft = LIGNE([0.0,0.0],[0.0,L]);
    BFront = [];
    BBack = [];

    SInit = final(SInit,'duplicate');
    ud = 0; % initial imposed displacement
    switch lower(loading)
        case 'tension'
            SInit = addcl(SInit,BU,{'UX','UY'},[0;ud]);
            SInit = addcl(SInit,BL,'UY');
        case 'shear'
            SInit = addcl(SInit,BU,{'UX','UY'},[ud;0]);
            SInit = addcl(SInit,BLeft,'UY');
            SInit = addcl(SInit,BRight,'UY');
            SInit = addcl(SInit,BL);
    end

%% Save problem
save(fullfile(pbPathname, [pbFileName '.mat']),'T','idSnap','S_phaseInit','SInit',...
    'D','C','BU','BL','BRight','BLeft','BFront','BBack');
else
    load(fullfile(pbPathname, [pbFileName '.mat']),'T','idSnap','S_phaseInit','SInit',...
    'D','C','BU','BL','BRight','BLeft','BFront','BBack');
end

% Mesh parameters, also used to reshape random fields
nbElemTot = getnbelem(SInit);
elem = SInit.groupelem{1};
% gauss = calc_gauss(elem,'mass');
gauss = calc_gauss(elem,'rigi');
nbGauss = gauss.nbgauss; % number of Gauss points per element

%% Solve Problem
if solveProblem
    %% Mean parameters
    % Phase field parameters
    %     [gc,l] = paramMatPhaseSingleEdgeCrack(symmetry,test); % mean phase field parameters
    %     gc = coeff_gc*gc;
    %     l = coeff_l*l;
    mat_phaseInit = FOUR_ISOT('k',gc*l,'r',gc/l);
    mat_phaseInit = setnumber(mat_phaseInit,1);
    S_phaseInit = setmaterial(S_phaseInit,mat_phaseInit);

    % Elasticity properties
    if Dim == 2, SInit = setoption(SInit,option); end
    switch lower(symmetry)
        case 'isotropic'
            [E,NU] = paramMatIsot(Dim,lambda,mu,'option',option); % mean elasticity properties
            matInit = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e);
        case 'anisotropic'
            matInit = ELAS_ANISOT();
    end

    %% Update material properties to take the initial crack into account
    d = calc_init_dirichlet(S_phaseInit);

    switch lower(symmetry)
        case 'isotropic'
            matInit = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel);
        case 'anisotropic'
            error('Not implemented yet')
    end
    matInit = setnumber(matInit,1);
    SInit = setmaterial(SInit,matInit);

    %% Compute the order for the spectral representation such that there is no periodicity in the domain
    if lcorr<Inf
        if test
            nu = 28;
        else
            % Node coordinates
            node = getnode(SInit);
            x = getcoord(node);
            domainExtent = max(x) - min(x);
            nu = ceil(max(domainExtent./lcorr')); % one-dimensional order nu
            % such that domainExtent(j) <= period(j)/2 = nu*lcorr(j) for all spatial dimensions j=1,...,dim
            nu = 2*floor((nu+1)/2); % ensure one-dimensional order nu is even
        end
    end

    %% Initialize statistical means and second-order moments
    sz_d = getnbddl(S_phaseInit);
    sz_u = getnbddl(SInit);
    dt_mean = zeros(sz_d,length(T));
    dt_moment2 = zeros(sz_d,length(T));
    ut_mean = zeros(sz_u,length(T));
    ut_moment2 = zeros(sz_u,length(T));

    % Initialize force-displacement curves
    f_sample = zeros(length(T),NMC);

    % Initialize array of examples to store
    dt_ex = zeros(sz_d,length(idSnap),nbEx);

    %% Monte-Carlo estimation
    delete(gcp('nocreate')); % delete the currently running parallel pool if it exists

    % Parallel computation
    runSequential = false;
    myparallel('start',nbWorkers); % start the pool with the correct number of workers
    parPoolObj = gcp();
    fprintf(['\nComputing ' num2str(NMC) ' parallel Monte Carlo iterations : '])
    progressMonitoring = ParforProgressbar(NMC,'title', ['Computing ' num2str(NMC) ' parallel Monte Carlo iterations : ']); % object for progress monitoring the 'parfor' loop
    timerMC = tic;
    ticBytes(parPoolObj);
    parfor kMC = 1:NMC % Monte-Carlo iterations

        %     % Sequential computation
        %     runSequential = true;
        %     fprintf(['\nComputing ' num2str(NMC) ' sequential Monte Carlo iterations : '])
        %     progressbar(['Computing ' num2str(NMC) ' sequential Monte Carlo iterations : '])
        %     timerMC = tic;
        %     for kMC = 1:NMC % Monte-Carlo iterations

        S_phase = S_phaseInit;
        S = SInit;

        %         % Recovery of homogeneous mean phase field parameters (necessary as
        %         % fracture toughness may vary from one computation to another)
        %         mats_phase = MATERIALS(S_phase);
        %         k = getparam(mats_phase{1},'k');
        %         r = getparam(mats_phase{1},'r');
        %         gc = sqrt(k*r);
        %         l = sqrt(k/r);

        %% Generation of random parameter/property fields
        substream = RandStream.create('mlfg6331_64','NumStreams',NMC,'StreamIndices',kMC);

        if lcorr==Inf % all parameters are homogeneous
            V = randn(substream,1,nbGermsTot);
        else % there are heterogeneous parameters
            %             xgauss = calc_gausscoord(S,'mass'); % gauss points coordinates
            xgauss = calc_gausscoord(S,'rigi'); % gauss points coordinates
            V = shinozukaSample(substream,xgauss,lcorr,nbGermsTot,'order',nu);
        end

        % Phase field parameters
        if (deltaP1==0)&&(deltaP2==0) % deterministic phase field parameters
            % Parameters remain unchanged. Nothing to do.

        else % random phase field parameters (homogeneous XOR heterogeneous)
            % Shape and scale parameters of the gamma distributions followed by phase field parameters
            [aP1,bP1,aP2,bP2] = hyperMatPhaseSingleEdgeCrack(gc,l,deltaP1,deltaP2);
            gc_k = gc; % fracture toughness sample [N/m]
            l_k = l; % regularization length sample [m]
            if deltaP2==0 % deterministic regularization length but random fracture toughness
                gc_k = gaminv(normcdf(V(:,1)),aP1,bP1);
            elseif deltaP1==0 % deterministic fracture toughness but random regularization length
                l_k = gaminv(normcdf(V(:,1)),aP2,bP2);
            else % random phase field parameters with a possible statistical correlation
                gc_k = gaminv(normcdf(V(:,1)),aP1,bP1);
                l_k = gaminv(normcdf(rhoPF*V(:,1) + sqrt(1-rhoPF^2)*V(:,2)),aP2,bP2);
            end

            if lcorr<Inf % heterogeneous properties
                if deltaP1~=0
                    gc_k = reshape(gc_k,1,1,nbElemTot,nbGauss);
                    gc_k = MYDOUBLEND(gc_k);
                    gc_k = FEELEMFIELD(cell(gc_k),'type','scalar','storage','gauss');
                end
                if deltaP2~=0
                    l_k = reshape(l_k,1,1,nbElemTot,nbGauss);
                    l_k = MYDOUBLEND(l_k);
                    l_k = FEELEMFIELD(cell(l_k),'type','scalar','storage','gauss');
                end
            end

            % Update phase field parameters
            mat_phase = FOUR_ISOT('k',gc_k.*l_k,'r',gc_k./l_k);
            mat_phase = setnumber(mat_phase,1);
            S_phase = setmaterial(S_phase,mat_phase);
        end

        % Elasticity properties
        if deltaC==0 % deterministic elasticity properties
            % Properties remain unchanged. Nothing to do.

        else % random elasticity properties
            % Initialize the material to prevent the 'Uninitialized Temporaries' warning display by Matlab
            mat = MATERIAL();
            switch lower(symmetry)
                case 'isotropic'
                    %                     [E,NU] = paramMatIsot(Dim,lambda,mu,'option',option);
                    % Shape and scale parameters of the independent gamma distributions followed by the bulk and the shear modulus
                    [aC1,bC1,aC2,bC2] = hyperMatIsot(Dim,lambda,mu,deltaC,'option',option);
                    C1_k = gaminv(normcdf(V(:,nbGermsPF + 1)),aC1,bC1); % samples for bulk modulus [Pa]
                    C2_k = gaminv(normcdf(V(:,nbGermsPF + 2)),aC2,bC2); % samples for shear modulus [Pa]
                    E_k = (9*C1_k.*C2_k)./(3*C1_k+C2_k); % [Pa]
                    NU_k = (3*C1_k-2*C2_k)./(6*C1_k+2*C2_k);

                    if lcorr<Inf % heterogeneous properties
                        E_k = reshape(E_k,1,1,nbElemTot,nbGauss);
                        E_k = MYDOUBLEND(E_k);
                        E_k = FEELEMFIELD(cell(E_k),'type','scalar','storage','gauss');
                        NU_k = reshape(NU_k,1,1,nbElemTot,nbGauss);
                        NU_k = MYDOUBLEND(NU_k);
                        NU_k = FEELEMFIELD(cell(NU_k),'type','scalar','storage','gauss');
                    end

                    mat = ELAS_ISOT('E',E_k,'NU',NU_k,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel);
                    % S = setoption(S,option);

                case 'anisotropic'
                    error('Not implemented yet')
            end
            % Update material properties
            mat = setnumber(mat,1);
            S = setmaterial(S,mat);
        end

        %% Solve deterministic problem
        [dt,ut,ft] = solvePFDetLinElasSingleEdgeCrack(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading);

        %% Compute second-order statistics
        dt_val = getvalue(dt);
        dt_mean = dt_mean + dt_val/NMC;
        dt_moment2 = dt_moment2 + dt_val.^2/NMC;
        ut_val = getvalue(ut);
        ut_mean = ut_mean + ut_val/NMC;
        ut_moment2 = ut_moment2 + ut_val.^2/NMC;

        %% Store force-displacement curve and example realizations
        f_sample(:,kMC) = ft;
        dt_ex(:,:,kMC) = dt_val(:,idSnap);

        %% Increment progress
        if runSequential % sequential loop case
            progressbar(kMC/NMC);
        else % parallel loop case
            progressMonitoring.increment();
        end
    end
    fprintf('done\n')
    if ~isempty(gcp('nocreate'))
        tocBytes(parPoolObj);
        delete(progressMonitoring); % delete the progress handle when the parfor loop is done
    end
    timeMC = toc(timerMC);
    myparallel('stop')

    %% Post-treatment
    % Compute unbiased variances
    if NMC>1
        dt_var = (NMC/(NMC-1))*(dt_moment2 - dt_mean.^2);
        ut_var = (NMC/(NMC-1))*(ut_moment2 - ut_mean.^2);
    else
        dt_var = zeros(sz_d,length(T));
        ut_var = zeros(sz_u,length(T));
    end

    % Statistical outputs of solution
    probs = [0.025 0.975];

    ft_mean = mean(f_sample,2);
    ft_std = std(f_sample,0,2);
    ft_ci = quantile(f_sample,probs,2);

    fmax = max(f_sample);
    fmax_mean = mean(fmax);
    fmax_std = std(fmax);
    fmax_ci = quantile(fmax,probs);

    npts = 100;
    [fmax_f,fmax_xi,fmax_bw] = ksdensity(fmax,'npoints',npts);

    %% Save solution
    save(fullfile(pathname, [solFileName '.mat']),'timeMC','nbWorkers','NMC',...
        'nu','runSequential','dt_ex',...
        'dt_mean','ut_mean','dt_var','ut_var','ft_mean','ft_std','ft_ci',...
        'fmax','fmax_mean','fmax_std','fmax_ci','probs','fmax_f','fmax_xi','fmax_bw');
else
    load(fullfile(pathname, [solFileName '.mat']),'timeMC','nbWorkers','NMC',...
        'nu','runSequential','dt_ex',...
        'dt_mean','ut_mean','dt_var','ut_var','ft_mean','ft_std','ft_ci',...
        'fmax','fmax_mean','fmax_std','fmax_ci','probs','fmax_f','fmax_xi','fmax_bw');
end

%% Outputs
if printInfos
    % Phase field element model
    fprintf('\nPhase field model :')
    fprintf('\ndim                  = %d',Dim)
    fprintf(['\nelement type         = ', elemtype])
    fprintf('\nloading              = %s',loading)
    fprintf('\nmat sym              = %s',symmetry)
    if strcmpi(symmetry,'anisotropic')
        fprintf('\nangle                = %g deg',ang)
    end
    fprintf('\nPF model             = %s',PFmodel);
    fprintf('\nnb elements          = %g',getnbelem(SInit))
    fprintf('\nnb nodes             = %g',getnbnode(SInit))
    fprintf('\nnb displacement dofs = %g',getnbddl(SInit))
    fprintf('\nnb time dofs         = %g',getnbtimedof(T))
    fprintf('\n')

    % Computation infos
    if runSequential
        fprintf('\nSequential computation :')
    else
        fprintf('\nParallel computation :')
        fprintf('\nnb workers       = %g',nbWorkers)
    end
    fprintf('\nnb samples       = %g',NMC)
    fprintf('\ncomputation time = %f s',timeMC)
    fprintf('\n')

    % Gaussian fields
    order = nu^Dim; % d-dimensional order of the spectral representation for all spatial dimensions
    xgauss = calc_gausscoord(SInit,'mass'); % gauss points coordinates
    nx = size(xgauss,1); % number of points
    nV = nbGermsTot*NMC; % number of independent realizations for all Gaussian random fields
    fprintf('\nGaussian fields :')
    fprintf('\nnumber of points  = %d',nx)
    fprintf('\nnumber of fields  = %d',nbGermsTot)
    fprintf('\nnumber of samples =%s for each Gaussian random field',sprintf(' %g',NMC))
    fprintf('\nnumber of samples =%s for all Gaussian random fields',sprintf(' %g',nV))
    fprintf('\nnumber of terms   =%s in the spectral representation',sprintf(' %g',order))
    fprintf('\n')

    % Maximum reaction force
    fprintf('\nMaximum reaction force :')
    if Dim==2
        fprintf('\nmean(fmax)    = %g kN/mm',fmax_mean*1e-6)
        fprintf('\nstd(fmax)     = %g kN/mm',fmax_std*1e-6)
    elseif Dim==3
        fprintf('\nmean(fmax)    = %g kN',fmax_mean*1e-3)
        fprintf('\nstd(fmax)     = %g kN',fmax_std*1e-3)
    end
    fprintf('\ndisp(fmax)    = %g',fmax_std/fmax_mean)
    if Dim==2
        fprintf('\n%d%% ci(fmax)  = [%g,%g] kN/mm',(probs(2)-probs(1))*100,fmax_ci(1)*1e-6,fmax_ci(2)*1e-6)
    elseif Dim==3
        fprintf('\n%d%% ci(fmax)  = [%g,%g] kN',(probs(2)-probs(1))*100,fmax_ci(1)*1e-3,fmax_ci(2)*1e-3)
    end
    fprintf('\n')
end

if (~exist(figPathname,'dir'))&&(displayModel||displaySolution||makeMovie||displayExamples)
    mkdir(figPathname);
end

%% Display domains, boundary conditions and meshes
if displayModel
    plotDomain({D,C},'legend',false);
    mysaveas(figPathname,'domain',formats,renderer);
    %     mymatlab2tikz(figPathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(SInit,'legend',false);
    ampl = 0.5;
    v = calc_init_dirichlet(SInit);
    [hN,legN] = vectorplot(SInit,'U',v,ampl,'r','LineWidth',1);
    %     legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(figPathname,'boundary_conditions_displacement',formats,renderer);
    
    [hD_phase,legD_phase] = plotBoundaryConditions(S_phaseInit,'legend',false);
    %     legend([hD_phase,hN_phase],[legD_phase,legN_phase],'Location','NorthEastOutside')
    mysaveas(figPathname,'boundary_conditions_damage',formats,renderer);
    
    %     plotModel(SInit,'legend',false);
    %     mysaveas(figPathname,'mesh',formats,renderer);
    
    plotModel(SInit,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(figPathname,'mesh',formats,renderer);
    
    u = ut_mean(:,:,end);
    for k=idSnap
        ampl = getsize(SInit)/max(abs(u(:)))/20;
        plotModelDeflection(SInit,unfreevector(SInit,u(:,k)),'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
        mysaveas(figPathname,['mesh_deflected_sample' solFileAppend '_t' num2str(k)],formats,renderer);

        figure('Name','Meshes')
        clf
        plot(SInit,'Color','k','FaceColor','k','FaceAlpha',0.1);
        plot(SInit+ampl*unfreevector(SInit,u(:,k)),'Color','b','FaceColor','b','FaceAlpha',0.1);
        mysaveas(figPathname,['meshes_deflected' solFileAppend '_t' num2str(k)],formats,renderer);
    end
end

%% Display statistics of solutions
if displaySolution    
        %% Display force-displacement curve
        [t,~] = gettevol(T);
        figure('Name','Force-displacement')
        clf
        ciplot(ft_ci(:,1)*((Dim==2)*1e-6+(Dim==3)*1e-3),ft_ci(:,2)*((Dim==2)*1e-6+(Dim==3)*1e-3),t*1e3,'b');
        alpha(0.2)
        hold on
        plot(t*1e3,ft_mean*((Dim==2)*1e-6+(Dim==3)*1e-3),'-b','LineWidth',linewidth)
        grid on
        box on
        set(gca,'FontSize',fontsize)
        xlabel('Displacement [mm]','Interpreter',interpreter)
        ylabel('Force [kN]','Interpreter',interpreter)
        leg = legend({['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
            'mean value'},'Location','NorthWest');
        set(leg,'Interpreter','latex')
        mysaveas(figPathname,['force_displacement' solFileAppend],formats);
        % mymatlab2tikz(figPathname,['force_displacement' solFileAppend '.tex']);
    
        %% Display pdf of critical force
        figure('Name','Probability Density Estimate: Critical force')
        clf
        plot(fmax_xi*((Dim==2)*1e-6+(Dim==3)*1e-3),fmax_f,'-b','LineWidth',linewidth)
        hold on
        ind_fmax = find(fmax_xi>=fmax_ci(1) & fmax_xi<fmax_ci(2));
        area(fmax_xi(ind_fmax)*((Dim==2)*1e-6+(Dim==3)*1e-3),fmax_f(ind_fmax),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
        scatter(fmax_mean*((Dim==2)*1e-6+(Dim==3)*1e-3),0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b')
        hold off
        grid on
        box on
        set(gca,'FontSize',fontsize)
        xlabel('$f$ [kN]','Interpreter',interpreter)
        ylabel('$p_{F_c}(f)$','Interpreter',interpreter)
        leg = legend('pdf',...
            ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
            'mean value');
        set(leg,'Interpreter',interpreter)
        mysaveas(figPathname,['pdf_fmax' solFileAppend],formats,renderer);
        % mymatlab2tikz(figPathname,['pdf_fmax' solFileAppend '.tex']);
    
    %% Display means and variances of solutions at different instants
    for j=1:length(idSnap)
        dj = dt_mean(:,idSnap(j));
        dj_var = dt_var(:,idSnap(j));

        plotSolution(S_phaseInit,dj);
        set(gcf,'Name','d_mean')
        mysaveas(figPathname,['damage_mean' solFileAppend '_t' num2str(idSnap(j))],formats,renderer);
        plotSolution(S_phaseInit,dj_var);
        set(gcf,'Name','d_variance')
        mysaveas(figPathname,['damage_var' solFileAppend '_t' num2str(idSnap(j))],formats,renderer);

        %         ampl = 0;
        %         uj = ut_mean(:,idSnap(j));
        %         uj_var = ut_var(:,idSnap(j));
        %         for i=1:Dim
        %             plotSolution(SInit,uj,'displ',i,'ampl',ampl);
        %             set(gcf,'Name',['u_' num2str(i) '_mean'])
        %             mysaveas(figPathname,['displacement' num2str(i) '_mean' solFileAppend '_t' num2str(idSnap(j))],formats,renderer);
        %             plotSolution(SInit,uj_var,'displ_var',i,'ampl',ampl);
        %             set(gcf,'Name',['u_' num2str(i) '_variance'])
        %             mysaveas(figPathname,['displacement' num2str(i) '_var'  solFileAppend '_t' num2str(idSnap(j))],formats,renderer);
        %         end
        %
        %         for i=1:(Dim*(Dim+1)/2)
        %             plotSolution(SInit,uj,'epsilon',i,'ampl',ampl);
        %             mysaveas(figPathname,['epsilon' num2str(i) '_mean'  solFileAppend '_t' num2str(idSnap(j))],formats,renderer);
        %
        %             % Stress and energy cannot be computed since the material
        %             % properties are not saved
        %             % plotSolution(SInit,uj,'sigma',i,'ampl',ampl);
        %             % mysaveas(figPathname,['sigma' num2str(i) '_mean' solFileAppend '_t' num2str(idSnap(j))],formats,renderer);
        %         end
        %
        %         plotSolution(SInit,uj,'epsilon','mises','ampl',ampl);
        %         mysaveas(figPathname,['epsilonVonMises_mean' solFileAppend '_t' num2str(idSnap(j))],formats,renderer);
        %
        %         % Stress and energy cannot be computed since the material
        %         % properties are not saved
        %         % plotSolution(SInit,uj,'sigma','mises','ampl',ampl);
        %         % mysaveas(figPathname,['sigma_von_mises_mean' solFileAppend '_t' num2str(idSnap(j))],formats,renderer);
        %
        %         % plotSolution(SInit,uj,'energyint','','ampl',ampl);
        %         % mysaveas(figPathname,['internal_energy_mean' solFileAppend '_t' num2str(idSnap(j))],formats,renderer);
    end
    
end

%% Display evolution of means and variances of solutions
if makeMovie
    options = {'plotiter',true,'plottime',false};
    framerate = 80;
    
    sz_d = [getnbddl(S_phaseInit),getnbtimedof(T)];
    dk = TIMEMATRIX(reshape(dt_mean(:,:),sz_d),T);
    dk_var = TIMEMATRIX(reshape(dt_var(:,:),sz_d),T);

    evolSolution(S_phaseInit,dk,'FrameRate',framerate,'filename',['damage_mean' solFileAppend],'pathname',figPathname,options{:});
    evolSolution(S_phaseInit,dk_var,'FrameRate',framerate,'filename',['damage_var' solFileAppend],'pathname',figPathname,options{:});
    
    %     % ampl = getsize(SInit)/max(abs(ut_mean(:)))/20;
    %     ampl = 0;
    %     uk = TIMEMATRIX(reshape(ut_mean(:,:),sz_u),T);
    %     uk_var = TIMEMATRIX(reshape(ut_var(:,:),sz_u),T);
    %     sz_u = [getnbddl(SInit),getnbtimedof(T)];
    %     for i=1:Dim
    %         evolSolution(SInit,uk,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement' num2str(i) '_mean' solFileAppend],'pathname',figPathname,options{:});
    %         evolSolution(SInit,uk_var,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement' num2str(i) '_var' solFileAppend],'pathname',figPathname,options{:});
    %     end
    %
    %     % Stress and energy cannot be computed since the material
    %     % properties are not saved
    %     for i=1:(Dim*(Dim+1)/2)
    %         evolSolution(SInit,uk,'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['epsilon' num2str(i) '_mean' solFileAppend],'pathname',figPathname,options{:});
    %         % evolSolution(SInit,uk,'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['sigma' num2str(i) 'mean' solFileAppend],'pathname',figPathname,options{:});
    %     end
    %
    %     evolSolution(SInit,uk,'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename',['epsilonVonMises_mean' solFileAppend],'pathname',figPathname,options{:});
    %     % evolSolution(SInit,uk,'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename',['sigmaVonMises_mean' solFileAppend],'pathname',figPathname,options{:});
    %     % evolSolution(SInit,uk,'energyint','','ampl',ampl,'FrameRate',framerate,'filename',['internal_energy_mean' solFileAppend],'pathname',figPathname,options{:});
end

%% Save means and variances of solutions
if saveParaview
    if ~exist(visualPathname,'dir')
        mkdir(visualPathname);
    end
    [~,rep] = gettevol(T);
    for i=1:length(T)
        di = dt_mean(:,rep(i))';
        ui = ut_mean(:,rep(i))';
        dvi = dt_var(:,rep(i))';
        uvi = ut_var(:,rep(i))';
        
        write_vtk_mesh(SInit,{di,ui,dvi,uvi},[],...
            {'damage_mean','displacement_mean','damage_variance','displacement_variance'},[],...
            visualPathname,'solution',1,i-1);
    end
    make_pvd_file(visualPathname,'solution',1,length(T));
end

%% Display some realizations as examples
if displayExamples
    for i=1:min(nbEx,NMC)
        for j=1:length(idSnap)
            dj = dt_ex(:,j,i);
            plotSolution(S_phaseInit,dj);
            set(gcf,'Name',['d_real' num2str(i) '_t' num2str(idSnap(j))])
            mysaveas(figPathname,['damage' exFileAppend num2str(i) '_t' num2str(idSnap(j))],formats,renderer);
        end
    end
end