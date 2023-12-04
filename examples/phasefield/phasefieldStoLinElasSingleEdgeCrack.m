%% Phase field fracture model - stochastic linear elasticity problem with single edge crack %%
%%------------------------------------------------------------------------------------------%%
% [Bourdin, Francfort, Marigo, 2000, JMPS] (isotropic phase field model with no split of Bourdin et al.)
% [Miehe, Welschinger, Hofacker, 2010 IJNME] (anisotropic phase field model of Miehe et al.)
% [Miehe, Hofacker, Welschinger, 2010, CMAME] (anisotropic phase field model of Miehe et al.)
% [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME] (anisotropic phase field model of Miehe et al.)
% [Hesch, Weinberg, 2014, IJNME] (anisotropic phase field model of Miehe et al.)
% [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM] (anisotropic phase field model of Miehe et al.)
% [Ambati, Gerasimov, De Lorenzis, 2015, CM] (hybrid isotropic-anisotropic phase field model of Ambati et al. with Miehe et al. decomposition compared with the isotropic one of Bourdin et al. and the anisotropic ones of Amor et al. and Miehe et al.)
% [Liu, Li, Msekh, Zuo, 2016, CMS] (anisotropic phase field model of Miehe et al.)
% [Zhou, Rabczuk, Zhuang, 2018, AES] (anisotropic phase field model of Miehe et al.)
% [Wu, Nguyen, 2018, JMPS] (hybrid isotropic-anisotropic phase field model of Wu et al.)
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM] (anisotropic phase field model of Wu et al.)
% [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US] (anisotropic phase field model of Amor et al.)
% [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME] (anisotropic phase field model of Miehe et al.)
% [Wu, Nguyen, Zhou, Huang, 2020, CMAME] (anisotropic phase field model of Wu et al.)
% [Kristensen, Martinez-Paneda, 2020, TAFM] (hybrid isotropic-anisotropic phase field model of Ambati et al. with Amor et al. decomposition)
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME] (anisotropic phase field model of Nguyen et al.)
% [Hu, Guilleminot, Dolbow, 2020, CMAME] (anisotropic phase field model of Hu et al.)
% [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME] (anisotropic phase field model of Miehe et al.)

% clc
clearvars
close all
% rng('default');

%% Input data
setProblem = true;
solveProblem = true;
displayModel = false;
displaySolution = false;
makeMovie = false;
saveParaview = false;

test = true; % coarse mesh
% test = false; % fine mesh

numWorkers = 4;
% numWorkers = 1; maxNumCompThreads(1); % mono-thread computation

% Deterministic model parameters
Dim = 2; % space dimension Dim = 2, 3
symmetry = 'Isot'; % 'Isot', 'MeanIsot', 'Anisot'. Material symmetry
ang = 45; % clockwise material orientation angle around z-axis for anisotopic material [deg]
loading = 'Shear'; % 'Tension' or 'Shear'
PFmodel = 'Miehe'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 1; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
critConv = 'Energy'; % 'Solution', 'Residual', 'Energy'
initialCrack = 'GeometricCrack'; % 'GeometricCrack', 'GeometricNotch', 'InitialPhaseField'
FEmesh = 'Optim'; % 'Unif' or 'Optim'

% Random model parameters
% N = 500; % number of samples
N = numWorkers;
randMat = struct('delta',0.2,'lcorr',1e-4); % random material parameters model
aGc = 0;
bGc = 0;
% gc = 2.7e3;
% aGc = 0.6*gc;
% bGc = 1.4*gc;
% aGc = [0.7,1.2]*gc;
% bGc = [0.8,1.3]*gc;
randPF = struct('aGc',aGc,'bGc',bGc,'lcorr',Inf); % random phase field parameters model

suffix = '';

foldername = ['singleEdgeCrack' loading '_' num2str(Dim) 'D'];
filename = ['linElas' symmetry];
if strcmpi(symmetry,'anisot') % anisotropic material
    filename = [filename num2str(ang) 'deg'];
end
filename = [filename PFmodel PFsplit PFregularization PFsolver initialCrack...
    'MaxIter' num2str(maxIter)];
if maxIter>1
    filename = [filename 'Tol' num2str(tolConv) num2str(critConv)];
end
filename = [filename 'Mesh' FEmesh '_' num2str(N) 'samples'];
if any(randMat.delta)
    filename = [filename '_RandMat_Delta' num2str(randMat.delta,'_%g') '_Lcorr' num2str(randMat.lcorr,'_%g')];
end
if any(randPF.aGc) && any(randPF.bGc)
    gcbounds = [randPF.aGc(:),randPF.bGc(:)]';
    filename = [filename '_RandPF_Gc' num2str(gcbounds(:)','_%g') '_Lcorr' num2str(randPF.lcorr,'_%g')];
end
filename = [filename suffix];

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','phasefieldSto',foldername,filename);
if test
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','phasefieldSto_test',foldername,filename);
end
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    L = 1e-3;
    a = L/2;
    b = L/2;
    if Dim==2
        e = 1;
        D = DOMAIN(2,[0.0,0.0],[L,L]);
        C = LIGNE([0.0,b],[a,b]);
    elseif Dim==3
        e = 0.1e-3;
        D = DOMAIN(3,[0.0,0.0,0.0],[L,L,e]);
        C = QUADRANGLE([0.0,b,0.0],[a,b,0.0],[a,b,e],[0.0,b,e]);
    end
    
    % clD = 6.25e-5; % [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME]
    % clD = 5e-5; % [Hu, Guilleminot, Dolbow, 2020, CMAME]
    % clD = 3e-5; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
    % clD = 2e-5; % [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
    % clD = 3.96e-6; % [Zhou, Rabczuk, Zhuang, 2018, AES]
    % clD = 3.9e-6; % [Hesch, Weinberg, 2014, IJNME], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % clD = 2e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
    % clD = 1e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
    
    % clC = 5e-6; % [Hu, Guilleminot, Dolbow, 2020, CMAME]
    % clC = 3.96e-6; % [Zhou, Rabczuk, Zhuang, 2018, AES]
    % clC = 3.906e-6; % [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME]
    % clC = 3.9e-6; % [Hesch, Weinberg, 2014, IJNME], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % clC = 3.75e-6; % (shear test) [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
    % clC = 2.5e-6; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
    % clC = 2e-6; % (shear test) [Miehe, Hofacker, Welschinger, 2010, CMAME]
    % clC = 1e-6; % (tension test) [Miehe, Welschinger, Hofacker, 2010 IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
    % clC = 2e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
    % clC = 1e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
    % clC = 6e-7; % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
    switch lower(FEmesh)
        case 'unif' % uniform mesh
            switch lower(symmetry)
                case 'isot' % isotropic material
                    cl = 5e-6;
                case 'anisot' % anisotropic material
                    cl = 4.25e-6;
                otherwise
                    error('Wrong material symmetry class');
            end
            if test
                if Dim==2
                    cl = 1e-5;
                elseif Dim==3
                    cl = 2e-5;
                end
            end
            clD = cl;
            clC = cl;
            B = [];
        case 'optim' % optimized mesh
            if Dim==2
                clD = 2.5e-5;
                clC = 2.5e-6;
            elseif Dim==3
                clD = 4e-5;
                clC = 5e-6;
            end
            if test
                clD = 4e-5;
                clC = 1e-5;
            end
            VIn = clC;
            VOut = clD;
            XMin = a; XMax = L;
            switch lower(loading)
                case 'tension'
                    if strcmpi(symmetry,'isot')
                        YMin = b-L/8; YMax = b+L/8;
                    else
                        YMin = 0; YMax = L;
                    end
                case 'shear'
                    if strcmpi(PFmodel,'bourdin')
                        YMin = 0; YMax = L;
                    else
                        YMin = 0; YMax = b;
                    end
                otherwise
                    error('Wrong loading case');
            end
            ZMin = 0;
            ZMax = e;
            Thickness = a;
            % Thickness = 0;
            B = struct('VIn',VIn,'VOut',VOut,'XMin',XMin,'XMax',XMax,'YMin',YMin,'YMax',YMax,'ZMin',ZMin,'ZMax',ZMax,'Thickness',Thickness);
        otherwise
            error('Wrong FE mesh')
    end
    switch lower(initialCrack)
        case 'geometriccrack'
            S_phase = gmshdomainwithedgecrack(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,'Box',B);
        case 'geometricnotch'
            c = 1e-5; % crack width
            S_phase = gmshdomainwithedgenotch(D,C,c,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,'Box',B);
        case 'initialphasefield'
            S_phase = gmshdomainwithedgecrack(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,'noduplicate','Box',B);
        otherwise
            error('Wrong model for initial crack');
    end
    S = S_phase;
    
    %% Phase field problem
    %% Material
    switch lower(symmetry)
        case {'isot','meanisot'} % almost surely or mean isotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 2.7e3; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
            % Regularization parameter (width of the smeared crack)
            % l = 3.75e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME]
            % l = 3e-5; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
            % l = 1.5e-5; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Hesch, Weinberg, 2014, IJNME], [Liu, Li, Msekh, Zuo, 2016, CMS], [Zhou, Rabczuk, Zhuang, 2018, AES], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM], [Hu, Guilleminot, Dolbow, 2020, CMAME]
            l = 1e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
            % l = 7.5e-6; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Hesch, Weinberg, 2014, IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Zhou, Rabczuk, Zhuang, 2018, AES], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
            % l = 5e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
            % l = 4e-6; % [Ambati, Gerasimov, De Lorenzis, 2015, CM]
            % eta = 0.052; w0 = 75.94; l = eta/sqrt(w0)*1e-3; % l = 6e-7; % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
        case 'anisot' % anisotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 1e3; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            % Regularization parameter (width of the smeared crack)
            l = 8.5e-6; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        otherwise
            error('Wrong material symmetry class');
    end
    % Small artificial residual stiffness
    % k = 1e-12;
    k = 0;
    
    % Material
    switch lower(PFregularization)
        case 'at1'
            % c0 = 8/3;
            K = 3/4*gc*l; % K = 2*(gc*l)/c0;
            R = 0;
            Qn = -3/8*gc/l; % Qn = -(gc/l)/c0;
        case 'at2'
            % c0 = 2;
            K = gc*l; % K = 2*(gc*l)/c0;
            R = gc/l; % R = 2*(gc/l)/c0;
            Qn = 0;
        otherwise
            error('Wrong regularization model');
    end
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'PFregularization',PFregularization,'aGc',randPF.aGc,'bGc',randPF.bGc,'lcorr',randPF.lcorr);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BRight = LIGNE([L,0.0],[L,L]);
    elseif Dim==3
        BRight = PLAN([L,0.0,0.0],[L,L,0.0],[L,0.0,e]);
    end
    
    findddlboundary = @(S_phase) findddl(S_phase,'T',BRight);
    
    if strcmpi(initialCrack,'geometriccrack')
        S_phase = final(S_phase,'duplicate');
    else
        S_phase = final(S_phase);
    end
    
    if strcmpi(initialCrack,'initialphasefield')
        S_phase = addcl(S_phase,C,'T',1);
    end
    
    %% Stiffness matrices and sollicitation vectors
    % a_phase = BILINFORM(1,1,K); % uniform values
    % % a_phase = DIFFUSIONFORM(K);
    % a_phase = setfree(a_phase,0);
    % K_phase = calc_matrix(a_phase,S_phase); % quadorder=0, nbgauss=1
    % % K_phase = calc_matrix(a_phase,S_phase,[],[],'quadorder',2);
    % b_phase = calc_nonhomogeneous_vector(S_phase,K_phase);
    % K_phase = freematrix(S_phase,K_phase);
    
    % r_phase = BILINFORM(0,0,R); % uniform values
    % R_phase = calc_matrix(r_phase,S_phase); % quadorder=2, nbgauss=3
    % A_phase = K_phase + R_phase;
    
    % l_phase = LINFORM(0,Qn); % uniform values
    % l_phase = setfree(l_phase,1);
    % b_phase = -b_phase + calc_vector(l_phase,S_phase);
    
    % [A_phase,b_phase] = calc_rigi(S_phase);
    % b_phase = -b_phase + bodyload(S_phase,[],'QN',Qn);
    
    %% Linear elastic displacement field problem
    %% Materials
    % Option
    option = 'DEFO'; % plane strain [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Hesch, Weinberg, 2014, IJNME], [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Zhou, Rabczuk, Zhuang, 2018, AES], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
    % option = 'CONT'; % plane stress [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
    switch lower(symmetry)
        case {'isot','meanisot'} % almost surely or mean isotropic material
            % Lame coefficients
            % lambda = 121.1538e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
            % mu = 80.7692e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
            % lambda = 121.154e9; % [Hesch, Weinberg, 2014, IJNME]
            % mu = 80.769e9; % [Hesch, Weinberg, 2014, IJNME]
            lambda = 121.15e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME], [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
            mu = 80.77e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME], [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
            if strcmpi(symmetry,'isot')
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
                    % E = 210e9; NU = 0.2; % [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
                    % E = 210e9; NU = 0.3; % [Zhou, Rabczuk, Zhuang, 2018, AES], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME], [Kristensen, Martinez-Paneda, 2020, TAFM]
                    % kappa = 121030e6; NU=0.227; lambda=3*kappa*NU/(1+NU); mu = 3*kappa*(1-2*NU)/(2*(1+NU)); E = 3*kappa*(1-2*NU); % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
                elseif Dim==3
                    E = mu*(3*lambda+2*mu)/(lambda+mu);
                    NU = lambda/(lambda+mu)/2;
                end
            elseif strcmpi(symmetry,'meanisot')
                % Elasticity matrix
                if Dim==2
                    Cmat = e*...
                        [lambda+2*mu,lambda,0;...
                        lambda,lambda+2*mu,0;...
                        0,0,mu];
                elseif Dim==3
                    Cmat = [lambda+2*mu,lambda,lambda,0,0,0;...
                        lambda,lambda+2*mu,lambda,0,0,0;...
                        lambda,lambda,lambda+2*mu,0,0,0;...
                        0,0,0,mu,0,0;...
                        0,0,0,0,mu,0;...
                        0,0,0,0,0,mu];
                end
            end
            
        case 'anisot' % anisotropic material
            if Dim==2
                switch lower(option)
                    case 'defo'
                        % [Nguyen, Yvonnet, Waldmann, He, 2020,IJNME]
                        % Elasticity matrix in reference material coordinate system [Pa]
                        Cmat = e*...
                            [65 20 0;
                            20 260 0;
                            0 0 30]*1e9;
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
                        error('Not implemented yet');
                end
            elseif Dim==3
                error('Not implemented yet');
            end
        otherwise
            error('Wrong material symmetry class');
    end
    % Energetic degradation function
    g = @(d) (1-d).^2;
    % Density
    RHO = 1;
    
    % Material
    d = calc_init_dirichlet(S_phase);
    switch lower(symmetry)
        case 'isot' % almost surely isotropic material
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit,'delta',randMat.delta,'lcorr',randMat.lcorr);
        case {'meanisot','anisot'} % mean isotropic or anisotropic material
            mat = ELAS_ANISOT('C',Cmat,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit,'delta',randMat.delta,'lcorr',randMat.lcorr);
        otherwise
            error('Wrong material symmetry class');
    end
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BU = LIGNE([0.0,L],[L,L]);
        BL = LIGNE([0.0,0.0],[L,0.0]);
        % BRight = LIGNE([L,0.0],[L,L]);
        BLeft = LIGNE([0.0,0.0],[0.0,L]);
        BFront = [];
        BBack = [];
    elseif Dim==3
        BU = PLAN([0.0,L,0.0],[L,L,0.0],[0.0,L,e]);
        BL = PLAN([0.0,0.0,0.0],[L,0.0,0.0],[0.0,0.0,e]);
        % BRight = PLAN([L,0.0,0.0],[L,L,0.0],[L,0.0,e]);
        BLeft = PLAN([0.0,0.0,0.0],[0.0,L,0.0],[0.0,0.0,e]);
        BFront = PLAN([0.0,0.0,e],[L,0.0,e],[0.0,L,e]);
        BBack = PLAN([0.0,0.0,0.0],[L,0.0,0.0],[0.0,L,0.0]);
    end
    
    addbc = @(S,ud) addbcSingleEdgeCrack(S,ud,BU,BL,BLeft,BRight,BFront,BBack,loading);
    switch lower(loading)
        case 'tension'
            findddlforce = @(S) findddl(S,'UY',BU);
        case 'shear'
            findddlforce = @(S) findddl(S,'UX',BU);
        otherwise
            error('Wrong loading case');
    end
    
    if strcmpi(initialCrack,'geometriccrack')
        S = final(S,'duplicate');
    else
        S = final(S);
    end
    
    ud = 0;
    S = addbc(S,ud);
    
    u = calc_init_dirichlet(S);
    mats = MATERIALS(S);
    for m=1:length(mats)
        mats{m} = setparam(mats{m},'u',u);
    end
    S = actualisematerials(S,mats);
    
    %% Stiffness matrices and sollicitation vectors
    % [A,b] = calc_rigi(S);
    % b = -b;
    
    %% Time scheme
    switch lower(symmetry)
        case {'isot','meanisot'} % almost surely or mean isotropic material
            if Dim==2
                switch lower(loading)
                    case 'tension'
                        % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM]
                        % du = 1e-5 mm during the first 500 time steps (up to u = 5e-3 mm)
                        % du = 1e-6 mm during the last 1300 time steps (up to u = 6.3e-3 mm)
                        % dt0 = 1e-8;
                        % nt0 = 500;
                        % dt1 = 1e-9;
                        % nt1 = 1300;
                        % t0 = linspace(dt0,nt0*dt0,nt0);
                        % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        % t = [t0,t1];
                        
                        % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Hesch, Weinberg, 2014, IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
                        % du = 1e-5 mm during 630 time steps (up to u = 6.3e-3 mm)
                        % dt = 1e-8;
                        % nt = 630;
                        % t = linspace(dt,nt*dt,nt);
                        
                        % [Zhou, Rabczuk, Zhuang, 2018, AES]
                        % du = 1e-5 mm during the first 450 time steps (up to u = 4.5e-3 mm)
                        % du = 1e-6 mm during the last 1800 time steps (up to u = 6.3e-3 mm)
                        % dt0 = 1e-8;
                        % nt0 = 450;
                        % dt1 = 1e-9;
                        % nt1 = 1800;
                        % t0 = linspace(dt0,nt0*dt0,nt0);
                        % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        % t = [t0,t1];

                        % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
                        % du = 1e-4 mm during 63 time steps (up to u = 6.3e-3 mm)
                        % dt = 1e-7;
                        % nt = 63;
                        % t = linspace(dt,nt*dt,nt);
                        
                        % [Liu, Li, Msekh, Zuo, 2016, CMS]
                        % du = 1e-4 mm during the first 50 time steps (up to u = 5e-3 mm)
                        % du = 1e-6 mm during the last 1300 time steps (up to u = 6.3e-3 mm)
                        % dt0 = 1e-7;
                        % nt0 = 50;
                        % dt1 = 1e-9;
                        % nt1 = 1300;
                        % t0 = linspace(dt0,nt0*dt0,nt0);
                        % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        % t = [t0,t1];
                        
                        % [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
                        % du = 2e-4 mm during 32 time steps (up to u = 6.4e-3 mm)
                        % dt = 2e-7;
                        % nt = 32;
                        % t = linspace(dt,nt*dt,nt);
                        
                        % du = 1e-5 mm during the first 400 time steps (up to u = 4e-3 mm)
                        % du = 1e-6 mm during the last 4000 time steps (up to u = 8e-3 mm)
                        dt0 = 1e-8;
                        nt0 = 400;
                        dt1 = 1e-9;
                        nt1 = 4000;
                        if test
                            dt0 = 1e-7;
                            nt0 = 40;
                            dt1 = 1e-8;
                            nt1 = 400;
                        end
                        t0 = linspace(dt0,nt0*dt0,nt0);
                        t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        t = [t0,t1];
                    case 'shear'
                        % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Hesch, Weinberg, 2014, IJNME]
                        % du = 1e-4 mm during the first 100 time steps (up to u = 10e-3 mm)
                        % du = 1e-6 mm during the last 10 000 time steps (up to u = 20e-3 mm)
                        % dt0 = 1e-7;
                        % nt0 = 100;
                        % dt1 = 1e-9;
                        % nt1 = 10000;
                        % t0 = linspace(dt0,nt0*dt0,nt0);
                        % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        % t = [t0,t1];
                        
                        % [Liu, Li, Msekh, Zuo, 2016, CMS]
                        % du = 1e-4 mm during the first 50 time steps (up to u = 5e-3 mm)
                        % du = 1e-5 mm during the last 1500 time steps (up to u = 20e-3 mm)
                        % dt0 = 1e-7;
                        % nt0 = 50;
                        % dt1 = 1e-8;
                        % nt1 = 1500;
                        % t0 = linspace(dt0,nt0*dt0,nt0);
                        % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        % t = [t0,t1];
                        
                        % [Zhou, Rabczuk, Zhuang, 2018, AES]
                        % du = 1e-4 mm during the first 80 time steps (up to u = 8e-3 mm)
                        % du = 1e-5 mm during the last 1200 time steps (up to u = 20e-3 mm)
                        % dt0 = 1e-7;
                        % nt0 = 80;
                        % dt1 = 1e-8;
                        % nt1 = 1200;
                        % t0 = linspace(dt0,nt0*dt0,nt0);
                        % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        % t = [t0,t1];
                        
                        % [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
                        % du = 1e-4 mm during 200 time steps (up to u = 20e-3 mm)
                        % dt = 1e-7;
                        % nt = 200;
                        % t = linspace(dt,nt*dt,nt);
                        
                        % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
                        % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM],
                        % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
                        % du = 1e-5 mm during 1500 time steps (up to u = 15e-3 mm)
                        % dt = 1e-8;
                        % nt = 1500;
                        
                        % du = 1e-5 mm during 2000 time steps (up to u = 20e-3 mm)
                        dt = 1e-8;
                        nt = 2000;
                        if test
                            dt = 5e-8;
                            nt = 400;
                        end
                        t = linspace(dt,nt*dt,nt);
                end
            elseif Dim==3
                % du = 1e-5 mm during 2500 time steps (up to u = 25e-3 mm)
                dt = 1e-8;
                nt = 2500;
                if test
                    dt = 1e-7;
                    nt = 250;
                end
                t = linspace(dt,nt*dt,nt);
            end
            T = TIMEMODEL(t);
            
        case 'anisot' % anisotropic material
            if Dim==2
                switch lower(loading)
                    case 'tension'
                        % du = 1e-5 mm (up to u = 10e-3 mm)
                        dt = 1e-8;
                        nt = 1000;
                        if test
                            dt1 = 4e-8;
                            nt = 250;
                        end
                    case 'shear'
                        % du = 1e-5 mm (up to u = 20e-3 mm)
                        dt = 1e-8;
                        nt = 2000;
                        if test
                            dt = 4e-8;
                            nt = 500;
                        end
                end
                t = linspace(dt,nt*dt,nt);
            elseif Dim==3
                % du = 1e-5 mm (up to u = 20e-3 mm)
                dt = 1e-8;
                nt = 2000;
                if test
                    dt = 2e-7;
                    nt = 100;
                end
                t = linspace(dt,nt*dt,nt);
            end
            T = TIMEMODEL(t);
            
        otherwise
            error('Wrong material symmetry class');
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','addbc','findddlforce','findddlboundary');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','addbc','findddlforce','findddlboundary');
end

%% Solution
if solveProblem
    myparallel('start',numWorkers);
    
    %% Solution
    tTotal = tic;
    
    nbSamples = 1;
    fun = @(S_phase,S) solvePFDetLinElas(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,'maxiter',maxIter,'tol',tolConv,'crit',critConv);
    % fun = @(S_phase,S) solvePFDetLinElasSingleEdgeCrack(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,'maxiter',maxIter,'tol',tolConv,'crit',critConv);
    [ft,dmaxt,dt_mean,ut_mean,dt_var,ut_var,dt_sample,ut_sample] = solvePFStoLinElas(S_phase,S,T,fun,N,'nbsamples',nbSamples);
    t = gettevol(T);
    idc = arrayfun(@(i) find(dmaxt(i,:)>=min(0.75,max(dmaxt(i,:))),1),1:N)';
    fc = arrayfun(@(i) ft(i,idc(i)),1:N)';
    udc = t(idc);
    [fmax,idmax] = max(ft,[],2);
    udmax = t(idmax);
    
    time = toc(tTotal);
    
    myparallel('stop');
    
    %% Statistical outputs of solution
    probs = [0.025 0.975];
    
    ft_mean = mean(ft);
    ft_std = std(ft);
    ft_ci = quantile(ft,probs);
    
    fmax_mean = mean(fmax);
    fmax_std = std(fmax);
    fmax_ci = quantile(fmax,probs);
    
    udmax_mean = mean(udmax);
    udmax_std = std(udmax);
    udmax_ci = quantile(udmax,probs);
    
    fc_mean = mean(fc);
    fc_std = std(fc);
    fc_ci = quantile(fc,probs);
    
    udc_mean = mean(udc);
    udc_std = std(udc);
    udc_ci = quantile(udc,probs);
    
    npts = 100;
    [fmax_f,fmax_xi,fmax_bw] = ksdensity(fmax,'npoints',npts);
    [udmax_f,udmax_xi,udmax_bw] = ksdensity(udmax,'npoints',npts);
    [fc_f,fc_xi,fc_bw] = ksdensity(fc,'npoints',npts);
    [udc_f,udc_xi,udc_bw] = ksdensity(udc,'npoints',npts);
    
    save(fullfile(pathname,'solution.mat'),'N','ft','dmaxt','dt_mean','ut_mean',...
        'dt_var','ut_var','dt_sample','ut_sample',...
        'ft_mean','ft_std','ft_ci','probs','time',...
        'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
        'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw',...
        'fc','fc_mean','fc_std','fc_ci','fc_f','fc_xi','fc_bw',...
        'udc','udc_mean','udc_std','udc_ci','udc_f','udc_xi','udc_bw');
else
    load(fullfile(pathname,'solution.mat'),'N','ft','dmaxt','dt_mean','ut_mean',...
        'dt_var','ut_var','dt_sample','ut_sample',...
        'ft_mean','ft_std','ft_ci','probs','time',...
        'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
        'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw',...
        'fc','fc_mean','fc_std','fc_ci','fc_f','fc_xi','fc_bw',...
        'udc','udc_mean','udc_std','udc_ci','udc_f','udc_xi','udc_bw');
end

%% Outputs
fid = fopen(fullfile(pathname,'results.txt'),'w');
fprintf(fid,'Single edge crack\n');
fprintf(fid,'\n');
fprintf(fid,'dim      = %d\n',Dim);
fprintf(fid,'loading  = %s\n',loading);
fprintf(fid,'mat sym  = %s\n',symmetry);
if strcmpi(symmetry,'anisot')
    fprintf(fid,'angle    = %g deg\n',ang);
end
fprintf(fid,'PF model = %s\n',PFmodel);
fprintf(fid,'PF split = %s\n',PFsplit);
fprintf(fid,'PF regularization = %s\n',PFregularization);
fprintf(fid,'PF solver = %s\n',PFsolver);
fprintf(fid,'nb elements = %g\n',getnbelem(S));
fprintf(fid,'nb nodes    = %g\n',getnbnode(S));
fprintf(fid,'nb dofs     = %g\n',getnbddl(S));
fprintf(fid,'nb time dofs = %g\n',getnbtimedof(T));
fprintf(fid,'nb samples = %g\n',N);
fprintf(fid,'elapsed time = %f s\n',time);
fprintf(fid,'\n');

if Dim==2
    fprintf(fid,'mean(fmax)   = %g kN/mm\n',fmax_mean*1e-6);
    fprintf(fid,'std(fmax)    = %g kN/mm\n',fmax_std*1e-6);
elseif Dim==3
    fprintf(fid,'mean(fmax)   = %g kN\n',fmax_mean*1e-3);
    fprintf(fid,'std(fmax)    = %g kN\n',fmax_std*1e-3);
end
fprintf(fid,'disp(fmax)   = %g\n',fmax_std/fmax_mean);
if Dim==2
    fprintf(fid,'%d%% ci(fmax) = [%g,%g] kN/mm\n',(probs(2)-probs(1))*100,fmax_ci(1)*1e-6,fmax_ci(2)*1e-6);
elseif Dim==3
    fprintf(fid,'%d%% ci(fmax) = [%g,%g] kN\n',(probs(2)-probs(1))*100,fmax_ci(1)*1e-3,fmax_ci(2)*1e-3);
end
fprintf(fid,'\n');

if Dim==2
    fprintf(fid,'mean(fc)   = %g kN/mm\n',fc_mean*1e-6);
    fprintf(fid,'std(fc)    = %g kN/mm\n',fc_std*1e-6);
elseif Dim==3
    fprintf(fid,'mean(fc)   = %g kN\n',fc_mean*1e-3);
    fprintf(fid,'std(fc)    = %g kN\n',fc_std*1e-3);
end
fprintf(fid,'disp(fc)   = %g\n',fc_std/fc_mean);
if Dim==2
    fprintf(fid,'%d%% ci(fc) = [%g,%g] kN/mm\n',(probs(2)-probs(1))*100,fc_ci(1)*1e-6,fc_ci(2)*1e-6);
elseif Dim==3
    fprintf(fid,'%d%% ci(fc) = [%g,%g] kN\n',(probs(2)-probs(1))*100,fc_ci(1)*1e-3,fc_ci(2)*1e-3);
end
fprintf(fid,'\n');

fprintf(fid,'mean(udmax)   = %g mm\n',udmax_mean*1e3);
fprintf(fid,'std(udmax)    = %g mm\n',udmax_std*1e3);
fprintf(fid,'disp(udmax)   = %g\n',udmax_std/udmax_mean);
fprintf(fid,'%d%% ci(udmax) = [%g,%g] mm\n',(probs(2)-probs(1))*100,udmax_ci(1)*1e3,udmax_ci(2)*1e3);
fprintf(fid,'\n');

fprintf(fid,'mean(udc)   = %g mm\n',udc_mean*1e3);
fprintf(fid,'std(udc)    = %g mm\n',udc_std*1e3);
fprintf(fid,'disp(udc)   = %g\n',udc_std/udc_mean);
fprintf(fid,'%d%% ci(udc) = [%g,%g] mm\n',(probs(2)-probs(1))*100,udc_ci(1)*1e3,udc_ci(2)*1e3);
fclose(fid);

%% Display
if displayModel
    [t,rep] = gettevol(T);
    
    %% Display domains, boundary conditions and meshes
%     plotDomain({D,C},'legend',false);
%     mysaveas(pathname,'domain',formats,renderer);
%     mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 0.5;
    v = calc_init_dirichlet(S);
    [hN,legN] = vectorplot(S,'U',v,ampl,'r','LineWidth',1);
    % legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions_displacement',formats,renderer);
    
    [hD_phase,legD_phase] = plotBoundaryConditions(S_phase,'legend',false);
    % legend([hD_phase,hN_phase],[legD_phase,legN_phase],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions_damage',formats,renderer);
    
    % plotModel(S,'legend',false);
    % mysaveas(pathname,'mesh',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    % u = ut(:,:,end);
    % for k=1:size(u,1)
    %     ampl = getsize(S)/max(abs(u(k,:)))/20;
    %     plotModelDeflection(S,u(k,:)','ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    %     mysaveas(pathname,['mesh_deflected_sample_' num2str(k)],formats,renderer);
    %
    %     figure('Name','Meshes')
    %     clf
    %     plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    %     plot(S+ampl*unfreevector(S,u(k,:)'),'Color','b','FaceColor','b','FaceAlpha',0.1);
    %     mysaveas(pathname,['meshes_deflected_' num2str(k)],formats,renderer);
    % end
end

%% Display statistics of solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display force-displacement curve
    figure('Name','Force vs displacement')
    clf
    plot(t*1e3,ft_mean*((Dim==2)*1e-6+(Dim==3)*1e-3),'-b','Linewidth',linewidth)
    hold on
    ciplot(ft_ci(1,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),ft_ci(2,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),t*1e3,'b');
    alpha(0.2)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    l = legend('mean function',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'Location','NorthWest');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    figure('Name','Forces vs displacement')
    clf
    color = distinguishable_colors(N);
    for i=1:N
        plot(t*1e3,ft(i,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),'LineStyle','-','Color',color(i,:),'Linewidth',linewidth)
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    mysaveas(pathname,'forces_displacement',formats);
    mymatlab2tikz(pathname,'forces_displacement.tex');
    
    %% Display pdf of maximum force
    figure('Name','Probability Density Estimate: Maximum force')
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
    ylabel('$p_{F_{\mathrm{max}}}(f)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_fmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_fmax.tex');
    
    %% Display pdf of critical force
    figure('Name','Probability Density Estimate: Critical force')
    clf
    plot(fc_xi*((Dim==2)*1e-6+(Dim==3)*1e-3),fc_f,'-b','LineWidth',linewidth)
    hold on
    ind_fc = find(fc_xi>=fc_ci(1) & fc_xi<fc_ci(2));
    area(fc_xi(ind_fc)*((Dim==2)*1e-6+(Dim==3)*1e-3),fc_f(ind_fc),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(fc_mean*((Dim==2)*1e-6+(Dim==3)*1e-3),0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$f$ [kN]','Interpreter',interpreter)
    ylabel('$p_{F_c}(f)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_fc',formats,renderer);
    mymatlab2tikz(pathname,'pdf_fc.tex');
    
    %% Display pdf of maximum displacement
    figure('Name','Probability Density Estimate: Maximum displacement')
    clf
    plot(udmax_xi*1e3,udmax_f,'-b','LineWidth',linewidth)
    hold on
    ind_udmax = find(udmax_xi>=udmax_ci(1) & udmax_xi<udmax_ci(2));
    area(udmax_xi(ind_udmax)*1e3,udmax_f(ind_udmax),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(udmax_mean*1e3,0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$u$ [mm]','Interpreter',interpreter)
    ylabel('$p_{U_{D,\mathrm{max}}}(u)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_udmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_udmax.tex');
    
    %% Display pdf of critical displacement
    figure('Name','Probability Density Estimate: Critical displacement')
    clf
    plot(udc_xi*1e3,udc_f,'-b','LineWidth',linewidth)
    hold on
    ind_udc = find(udc_xi>=udc_ci(1) & udc_xi<udc_ci(2));
    area(udc_xi(ind_udc)*1e3,udc_f(ind_udc),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(udc_mean*1e3,0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$u$ [mm]','Interpreter',interpreter)
    ylabel('$p_{U_{D,c}}(u)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_udc',formats,renderer);
    mymatlab2tikz(pathname,'pdf_udc.tex');
    
    %% Display means, variances and samples of solutions at different instants
    ampl = 0;
    switch lower(symmetry)
        case 'isot' % isotropic material
            switch lower(loading)
                case 'tension'
                    tSnapshots = [5.5 5.75 6 6.15 6.25 6.30 6.45 6.5]*1e-6;
                case 'shear'
                    tSnapshots = [1 1.25 1.35 1.5 1.75]*1e-5;
                otherwise
                    error('Wrong loading case');
            end
        case 'anisot' % anisotropic material
            switch lower(loading)
                case 'tension'
                    tSnapshots = [5 6 7 8 9]*1e-6;
                case 'shear'
                    tSnapshots = [1 1.25 1.35 1.5 1.75]*1e-5;
                otherwise
                    error('Wrong loading case');
            end
        otherwise
            error('Wrong material symmetry class');
    end
    rep = arrayfun(@(x) find(t<x+eps,1,'last'),tSnapshots);
    rep = [rep,length(T)];
    % tSnapshots = [tSnapshots,gett1(T)];
    % rep = arrayfun(@(x) find(t<x+eps,1,'last'),tSnapshots);
    
    for j=1:length(rep)
        dj = dt_mean(:,rep(j));
        dj_var = dt_var(:,rep(j));
        uj = ut_mean(:,rep(j));
        uj_var = ut_var(:,rep(j));
        
        plotSolution(S_phase,dj);
        mysaveas(pathname,['damage_mean_t' num2str(rep(j))],formats,renderer);
        plotSolution(S_phase,dj_var);
        mysaveas(pathname,['damage_var_t' num2str(rep(j))],formats,renderer);
        
        for i=1:Dim
            plotSolution(S,uj,'displ',i,'ampl',ampl);
            mysaveas(pathname,['displacement_' num2str(i) '_mean_t' num2str(rep(j))],formats,renderer);
            plotSolution(S,uj_var,'displ',i,'ampl',ampl);
            mysaveas(pathname,['displacement_' num2str(i) '_var_t' num2str(rep(j))],formats,renderer);
        end
        
        % for i=1:(Dim*(Dim+1)/2)
        %     plotSolution(S,uj,'epsilon',i,'ampl',ampl);
        %     mysaveas(pathname,['epsilon_' num2str(i) '_mean_t' num2str(rep(j))],formats,renderer);
        %
        %     plotSolution(S,uj,'sigma',i,'ampl',ampl);
        %     mysaveas(pathname,['sigma_' num2str(i) '_mean_t' num2str(rep(j))],formats,renderer);
        % end
        %
        % plotSolution(S,uj,'epsilon','mises','ampl',ampl);
        % mysaveas(pathname,['epsilon_von_mises_mean_t' num2str(rep(j))],formats,renderer);
        %
        % plotSolution(S,uj,'sigma','mises','ampl',ampl);
        % mysaveas(pathname,['sigma_von_mises_mean_t' num2str(rep(j))],formats,renderer);
        %
        % plotSolution(S,uj,'energyint','','ampl',ampl);
        % mysaveas(pathname,['internal_energy_mean_t' num2str(rep(j))],formats,renderer);
    end

    for k=1:size(dt_sample,1)
    for j=1:length(rep)
        dj = dt_sample(k,:,rep(j))';
        uj = ut_sample(k,:,rep(j))';
        
        plotSolution(S_phase,dj);
        mysaveas(pathname,['damage_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        
        for i=1:Dim
            plotSolution(S,uj,'displ',i,'ampl',ampl);
            mysaveas(pathname,['displacement_' num2str(i) '_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        end
        
        % for i=1:(Dim*(Dim+1)/2)
        %     plotSolution(S,uj,'epsilon',i,'ampl',ampl);
        %     mysaveas(pathname,['epsilon_' num2str(i) '_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        %
        %     plotSolution(S,uj,'sigma',i,'ampl',ampl);
        %     mysaveas(pathname,['sigma_' num2str(i) '_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        % end
        %
        % plotSolution(S,uj,'epsilon','mises','ampl',ampl);
        % mysaveas(pathname,['epsilon_von_mises_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        %
        % plotSolution(S,uj,'sigma','mises','ampl',ampl);
        % mysaveas(pathname,['sigma_von_mises_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        %
        % plotSolution(S,uj,'energyint','','ampl',ampl);
        % mysaveas(pathname,['internal_energy_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
    end
    end
    
end

%% Display evolution of means, variances and samples of solutions
if makeMovie
    sz_d = [getnbddl(S_phase),getnbtimedof(T)];
    sz_u = [getnbddl(S),getnbtimedof(T)];
    ampl = 0;
    % ampl = getsize(S)/max(max(max(abs(ut))))/20;
    
    options = {'plotiter',true,'plottime',false};
    % framerate = 80;
    switch lower(loading)
        case 'tension'
            framerate = 400;
        case 'shear'
            framerate = 200;
        otherwise
            error('Wrong loading case');
    end
    
    dk = TIMEMATRIX(reshape(dt_mean(:,:),sz_d),T);
    dk_var = TIMEMATRIX(reshape(dt_var(:,:),sz_d),T);
    % uk = TIMEMATRIX(reshape(ut_mean(:,:),sz_u),T);
    % uk_var = TIMEMATRIX(reshape(ut_var(:,:),sz_u),T);
    
    evolSolution(S_phase,dk,'FrameRate',framerate,'filename','damage_mean','pathname',pathname,options{:});
    evolSolution(S_phase,dk_var,'FrameRate',framerate,'filename','damage_var','pathname',pathname,options{:});
    % for i=1:Dim
    %     evolSolution(S,uk,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i) '_mean'],'pathname',pathname,options{:});
    %     evolSolution(S,uk_var,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i) '_var'],'pathname',pathname,options{:});
    % end
    %
    % for i=1:(Dim*(Dim+1)/2)
    %     evolSolution(S,uk,'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['epsilon_' num2str(i) '_mean'],'pathname',pathname,options{:});
    %     evolSolution(S,uk,'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['sigma_' num2str(i) '_mean'],'pathname',pathname,options{:});
    % end
    %
    % evolSolution(S,uk,'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename','epsilon_von_mises_mean','pathname',pathname,options{:});
    % evolSolution(S,uk,'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename','sigma_von_mises_mean','pathname',pathname,options{:});
    % evolSolution(S,uk,'energyint','','ampl',ampl,'FrameRate',framerate,'filename','internal_energy_mean','pathname',pathname,options{:});

    for k=1:size(dt_sample,1)
        dk = TIMEMATRIX(reshape(dt_sample(k,:,:),sz_d),T);
        % uk = TIMEMATRIX(reshape(ut_sample(k,:,:),sz_u),T);
        
        evolSolution(S_phase,dk,'FrameRate',framerate,'filename',['damage_sample_' num2str(k)],'pathname',pathname,options{:});
        % for i=1:Dim
        %     evolSolution(S,uk,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i) '_sample_' num2str(k)],'pathname',pathname,options{:});
        % end
        %
        % for i=1:(Dim*(Dim+1)/2)
        %     evolSolution(S,uk,'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['epsilon_' num2str(i) '_sample_' num2str(k)],'pathname',pathname,options{:});
        %     evolSolution(S,uk,'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['sigma_' num2str(i) '_sample_' num2str(k)],'pathname',pathname,options{:});
        % end
        %
        % evolSolution(S,uk,'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename',['epsilon_von_mises_sample_' num2str(k)],'pathname',pathname,options{:});
        % evolSolution(S,uk,'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename',['sigma_von_mises_sample_' num2str(k)],'pathname',pathname,options{:});
        % evolSolution(S,uk,'energyint','','ampl',ampl,'FrameRate',framerate,'filename',['internal_energy_sample_' num2str(k)],'pathname',pathname,options{:});
    end
end

%% Save means, variances and samples of solutions
if saveParaview
    [t,rep] = gettevol(T);
    for i=1:length(T)
        di = dt_mean(:,rep(i))';
        ui = ut_mean(:,rep(i))';
        dvi = dt_var(:,rep(i))';
        uvi = ut_var(:,rep(i))';
        
        write_vtk_mesh(S,{di,ui,dvi,uvi},[],...
            {'damage_mean','displacement_mean','damage_variance','displacement_variance'},[],...
            pathname,'solution_mean_variance',1,i-1);
    end
    make_pvd_file(pathname,'solution_mean_variance',1,length(T));

    for k=1:size(dt_sample,1)
        for i=1:length(T)
            di = dt_sample(k,:,rep(i))';
            ui = ut_sample(k,:,rep(i))';
            
            write_vtk_mesh(S,{di,ui},[],...
                {'damage','displacement'},[],...
                pathname,['solution_sample_' num2str(k)],1,i-1);
        end
        make_pvd_file(pathname,['solution_sample_' num2str(k)],1,length(T));
    end
end
