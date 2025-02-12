%% Phase field fracture model - deterministic linear elasticity problem with single edge crack %%
%%---------------------------------------------------------------------------------------------%%
% [Bourdin, Francfort, Marigo, 2000, JMPS] (isotropic phase field model with no split of Bourdin et al.)
% [Miehe, Welschinger, Hofacker, 2010 IJNME] (anisotropic phase field model of Miehe et al.)
% [Miehe, Hofacker, Welschinger, 2010, CMAME] (anisotropic phase field model of Miehe et al.)
% [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME] (anisotropic phase field model of Miehe et al.)
% [Hesch, Weinberg, 2014, IJNME] (anisotropic phase field model of Miehe et al.)
% [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM] (anisotropic phase field model of Miehe et al.)
% [Ambati, Gerasimov, De Lorenzis, 2015, CM] (hybrid isotropic-anisotropic phase field model of Ambati et al. with Miehe et al. decomposition compared with the isotropic one of Bourdin et al. and the anisotropic ones of Amor et al. and Miehe et al.)
% [Gerasimov, De Lorenzis, 2016, CMAME] (anisotropic phase field model of Miehe et al.)
% [Liu, Li, Msekh, Zuo, 2016, CMS] (anisotropic phase field model of Miehe et al.)
% [Zhou, Rabczuk, Zhuang, 2018, AES] (anisotropic phase field model of Miehe et al.)
% [Wu, Nguyen, 2018, JMPS] (hybrid isotropic-anisotropic phase field model of Wu et al.)
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM] (anisotropic phase field model of Wu et al.)
% [Gerasimov, De Lorenzis, 2019, CMAME] (anisotropic phase field model of Amor et al. and Miehe et al.)
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
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displayModel = false;
displaySolution = false;
makeMovie = false;
saveParaview = false;

% test = true; % coarse mesh
test = false; % fine mesh

Dim = 2; % space dimension Dim = 2, 3
symmetry = 'Isot'; % 'Isot' or 'Anisot'. Material symmetry
ang = 45; % clockwise material orientation angle around z-axis for anisotopic material [deg]
loading = 'Shear'; % 'Tension' or 'Shear'
PFmodel = 'Miehe'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 1; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
initialCrack = 'GeometricCrack'; % 'GeometricCrack', 'GeometricNotch', 'InitialPhaseField'
FEmesh = 'Optim'; % 'Unif' or 'Optim'
coeff_gc = 1.0;

% PFmodels = {'Bourdin','Amor','Miehe','HeAmor','HeFreddi','Zhang'};
% PFsplits = {'Strain','Stress'};
% PFregularizations = {'AT1','AT2'};
% PFsolvers = {'HistoryFieldElem','BoundConstrainedOptim'};
% initialCracks = {'GeometricCrack','InitialPhaseField'};
% maxIters = [1,100];
% coeffs_gc = [0.6,0.8,1.0,1.2,1.4];

% for iPFmodel=1:length(PFmodels)
% PFmodel = PFmodels{iPFmodel};
% for iPFsplit=1:length(PFsplits)
% PFsplit = PFsplits{iPFsplit};
% for iPFRegularization=1:length(PFregularizations)
% PFregularization = PFregularizations{iPFRegularization};
% for iPFsolver=1:length(PFsolvers)
% PFsolver = PFsolvers{iPFsolver};
% for iinitialCrack=1:length(initialCracks)
% initialCrack = initialCracks{iinitialCrack};
% for imaxIter=1:length(maxIters)
% maxIter = maxIters(imaxIter);
% for icoeff_gc=1:length(coeffs_gc)
% coeff_gc = coeffs_gc(icoeff_gc);
% close all

suffix = ['_coeffgc' num2str(coeff_gc,'_%g')];

foldername = ['singleEdgeCrack' loading '_' num2str(Dim) 'D'];
filename = ['linElas' symmetry];
if strcmpi(symmetry,'anisot') % anisotropic material
    filename = [filename num2str(ang) 'deg'];
end
filename = [filename PFmodel PFsplit PFregularization PFsolver initialCrack...% 'MaxIter' num2str(maxIter) 'Tol' num2str(tolConv)
    'Mesh' FEmesh];
filename = [filename suffix];

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','phasefieldDet',foldername,filename);
if test
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','phasefieldDet_test',foldername,filename);
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
    % clC = 2.5e-6; % [Gerasimov, De Lorenzis, 2019, CMAME], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
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
            ZMin = 0; ZMax = e;
            Thickness = a;
            % Thickness = 0;
            B = struct('VIn',VIn,'VOut',VOut,'XMin',XMin,'XMax',XMax,'YMin',YMin,'YMax',YMax,'ZMin',ZMin,'ZMax',ZMax,'Thickness',Thickness);
        otherwise
            error('Wrong FE mesh')
    end
    % switch lower(initialCrack)
    %     case 'geometriccrack'
    %         S_phase = gmshdomainwithedgecrack(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,'Box',B);
    %     case 'geometricnotch'
    %         c = 1e-5; % crack width
    %         S_phase = gmshdomainwithedgenotch(D,C,c,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,'Box',B);
    %     case 'initialphasefield'
    %         S_phase = gmshdomainwithedgecrack(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,'noduplicate','Box',B);
    %     otherwise
    %         error('Wrong model for initial crack');
    % end
    S_phase = gmsh2femobject(Dim,fullfile(pathname,'gmsh_domain_single_edge_crack.msh'),Dim);
    S = S_phase;
    
    %% Phase field problem
    %% Material
    switch lower(symmetry)
        case 'isot' % isotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 2.7e3; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
            % Regularization parameter (width of the smeared crack)
            % l = 3.75e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME]
            % l = 3e-5; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
            % l = 1.5e-5; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Hesch, Weinberg, 2014, IJNME], [Liu, Li, Msekh, Zuo, 2016, CMS], [Zhou, Rabczuk, Zhuang, 2018, AES], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM], [Hu, Guilleminot, Dolbow, 2020, CMAME]
            l = 1e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Gerasimov, De Lorenzis, 2016, CMAME], [Wu, Nguyen, 2018, JMPS], [Gerasimov, De Lorenzis, 2019, CMAME], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
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
    gc = gc*coeff_gc;
    % Small artificial residual stiffness
    % k = 1e-12;
    k = 0;
    
    % Material
    [K,R,Qn] = setphasefieldparam(gc,l,PFregularization);
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'PFregularization',PFregularization);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    switch lower(initialCrack)
        case 'geometriccrack'
            S_phase = final(S_phase,'duplicate');
        case 'geometricnotch'
            S_phase = final(S_phase);
        case 'initialphasefield'
            S_phase = final(S_phase);
            S_phase = addcl(S_phase,C,'T',1);
        otherwise
            error('Wrong model for initial crack');
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
    option = 'DEFO'; % plane strain [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Hesch, Weinberg, 2014, IJNME], [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Gerasimov, De Lorenzis, 2016, CMAME], [Zhou, Rabczuk, Zhuang, 2018, AES], [Gerasimov, De Lorenzis, 2019, CMAME], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
    % option = 'CONT'; % plane stress [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
    switch lower(symmetry)
        case 'isot' % isotropic material
            % Lame coefficients
            % lambda = 121.1538e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
            % mu = 80.7692e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
            % lambda = 121.154e9; % [Hesch, Weinberg, 2014, IJNME]
            % mu = 80.769e9; % [Hesch, Weinberg, 2014, IJNME]
            lambda = 121.15e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME], [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
            mu = 80.77e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME], [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
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
                % E = 210e9; NU = 0.3; % [Gerasimov, De Lorenzis, 2016, CMAME], [Zhou, Rabczuk, Zhuang, 2018, AES], [Wu, Nguyen, 2018, JMPS], [Gerasimov, De Lorenzis, 2019, CMAME], [Wu, Nguyen, Zhou, Huang, 2020, CMAME], [Kristensen, Martinez-Paneda, 2020, TAFM]
                % kappa = 121030e6; NU=0.227; lambda=3*kappa*NU/(1+NU); mu = 3*kappa*(1-2*NU)/(2*(1+NU)); E = 3*kappa*(1-2*NU); % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
            elseif Dim==3
                E = mu*(3*lambda+2*mu)/(lambda+mu);
                NU = lambda/(lambda+mu)/2;
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
        case 'isot' % isotropic material
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
        case 'anisot' % anisotropic material
            mat = ELAS_ANISOT('C',Cmat,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
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
        BRight = LIGNE([L,0.0],[L,L]);
        BLeft = LIGNE([0.0,0.0],[0.0,L]);
        BFront = [];
        BBack = [];
    elseif Dim==3
        BU = PLAN([0.0,L,0.0],[L,L,0.0],[0.0,L,e]);
        BL = PLAN([0.0,0.0,0.0],[L,0.0,0.0],[0.0,0.0,e]);
        BRight = PLAN([L,0.0,0.0],[L,L,0.0],[L,0.0,e]);
        BLeft = PLAN([0.0,0.0,0.0],[0.0,L,0.0],[0.0,0.0,e]);
        BFront = PLAN([0.0,0.0,e],[L,0.0,e],[0.0,L,e]);
        BBack = PLAN([0.0,0.0,0.0],[L,0.0,0.0],[0.0,L,0.0]);
    end
    
    switch lower(initialCrack)
        case 'geometriccrack'
            S = final(S,'duplicate');
        otherwise
            S = final(S);
    end
    
    ud = 0;
    switch lower(loading)
        case 'tension'
            if Dim==2
                S = addcl(S,BU,{'UX','UY'},[0;ud]);
            elseif Dim==3
                S = addcl(S,BU,{'UX','UY','UZ'},[0;ud;0]);
            end
            S = addcl(S,BL,'UY');
        case 'shear'
            if Dim==2
                S = addcl(S,BU,{'UX','UY'},[ud;0]);
                S = addcl(S,BLeft,'UY');
                S = addcl(S,BRight,'UY');
            elseif Dim==3
                S = addcl(S,BU,{'UX','UY','UZ'},[ud;0;0]);
                S = addcl(S,BLeft,{'UY','UZ'});
                S = addcl(S,BRight,{'UY','UZ'});
                S = addcl(S,BFront,{'UY','UZ'});
                S = addcl(S,BBack,{'UY','UZ'});
            end
            S = addcl(S,BL);
        otherwise
            error('Wrong loading case');
    end
    
    %% Stiffness matrices and sollicitation vectors
    % [A,b] = calc_rigi(S);
    % b = -b;
    
    %% Time scheme
    switch lower(symmetry)
        case 'isot' % isotropic material
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
                        % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
                        % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM],
                        % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
                        % du = 1e-5 mm during 1500 time steps (up to u = 15e-3 mm)
                        % dt = 1e-8;
                        % nt = 1500;
                        
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
                        
                        % [Gerasimov, De Lorenzis, 2016, CMAME]
                        % du = 6e-3 mm at the first time step
                        % du = 3e-4 mm during the last 46 time steps (up to u = 20e-3 mm)
                        % t0 = 6e-6;
                        % dt1 = 3e-7;
                        % nt1 = 47;
                        % t1 = linspace(t0+dt1,t0+nt1*dt1,nt1);
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
                        % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
                        % du = 6e-5 mm during the first stage (until the phase field reaches the threshold value)
                        % du = 2e-5 mm during the last stage (as soon as the phase field exceeds the threshold value, up to u = 10e-3 mm)
                        % dt0 = 6e-8;
                        % dt1 = 2e-8;
                        % if test
                        %     dt0 = 12e-8;
                        %     dt1 = 4e-8;
                        % end
                        % tf = 10e-6;
                        % dthreshold = 0.6;
                        
                        % du = 1e-5 mm during the first stage (until the phase field reaches the threshold value)
                        % du = 1e-6 mm during the last stage (as soon as the phase field exceeds the threshold value, up to u = 10e-3 mm)
                        dt0 = 1e-8;
                        dt1 = 1e-9;
                        if test
                            dt0 = 1e-7;
                            dt1 = 1e-8;
                        end
                        tf = 10e-6;
                        dthreshold = 0.6;
                    case 'shear'
                        % du = 6e-5 mm during the first stage (until the phase field reaches the threshold value)
                        % du = 2e-5 mm during the last stage (as soon as the phase field exceeds the threshold value, up to u = 20e-3 mm)
                        % dt0 = 6e-8;
                        % dt1 = 2e-8;
                        % if test
                        %     dt0 = 12e-8;
                        %     dt1 = 4e-8;
                        % end
                        % tf = 20e-6;
                        % dthreshold = 0.6;
                        
                        % du = 1e-5 mm during the first stage (until the phase field reaches the threshold value)
                        % du = 1e-5 mm during the last stage (as soon as the phase field exceeds the threshold value, up to u = 20e-3 mm)
                        dt0 = 1e-8;
                        dt1 = 1e-8;
                        if test
                            dt0 = 1e-7;
                            dt1 = 1e-7;
                        end
                        tf = 20e-6;
                        dthreshold = 0.6;
                end
            elseif Dim==3
                % du = 2e-5 mm (up to u = 20e-3 mm)
                dt0 = 2e-8;
                dt1 = 2e-8;
                if test
                    dt0 = 2e-7;
                    dt1 = 2e-7;
                end
                tf = 20e-6;
                dthreshold = 0.6;
            end
            T = struct('dt0',dt0,'dt1',dt1,'tf',tf,'dthreshold',dthreshold);
            
        otherwise
            error('Wrong material symmetry class');
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','BU','BL','BRight','BLeft','BFront','BBack','loading','symmetry','ang');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','BU','BL','BRight','BLeft','BFront','BBack','loading','symmetry','ang');
end

%% Solution
if solveProblem
    tTotal = tic;
    
    switch lower(symmetry)
        case 'isot' % isotropic material
            fun = @solvePFDetLinElasSingleEdgeCrack;
        case 'anisot' % anisotropic material
            fun = @solvePFDetLinElasSingleEdgeCrackThreshold;
        otherwise
            error('Wrong material symmetry class');
    end
    switch lower(PFsolver)
        case {'historyfieldelem','historyfieldnode'}
            [dt,ut,ft,Ht] = fun(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,'maxiter',maxIter,'tol',tolConv);
        otherwise
            [dt,ut,ft] = fun(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,'maxiter',maxIter,'tol',tolConv);
    end
    
    if strcmpi(symmetry,'anisot')
        T = gettimemodel(dt);
    end
    t = gettevol(T);
    dt_val = getvalue(dt);
    dmaxt = max(dt_val);
    idc = find(dmaxt>=min(0.75,max(dmaxt)),1);
    fc = ft(idc);
    udc = t(idc);
    [fmax,idmax] = max(ft,[],2);
    udmax = t(idmax);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'dt','ut','ft','dmaxt','fmax','udmax','fc','udc','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        save(fullfile(pathname,'solution.mat'),'Ht','-append');
    end
    if strcmpi(symmetry,'anisot')
        save(fullfile(pathname,'solution.mat'),'T','-append');
    end
else
    load(fullfile(pathname,'solution.mat'),'dt','ut','ft','dmaxt','fmax','udmax','fc','udc','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        load(fullfile(pathname,'solution.mat'),'Ht');
    end
    if strcmpi(symmetry,'anisot')
        load(fullfile(pathname,'solution.mat'),'T');
    end
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
fprintf(fid,'elapsed time = %f s\n',time);
fprintf(fid,'\n');

if Dim==2
    fprintf(fid,'fmax  = %g kN/mm\n',fmax*1e-6);
    fprintf(fid,'fc    = %g kN/mm\n',fc*1e-6);
elseif Dim==3
    fprintf(fid,'fmax  = %g kN\n',fmax*1e-3);
    fprintf(fid,'fc    = %g kN\n',fc*1e-3);
end
fprintf(fid,'udmax = %g mm\n',udmax*1e3);
fprintf(fid,'udc   = %g mm\n',udc*1e3);
fclose(fid);

%% Display
if displayModel
    [t,rep] = gettevol(T);
    
    %% Display domains, boundary conditions and meshes
    % plotDomain({D,C},'legend',false);
    % mysaveas(pathname,'domain',formats,renderer);
    % mymatlab2tikz(pathname,'domain.tex');
    
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
    
    u = getmatrixatstep(ut,rep(end));
    ampl = getsize(S)/max(abs(u))/20;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*unfreevector(S,u),'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
end

%% Display solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display force-displacement curve
    figure('Name','Force vs displacement')
    clf
    plot(t*1e3,ft*((Dim==2)*1e-6+(Dim==3)*1e-3),'-b','LineWidth',linewidth)
    % hold on
    % scatter(udmax*1e3,fmax*((Dim==2)*1e-6+(Dim==3)*1e-3),'Marker','+','MarkerEdgeColor','b','LineWidth',linewidth)
    % scatter(udc*1e3,fc*((Dim==2)*1e-6+(Dim==3)*1e-3),'Marker','+','MarkerEdgeColor','r','LineWidth',linewidth)
    % hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    %% Display maximum damage-displacement curve
    figure('Name','Maximum damage vs displacement')
    clf
    plot(t*1e3,dmaxt,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Maximum damage','Interpreter',interpreter)
    mysaveas(pathname,'max_damage_displacement',formats);
    mymatlab2tikz(pathname,'max_damage_displacement.tex');
    
    %% Display solutions at different instants
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
        dj = getmatrixatstep(dt,rep(j));
        uj = getmatrixatstep(ut,rep(j));
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            Hj = getmatrixatstep(Ht,rep(j));
        end
        
        plotSolution(S_phase,dj);
        mysaveas(pathname,['damage_t' num2str(rep(j))],formats,renderer);
        
        for i=1:Dim
            plotSolution(S,uj,'displ',i,'ampl',ampl);
            mysaveas(pathname,['displacement_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        end
        
        % for i=1:(Dim*(Dim+1)/2)
        %     plotSolution(S,uj,'epsilon',i,'ampl',ampl);
        %     mysaveas(pathname,['epsilon_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        %
        %     plotSolution(S,uj,'sigma',i,'ampl',ampl);
        %     mysaveas(pathname,['sigma_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        % end
        %
        % plotSolution(S,uj,'epsilon','mises','ampl',ampl);
        % mysaveas(pathname,['epsilon_von_mises_t' num2str(rep(j))],formats,renderer);
        %
        % plotSolution(S,uj,'sigma','mises','ampl',ampl);
        % mysaveas(pathname,['sigma_von_mises_t' num2str(rep(j))],formats,renderer);
        %
        % plotSolution(S,uj,'energyint','','ampl',ampl);
        % mysaveas(pathname,['internal_energy_density_t' num2str(rep(j))],formats,renderer);
        %
        % if strcmpi(PFsolver,'historyfieldelem')
        %     figure('Name','Solution H')
        %     clf
        %     plot(Hj,S_phase);
        %     colorbar
        %     set(gca,'FontSize',fontsize)
        %     mysaveas(pathname,['internal_energy_density_history_t' num2str(rep(j))],formats,renderer);
        % elseif strcmpi(PFsolver,'historyfieldnode')
        %     plotSolution(S_phase,Hj,'ampl',ampl);
        %     mysaveas(pathname,['internal_energy_density_history_t' num2str(rep(j))],formats,renderer);
        % end
    end
end

%% Display evolution of solutions
if makeMovie
    ampl = 0;
    % ampl = getsize(S)/max(max(abs(getvalue(ut))))/20;
    
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
    
    evolSolution(S_phase,dt,'FrameRate',framerate,'filename','damage','pathname',pathname,options{:});
    % for i=1:Dim
    %     evolSolution(S,ut,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i)],'pathname',pathname,options{:});
    % end
    %
    % for i=1:(Dim*(Dim+1)/2)
    %     evolSolution(S,ut,'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['epsilon_' num2str(i)],'pathname',pathname,options{:});
    %     evolSolution(S,ut,'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['sigma_' num2str(i)],'pathname',pathname,options{:});
    % end
    %
    % evolSolution(S,ut,'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename','epsilon_von_mises','pathname',pathname,options{:});
    % evolSolution(S,ut,'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename','sigma_von_mises','pathname',pathname,options{:});
    % evolSolution(S,ut,'energyint','','ampl',ampl,'FrameRate',framerate,'filename','internal_energy_density','pathname',pathname,options{:});
    % if strcmpi(PFsolver,'historyfieldelem')
    %     figure('Name','Solution H')
    %     clf
    %     T = setevolparam(T,'colorbar',true,'FontSize',fontsize,options{:});
    %     frame = evol(T,Ht,S_phase,'rescale',true);
    %     saveMovie(frame,'FrameRate',framerate,'filename','internal_energy_density_history','pathname',pathname);
    % elseif strcmpi(PFsolver,'historyfieldnode')
    %     evolSolution(S_phase,Ht,'ampl',ampl,'FrameRate',framerate,'filename','internal_energy_density_history','pathname',pathname,options{:});
    % end
end

%% Save solutions
if saveParaview
    [t,rep] = gettevol(T);
    for i=1:length(T)
        di = getmatrixatstep(dt,rep(i));
        ui = getmatrixatstep(ut,rep(i));
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            Hi = getmatrixatstep(Ht,rep(i));
        end
        
        switch lower(PFsolver)
            case 'historyfieldelem'
                write_vtk_mesh(S,{di,ui},{Hi},...
                    {'damage','displacement'},{'internal energy density history'},...
                    pathname,'solution',1,i-1);
            case 'historyfieldnode'
                write_vtk_mesh(S,{di,ui,Hi},[],...
                    {'damage','displacement','internal energy density history'},[],...
                    pathname,'solution',1,i-1);
            otherwise
                write_vtk_mesh(S,{di,ui},[],...
                    {'damage','displacement'},[],...
                    pathname,'solution',1,i-1);
        end
    end
    make_pvd_file(pathname,'solution',1,length(T));
end

% myparallel('stop');

% end
% end
% end
% end
% end
% end
% end
