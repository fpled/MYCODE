%% Phase-field fracture model - deterministic linear elasticity problem %%
%  Asymmetric notched plate with three holes under three-point bending  %%
%%----------------------------------------------------------------------%%
% [Ingraffea, Grigoriu, 1990] (experimental tests and LEFM-based FEM)
% [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM] (LEFM SIF-based method with local remeshing and special FE)
% [Ventura, Xu, Belytschko, 2002, IJNME] (vector level set method with discontinuous enrichment in meshless method)
% [Miehe, Gurses, 2007, IJNME] (FEM + R-adaptive mesh alignment)
% [Guidault, Allix, Champaney, Cornuault, 2008, CMAME] (MsXFEM)
% [Miehe, Welschinger, Hofacker, 2010, IJNME] (anisotropic phase-field model of Miehe et al.)
% [Miehe, Hofacker, Welschinger, 2010, CMAME] (anisotropic phase-field model of Miehe et al.)
% [Häusler, Lindhorst, Horst, 2011, IJNME] (XFEM)
% [Geniaut, Galenne, 2012, IJSS] (XFEM)
% [Passieux, Rethore, Gravouil, Baietto, 2013, CM] (XFEM)
% [Ambati, Gerasimov, De Lorenzis, 2015, CM] (hybrid isotropic-anisotropic phase-field model of Ambati et al. compared with the anisotropic one of Miehe et al.)
% [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME] (isotropic phase-field model with no split of Bourdin et al. compared to experimental data of [Winkler, 2001, PhD thesis])
% [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS] (isotropic phase-field model with no split of Bourdin et al.)
% [Molnar, Gravouil, 2017, FEAD] (isotropic phase-field model with no split of Bourdin et al.)
% [Cervera, Barbat, Chiumenti, 2017, CM] (LEFM Mixed FEM)
% [Khisamitov, Meschke, 2018, CMAME] (anisotropic phase-field model for interfacial elastic energy)
% [Wu, 2018, CMAME] (PF-CZM, hybrid isotropic-anisotropic phase-field model of Wu et al.)
% [Wu, Nguyen, 2018, JMPS] (PF-CZM, hybrid isotropic-anisotropic phase-field model of Wu et al.)
% [Bhowmick, Liu, 2018, EFM] (anisotropic phase-field model of Miehe et al. + CS-FEM)
% [Patil, Mishra, Singh, 2018, CMAME] (AMsPFM = anisotropic PFM of Miehe et.al + MsFEM)
% [Patil, Mishra, Singh, 2018, CMAME] (LMXPFM = anisotropic PFM of Miehe et.al + XFEM + MsFEM)
% [Mandal, Nguyen, Wu, 2019, EFM] (hybrid AT1, AT2 and PF-CZM)
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM] (PF-CZM, anisotropic phase-field model of Wu et al.)
% [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME] (anisotropic phase-field model of Miehe et al.)
% [Fu, Yi, Chen, Bui, Hu, Yao, 2020, TAFM] (LEFM, FNM + SASE)
% [Min, Hu, Yao, Bui, Zhang, 2022, CMAME] (isotropic phase-field model with no split of Bourdin et al.)
% [Rahimi, Moutsanidis, 2022, CMAME] (anisotropic phase-field model with Total Lagrangian SPH approximation)
% [Lee, Meng, Ashkpour, Rahmaninezhad, Iqbal, Mishra, Hubler, Sales, Farnam, Najafi, 2024, CBM] (anisotropic phase-field model of Miehe et al. with healing and re-damage)

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

test = true; % coarse mesh
% test = false; % fine mesh

Dim = 2; % space dimension Dim = 2
setup = 1; % notch geometry setup = 1, 2, 3, 4, 5, 6
PFmodel = 'Miehe'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 1; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
critConv = 'Energy'; % 'Solution', 'Residual', 'Energy'
initialCrack = 'GeometricCrack'; % 'GeometricCrack', 'GeometricNotch', 'InitialPhaseField'
FEmesh = 'Optim'; % 'Unif' or 'Optim'

suffix = '';

foldername = ['asymmetricNotchedPlateSetup' num2str(setup) '_' num2str(Dim) 'D'];
filename = ['linElas' PFmodel PFsplit PFregularization PFsolver initialCrack...
    'MaxIter' num2str(maxIter)];
if maxIter>1
    filename = [filename 'Tol' num2str(tolConv) num2str(critConv)];
end
filename = [filename 'Mesh' FEmesh suffix];

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
    unit = 1e-3; % in [mm] [Guidault, Allix, Champaney, Cornuault, 2008, CMAME], [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Passieux, Rethore, Gravouil, Baietto, 2013, CM],
    % [Molnar, Gravouil, 2017, FEAD], [Khisamitov, Meschke, 2018, CMAME], [Bhowmick, Liu, 2018, EFM], [Patil, Mishra, Singh, 2018, CMAME] (AMsPFM), [Patil, Mishra, Singh, 2018, CMAME] (LMXPFM), [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME], [Lee et al., 2024, CBM]
    % unit = 25.4e-3; % in [inch] [Ingraffea, Grigoriu, 1990], [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM], [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME],
    % [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS], [Cervera, Barbat, Chiumenti, 2017, CM], [Wu, Nguyen, 2018, JMPS], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
    switch setup
        case 1 % [Ingraffea, Grigoriu, 1990], [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM], [Ventura, Xu, Belytschko, 2002, IJNME],
            % [Miehe, Gurses, 2007, IJNME], [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Häusler, Lindhorst, Horst, 2011, IJNME],
            % [Geniaut, Galenne, 2012, IJSS], [Passieux, Rethore, Gravouil, Baietto, 2013, CM], [Ambati, Gerasimov, De Lorenzis, 2015, CM],
            % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS], [Molnar, Gravouil, 2017, FEAD], [Cervera, Barbat, Chiumenti, 2017, CM],
            % [Khisamitov, Meschke, 2018, CMAME], [Wu, Nguyen, 2018, JMPS], [Patil, Mishra, Singh, 2018, CMAME] (AMsPFM), [Patil, Mishra, Singh, 2018, CMAME] (LMXPFM), [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM],
            % [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME], [Fu, Yi, Chen, Bui, Hu, Yao, 2020, TAFM],
            % [Min, Hu, Yao, Bui, Zhang, 2022, CMAME], [Rahimi, Moutsanidis, 2022, CMAME], [Lee et al., 2024, CBM]
            a = 1*unit; % crack length
            b = 6*unit; % crack offset from the centerline
        case 2 % [Ingraffea, Grigoriu, 1990], [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
            a = 2.5*unit; % crack length
            b = 6*unit; % crack offset from the centerline
        case 3 % [Wu, Nguyen, 2018, JMPS], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
            a = 1.5*unit; % crack length
            b = 5.15*unit; % crack offset from the centerline
        case 4 % [Ingraffea, Grigoriu, 1990], [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM], [Ventura, Xu, Belytschko, 2002, IJNME],
            % [Miehe, Gurses, 2007, IJNME], [Guidault, Allix, Champaney, Cornuault, 2008, CMAME], [Häusler, Lindhorst, Horst, 2011, IJNME], [Geniaut, Galenne,2012, IJSS],
            % [Passieux, Rethore, Gravouil, Baietto, 2013, CM], [Molnar, Gravouil, 2017, FEAD], [Khisamitov, Meschke, 2018, CMAME]
            a = 1.5*unit; % crack length
            b = 5*unit; % crack offset from the centerline
        case 5 % [Ingraffea, Grigoriu, 1990], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
            a = 1.5*unit; % crack length
            b = 4.75*unit; % crack offset from the centerline
        case 6 % [Ingraffea, Grigoriu, 1990], [Cervera, Barbat, Chiumenti, 2017, CM], [Wu, 2018, CMAME]
            a = 1.5*unit; % crack length
            b = 4.5*unit; % crack offset from the centerline
        otherwise
            error('Wrong setup');
    end
    L = 10*unit; % half-length
    h = 4*unit; % half-height
    ls = 9*unit; % location of the support from the centerline
    lh = 4*unit; % location of the holes from the centerline
    dh = 2*unit; % distance between the holes
    ph = 1.25*unit; % location of the top hole from the top
    r = 0.25*unit; % radius of the holes
    e = 1; % thickness
    % e = 0.5*unit; % [Cervera, Barbat, Chiumenti, 2017, CM], [Wu, Nguyen, 2018, JMPS], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
    
    % cl = 1e-3; % [mm] [Cervera, Barbat, Chiumenti, 2017, CM]
    % cl = 0.135*unit; % [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME]
    % cl = 0.066*unit; % [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME]
    % cl = 0.05e-3; % [mm] [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS]
    % cl = 0.033*unit; % [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME]
    % cl = 0.02*unit; % [Khisamitov, Meschke, 2018, CMAME] (setup 1), [Rahimi, Moutsanidis, 2022, CMAME]
    cl = 0.0125*unit; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME]
    % cl = 0.01*unit; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Molnar, Gravouil, 2017, FEAD], [Mandal, Nguyen, Wu, 2019, EFM]
    % cl = 0.005*unit; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Khisamitov, Meschke, 2018, CMAME] (setup 4), [Wu, 2018, CMAME], [Wu, Nguyen, 2018, JMPS], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
    % cl = 0.003*unit; % [Mandal, Nguyen, Wu, 2019, EFM]
    % cl = 0.0008*unit; % [Lee et al., 2024, CBM]
    if test
        cl = 0.05*unit;
    end
    switch lower(FEmesh)
        case 'unif' % uniform mesh
            clD = cl; % characteristic length for domain
            B = [];
        case 'optim' % optimized mesh
            clD = 0.1*unit; % characteristic length for domain
            % clD = 0.01*unit; % [Khisamitov, Meschke, 2018, CMAME] (setup 4)
            % clD = 2.54e-3; % [mm] [Cervera, Barbat, Chiumenti, 2017, CM]
            if test
                clD = 0.5*unit;
            end
            VIn = cl; VOut = clD;
            XMin = -b-0.5*unit; XMax = -lh+1*unit;
            YMin = -h; YMax = h;
            Thickness = ls-b-0.5*unit;
            % Thickness = 0;
            if Dim==2
                B = struct('VIn',VIn,'VOut',VOut,'XMin',XMin,'XMax',XMax,'YMin',YMin,'YMax',YMax,'Thickness',Thickness);
            elseif Dim==3
                ZMin = 0; ZMax = e;
                B = struct('VIn',VIn,'VOut',VOut,'XMin',XMin,'XMax',XMax,'YMin',YMin,'YMax',YMax,'ZMin',ZMin,'ZMax',ZMax,'Thickness',Thickness);
            end
        otherwise
            error('Wrong FE mesh')
    end
    clC = cl; % characteristic length for edge crack/notch
    clH = cl; % characteristic length for circular holes
    switch lower(initialCrack)
        case 'geometriccrack'
            S_phase = gmshAsymmetricPlateWithSingleEdgeCrackThreeHoles(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_asymmetric_notched_plate'),Dim,'Box',B);
        case 'geometricnotch'
            c = 0.05*unit; % crack width
            S_phase = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,fullfile(pathname,'gmsh_asymmetric_notched_plate'),Dim,'Box',B);
        case 'initialphasefield'
            S_phase = gmshAsymmetricPlateWithSingleEdgeCrackThreeHoles(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_asymmetric_notched_plate'),Dim,'noduplicate','refinecrack','Box',B);
        otherwise
            error('Wrong model for initial crack');
    end
    S = S_phase;
    
    %% Phase-field problem
    %% Material
    % Critical energy release rate (or fracture toughness)
    gc = 1e3;
    % gc = 2.7e3; % [Bhowmick, Liu, 2018, EFM]
    % gc = 500; % [Cervera, Barbat, Chiumenti, 2017, CM], [Lee et al., 2024, CBM]
    % gc = 315; % [Wu, 2018, CMAME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
    % gc = 304.321; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % Regularization parameter (width of the smeared crack)
    % l = 0.27*unit; % [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME]
    % l = 0.2e-3; % [mm] [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS]
    % l = 0.15e-3; % [mm] [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS]
    % l = 0.132*unit; % [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME]
    % l = 0.075*unit; % [Bhowmick, Liu, 2018, EFM], [Mandal, Nguyen, Wu, 2019, EFM]
    % l = 0.066*unit; % [Kirkesaether Brun, Wick, Berre, Nordbotten, Radu, 2020, CMAME]
    % l = 0.05*unit; % [Wu, Nguyen, 2018, JMPS], [Mandal, Nguyen, Wu, 2019, EFM]
    % l = 0.043*unit; % [Patil, Mishra, Singh, 2018, CMAME] (LMXPFM)
    % l = 0.04*unit; % [Patil, Mishra, Singh, 2018, CMAME] (AMsPFM)
    % l = 0.0375*unit; % [Mandal, Nguyen, Wu, 2019, EFM]
    % l = 0.03*unit; % [Lee et al., 2024, CBM]
    l = 0.025*unit; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Molnar, Gravouil, 2017, FEAD], [Wu, 2018, CMAME], [Wu, Nguyen, 2018, JMPS], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Min, Hu, Yao, Bui, Zhang, 2022, CMAME]
    % l = 0.02*unit; % [Rahimi, Moutsanidis, 2022, CMAME]
    % l = 0.0125*unit; % [Mandal, Nguyen, Wu, 2019, EFM]
    % l = 0.01*unit; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % Small artificial residual stiffness
    % k = 1e-12;
    k = 0;
    
    % Material
    [K,R,Qn] = setphasefieldparam(gc,l,PFregularization);
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'DIM3',e,'PFregularization',PFregularization);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    switch lower(initialCrack)
        case 'geometriccrack'
            C = POINT([-b,-h+a]); % crack tip
        case 'geometricnotch'
            C = CIRCLE(-b,-h+a-c/2,c/2); % circular notch
            % C = LINE([-b-c/2,-h+a],[-b+c/2,-h+a]); % rectangular notch
            % C = POINT([-b,-h+a]); % V notch
        case 'initialphasefield'
            C = LINE([-b,-h],[-b,-h+a]); % crack line
        otherwise
            error('Wrong model for initial crack');
    end
    R0 = 2*unit;
    BU = CIRCLE(0.0,h,R0);
    BL = CIRCLE(-ls,-h,R0);
    BR = CIRCLE(ls,-h,R0);
    % LU = LINE([-L,h],[L,h]);
    
    % findddlboundary = @(S_phase) findddl(S_phase,'T',LU);
    findddlboundary = @(S_phase) [];
    
    if strcmpi(initialCrack,'geometriccrack')
        S_phase = final(S_phase,'duplicate');
    else
        S_phase = final(S_phase);
    end
    
    S_phase = addbcdamageAsymmetricNotchedPlate(S_phase,C,BU,BL,BR,initialCrack);
    
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
    option = 'DEFO'; % plane strain [Miehe, Gurses, 2007, IJNME], [Guidault, Allix, Champaney, Cornuault, 2008, CMAME], [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Molnar, Gravouil, 2017, FEAD], [Khisamitov, Meschke, 2018, CMAME], [Bhowmick, Liu, 2018, EFM], [Patil, Mishra, Singh, 2018, CMAME] (AMsPFM), [Rahimi, Moutsanidis, 2022, CMAME]
    % option = 'CONT'; % plane stress [Passieux, Rethore, Gravouil, Baietto, 2013, CM], [Cervera, Barbat, Chiumenti, 2017, CM], [Wu, 2018, CMAME], [Wu, Nguyen, 2018, JMPS], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
    % Lame coefficients
    lambda = 12e9;
    mu = 8e9;
    % Young modulus and Poisson ratio
    switch lower(option)
        case 'defo'
            E = mu*(3*lambda+2*mu)/(lambda+mu); % E = 20.8e9;
            NU = lambda/(lambda+mu)/2; % NU = 0.3;
        case 'cont'
            E = 4*mu*(lambda+mu)/(lambda+2*mu);
            NU = lambda/(lambda+2*mu);
    end
    % E = 20e9; NU = 0.3; % [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS]
    % E = 3.102e9; NU = 0.35; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Cervera, Barbat, Chiumenti, 2017, CM]
    % E = 3.275e9; NU = 0.35; % [Ingraffea, Grigoriu, 1990], [Wu, 2018, CMAME], [Wu, Nguyen, 2018, JMPS], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
    % E = 200e9; NU = 0.3; % [Guidault, Allix, Champaney, Cornuault, 2008, CMAME]
    % E = 210e9; NU = 0.2; % [Bhowmick, Liu, 2018, EFM]
    % E = 12e9; NU = 0.2; % [Lee et al., 2024, CBM]
    % Energetic degradation function
    g = @(d) (1-d).^2;
    % Density
    RHO = 1;
    % RHO = 1190; % [Rahimi, Moutsanidis, 2022, CMAME]
    
    % Material
    d = calc_init_dirichlet(S_phase);
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    PU = POINT([0.0,h]);
    PL = POINT([-ls,-h]);
    PR = POINT([ls,-h]);
    
    addbc = @(S,ud) addbcAsymmetricNotchedPlate(S,ud,PU,PL,PR);
    findddlforce = @(S) findddl(S,'UY',PU);
    
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
    % [Wu, 2018, CMAME]
    % du = 1e-2 inches during 100 time steps (up to u = 1.0 inches)
    % dt = 1e-2*unit;
    % nt = 100;
    % t = linspace(dt,nt*dt,nt);
    
    % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
    % du = 1e-4 inches during 2500 time steps (up to u = 0.25 inches)
    % dt = 1e-4*unit;
    % nt = 2500;
    % t = linspace(dt,nt*dt,nt);
    
    % [Molnar, Gravouil, 2017, FEAD], [Bhowmick, Liu, 2018, EFM]
    % du = 1e-3 mm during the first 150 time steps (up to u = 0.15 mm)
    % du = 1e-4 mm during the last  1000 time steps (up to u = 0.25 mm)
    % dt0 = 1e-3*unit;
    % nt0 = 150;
    % dt1 = 1e-4*unit;
    % nt1 = 1000;
    % t = linspace(dt,nt*dt,nt);
    
    % [Ambati, Gerasimov, De Lorenzis, 2015, CM]
    % du = 1e-3 mm during the first 200 time steps (up to u = 0.2 mm)
    % du = 1e-4 mm during the last  500 time steps (up to u = 0.25 mm)
    dt0 = 1e-3*unit;
    nt0 = 200;
    dt1 = 1e-4*unit;
    nt1 = 500;
    if test
        dt0 = 2e-3*unit;
        nt0 = 100;
        dt1 = 2e-4*unit;
        nt1 = 250;
    end
    t0 = linspace(dt0,nt0*dt0,nt0);
    t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
    t = [t0,t1];
    
    % [Patil, Mishra, Singh, 2018, CMAME] (AMsPFM)
    % du = 1e-3 mm during the first 180 time steps (up to u = 0.18 mm)
    % du = 1e-4 mm during the last  100 time steps (up to u = 0.19 mm)
    % dt0 = 1e-3*unit;
    % nt0 = 180;
    % dt1 = 1e-4*unit;
    % nt1 = 100;
    % t0 = linspace(dt0,nt0*dt0,nt0);
    % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
    % t = [t0,t1];
    
    % [Patil, Mishra, Singh, 2018, CMAME] (LMXPFM)
    % du = 1e-3 mm during the first 170 time steps (up to u = 0.17 mm)
    % du = 1e-5 mm during the last 2000 time steps (up to u = 0.19 mm)
    % dt0 = 1e-3*unit;
    % nt0 = 170;
    % dt1 = 1e-5*unit;
    % nt1 = 2000;
    % t0 = linspace(dt0,nt0*dt0,nt0);
    % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
    % t = [t0,t1];
    
    T = TIMEMODEL(t);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'unit','T','S_phase','S','addbc','findddlforce','findddlboundary');
else
    load(fullfile(pathname,'problem.mat'),'unit','T','S_phase','S','addbc','findddlforce','findddlboundary');
end

%% Solution
if solveProblem
    tTotal = tic;
    
    switch lower(PFsolver)
        case {'historyfieldelem','historyfieldnode'}
            [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElas(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true);
        otherwise
            [dt,ut,ft,~,Edt,Eut,output] = solvePFDetLinElas(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true);
    end
    % switch lower(PFsolver)
    %     case {'historyfieldelem','historyfieldnode'}
    %         [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElasAsymmetricNotchedPlate(S_phase,S,T,PFsolver,PU,PL,PR,'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true);
    %     otherwise
    %         [dt,ut,ft,~,Edt,Eut,output] = solvePFDetLinElasAsymmetricNotchedPlate(S_phase,S,T,PFsolver,PU,PL,PR,'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true);
    % end
    
    t = gettevol(T);
    dt_val = getvalue(dt);
    dmaxt = max(dt_val);
    idc = find(dmaxt>=min(0.75,max(dmaxt)),1);
    fc = ft(idc);
    udc = t(idc);
    [fmax,idmax] = max(ft,[],2);
    udmax = t(idmax);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'dt','ut','ft','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        save(fullfile(pathname,'solution.mat'),'Ht','-append');
    end
else
    load(fullfile(pathname,'solution.mat'),'dt','ut','ft','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        load(fullfile(pathname,'solution.mat'),'Ht');
    end
end

%% Outputs
if solveProblem
    fid = fopen(fullfile(pathname,'results.txt'),'w');
    fprintf(fid,'Asymmetric notched plate\n');
    fprintf(fid,'\n');
    fprintf(fid,'setup    = %d\n',setup);
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
    
    fprintf(fid,'fmax  = %g kN/mm\n',fmax*1e-6);
    fprintf(fid,'fc    = %g kN/mm\n',fc*1e-6);
    fprintf(fid,'udmax = %g mm\n',udmax*1e3);
    fprintf(fid,'udc   = %g mm\n',udc*1e3);
    fclose(fid);
end

%% Display
if Dim==2
    facealpha = 0.1;
    facecolor = 'k';
    facecolordef = 'b';
elseif Dim==3
    facealpha = 1;
    facecolor = 'w';
    facecolordef = 'w';
end

if displayModel
    [t,rep] = gettevol(T);
    
    %% Display domains, boundary conditions and meshes
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 0.5;
    v = calc_init_dirichlet(S);
    [hN,legN] = vectorplot(S,'U',v,ampl,'r','LineWidth',1);
    % legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions_displacement',formats,renderer);
    
    [hD_phase,legD_phase] = plotBoundaryConditions(S_phase,'legend',false);
    % legend(hD_phase,legD_phase,'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions_damage',formats,renderer);
    
    % plotModel(S,'legend',false);
    % mysaveas(pathname,'mesh',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor',facecolor,'FaceAlpha',facealpha,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    u = getmatrixatstep(ut,rep(end));
    ampl = getsize(S)/max(abs(u))/20;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor',facecolordef,'FaceAlpha',facealpha,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor',facecolor,'FaceAlpha',facealpha);
    plot(S+ampl*unfreevector(S,u),'Color','b','FaceColor',facecolordef,'FaceAlpha',facealpha);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
end

%% Display solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display force-displacement curve
    figure('Name','Force vs displacement')
    clf
    plot(t*1e3,ft*1e-6,'-b','LineWidth',linewidth)
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
    
    %% Display energy-displacement curves
    figure('Name','Energies vs displacement')
    clf
    plot(t*1e3,Eut,'-b','LineWidth',linewidth)
    hold on
    plot(t*1e3,Edt,'-r','LineWidth',linewidth)
    plot(t*1e3,Eut+Edt,'-k','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Energy [J]','Interpreter',interpreter)
    legend('$\Psi_u$','$\Psi_c$','$\Psi_{\mathrm{tot}}$',...
        'Location','NorthWest','Interpreter','latex')
    mysaveas(pathname,'energies_displacement',formats);
    mymatlab2tikz(pathname,'energies_displacement.tex');
    
    %% Display outputs of iterative resolution
    figure('Name','Number of iterations vs displacement')
    clf
    plot(t*1e3,output.iteration,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Number of iterations','Interpreter',interpreter)
    mysaveas(pathname,'nb_iterations_displacement',formats);
    mymatlab2tikz(pathname,'nb_iterations_displacement.tex');
    
    figure('Name','Computing time vs displacement')
    clf
    plot(t*1e3,output.time,'-r','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Computing time [s]','Interpreter',interpreter)
    mysaveas(pathname,'cpu_time_displacement',formats);
    mymatlab2tikz(pathname,'cpu_time_displacement.tex');
    
    figure('Name','Error vs displacement')
    clf
    plot(t*1e3,output.error,'-k','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Error','Interpreter',interpreter)
    mysaveas(pathname,'error_displacement',formats);
    mymatlab2tikz(pathname,'error_displacement.tex');
    
    %% Display solutions at different instants
    ampl = 0;
    switch setup
        case {1,2}
            tSnapshots = [0.210 0.215 0.218 0.219 0.220 0.222]*unit;
        case {3,4,5,6}
            tSnapshots = [0.190 0.201 0.203 0.204 0.205 0.207]*unit;
        otherwise
            error('Wrong setup');
    end
    rep = arrayfun(@(x) find(t>x-eps,1),tSnapshots);
    rep = [rep,length(T)];
    % tSnapshots = [tSnapshots,gett1(T)];
    % rep = arrayfun(@(x) find(t>x-eps,1),tSnapshots);
    
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
        % plotSolution(S,uj,'energyint','local','ampl',ampl);
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
    framerate = 80;
    
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
    % evolSolution(S,ut,'energyint','local','ampl',ampl,'FrameRate',framerate,'filename','internal_energy_density','pathname',pathname,options{:});
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
