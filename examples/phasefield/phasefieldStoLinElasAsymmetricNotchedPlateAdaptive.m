%% Phase-field fracture model - stochastic linear elasticity problem   %%
%  Asymmetric notched plate with three holes under three-point bending %%
%%---------------------------------------------------------------------%%
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
% [Wu, 2018, CMAME] (hybrid isotropic-anisotropic phase-field model of Wu et al.)
% [Wu, Nguyen, 2018, JMPS] (hybrid isotropic-anisotropic phase-field model of Wu et al.)
% [Bhowmick, Liu, 2018, EFM] (anisotropic phase-field model of Miehe et al. + CS-FEM)
% [Mandal, Nguyen, Wu, 2019, EFM] (hybrid AT1, AT2 and PF-CZM)
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM] (anisotropic phase-field model of Wu et al.)
% [Fu, Yi, Chen, Bui, Hu, Yao, 2020, TAFM] (LEFM, FNM + SASE)
% [Min, Hu, Yao, Bui, Zhang, 2022, CMAME] (isotropic phase-field model with no split of Bourdin et al.)
% [Rahimi, Moutsanidis, 2022, CMAME] (anisotropic phase-field model with Total Lagrangian SPH approximation)
% [Lee, Meng, Ashkpour, Rahmaninezhad, Iqbal, Mishra, Hubler, Sales, Farnam, Najafi, 2024, CBM] (anisotropic phase-field model of Miehe et al. with healing and re-damage)

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
Dim = 2; % space dimension Dim = 2
setup = 1; % notch geometry setup = 1, 2, 3, 4, 5, 6
PFmodel = 'Miehe'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 1; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
critConv = 'Energy'; % 'Solution', 'Residual', 'Energy'
initialCrack = 'GeometricNotch'; % 'GeometricCrack', 'GeometricNotch', 'InitialPhaseField'

% Random model parameters
% N = 500; % number of samples
N = numWorkers;
randMat = struct('delta',0.2,'lcorr',1e-4); % random material parameters model
randPF = struct('aGc',0,'bGc',0,'lcorr',Inf); % random phase-field parameters model

suffix = '';

foldername = ['asymmetricNotchedPlateSetup' num2str(setup) '_' num2str(Dim) 'D'];
filename = ['linElas' PFmodel PFsplit PFregularization PFsolver initialCrack...
    'MaxIter' num2str(maxIter)];
if maxIter>1
    filename = [filename 'Tol' num2str(tolConv) num2str(critConv)];
end
filename = [filename 'MeshAdapt_' num2str(N) 'samples'];
% if any(randMat.delta)
%     filename = [filename '_RandMat_Delta' num2str(randMat.delta,'_%g') '_Lcorr' num2str(randMat.lcorr,'_%g')];
% end
% if any(randPF.aGc) && any(randPF.bGc)
%     gcbounds = [randPF.aGc(:),randPF.bGc(:)]';
%     filename = [filename '_RandPF_Gc' num2str(gcbounds(:)','_%g') '_Lcorr' num2str(randPF.lcorr,'_%g')];
% end
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

gmshoptions = '-v 0';
mmgoptions = '-nomove -hausd 0.000005 -hgrad 1.1 -v -1';
% gmshoptions = '-v 5';
% mmgoptions = '-nomove -hausd 0.01 -hgrad 1.3 -v 1';

%% Problem
if setProblem
    %% Domains and meshes
    unit = 1e-3; % in [mm] % [Guidault, Allix, Champaney, Cornuault, 2008, CMAME], [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME],
    % [Passieux, Rethore, Gravouil, Baietto, 2013, CM], [Molnar, Gravouil, 2017, FEAD], [Khisamitov, Meschke, 2018, CMAME], [Bhowmick, Liu, 2018, EFM], [Lee et al., 2024, CBM]
    % unit = 25.4e-3; % in [inch] % [Ingraffea, Grigoriu, 1990], [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM], [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME],
    % [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS], [Cervera, Barbat, Chiumenti, 2017, CM], [Wu, Nguyen, 2018, JMPS], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
    switch setup
        case 1 % [Ingraffea, Grigoriu, 1990], [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM], [Ventura, Xu, Belytschko, 2002, IJNME],
            % [Miehe, Gurses, 2007, IJNME], [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Häusler, Lindhorst, Horst, 2011, IJNME],
            % [Geniaut, Galenne, 2012, IJSS], [Passieux, Rethore, Gravouil, Baietto, 2013, CM], [Ambati, Gerasimov, De Lorenzis, 2015, CM],
            % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS], [Molnar, Gravouil, 2017, FEAD], [Cervera, Barbat, Chiumenti, 2017, CM],
            % [Khisamitov, Meschke, 2018, CMAME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Fu, Yi, Chen, Bui, Hu, Yao, 2020, TAFM],
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
    
    clD = 0.1*unit; % characteristic length for domain
    % clD = 0.01*unit; % [Khisamitov, Meschke, 2018, CMAME] (setup 4)
    % clD = 2.54e-3; % [mm] % [Cervera, Barbat, Chiumenti, 2017, CM]
    % cl = 1e-3; % [mm] % [Cervera, Barbat, Chiumenti, 2017, CM]
    % cl = 0.05e-3; % [mm] % [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS]
    % cl = 0.02*unit; % [Khisamitov, Meschke, 2018, CMAME] (setup 1), [Rahimi, Moutsanidis, 2022, CMAME]
    cl = 0.025*unit/2; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME]
    % cl = 0.01*unit; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Molnar, Gravouil, 2017, FEAD], [Mandal, Nguyen, Wu, 2019, EFM]
    % cl = 0.005*unit; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Khisamitov, Meschke, 2018, CMAME] (setup 4), [Wu, 2018, CMAME], [Wu, Nguyen, 2018, JMPS], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
    % cl = 0.003*nuit; % [Mandal, Nguyen, Wu, 2019, EFM]
    % cl = 0.0008*unit; % [Lee et al., 2024, CBM]
    if test
        clD = 0.25*unit;
        cl = 0.025*unit;
    % else
    %     clD = min(min(min(randMat.lcorr),min(randPF.lcorr))/4,clD);
    %     cl = min(min(min(randMat.lcorr),min(randPF.lcorr))/4,cl);
    end
    clC = cl; % characteristic length for edge crack/notch
    clH = cl; % characteristic length for circular holes
    switch lower(initialCrack)
        case 'geometriccrack'
            S_phase = gmshAsymmetricPlateWithSingleEdgeCrackThreeHoles(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_asymmetric_notched_plate'));
        case 'geometricnotch'
            c = 0.025*unit; % crack width
            S_phase = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,fullfile(pathname,'gmsh_asymmetric_notched_plate'));
        case 'initialphasefield'
            S_phase = gmshAsymmetricPlateWithSingleEdgeCrackThreeHoles(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_asymmetric_notched_plate'),Dim,'noduplicate');
        otherwise
            error('Wrong model for initial crack');
    end
    
    sizemap = @(d) (clC-clD)*d+clD;
    % sizemap = @(d) clD*clC./((clD-clC)*d+clC);
    
    %% Phase-field problem
    %% Material
    % Critical energy release rate (or fracture toughness)
    gc = 1e3;
    % gc = 2.7e3; % [Bhowmick, Liu, 2018, EFM]
    % gc = 500; % [Cervera, Barbat, Chiumenti, 2017, CM], [Lee et al., 2024, CBM]
    % gc = 315; % [Wu, 2018, CMAME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
    % gc = 304.321; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % Regularization parameter (width of the smeared crack)
    % l = 0.2e-3; % [mm] % [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS]
    % l = 0.15e-3; % [mm] % [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS]
    % l = 0.075*unit; % [Bhowmick, Liu, 2018, EFM], [Mandal, Nguyen, Wu, 2019, EFM]
    % l = 0.05*unit; % [Wu, Nguyen, 2018, JMPS], [Mandal, Nguyen, Wu, 2019, EFM]
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
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'DIM3',e,'PFregularization',PFregularization,'aGc',randPF.aGc,'bGc',randPF.bGc,'lcorr',randPF.lcorr);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    switch lower(initialCrack)
        case 'geometriccrack'
            C = POINT([-b,-h+a]);
        case 'geometricnotch'
            C = CIRCLE(-b,-h+a-c/2,c/2); % circular notch
            % C = LIGNE([-b-c/2,-h+a],[-b+c/2,-h+a]); % rectangular notch
            % C = POINT([-b,-h+a]); % V notch
        case 'initialphasefield'
            C = LIGNE([-b,-h],[-b,-h+a]);
        otherwise
            error('Wrong model for initial crack');
    end
    R0 = 2*unit;
    BU = CIRCLE(0.0,h,R0);
    BL = CIRCLE(-ls,-h,R0);
    BR = CIRCLE(ls,-h,R0);
    H1 = CIRCLE(-lh,h-ph-2*dh,r);
    H2 = CIRCLE(-lh,h-ph-dh,r);
    H3 = CIRCLE(-lh,h-ph,r);
    % LU = LIGNE([-L,h],[L,h]);
    
    addbcdamage = @(S_phase) addbcdamageAsymmetricNotchedPlate(S_phase,C,BU,BL,BR,initialCrack);
    addbcdamageadapt = @(S_phase) addbcdamageAsymmetricNotchedPlateAdaptive(S_phase,C,H1,H2,H3);
    % findddlboundary = @(S_phase) findddl(S_phase,'T',LU);
    findddlboundary = @(S_phase) [];
    
    if strcmpi(initialCrack,'geometriccrack')
        S_phase = final(S_phase,'duplicate');
    else
        S_phase = final(S_phase);
    end
    
    S_phase = addbcdamageadapt(S_phase);
    
    % [A_phase,b_phase] = calc_rigi(S_phase);
    % b_phase = -b_phase;
    % d = A_phase\b_phase;
    % d = unfreevector(S_phase,d);
    d = calc_init_dirichlet(S_phase);
    cl = sizemap(d);
    S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_asymmetric_notched_plate'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    S = S_phase;
    
    S_phase = setmaterial(S_phase,mat_phase);
    
    if strcmpi(initialCrack,'geometriccrack')
        S_phase = final(S_phase,'duplicate');
    else
        S_phase = final(S_phase);
    end
    
    S_phase = addbcdamage(S_phase);
    
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
    option = 'DEFO'; % plane strain [Miehe, Gurses, 2007, IJNME], [Guidault, Allix, Champaney, Cornuault, 2008, CMAME], [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Molnar, Gravouil, 2017, FEAD], [Khisamitov, Meschke, 2018, CMAME], [Bhowmick, Liu, 2018, EFM], [Rahimi, Moutsanidis, 2022, CMAME]
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
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit,'delta',randMat.delta,'lcorr',randMat.lcorr);
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
    
    T = TIMEMODEL(t);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'unit','T','S_phase','S','sizemap','addbc','addbcdamage','addbcdamageadapt','findddlforce','findddlboundary');
else
    load(fullfile(pathname,'problem.mat'),'unit','T','S_phase','S','sizemap','addbc','addbcdamage','addbcdamageadapt','findddlforce','findddlboundary');
end

%% Solution
if solveProblem
    myparallel('start',numWorkers);
    
    %% Solution
    tTotal = tic;
    
    nbSamples = 1;
    fun = @(S_phase,S,filename) solvePFDetLinElasAdaptive(S_phase,S,T,PFsolver,addbc,addbcdamage,addbcdamageadapt,findddlforce,findddlboundary,sizemap,...
        'maxiter',maxIter,'tol',tolConv,'crit',critConv,'filename',filename,'pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    % fun = @(S_phase,S,filename) solvePFDetLinElasAsymmetricNotchedPlateAdaptive(S_phase,S,T,PFsolver,C,BU,BL,BR,H1,H2,H3,PU,PL,PR,initialCrack,sizemap,...
    %     'maxiter',maxIter,'tol',tolConv,'crit',critConv,'filename',filename,'pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    [ft,dmaxt,dt,ut,St_phase,St] = solvePFStoLinElasAdaptive(S_phase,S,T,fun,N,'filename','gmsh_asymmetric_notched_plate','pathname',pathname,'nbsamples',nbSamples);
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
    [fmax_f,fmax_xi,fmax_bw] = ksdensity(fmax,'NumPoints',npts);
    [udmax_f,udmax_xi,udmax_bw] = ksdensity(udmax,'NumPoints',npts);
    [fc_f,fc_xi,fc_bw] = ksdensity(fc,'NumPoints',npts);
    [udc_f,udc_xi,udc_bw] = ksdensity(udc,'NumPoints',npts);
    
    save(fullfile(pathname,'solution.mat'),'N','ft','dmaxt','dt','ut','St_phase','St',...
        'ft_mean','ft_std','ft_ci','probs','time',...
        'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
        'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw',...
        'fc','fc_mean','fc_std','fc_ci','fc_f','fc_xi','fc_bw',...
        'udc','udc_mean','udc_std','udc_ci','udc_f','udc_xi','udc_bw');
else
    load(fullfile(pathname,'solution.mat'),'N','ft','dmaxt','dt','ut','St_phase','St',...
        'ft_mean','ft_std','ft_ci','probs','time',...
        'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
        'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw',...
        'fc','fc_mean','fc_std','fc_ci','fc_f','fc_xi','fc_bw',...
        'udc','udc_mean','udc_std','udc_ci','udc_f','udc_xi','udc_bw');
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
    fprintf(fid,'nb elements = %g (initial) - %g (final)\n',getnbelem(S),getnbelem(St{end}));
    fprintf(fid,'nb nodes    = %g (initial) - %g (final)\n',getnbnode(S),getnbnode(St{end}));
    fprintf(fid,'nb dofs     = %g (initial) - %g (final)\n',getnbddl(S),getnbddl(St{end}));
    fprintf(fid,'nb time dofs = %g\n',getnbtimedof(T));
    fprintf(fid,'nb samples = %g\n',N);
    fprintf(fid,'elapsed time = %f s\n',time);
    fprintf(fid,'\n');
    
    fprintf(fid,'mean(fmax)   = %g kN/mm\n',fmax_mean*1e-6);
    fprintf(fid,'std(fmax)    = %g kN/mm\n',fmax_std*1e-6);
    fprintf(fid,'disp(fmax)   = %g\n',fmax_std/fmax_mean);
    fprintf(fid,'%d%% ci(fmax) = [%g,%g] kN/mm\n',(probs(2)-probs(1))*100,fmax_ci(1)*1e-6,fmax_ci(2)*1e-6);
    fprintf(fid,'\n');
    
    fprintf(fid,'mean(fc)   = %g kN/mm\n',fc_mean*1e-6);
    fprintf(fid,'std(fc)    = %g kN/mm\n',fc_std*1e-6);
    fprintf(fid,'disp(fc)   = %g\n',fc_std/fc_mean);
    fprintf(fid,'%d%% ci(fc) = [%g,%g] kN/mm\n',(probs(2)-probs(1))*100,fc_ci(1)*1e-6,fc_ci(2)*1e-6);
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
end

%% Display
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
    % legend([hD_phase,hN_phase],[legD_phase,legN_phase],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions_damage',formats,renderer);
    
    % plotModel(S,'legend',false);
    % mysaveas(pathname,'mesh_init',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_init',formats,renderer);
    
    % u = ut(:,rep(end));
    u = ut(:,end);
    S_final = St(:,end);
    for k=1:numel(S_final)
        plotModel(S_final{k},'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
        mysaveas(pathname,['mesh_final_sample_' num2str(k)],formats,renderer);
        
        ampl = getsize(S_final{k})/max(abs(u{k}))/20;
        plotModelDeflection(S_final{k},u{k},'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
        mysaveas(pathname,['mesh_deflected_sample_' num2str(k)],formats,renderer);
        
        figure('Name','Meshes')
        clf
        plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
        plot(S_final{k}+ampl*unfreevector(S_final{k},u{k}),'Color','b','FaceColor','b','FaceAlpha',0.1);
        mysaveas(pathname,['meshes_deflected_sample_' num2str(k)],formats,renderer);
    end
end

%% Display statistics of solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display force-displacement curve
    figure('Name','Force vs displacement')
    clf
    plot(t*1e3,ft_mean*1e-6,'-b','LineWidth',linewidth)
    hold on
    ciplot(ft_ci(1,:)*1e-6,ft_ci(2,:)*1e-6,t*1e3,'b');
    alpha(0.2)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    legend('mean function',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'Location','NorthWest','Interpreter',interpreter)
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    colors = distinguishable_colors(N);
    figure('Name','Forces vs displacement')
    clf
    for i=1:N
        plot(t*1e3,ft(i,:)*1e-6,'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]'...,'Interpreter',interpreter...
        )
    ylabel('Force [kN]'...,'Interpreter',interpreter...
        )
    mysaveas(pathname,'forces_displacement',formats);
    mymatlab2tikz(pathname,'forces_displacement.tex');
    
    %% Display pdf of maximum force
    figure('Name','Probability Density Estimate: Maximum force')
    clf
    plot(fmax_xi*1e-6,fmax_f*1e6,'-r','LineWidth',linewidth)
    hold on
    ind_fmax = find(fmax_xi>=fmax_ci(1) & fmax_xi<fmax_ci(2));
    area(fmax_xi(ind_fmax)*1e-6,fmax_f(ind_fmax)*1e6,'FaceColor','r','EdgeColor','none','FaceAlpha',0.2)
    scatter(fmax_mean*1e-6,0,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$f$ [kN]','Interpreter',interpreter)
    ylabel('$p_{F_m}(f)$','Interpreter',interpreter)
    legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value',...
        'Direction','normal','Interpreter',interpreter)
    mysaveas(pathname,'pdf_fmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_fmax.tex');
    
    %% Display pdf of critical force
    figure('Name','Probability Density Estimate: Critical force')
    clf
    plot(fc_xi*1e-6,fc_f*1e6,'-b','LineWidth',linewidth)
    hold on
    ind_fc = find(fc_xi>=fc_ci(1) & fc_xi<fc_ci(2));
    area(fc_xi(ind_fc)*1e-6,fc_f(ind_fc)*1e6,'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(fc_mean*1e-6,0,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$f$ [kN]','Interpreter',interpreter)
    ylabel('$p_{F_c}(f)$','Interpreter',interpreter)
    legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value',...
        'Direction','normal','Interpreter',interpreter)
    mysaveas(pathname,'pdf_fc',formats,renderer);
    mymatlab2tikz(pathname,'pdf_fc.tex');
    
    %% Display pdf of maximum displacement
    figure('Name','Probability Density Estimate: Maximum displacement')
    clf
    plot(udmax_xi*1e3,udmax_f*1e-3,'-r','LineWidth',linewidth)
    hold on
    ind_udmax = find(udmax_xi>=udmax_ci(1) & udmax_xi<udmax_ci(2));
    area(udmax_xi(ind_udmax)*1e3,udmax_f(ind_udmax)*1e-3,'FaceColor','r','EdgeColor','none','FaceAlpha',0.2)
    scatter(udmax_mean*1e3,0,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$u$ [mm]','Interpreter',interpreter)
    ylabel('$p_{U_m}(u)$','Interpreter',interpreter)
    legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value',...
        'Direction','normal','Interpreter',interpreter)
    mysaveas(pathname,'pdf_udmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_udmax.tex');
    
    %% Display pdf of critical displacement
    figure('Name','Probability Density Estimate: Critical displacement')
    clf
    plot(udc_xi*1e3,udc_f*1e-3,'-b','LineWidth',linewidth)
    hold on
    ind_udc = find(udc_xi>=udc_ci(1) & udc_xi<udc_ci(2));
    area(udc_xi(ind_udc)*1e3,udc_f(ind_udc)*1e-3,'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(udc_mean*1e3,0,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$u$ [mm]','Interpreter',interpreter)
    ylabel('$p_{U_c}(u)$','Interpreter',interpreter)
    legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value',...
        'Direction','normal','Interpreter',interpreter)
    mysaveas(pathname,'pdf_udc',formats,renderer);
    mymatlab2tikz(pathname,'pdf_udc.tex');
    
    %% Display samples of solutions at different instants
    ampl = 0;
    switch setup
        case {1,2}
            tSnapshots = [0.210 0.215 0.218 0.219 0.220 0.222]*unit;
        case {3,4,5,6}
            tSnapshots = [0.190 0.201 0.203 0.204 0.205 0.207]*unit;
    end
    rep = arrayfun(@(x) find(t<x+eps,1,'last'),tSnapshots);
    rep = [rep,length(T)];
    % tSnapshots = [tSnapshots,gett1(T)];
    % rep = arrayfun(@(x) find(t<x+eps,1,'last'),tSnapshots);
    
    for k=1:size(St,1)
    for j=1:length(rep)
        dj = dt{k,rep(j)};
        uj = ut{k,rep(j)};
        Sj = St{k,rep(j)};
        Sj_phase = St_phase{k,rep(j)};
        
        plotModel(Sj,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
        mysaveas(pathname,['mesh_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        
        plotSolution(Sj_phase,dj);
        mysaveas(pathname,['damage_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        
        for i=1:Dim
            plotSolution(Sj,uj,'displ',i,'ampl',ampl);
            mysaveas(pathname,['displacement_' num2str(i) '_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        end
        
        % for i=1:(Dim*(Dim+1)/2)
        %     plotSolution(Sj,uj,'epsilon',i,'ampl',ampl);
        %     mysaveas(pathname,['epsilon_' num2str(i) '_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        %
        %     plotSolution(Sj,uj,'sigma',i,'ampl',ampl);
        %     mysaveas(pathname,['sigma_' num2str(i) '_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        % end
        %
        % plotSolution(Sj,uj,'epsilon','mises','ampl',ampl);
        % mysaveas(pathname,['epsilon_von_mises_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        %
        % plotSolution(Sj,uj,'sigma','mises','ampl',ampl);
        % mysaveas(pathname,['sigma_von_mises_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        %
        % plotSolution(Sj,uj,'energyint','','ampl',ampl);
        % mysaveas(pathname,['internal_energy_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
    end
    end
    
end

%% Display evolution of samples of solutions
if makeMovie
    ampl = 0;
    % umax = cellfun(@(u) max(abs(u)),ut,'UniformOutput',false);
    % ampl = getsize(S)/max([umax{:}])/20;
    
    options = {'plotiter',true,'plottime',false};
    framerate = 80;
    
    for k=1:size(St,1)
        evolModel(T,St(k,:),'FrameRate',framerate,'filename',['mesh_sample_' num2str(k)],'pathname',pathname,options{:});
        
        evolSolutionCell(T,St_phase(k,:),dt(k,:),'FrameRate',framerate,'filename',['damage_sample_' num2str(k)],'pathname',pathname,options{:});
        % for i=1:Dim
        %     evolSolutionCell(T,St(k,:),ut(k,:),'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i) '_sample_' num2str(k)],'pathname',pathname,options{:});
        % end
        %
        % for i=1:(Dim*(Dim+1)/2)
        %     evolSolutionCell(T,St(k,:),ut(k,:),'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['epsilon_' num2str(i) '_sample_' num2str(k)],'pathname',pathname,options{:});
        %     evolSolutionCell(T,St(k,:),ut(k,:),'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['sigma_' num2str(i) '_sample_' num2str(k)],'pathname',pathname,options{:});
        % end
        %
        % evolSolutionCell(T,St(k,:),ut(k,:),'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename',['epsilon_von_mises_sample_' num2str(k)],'pathname',pathname,options{:});
        % evolSolutionCell(T,St(k,:),ut(k,:),'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename',['sigma_von_mises_sample_' num2str(k)],'pathname',pathname,options{:});
        % evolSolutionCell(T,St(k,:),ut(k,:),'energyint','','ampl',ampl,'FrameRate',framerate,'filename',['internal_energy_sample_' num2str(k)],'pathname',pathname,options{:});
    end
end

%% Save samples of solutions
if saveParaview
    [t,rep] = gettevol(T);
    for k=1:size(St,1)
        for i=1:length(T)
            di = dt{k,rep(i)};
            ui = ut{k,rep(i)};
            Si = St{k,rep(i)};
            % Si_phase = St_phase{k,rep(i)};
            
            write_vtk_mesh(Si,{di,ui},[],...
                {'damage','displacement'},[],...
                pathname,['solution_sample_' num2str(k)],1,i-1);
        end
        make_pvd_file(pathname,['solution_sample_' num2str(k)],1,length(T));
    end
end
