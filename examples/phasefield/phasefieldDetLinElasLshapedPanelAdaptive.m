%% Phase-field fracture model - deterministic linear elasticity problem %%
%  L-shaped concrete notched panel under mixed-mode failure             %%
%%----------------------------------------------------------------------%%
% [Winkler, 2001, PhD thesis] (experimental tests)
% [Mosler, Meschke, 2004, CMAME]
% [Unger, Eckardt, Konke, 2007, CMAME] (XFEM)
% [Dumstorff, Meschke, 2007, IJNME] (XFEM)
% [Meschke, Dumstorff, 2007, CMAME] (XFEM)
% [Jager, Steinmann, Kuhl, 2008, IJNME] (CZM)
% [Zamani, Gracie, Eslami, 2012, IJNME] (XFEM, GFEM)
% [Ambati, Gerasimov, De Lorenzis, 2015, CM] (hybrid isotropic-anisotropic phase-field model of Ambati et al. compared with the anisotropic one of Miehe et al.)
% [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME] (isotropic phase-field model with no split of Bourdin et al. compared to experimental data of [Winkler, 2001, PhD thesis])
% [Gerasimov, De Lorenzis, 2016, CMAME] (anisotropic phase-field model of Miehe et al.)
% [Cervera, Barbat, Chiumenti, 2017, CM] (LEFM Mixed FEM)
% [Wick, 2017, SIAM JSC] (anisotropic phase-field model of Miehe et al.)
% [Zhang, Vignes, Sloan, Sheng, 2017, CM] (anisotropic phase-field model of Miehe et al.)
% [Wu, 2018, CMAME] (PF-CZM)
% [Kumar, Franckfort, Lopez-Pamies, 2018, JMPS] (elastomers undergoing large deformations)
% [Patil, Mishra, Singh, 2018, CMAME] (LMXPFM = anisotropic PFM of Miehe et.al + XFEM + MsFEM)
% [Gerasimov, De Lorenzis, 2019, CMAME] (anisotropic phase-field model of Amor et al. and Miehe et al.)
% [Hirshikesh, Jansari, Kannan, Annabattula, Natarajan, 2019, EFM] (hybrid isotropic-anisotropic phase-field model of Ambati et al.)
% [Kakouris, Triantafyllou, 2019, IJNME] (PF MPM)
% [Mandal, Nguyen, Wu, 2019, EFM] (hybrid AT1, AT2 and PF-CZM)
% [Geelen, 2020, PhD thesis] (anisotropic phase-field model of Miehe et al.)
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM] (anisotropic phase-field model of Wu et al.)
% [Mang, Wick, Wollner, 2020, CM] (anisotropic phase-field model of Miehe et al.)
% [Jodlbauer, Langer, Wick, 2020, CMAME] (anisotropic phase-field model of Miehe et al.)
% [Kirkesaether Brun, Wick, Inga Berre, Nordbotten, Radu, 2020, CMAME] (anisotropic phase-field model of Miehe et al.)
% [Wu, Huang, Nguyen, 2020, CMAME] (PF-CZM)
% [Yang, He, Liu, Deng, Huang, 2020, IJMS] (PD CZM)
% [Muixi, Marco, Rodriguez-Ferran, Fernandez-Mendez, 2021, CM] (PF-XFEM = XFEM + hybrid isotropic-anisotropic PFM of Ambati et al.)
% [Huang, Zhang, Li, Yang, Wu, Withers, 2021, EFM] (hybrid phase-field model of Wu and Cervera, 2018, IJSS]
% [Lampron, Therriault, Levesque, 2021, CMAME] (anisotropic phase-field model of Miehe et al.)
% [Si, Yu, Li, Natarajan, 2023, CMAME] (hybrid isotropic-anisotropic phase-field model of Ambati et al. with multi-patch adaptive isogeometric phase-field method based on Nitsche's method)
% [Bharali, Larsson, Janicke, 2024, CM] (anisotropic phase-field model of Wu et al.)
% [Hu, Tan, Xia, Min, Xu, Yao, Sun, Zhang, Quoc Bui, Zhuang, Rabczuk, 2023] (hybrid isotropic-anisotropic phase-field model of Ambati et al.)

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
loading = 'Monotonic'; % 'Monotonic' or 'Cyclic'
PFmodel = 'Miehe'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 1; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
critConv = 'Energy'; % 'Solution', 'Residual', 'Energy'

suffix = '';

foldername = ['LshapedPanel' num2str(loading) '_' num2str(Dim) 'D'];
filename = ['linElas' PFmodel PFsplit PFregularization PFsolver...
    'MaxIter' num2str(maxIter)];
if maxIter>1
    filename = [filename 'Tol' num2str(tolConv) num2str(critConv)];
end
filename = [filename 'MeshAdapt' suffix];

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

gmshoptions = '-v 0';
mmgoptions = '-nomove -hausd 0.01 -hgrad 1.1 -v -1';
% gmshoptions = '-v 5';
% mmgoptions = '-nomove -hausd 0.01 -hgrad 1.3 -v 1';

%% Problem
if setProblem
    %% Domains and meshes
    a = 250e-3; % half-length
    b = 30e-3; % distance of applied load from the right edge
    e = 100e-3; % thickness
    
    clD = 4*5e-3; % [Gerasimov, De Lorenzis, 2019, CMAME]
    % clD = 10e-3; % [Muixi, Marco, Rodriguez-Ferran, Fernandez-Mendez, 2021, CM]
    if Dim==2
        % cl = 14.577e-3; % [Mang, Wick, Wollner, 2020, CM], [Kirkesaether Brun, Wick, Inga Berre, Nordbotten, Radu, 2020, CMAME]
        % cl = 12e-3; % [Cervera, Barbat, Chiumenti, 2017, CM]
        % cl = 8.33e-3; % [Cervera, Barbat, Chiumenti, 2017, CM]
        % cl = 8e-3; % [Cervera, Barbat, Chiumenti, 2017, CM]
        % cl = 7.289e-3; % [Kirkesaether Brun, Wick, Inga Berre, Nordbotten, Radu, 2020, CMAME]
        % cl = 5e-3; % [Kakouris, Triantafyllou, 2019, IJNME]
        % cl = 4e-3; % [Cervera, Barbat, Chiumenti, 2017, CM]
        % cl = 2.5e-3; % [Kakouris, Triantafyllou, 2019, IJNME], [Gerasimov, De Lorenzis, 2019, CMAME]
        % cl = 2e-3; [Yang, He, Liu, Deng, Huang, 2020, IJMS]
        % cl = 1.822e-3; % [Wick, 2017, SIAM JSC]
        % cl = 1.776e-3; % [Zhang, Vignes, Sloan, Sheng, 2017, CM]
        % cl = 1/3*5e-3; % [Gerasimov, De Lorenzis, 2019, CMAME]
        % cl = 1.344e-3; % [Huang, Zhang, Li, Yang, Wu, Withers, 2021, EFM]
        % cl = 1.25e-3; % [Kakouris, Triantafyllou, 2019, IJNME], [Gerasimov, De Lorenzis, 2019, CMAME]
        cl = 1e-3; % [Wu, 2018, CMAME], [Wu, Huang, Nguyen, 2020, CMAME]
        % cl = 0.625e-3; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Kumar, Franckfort, Lopez-Pamies, 2018, JMPS], [Lampron, Therriault, Levesque, 2021, CMAME]
        % cl = 0.5e-3; % [Wu, 2018, CMAME], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Muixi, Marco, Rodriguez-Ferran, Fernandez-Mendez, 2021, CM], [Kakouris, Triantafyllou, 2019, IJNME], [Hu et al., 2023]
        % cl = 0.3125e-3; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
        % cl = 0.25e-3; % [Wu, 2018, CMAME], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
    elseif Dim==3
        cl = 2.5e-3; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
        % cl = 8e-3; % [Cervera, Barbat, Chiumenti, 2017, CM]
    end
    if test
        cl = 5e-3;
    end
    clC = cl; % characteristic length for crack
    S_phase = gmshLshapedpanel(a,b,e,clD,clC,fullfile(pathname,'gmsh_Lshaped_panel'));
    
    sizemap = @(d) (clC-clD)*d+clD;
    % sizemap = @(d) clD*clC./((clD-clC)*d+clC);
    
    %% Phase-field problem
    %% Material
    % Critical energy release rate (or fracture toughness)
    % gc = 890e3; % [Kakouris, Triantafyllou, 2019, IJNME]
    % gc = 160; % [Cervera, Barbat, Chiumenti, 2017, CM]
    % gc = 130; % [Unger, Eckardt, Konke, 2007, CMAME], [Zamani, Gracie, Eslami, 2012, IJNME], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Wu, Huang, Nguyen, 2020, CMAME], [Yang, He, Liu, Deng, Huang, 2020, IJMS], [Bharali, Larsson, Janicke, 2024, CM]
    % gc = 120; % [Huang, Zhang, Li, Yang, Wu, Withers, 2021, EFM]
    % gc = 95; % [Dumstorff, Meschke, 2007, IJNME], [Meschke, Dumstorff, 2007, CMAME], [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Zhang, Vignes, Sloan, Sheng, 2017, CM], [Kumar, Franckfort, Lopez-Pamies, 2018, JMPS], [Patil, Mishra, Singh, 2018, CMAME], [Hirshikesh, Jansari, Kannan, Annabattula, Natarajan, 2019, EFM], [Kirkesaether Brun, Wick, Inga Berre, Nordbotten, Radu, 2020, CMAME], [Lampron, Therriault, Levesque, 2021, CMAME], [Si, Yu, Li, Natarajan, 2023, CMAME]
    % gc = 90; % [Wu, 2018, CMAME], [Hu et al., 2023]
    gc = 89; % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Gerasimov, De Lorenzis, 2016, CMAME], [Wick, 2017, SIAM JSC], [Kakouris, Triantafyllou, 2019, IJNME], [Mang, Wick, Wollner, 2020, CM], [Jodlbauer, Langer, Wick, 2020, CMAME], [Geelen, 2020, PhD thesis]
    % gc = 65; % [Jager, Steinmann, Kuhl, 2008, IJNME]
    % Regularization parameter (width of the smeared crack)
    % l = 126.32e-3; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % l = 29.154e-3; % [Mang, Wick, Wollner, 2020, CM], [Kirkesaether Brun, Wick, Inga Berre, Nordbotten, Radu, 2020, CMAME]
    % l = 17.76e-3; % [Zhang, Vignes, Sloan, Sheng, 2017, CM]
    % l = 15e-3; % [Hu et al., 2023]
    % l = 14.577e-3; % [Kirkesaether Brun, Wick, Inga Berre, Nordbotten, Radu, 2020, CMAME]
    % l = 11e-3; % [Jodlbauer, Langer, Wick, 2020, CMAME]
    % l = 10e-3; % [Wu, 2018, CMAME], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Kakouris, Triantafyllou, 2019, IJNME], [Bharali, Larsson, Janicke, 2024, CM], [Hu et al., 2023]
    % l = 6.72e-3; % [Huang, Zhang, Li, Yang, Wu, Withers, 2021, EFM]
    l = 5e-3; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME] (3D), [Wu, 2018, CMAME], [Gerasimov, De Lorenzis, 2019, CMAME], [Kakouris, Triantafyllou, 2019, IJNME], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Wu, Huang, Nguyen, 2020, CMAME], [Hu et al., 2023]
    % l = 3.644e-3; % [Wick, 2017, SIAM JSC]
    % l = 3.125e-3; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Kumar, Franckfort, Lopez-Pamies, 2018, JMPS], [Lampron, Therriault, Levesque, 2021, CMAME]
    % l = 2.5e-3; % [Kakouris, Triantafyllou, 2019, IJNME], [Mandal, Nguyen, Wu, 2019, EFM], [Geelen, 2020, PhD thesis], [Muixi, Marco, Rodriguez-Ferran, Fernandez-Mendez, 2021, CM]
    % l = 2e-3; % [Zhang, Vignes, Sloan, Sheng, 2017, CM], [Patil, Mishra, Singh, 2018, CMAME]
    % l = 1.792e-3; % [Si, Yu, Li, Natarajan, 2023, CMAME]
    % l = 1.5625e-3; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % l = 1.1875e-3; % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Gerasimov, De Lorenzis, 2016, CMAME]
    % l = 1e-3; % [Kakouris, Triantafyllou, 2019, IJNME]
    % l = 0.2e-3; % [Hirshikesh, Jansari, Kannan, Annabattula, Natarajan, 2019, EFM]
    % Small artificial residual stiffness
    % k = 1e-12;
    k = 0;
    
    % Material
    [K,R,Qn] = setphasefieldparam(gc,l,PFregularization);
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'DIM3',e,'PFregularization',PFregularization);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BLeft = LIGNE([-a,-a],[-a,a]);
        BR = DOMAIN(2,[a-2*b,0.0],[a,a]);
    elseif Dim==3
        BLeft = PLAN([-a,-a,0.0],[-a,a,0.0],[-a,-a,e]);
        BR = DOMAIN(3,[a-2*b,0.0,0.0],[a,a,e]);
    end
    
    addbcdamage = @(S_phase) addcl(S_phase,BR,'T');
    addbcdamageadapt = @(S_phase) S_phase;
    findddlboundary = @(S_phase) findddl(S_phase,'T',BLeft);
    
    S_phase = final(S_phase);
    
    % [A_phase,b_phase] = calc_rigi(S_phase);
    % b_phase = -b_phase;
    % d = A_phase\b_phase;
    % d = unfreevector(S_phase,d);
    d = calc_init_dirichlet(S_phase);
    cl = sizemap(d);
    S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_Lshaped_panel'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    S = S_phase;
    
    S_phase = setmaterial(S_phase,mat_phase);
    S_phase = final(S_phase);
    
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
    % option = 'DEFO'; % plane strain [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Gerasimov, De Lorenzis, 2016, CMAME], [Kumar, Franckfort, Lopez-Pamies, 2018, JMPS], [Gerasimov, De Lorenzis, 2019, CMAME], [Geelen, 2020, PhD thesis], [Muixi, Marco, Rodriguez-Ferran, Fernandez-Mendez, 2021, CM], [Lampron, Therriault, Levesque, 2021, CMAME]
    option = 'CONT'; % plane stress [Dumstorff, Meschke, 2007, IJNME], [Meschke, Dumstorff, 2007, CMAME], [Zamani, Gracie, Eslami, 2012, IJNME], [Cervera, Barbat, Chiumenti, 2017, CM], [Hirshikesh, Jansari, Kannan, Annabattula, Natarajan, 2019, EFM], [Kakouris, Triantafyllou, 2019, IJNME], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Huang, Nguyen, 2020, CMAME], [Yang, He, Liu, Deng, Huang, 2020, IJMS], [Huang, Zhang, Li, Yang, Wu, Withers, 2021, EFM]
    % Lame coefficients
    % lambda = 6.161e9; % [Jager, Steinmann, Kuhl, 2008, IJNME]
    % mu = 10.953e9; % [Jager, Steinmann, Kuhl, 2008, IJNME]
    % lambda = 6.16e9; % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Gerasimov, De Lorenzis, 2016, CMAME], [Wick, 2017, SIAM JSC], [Kumar, Franckfort, Lopez-Pamies, 2018, JMPS], [Geelen, 2020, PhD thesis], [Mang, Wick, Wollner, 2020, CM], [Jodlbauer, Langer, Wick, 2020, CMAME]
    % mu = 10.95e9; % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Gerasimov, De Lorenzis, 2016, CMAME], [Wick, 2017, SIAM JSC], [Kumar, Franckfort, Lopez-Pamies, 2018, JMPS], [Geelen, 2020, PhD thesis], [Mang, Wick, Wollner, 2020, CM], [Jodlbauer, Langer, Wick, 2020, CMAME]
    % Young modulus and Poisson ratio
    % switch lower(option)
    %     case 'defo'
    %         E = mu*(3*lambda+2*mu)/(lambda+mu); % E = 25.85e9;
    %         NU = lambda/(lambda+mu)/2; % NU = 0.18;
    %     case 'cont'
    %         E = 4*mu*(lambda+mu)/(lambda+2*mu);
    %         NU = lambda/(lambda+2*mu);
    % end
    % E = 20e9; % [Unger, Eckardt, Konke, 2007, CMAME], [Zamani, Gracie, Eslami, 2012, IJNME], [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Wu, Huang, Nguyen, 2020, CMAME], [Yang, He, Liu, Deng, Huang, 2020, IJMS], [Huang, Zhang, Li, Yang, Wu, Withers, 2021, EFM], [Bharali, Larsson, Janicke, 2024, CM]
    % E = 25.84e9; % [Muixi, Marco, Rodriguez-Ferran, Fernandez-Mendez, 2021, CM]
    E = 25.85e9; % [Dumstorff, Meschke, 2007, IJNME], [Meschke, Dumstorff, 2007, CMAME], [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Cervera, Barbat, Chiumenti, 2017, CM], [Zhang, Vignes, Sloan, Sheng, 2017, CM], [Wu, 2018, CMAME], [Patil, Mishra, Singh, 2018, CMAME], [Gerasimov, De Lorenzis, 2019, CMAME], [Hirshikesh, Jansari, Kannan, Annabattula, Natarajan, 2019, EFM], [Kirkesaether Brun, Wick, Inga Berre, Nordbotten, Radu, 2020, CMAME], [Lampron, Therriault, Levesque, 2021, CMAME], [Si, Yu, Li, Natarajan, 2023, CMAME]
    NU = 0.18; % [Unger, Eckardt, Konke, 2007, CMAME], [Dumstorff, Meschke, 2007, IJNME], [Meschke, Dumstorff, 2007, CMAME], [Zamani, Gracie, Eslami, 2012, IJNME], [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Cervera, Barbat, Chiumenti, 2017, CM], [Zhang, Vignes, Sloan, Sheng, 2017, CM], [Wu, 2018, CMAME], [Patil, Mishra, Singh, 2018, CMAME], [Gerasimov, De Lorenzis, 2019, CMAME], [Hirshikesh, Jansari, Kannan, Annabattula, Natarajan, 2019, EFM],
    % [Mandal, Nguyen, Wu, 2019, EFM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Wu, Huang, Nguyen, 2020, CMAME], [Kirkesaether Brun, Wick, Inga Berre, Nordbotten, Radu, 2020, CMAME], [Mang, Wick, Wollner, 2020, CM], [Yang, He, Liu, Deng, Huang, 2020, IJMS], [Muixi, Marco, Rodriguez-Ferran, Fernandez-Mendez, 2021, CM], [Huang, Zhang, Li, Yang, Wu, Withers, 2021, EFM], [Lampron, Therriault, Levesque, 2021, CMAME], [Si, Yu, Li, Natarajan, 2023, CMAME], [Bharali, Larsson, Janicke, 2024, CM]
    % NU = 0.4999; % [Mang, Wick, Wollner, 2020, CM]
    % Energetic degradation function
    g = @(d) (1-d).^2;
    % Density
    RHO = 1;
    
    % Material
    d = calc_init_dirichlet(S_phase);
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BL = LIGNE([-a,-a],[0.0,-a]);
        BRight = LIGNE([a-b,0.0],[a,0.0]);
        % BRight = POINT([a-b,0.0]);
    elseif Dim==3
        BL = PLAN([-a,-a,0.0],[0.0,-a,0.0],[-a,-a,e]);
        BRight = QUADRANGLE([a-b,0.0,0.0],[a,0.0,0.0],[a,0.0,e],[a-b,0.0,e]);
        % BRight = LIGNE([a-b,0.0,0.0],[a-b,0.0,e]);
    end
    
    addbc = @(S,ud) addbcLshapedPanel(S,ud,BL,BRight);
    findddlforce = @(S) findddl(S,'UY',BRight);
    
    S = final(S);
    
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
    switch lower(loading)
        case 'monotonic'
            % [Patil, Mishra, Singh, 2018, CMAME]
            % du = 5e-4 mm during 2000 time steps (up to u = 1.0 mm)
            % dt = 5e-7;
            % nt = 2000;
            % t = linspace(dt,nt*dt,nt);
            
            % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Kirkesaether Brun, Wick, Inga Berre, Nordbotten, Radu, 2020, CMAME], [Bharali, Larsson, Janicke, 2024, CM]
            % du = 1e-3 mm during 800 time steps (up to u = 0.8 mm)
            % dt = 1e-6;
            % nt = 800;
            % t = linspace(dt,nt*dt,nt);
            
            % [Muixi, Marco, Rodriguez-Ferran, Fernandez-Mendez, 2021, CM], [Si, Yu, Li, Natarajan, 2023, CMAME]
            % du = 1e-3 mm during 600 time steps (up to u = 0.6 mm)
            % dt = 1e-6;
            % nt = 600;
            % t = linspace(dt,nt*dt,nt);
            
            % [Zhang, Vignes, Sloan, Sheng, 2017, CM]
            % du = 2e-3 mm during 400 time steps (up to u = 0.8 mm)
            % dt = 2e-6;
            % nt = 400;
            % t = linspace(dt,nt*dt,nt);
            
            % [Hirshikesh, Jansari, Kannan, Annabattula, Natarajan, 2019, EFM]
            % du = 2e-3 mm during 300 time steps (up to u = 0.6 mm)
            % dt = 2e-6;
            % nt = 300;
            % t = linspace(dt,nt*dt,nt);
            
            % [Wu, Huang, Nguyen, 2020, CMAME], [Lampron, Therriault, Levesque, 2021, CMAME]
            % du = 2e-2 mm during 50 time steps (up to u = 1 mm)
            % dt = 2e-5;
            % nt = 50;
            % t = linspace(dt,nt*dt,nt);
            
            % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Mang, Wick, Wollner, 2020, CM]
            % du = 1e-3 mm during 1000 time steps (up to u = 1 mm)
            dt = 1e-6;
            nt = 1000;
            if test
                dt = 1e-5;
                nt = 100;
            end
            t = linspace(dt,nt*dt,nt);
        case 'cyclic'
            % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Wick, 2017, SIAM JSC], [Kakouris, Triantafyllou, 2019, IJNME], [Mang, Wick, Wollner, 2020, CM], [Jodlbauer, Langer, Wick, 2020, CMAME]
            % du = 1e-3 mm during the first 300 time steps (up to u = 0.3 mm)
            % du = -1e-3 mm during the next 500 time steps (down to u = -0.2 mm)
            % du = 1e-3 mm during the last 1200 time steps (up to u = 1 mm)
            dt = 1e-6;
            nt0 = 300;
            nt1 = 500;
            nt2 = 1200;
            if test
                dt = 1e-5;
                nt0 = 30;
                nt1 = 50;
                nt2 = 120;
            end
            t0 = linspace(dt,nt0*dt,nt0);
            t1 = linspace(t0(end)-dt,t0(end)-nt1*dt,nt1);
            t2 = linspace(t1(end)+dt,t1(end)+nt2*dt,nt2);
            t = [t0,t1,t2];
            
            % [Gerasimov, De Lorenzis, 2019, CMAME]
            % du = 1e-2 mm during the first 36 time steps (up to u = 0.36 mm)
            % du = -3*1e-2 mm during the last 11 time steps (down to u = 0.03 mm)
            % dt0 = 1e-5;
            % nt0 = 36;
            % dt1 = -3*1e-5;
            % nt1 = 11;
            % t0 = linspace(dt0,nt0*dt0,nt0);
            % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
            % t = [t0,t1];
        otherwise
            error('Wrong loading case');
    end
    
    T = TIMEMODEL(t);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','addbc','addbcdamage','addbcdamageadapt','findddlforce','findddlboundary');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','addbc','addbcdamage','addbcdamageadapt','findddlforce','findddlboundary');
end

%% Solution
if solveProblem
    tTotal = tic;
    
    switch lower(PFsolver)
        case {'historyfieldelem','historyfieldnode'}
            [dt,ut,ft,St_phase,St,Ht,Edt,Eut,output] = solvePFDetLinElasAdaptive(S_phase,S,T,PFsolver,addbc,addbcdamage,addbcdamageadapt,findddlforce,findddlboundary,sizemap,...
                'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true,'filename','gmsh_Lshaped_panel','pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
        otherwise
            [dt,ut,ft,St_phase,St,~,Edt,Eut,output] = solvePFDetLinElasAdaptive(S_phase,S,T,PFsolver,addbc,addbcdamage,addbcdamageadapt,findddlforce,findddlboundary,sizemap,...
                'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true,'filename','gmsh_Lshaped_panel','pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    end
    % switch lower(PFsolver)
    %     case {'historyfieldelem','historyfieldnode'}
    %         [dt,ut,ft,St_phase,St,Ht,Edt,Eut,output] = solvePFDetLinElasLshapedPanelAdaptive(S_phase,S,T,PFsolver,BR,BL,BRight,BLeft,sizemap,...
    %             'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true,'filename','gmsh_Lshaped_panel','pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    %     otherwise
    %         [dt,ut,ft,St_phase,St,~,Edt,Eut,output] = solvePFDetLinElasLshapedPanelAdaptive(S_phase,S,T,PFsolver,BR,BL,BRight,BLeft,sizemap,...
    %             'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true,'filename','gmsh_Lshaped_panel','pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    % end
    
    t = gettevol(T);
    dmaxt = cellfun(@(d) max(d),dt);
    idc = find(dmaxt>=min(0.75,max(dmaxt)),1);
    fc = ft(idc);
    udc = t(idc);
    [fmax,idmax] = max(ft,[],2);
    udmax = t(idmax);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'dt','ut','ft','T','St_phase','St','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        save(fullfile(pathname,'solution.mat'),'Ht','-append');
    end
else
    load(fullfile(pathname,'solution.mat'),'dt','ut','ft','T','St_phase','St','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        load(fullfile(pathname,'solution.mat'),'Ht');
    end
end

%% Outputs
if solveProblem
    fid = fopen(fullfile(pathname,'results.txt'),'w');
    fprintf(fid,'L-shaped panel\n');
    fprintf(fid,'\n');
    fprintf(fid,'dim      = %d\n',Dim);
    fprintf(fid,'loading  = %s\n',loading);
    fprintf(fid,'PF model = %s\n',PFmodel);
    fprintf(fid,'PF split = %s\n',PFsplit);
    fprintf(fid,'PF regularization = %s\n',PFregularization);
    fprintf(fid,'PF solver = %s\n',PFsolver);
    fprintf(fid,'nb elements = %g (initial) - %g (final)\n',getnbelem(S),getnbelem(St{end}));
    fprintf(fid,'nb nodes    = %g (initial) - %g (final)\n',getnbnode(S),getnbnode(St{end}));
    fprintf(fid,'nb dofs     = %g (initial) - %g (final)\n',getnbddl(S),getnbddl(St{end}));
    fprintf(fid,'nb time dofs = %g\n',getnbtimedof(T));
    fprintf(fid,'elapsed time = %f s\n',time);
    fprintf(fid,'\n');
    
    if Dim==2
        fprintf(fid,'fmax  = %g kN/m\n',fmax*1e-3);
        fprintf(fid,'fc    = %g kN/m\n',fc*1e-3);
    elseif Dim==3
        fprintf(fid,'fmax  = %g kN\n',fmax*1e-3);
        fprintf(fid,'fc    = %g kN\n',fc*1e-3);
    end
    fprintf(fid,'udmax = %g mm\n',udmax*1e3);
    fprintf(fid,'udc   = %g mm\n',udc*1e3);
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
    
    % u = ut{rep(end)};
    u = ut{end};
    S_final = St{end};
    
    plotModel(S_final,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_final',formats,renderer);
    
    ampl = getsize(S_final)/max(abs(u))/20;
    plotModelDeflection(S_final,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S_final+ampl*unfreevector(S_final,u),'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
end

%% Display solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display displacement-loading step curve
    figure('Name','Displacement vs loading step')
    clf
    plot(1:length(t),t*1e3,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Loading step','Interpreter',interpreter)
    ylabel('Displacement [mm]','Interpreter',interpreter)
    mysaveas(pathname,'displacement_loading_step',formats);
    mymatlab2tikz(pathname,'displacement_loading_step.tex');
    
    %% Display force-displacement curve
    figure('Name','Force vs displacement')
    clf
    plot(t*1e3,ft*1e-3,'-b','LineWidth',linewidth)
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
    switch lower(loading)
        case 'monotonic'
            tSnapshots = [0.2 0.4 0.6 0.8 1.0]*1e-3;
        case 'cyclic'
            % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Wick, 2017, SIAM JSC], [Jodlbauer, Langer, Wick, 2020, CMAME]
            tSnapshots = [0.22 0.30 0.45 1.0]*1e-3;
            % [Kakouris, Triantafyllou, 2019, IJNME]
            % tSnapshots = [0.27 0.30 0.45 1.0]*1e-3;
    end
    rep = arrayfun(@(x) find(t<x+eps,1,'last'),tSnapshots);
    rep = [rep,length(T)];
    % tSnapshots = [tSnapshots,gett1(T)];
    % rep = arrayfun(@(x) find(t<x+eps,1,'last'),tSnapshots);
    
    for j=1:length(rep)
        dj = dt{rep(j)};
        uj = ut{rep(j)};
        Sj = St{rep(j)};
        Sj_phase = St_phase{rep(j)};
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            Hj = Ht{rep(j)};
        end
        
        plotModel(Sj,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
        mysaveas(pathname,['mesh_t' num2str(rep(j))],formats,renderer);
        
        plotSolution(Sj_phase,dj);
        mysaveas(pathname,['damage_t' num2str(rep(j))],formats,renderer);
        
        for i=1:Dim
            plotSolution(Sj,uj,'displ',i,'ampl',ampl);
            mysaveas(pathname,['displacement_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        end
        
        % for i=1:(Dim*(Dim+1)/2)
        %     plotSolution(Sj,uj,'epsilon',i,'ampl',ampl);
        %     mysaveas(pathname,['epsilon_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        %
        %     plotSolution(Sj,uj,'sigma',i,'ampl',ampl);
        %     mysaveas(pathname,['sigma_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        % end
        %
        % plotSolution(Sj,uj,'epsilon','mises','ampl',ampl);
        % mysaveas(pathname,['epsilon_von_mises_t' num2str(rep(j))],formats,renderer);
        %
        % plotSolution(Sj,uj,'sigma','mises','ampl',ampl);
        % mysaveas(pathname,['sigma_von_mises_t' num2str(rep(j))],formats,renderer);
        %
        % plotSolution(Sj,uj,'energyint','','ampl',ampl);
        % mysaveas(pathname,['internal_energy_density_t' num2str(rep(j))],formats,renderer);
        %
        % if strcmpi(PFsolver,'historyfieldelem')
        %     figure('Name','Solution H')
        %     clf
        %     plot(Hj,Sj_phase);
        %     colorbar
        %     set(gca,'FontSize',fontsize)
        %     mysaveas(pathname,['internal_energy_density_history_t' num2str(rep(j))],formats,renderer);
        % elseif strcmpi(PFsolver,'historyfieldnode')
        %     plotSolution(Sj_phase,Hj,'ampl',ampl);
        %     mysaveas(pathname,['internal_energy_density_history_t' num2str(rep(j))],formats,renderer);
        % end
    end
end

%% Display evolution of solutions
if makeMovie
    ampl = 0;
    % umax = cellfun(@(u) max(abs(u)),ut,'UniformOutput',false);
    % ampl = getsize(S)/max([umax{:}])/20;
    
    options = {'plotiter',true,'plottime',false};
    framerate = 80;
    
    evolModel(T,St,'FrameRate',framerate,'filename','mesh','pathname',pathname,options{:});
    
    evolSolutionCell(T,St_phase,dt,'FrameRate',framerate,'filename','damage','pathname',pathname,options{:});
    % for i=1:Dim
    %     evolSolutionCell(T,St,ut,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i)],'pathname',pathname,options{:});
    % end
    %
    % for i=1:(Dim*(Dim+1)/2)
    %     evolSolutionCell(T,St,ut,'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['epsilon_' num2str(i)],'pathname',pathname,options{:});
    %     evolSolutionCell(T,St,ut,'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['sigma_' num2str(i)],'pathname',pathname,options{:});
    % end
    %
    % evolSolutionCell(T,St,ut,'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename','epsilon_von_mises','pathname',pathname,options{:});
    % evolSolutionCell(T,St,ut,'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename','sigma_von_mises','pathname',pathname,options{:});
    % evolSolutionCell(T,St,ut,'energyint','','ampl',ampl,'FrameRate',framerate,'filename','internal_energy_density','pathname',pathname,options{:});
    % if strcmpi(PFsolver,'historyfieldnode')
    %     evolSolutionCell(T,St_phase,Ht,'ampl',ampl,'FrameRate',framerate,'filename','internal_energy_density_history','pathname',pathname,options{:});
    % end
end

%% Save solutions
if saveParaview
    [t,rep] = gettevol(T);
    for i=1:length(T)
        di = dt{rep(i)};
        ui = ut{rep(i)};
        Si = St{rep(i)};
        % Si_phase = St_phase{rep(i)};
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            Hi = Ht{rep(i)};
        end
        
        switch lower(PFsolver)
            case 'historyfieldelem'
                write_vtk_mesh(Si,{di,ui},{Hi},...
                    {'damage','displacement'},{'internal energy density history'},...
                    pathname,'solution',1,i-1);
            case 'historyfieldnode'
                write_vtk_mesh(Si,{di,ui,Hi},[],...
                    {'damage','displacement','internal energy density history'},[],...
                    pathname,'solution',1,i-1);
            otherwise
                write_vtk_mesh(Si,{di,ui},[],...
                    {'damage','displacement'},[],...
                    pathname,'solution',1,i-1);
        end
    end
    make_pvd_file(pathname,'solution',1,length(T));
end

% myparallel('stop');
