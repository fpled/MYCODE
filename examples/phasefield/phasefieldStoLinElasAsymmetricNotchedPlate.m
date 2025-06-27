%% Phase-field fracture model - stochastic linear elasticity problem   %%
%  Asymmetric notched plate with three holes under three-point bending %%
%%---------------------------------------------------------------------%%
% [Ingraffea, Grigoriu, 1990] (experimental tests)
% [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM] (SIF-based method with local remeshing and special FE)
% [Ventura, Xu, Belytschko, 2002, IJNME] (vector level set method with discontinuous enrichment in meshless method)
% [Guidault, Allix, Champaney, Cornuault, 2008, CMAME] (MsXFEM)
% [Miehe, Welschinger, Hofacker, 2010, IJNME] (anisotropic phase-field model of Miehe et al.)
% [Miehe, Hofacker, Welschinger, 2010, CMAME] (anisotropic phase-field model of Miehe et al.)
% [Häusler, Lindhorst, Horst, 2011, IJNME] (XFEM)
% [Geniaut, Galenne, 2012, IJSS] (XFEM)
% [Passieux, Rethore, Gravouil, Baietto, 2013, CM] (XFEM)
% [Ambati, Gerasimov, De Lorenzis, 2015, CM] (hybrid isotropic-anisotropic phase-field model of Ambati et al. compared with the anisotropic one of Miehe et al.)
% [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME] (isotropic phase-field model with no split of Bourdin et al. compared to experimental data of [Winkler, 2001, PhD thesis])
% [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS] (isotropic phase-field model with no split of Bourdin et al.)
% [Molnar, Gravouil, 2015, FEAD] (isotropic phase-field model with no split of Bourdin et al.)
% [Khisamitov, Meschke, 2018, CMAME] (anisotropic phase-field model for interfacial elastic energy)
% [Wu, Nguyen, 2018, JMPS] (hybrid isotropic-anisotropic phase-field model of Wu et al.)
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM] (anisotropic phase-field model of Wu et al.)
% [Rahimi, Moutsanidis, 2022, CMAME] (anisotropic phase-field model with Total Lagrangian SPH approximation)

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
setup = 2; % notch geometry setup = 1, 2, 3, 4, 5
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
randPF = struct('aGc',0,'bGc',0,'lcorr',Inf); % random phase-field parameters model

suffix = '';

foldername = ['asymmetricNotchedPlateSetup' num2str(setup) '_' num2str(Dim) 'D'];
filename = ['linElas' PFmodel PFsplit PFregularization PFsolver initialCrack...
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
    unit = 1e-3; % for mm % [Guidault, Allix, Champaney, Cornuault, 2008, CMAME], [Miehe, Welschinger, Hofacker, 2010, IJNME],
    % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Passieux, Rethore, Gravouil, Baietto, 2013, CM], [Molnar, Gravouil, 2015, FEAD], [Khisamitov, Meschke, 2018, CMAME]
    % unit = 25.4e-3; % for inch % [Ingraffea, Grigoriu, 1990], [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM],
    % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    switch setup
        case 1 % [Ingraffea, Grigoriu, 1990], [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM], [Ventura, Xu, Belytschko, 2002, IJNME],
            % [Guidault, Allix, Champaney, Cornuault, 2008, CMAME], [Häusler, Lindhorst, Horst, 2011, IJNME], [Geniaut, Galenne,2012, IJSS],
            % [Passieux, Rethore, Gravouil, Baietto, 2013, CM], [Molnar, Gravouil, 2015, FEAD], [Khisamitov, Meschke, 2018, CMAME]
            a = 1.5*unit; % crack length
            b = 5*unit; % crack offset from the centerline
        case 2 % [Ingraffea, Grigoriu, 1990], [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM], [Ventura, Xu, Belytschko, 2002, IJNME],
            % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Häusler, Lindhorst, Horst, 2011, IJNME],
            % [Geniaut, Galenne, 2012, IJSS], [Passieux, Rethore, Gravouil, Baietto, 2013, CM], [Molnar, Gravouil, 2015, FEAD], [Ambati, Gerasimov, De Lorenzis, 2015, CM],
            % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS], [Khisamitov, Meschke, 2018, CMAME],
            % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM], [Rahimi, Moutsanidis, 2022, CMAME]
            a = 1*unit; % crack length
            b = 6*unit; % crack offset from the centerline
        case 3 % [Ingraffea, Grigoriu, 1990], [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
            a = 2.5*unit; % crack length
            b = 6*unit; % crack offset from the centerline
        case 4 % [Ingraffea, Grigoriu, 1990], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
            a = 1.5*unit; % crack length
            b = 4.75*unit; % crack offset from the centerline
        case 5 % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
            a = 1.5*unit; % crack length
            b = 5.15*unit; % crack offset from the centerline
    end
    L = 10*unit; % half-length
    h = 4*unit; % half-height
    ls = 9*unit; % location of the support from the centerline
    lh = 4*unit; % location of the holes from the centerline
    dh = 2*unit; % distance between the holes
    ph = 1.25*unit; % location of the top hole from the top
    r = 0.25*unit; % radius of the holes
    e = 1; % thickness
    % e = 0.5*unit; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    
    % cl = 0.05e-3; % [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS]
    % cl = 0.02*unit; % [Khisamitov, Meschke, 2018, CMAME] (setup 1)
    cl = 0.025*unit/2; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME]
    % cl = 0.01*unit; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Molnar, Gravouil, 2015, FEAD]
    % cl = 0.01*unit/2; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Khisamitov, Meschke, 2018, CMAME] (setup 2), [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    if test
        cl = 0.025*unit;
    end
    switch lower(FEmesh)
        case 'unif'
            clD = cl; % characteristic length for domain
        case 'optim'
            clD = 0.1*unit; % characteristic length for domain
            % clD = 0.01*unit; % [Khisamitov, Meschke, 2018, CMAME] (setup 2)
            if test
                clD = 0.2*unit;
            end
        otherwise
            error('Wrong FE mesh')
    end
    clC = cl; % characteristic length for edge crack/notch
    clH = cl; % characteristic length for circular holes
    switch lower(initialCrack)
        case 'geometriccrack'
            S_phase = gmshasymmetricnotchedplatewithedgecrack(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_asymmetric_notched_plate'));
        case 'geometricnotch'
            c = 0.025*unit; % crack width
            S_phase = gmshasymmetricnotchedplatewithedgenotch(a,b,c,clD,clC,clH,unit,fullfile(pathname,'gmsh_asymmetric_notched_plate'));
        case 'initialphasefield'
            S_phase = gmshasymmetricnotchedplatewithedgecrack(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_asymmetric_notched_plate'),Dim,'noduplicate');
        otherwise
            error('Wrong model for initial crack');
    end
    S = S_phase;
    
    %% Phase-field problem
    %% Material
    % Critical energy release rate (or fracture toughness)
    gc = 1e3;
    % gc = 304.321; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % gc = 315; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % Regularization parameter (width of the smeared crack)
    % l = 0.2e-3; % [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS]
    % l = 0.15e-3; % [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS]
    % l = 0.05*unit; % [Wu, Nguyen, 2018, JMPS]
    l = 0.025*unit; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Molnar, Gravouil, 2015, FEAD], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
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
    C = LIGNE([-b,-h],[-b,-h+a]);
    R0 = 2*unit;
    BU = CIRCLE(0.0,h,R0);
    BL = CIRCLE(-ls,-h,R0);
    BR = CIRCLE(ls,-h,R0);
    % LU = LIGNE([-L,h],[L,h]);
    
    % findddlboundary = @(S_phase) findddl(S_phase,'T',LU);
    findddlboundary = @(S_phase) [];
    
    if strcmpi(initialCrack,'geometriccrack')
        S_phase = final(S_phase,'duplicate');
    else
        S_phase = final(S_phase);
    end
    
    if strcmpi(initialCrack,'initialphasefield')
        S_phase = addcl(S_phase,C,'T',1);
    end
    S_phase = addcl(S_phase,BU,'T');
    S_phase = addcl(S_phase,BL,'T');
    S_phase = addcl(S_phase,BR,'T');
    
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
    option = 'DEFO'; % plane strain [Guidault, Allix, Champaney, Cornuault, 2008, CMAME], [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Molnar, Gravouil, 2015, FEAD], [Khisamitov, Meschke, 2018, CMAME], [Rahimi, Moutsanidis, 2022, CMAME]
    % option = 'CONT'; % plane stress [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
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
    % E = 20e9; % [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS]
    % NU = 0.3; % [Msekh, Sargado, Jamshidian, Areias, Rabczuk, 2015, CMS]
    % E = 3.102e9; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % E = 3.275e9; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % NU = 0.35; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % E = 200e9; % [Guidault, Allix, Champaney, Cornuault, 2008, CMAME]
    % NU = 0.3; % [Guidault, Allix, Champaney, Cornuault, 2008, CMAME]
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
    % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % du = 1e-4 mm during 2500 time steps (up to u = 0.25 mm)
    % dt = 1e-4*unit;
    % nt = 2500;
    % t = linspace(dt,nt*dt,nt);
    
    % [Molnar, Gravouil, 2015, FEAD]
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
    save(fullfile(pathname,'problem.mat'),'unit','T','S_phase','S','addbc','findddlforce','findddlboundary');
else
    load(fullfile(pathname,'problem.mat'),'unit','T','S_phase','S','addbc','findddlforce','findddlboundary');
end

%% Solution
if solveProblem
    myparallel('start',numWorkers);
    
    %% Solution
    tTotal = tic;
    
    nbSamples = 1;
    fun = @(S_phase,S) solvePFDetLinElas(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,'maxiter',maxIter,'tol',tolConv,'crit',critConv);
    % fun = @(S_phase,S) solvePFDetLinElasAsymmetricNotchedPlate(S_phase,S,T,PFsolver,PU,PL,PR,'maxiter',maxIter,'tol',tolConv,'crit',critConv);
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
    [fmax_f,fmax_xi,fmax_bw] = ksdensity(fmax,'NumPoints',npts);
    [udmax_f,udmax_xi,udmax_bw] = ksdensity(udmax,'NumPoints',npts);
    [fc_f,fc_xi,fc_bw] = ksdensity(fc,'NumPoints',npts);
    [udc_f,udc_xi,udc_bw] = ksdensity(udc,'NumPoints',npts);
    
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
    
    %% Display means, variances and samples of solutions at different instants
    ampl = 0;
    switch setup
        case {1,4,5}
            tSnapshots = [0.190 0.201 0.203 0.204 0.205 0.207]*unit;
        case {2,3}
            tSnapshots = [0.210 0.215 0.218 0.219 0.220 0.222]*unit;
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
    framerate = 80;
    
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
