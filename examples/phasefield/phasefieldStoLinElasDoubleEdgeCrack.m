%% Phase-field fracture model - stochastic linear elasticity problem with double edge crack (two symmetric or asymmetric edge cracks)    %%
%%---------------------------------------------------------------------------------------------------------------------------------------%%
%% Setup 1: square domain [bxb]=[40mmx40mm] with crack length a=b/5=0.2*b=8mm and eccentricity e=b/40=0.025*b=1mm
% [Sumi, Wang, 1998, MM] (FE, LEFM)
% [Moes, Stolz, Bernard, Chevaugeon, 2011, IJNME] (TLS)
% [Molnar, Gravouil, 2017, FEAD] (isotropic phase-field model with no split of Bourdin et al.)
% [Bhowmick, Liu, 2018, EFM] (isotropic phase-field model with no split of Bourdin et al. and anisotropic tension-compression split for the phase-field driving force + CS-FEM)
% [Chen, Vasiukov, Gelebart, Park, 2019, CMAME] (anisotropic phase-field model of Miehe et al.)
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME] (anisotropic phase-field model of He et al.)
% [Si, Yu, Li, Natarajan, 2023, CMAME] (hybrid isotropic-anisotropic phase-field model of Ambati et al. with multi-patch adaptive isogeometric phase-field method based on Nitsche's method)
%% Setup 2: square domain [bxb]=[1mmx1mm] with crack length a=b/5=0.2*b=0.2mm and eccentricity e=b/20=0.05*b=0.05mm
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM] (anisotropic phase-field model of Amor et al.)
%% Setup 3: rectangular domain [WxH]=[60mmx120mm] with crack length a=H/12=10mm and eccentricity e=0mm, e=H/48=2.5mm, e=H/24=5mm, e=H/16=7.5mm
% [Shi, van Dam, van Mier, Sluys, 2000, MBS] (experimental tests)
% [Alfaiate, Wells, Sluys, 2002, EFM] (Embedded Discontinuity Elements)
% [Nguyen, 2005, PhD thesis] (Non local Damage model)
% [Nguyen, Houlsby, 2007, IJNAMG] (Non local Damage model)
% [Nguyen, Houlsby, 2008, IJNAMG] (Non local Damage model)
% [Nguyen, 2008, IJSS] (Non local Damage model)
% [Nguyen, Korsunsky, 2008, IJSS] (Non local Damage model)
% [Nguyen, 2011, IJSS] (Non local Damage model)
% [Galvez, Planas, Sancho, Reyes, Cendon, Casati, 2013, EFM] (Embedded Cohesive Crack model)
% [Stefanou, Georgioudakis, Papadrakakis, 2014, MMUQMS] (SLA)
% [Badnava, Mashayekhi, Kadkhodaei, 2016, IJDM] (GED)
% [Le, Nguyen, Bui, Sheikh, Kotousov, 2018, IJES] (Cohesive-Frictional model)
% [Badnava, Msekh, Etemadi, Rabczuk, 2018, FEAD] (anisotropic phase-field model of Miehe et al.)
% [Tong, Shen, Shao, Chen, 2020, EFM] (PD)
% [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM] (isotropic coupled phase-field fracture and plasticity model with no split of Bourdin et al. in elasto-plasticity)
% [De Maio, Cendon, Greco, Leonetti, Blasi, Planas, 2021, TAFM] (ECM vs DIM = Embedded Crack Model vs Diffuse Interface Model)
% [Li, Lu, Huang, Yang, 2022, OE] (PD)
% [Han, Li, Yu, Li, Zhang, 2022, JMPS] (PD-CZM)
% [Ribeiro Nogueira, Rastiello, Giry, Gatuingt, Callari, 2023, AJCE] (Gradient-Enhanced Eikonal non-local damage model)
% [Liu, Chen, Yuan, 2024, AAM] (PD-FEM)
%% Other setups
% [Melin, 1993] (experimental tests)
% [Fender, Lechenault, Daniels, 2010, PRL] (LEFM)
% [Judt, Ricoeur, 2015, TAFM] (FE, LEFM) (rectangular domain [WxH]=[70mmx80mm] with crack length a=W/10=7mm and eccentricity e=W/14=H/16=5mm)
% [Koivisto, Dalbe, Alava, Santucci, 2016, SR] (experimental DIC in polycarbonate sheets)
% [Mesgarnejad, Imanian, Karma, 2019, TAFM] (phase-field model of Karma, Kessler and Levin (KKL) for fatigue fracture) (rectangular domain [WxH=2*W] with crack lentgh a=W/2.5=0.4*W=0.2*H and eccentricity e=W/40=0.025*W=0.0125*H)
% [Brighenti, Rabczuk, Zhuang, 2021, EJMAS] (isotropic phase-field model with no split of Bourdin et al. and anisotropic tension-compression split for the phase-field driving force in large deformations of rate-dependent polymers) (rectangular domain [WxH=2.4*W] with crack length a=0.3*W=0.125*H and eccentricity e=a/2=0.15*W or e=a/3=0.1*W)
% [Grossman-Ponemon, Mesgarnejad, Karma, 2022, EFM] (phase-field model of Karma, Kessler and Levin (KKL) for fatigue fracture) (rectangular domain [WxH=2*W] with crack lentgh a=W/2.5=0.4*W=0.2*H and eccentricity e=W/40=0.025*W=0.0125*H, e=3*W/80=0.0375*W=0.01875*H)

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

numWorkers = maxNumCompThreads;
% numWorkers = 1; maxNumCompThreads(1); % mono-thread computation

% Deterministic model parameters
Dim = 2; % space dimension Dim = 2, 3
setup = 1; % notch geometry setup = 1, 2, 3
symmetry = 'Isot'; % 'Isot', 'MeanIsot'. Material symmetry
PFmodel = 'Miehe'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 1; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
critConv = 'Energy'; % 'Solution', 'Residual', 'Energy'
initialCrack = 'GeometricCrack'; % 'GeometricCrack', 'GeometricNotch', 'InitialPhaseField'
if setup==3
    initialCrack = 'GeometricNotch';
end
FEmesh = 'Optim'; % 'Unif' or 'Optim'

% Random model parameters
% N = 500; % number of samples
N = numWorkers;
switch setup
    case 1
        H = 40e-3; w = H; % w = 40e-3; % domain width
    case 2
        H = 1e-3; w = H; % w = 1e-3; % domain width
    case 3
        H = 120e-3; w = 60e-3;
end
lcorr = min(H,w)/10; % lcorr = 40e-4; % setup 1
                     % lcorr = 1e-4;  % setup 2
                     % lcorr = 60e-4; % setup 3
randMat = struct('delta',0.2,'lcorr',lcorr); % random material parameters model
randPF = struct('aGc',0,'bGc',0,'lcorr',Inf); % random phase-field parameters model

suffix = '';

foldername = ['doubleEdgeCracksSetup' num2str(setup) '_' num2str(Dim) 'D'];
filename = ['linElas' symmetry PFmodel PFsplit PFregularization PFsolver initialCrack...
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
    switch setup
        case 1
            H = 40e-3; % domain height
            w = H; % w = 40e-3; % domain width
            a = 0.2*H; % a = 8e-3; % crack length
            e = 0.025*H; % e = 1e-3; % eccentricity from the middle line
            DIM3 = 1e-3; % thickness
        case 2
            H = 1e-3; % domain height
            w = H; % w = 1e-3; % domain width
            a = 0.2*H; % a = 0.2e-3; % crack length
            e = -0.05*H; % e = -0.05e-3; % eccentricity from the middle line
            DIM3 = 0.1e-3; % thickness
        case 3
            H = 120e-3; % domain height
            w = 60e-3; % domain width
            a = H/12; % a = 10e-3; % crack length
            % e = -H/48; % e = -2.5e-3; % eccentricity from the middle line
            % e = -H/24; % e = -5e-3; % eccentricity from the middle line
            e = -H/16; % e = -7.5e-3; % eccentricity from the middle line
            DIM3 = 10e-3; % thickness
    end
    
    if Dim==2
        D = DOMAIN(2,[0.0,0.0],[w,H]);
        Ca = LINE([0.0,H/2+e],[a,H/2+e]);
        Cb = LINE([w-a,H/2-e],[w,H/2-e]);
    elseif Dim==3
        D = DOMAIN(3,[0.0,0.0,0.0],[w,H,DIM3]);
        Ca = QUADRANGLE([0.0,H/2+e,0.0],[  a,H/2+e,0.0],[  a,H/2+e,DIM3],[0.0,H/2+e,DIM3]);
        Cb = QUADRANGLE([  w,H/2-e,0.0],[w-a,H/2-e,0.0],[w-a,H/2-e,DIM3],[  w,H/2-e,DIM3]);
    end
    
    switch setup
        case 1
            % clC = 1e-4; % [Molnar, Gravouil, 2017, FEAD], [Chen, Vasiukov, Gelebart, Park, 2019, CMAME], [Si, Yu, Li, Natarajan, 2023, CMAME]
            % clC = 5e-5; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            switch lower(FEmesh)
                case 'unif' % uniform mesh
                    cl = 1e-4;
                    if test
                        if Dim==2
                            cl = 1.75e-4;
                        elseif Dim==3
                            clC = 4e-4;
                        end
                    end
                    clD = cl;
                    clC = cl;
                case 'optim' % optimized mesh
                    clD = 1e-3;
                    clC = 1e-4;
                    if test
                        clD = 2e-3;
                        if Dim==2
                            clC = 1.75e-4;
                        elseif Dim==3
                            clC = 4e-4;
                        end
                    end
                otherwise
                    error('Wrong FE mesh')
            end
        case 2
            % cl = 3.3e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
            % cl = 2.5e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
            switch lower(FEmesh)
                case 'unif' % uniform mesh
                    cl = 2.5e-6;
                    if test
                        if Dim==2
                            cl = 5e-6;
                        elseif Dim==3
                            cl = 1e-5;
                        end
                    end
                    clD = cl;
                    clC = cl;
                case 'optim' % optimized mesh
                    clD = 2.5e-5;
                    clC = 2.5e-6;
                    if test
                        if Dim==2
                            clD = 5e-5;
                            clC = 5e-6;
                        elseif Dim==3
                            clD = 8e-5;
                            clC = 1.75e-5;
                        end
                    end
                otherwise
                    error('Wrong FE mesh')
            end
        case 3
            % clD = 4.3e-3; % [Badnava, Msekh, Etemadi, Rabczuk, 2018, FEAD]
            % clC = 3e-3; % [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM]
            % clC = 1.5e-3; % [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM]
            % clC = 1.08e-3; % [Badnava, Msekh, Etemadi, Rabczuk, 2018, FEAD]
            % clC = 1e-3; % [Tong, Shen, Shao, Chen, 2020, EFM]
            % clC = 0.6e-3; % [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM]
            % clC = 0.5e-3; % [Li, Lu, Huang, Yang, 2022, OE]
            % clC = 0.1e-3; % [Han, Li, Yu, Li, Zhang, 2022, JMPS]
            switch lower(FEmesh)
                case 'unif' % uniform mesh
                    cl = 0.5e-3;
                    if test
                        cl = 0.75e-3;
                    end
                    clD = cl;
                    clC = cl;
                case 'optim' % optimized mesh
                    clD = 5e-3;
                    clC = 0.5e-3;
                    if test
                        if Dim==2
                            clD = 5e-3;
                            clC = 0.75e-3;
                        elseif Dim==3
                            clD = 6e-3;
                            clC = 1.5e-3;
                            % clC = 1.75e-3;
                        end
                    end
                otherwise
                    error('Wrong FE mesh')
            end
    end
    switch lower(FEmesh)
        case 'unif' % uniform mesh
            B = [];
        case 'optim' % optimized mesh
            VIn = clC; VOut = clD;
            % XMin = a; XMax = w-a;
            XMin = 0; XMax = w;
            switch setup
                case {1,2}
                    YMin = H/2-abs(3*e); YMax = H/2+abs(3*e);
                    Thickness = H/2-abs(3*e);
                case 3
                    % YMin = H/2-3*H/16; YMax = H/2+3*H/16;
                    % Thickness = H/2-3*H/16;
                    YMin = H/2-H/8; YMax = H/2+H/8;
                    Thickness = H/2-H/8;
            end
            % Thickness = 0;
            if Dim==2
                B = struct('VIn',VIn,'VOut',VOut,'XMin',XMin,'XMax',XMax,'YMin',YMin,'YMax',YMax,'Thickness',Thickness);
            elseif Dim==3
                ZMin = 0; ZMax = DIM3;
                B = struct('VIn',VIn,'VOut',VOut,'XMin',XMin,'XMax',XMax,'YMin',YMin,'YMax',YMax,'ZMin',ZMin,'ZMax',ZMax,'Thickness',Thickness);
            end
    end
    switch lower(initialCrack)
        case 'geometriccrack'
            S_phase = gmshDomainWithDoubleEdgeCrack(D,Ca,Cb,clD,clC,fullfile(pathname,'gmsh_domain_double_edge_crack'),Dim,'Box',B);
        case 'geometricnotch'
            switch setup
                case 1
                    c = 4e-4; % crack width
                    clC = min(clC,c/2);
                    S_phase = gmshDomainWithDoubleEdgeNotch(D,Ca,Cb,c,clD,clC,fullfile(pathname,'gmsh_domain_double_edge_crack'),Dim,'Box',B);
                case 2
                    % c = 6.6e-6; % crack width
                    c = 1e-5; % crack width
                    clC = min(clC,c/2);
                    S_phase = gmshDomainWithDoubleEdgeNotch(D,Ca,Cb,c,clD,clC,fullfile(pathname,'gmsh_domain_double_edge_crack'),Dim,'Box',B);
                case 3
                    c = 2e-3; % notch width
                    S_phase = gmshDomainWithDoubleEdgeNotch(D,Ca,Cb,c,clD,clC,fullfile(pathname,'gmsh_domain_double_edge_crack'),Dim,'Box',B,'rectangular');
            end
        case 'initialphasefield'
            S_phase = gmshDomainWithDoubleEdgeCrack(D,Ca,Cb,clD,clC,fullfile(pathname,'gmsh_domain_double_edge_crack'),Dim,'noduplicate','refinecrack','Box',B);
        otherwise
            error('Wrong model for initial crack');
    end
    S = S_phase;
    
    %% Phase-field problem
    %% Material
    switch setup
        case 1
            % Critical energy release rate (or fracture toughness)
            % gc = 3.4e3; % [Bhowmick, Liu, 2018, EFM]
            gc = 2.7e3; % [Molnar, Gravouil, 2017, FEAD], [Bhowmick, Liu, 2018, EFM], [Chen, Vasiukov, Gelebart, Park, 2019, CMAME], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Si, Yu, Li, Natarajan, 2023, CMAME]
            % Regularization parameter (width of the smeared crack)
            l = 2e-4; % [Molnar, Gravouil, 2017, FEAD], [Chen, Vasiukov, Gelebart, Park, 2019, CMAME], [Si, Yu, Li, Natarajan, 2023, CMAME]
            % l = 1.35e-4; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            % l = 1e-4;
            % l = 7.5e-3; % [Bhowmick, Liu, 2018, EFM]
        case 2
            % Critical energy release rate (or fracture toughness)
            gc = 1e3; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
            % Regularization parameter (width of the smeared crack)
            % l = 6.6e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
            l = 5e-6;
        case 3
            % Critical energy release rate (or fracture toughness)
            % gc = 63; % [Le, Nguyen, Bui, Sheikh, Kotousov, 2018, IJES], [Tong, Shen, Shao, Chen, 2020, EFM]
            gc = 59; % [Nguyen, 2005, PhD thesis], [Nguyen, Houlsby, 2007, IJNAMG], [Nguyen, Houlsby, 2008, IJNAMG], [Nguyen, 2008, IJSS], [Nguyen, Korsunsky, 2008, IJSS], [Nguyen, 2011, IJSS], [Stefanou, Georgioudakis, Papadrakakis, 2014, MMUQMS], [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM], [Alfaiate, Wells, Sluys, 2002, EFM]
            % gc = 50; % [Galvez, Planas, Sancho, Reyes, Cendon, Casati, 2013, EFM]
            % gc = 25; % [Han, Li, Yu, Li, Zhang, 2022, JMPS]
            % gc = 4.5; % [Badnava, Msekh, Etemadi, Rabczuk, 2018, FEAD]
            % Regularization parameter (width of the smeared crack)
            % l = 10e-3; % [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM]
            % l = 5e-3; % [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM]
            % l = 3e-3; % [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM]
            % l = 2.16e-3; % [Badnava, Msekh, Etemadi, Rabczuk, 2018, FEAD]
            l = 1e-3;
    end
    % Small artificial residual stiffness
    % k = 1e-12;
    k = 0;
    
    % Material
    [K,R,Qn] = setphasefieldparam(gc,l,PFregularization);
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'DIM3',DIM3,'PFregularization',PFregularization,'aGc',randPF.aGc,'bGc',randPF.bGc,'lcorr',randPF.lcorr);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BRight = LINE([  w,0.0],[  w,H]);
        BLeft  = LINE([0.0,0.0],[0.0,H]);
    elseif Dim==3
        BRight = PLANE([  w,0.0,0.0],[  w,H,0.0],[  w,0.0,DIM3]);
        BLeft  = PLANE([0.0,0.0,0.0],[0.0,H,0.0],[0.0,0.0,DIM3]);
    end
    
    findddlboundary = @(S_phase) union(findddl(S_phase,'T',BRight),findddl(S_phase,'T',BLeft));
    
    if strcmpi(initialCrack,'geometriccrack')
        S_phase = final(S_phase,'duplicate');
    else
        S_phase = final(S_phase);
    end
    
    S_phase = addbcdamageDoubleEdgeCrack(S_phase,Ca,Cb,initialCrack);
    
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
    switch setup
        case 1
            % Option
            option = 'DEFO'; % plane strain [Moes, Stolz, Bernard, Chevaugeon, 2011, IJNME], [Molnar, Gravouil, 2017, FEAD], [Bhowmick, Liu, 2018, EFM], [Chen, Vasiukov, Gelebart, Park, 2019, CMAME], [De Maio, Cendon, Greco, Leonetti, Blasi, Planas, 2021, TAFM], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Si, Yu, Li, Natarajan, 2023, CMAME]
            % Lame coefficients
            lambda = 121.15e9; % [Chen, Vasiukov, Gelebart, Park, 2019, CMAME], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            mu = 80.77e9; % [Chen, Vasiukov, Gelebart, Park, 2019, CMAME], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            % Young modulus and Poisson ratio
            switch lower(option)
                case 'defo'
                    E = mu*(3*lambda+2*mu)/(lambda+mu); % E = 210e9;
                    NU = lambda/(lambda+mu)/2; % NU = 0.3;
                case 'cont'
                    E = 4*mu*(lambda+mu)/(lambda+2*mu);
                    NU = lambda/(lambda+2*mu);
            end
            % E = 200e9; NU = 0.3; % [Sumi, Wang, 1998, MM], [Moes, Stolz, Bernard, Chevaugeon, 2011, IJNME]
            % E = 210e9; NU = 0.3; % [Molnar, Gravouil, 2017, FEAD], [Si, Yu, Li, Natarajan, 2023, CMAME]
            % E = 210e9; NU = 0.2; % [Bhowmick, Liu, 2018, EFM]
        case 2
            % Option
            option = 'CONT'; % plane stress [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
            % Young modulus and Poisson ratio
            E = 1e9; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
            NU = 0.3; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
        case 3
            % Option
            % option = 'DEFO'; % plane strain [Ribeiro Nogueira, Rastiello, Giry, Gatuingt, Callari, 2023, AJCE]
            option = 'CONT'; % plane stress [Nguyen, 2011, IJSS], [Stefanou, Georgioudakis, Papadrakakis, 2014, MMUQMS], [Tong, Shen, Shao, Chen, 2020, EFM], [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM], [Li, Lu, Huang, Yang, 2022, OE], [Han, Li, Yu, Li, Zhang, 2022, JMPS]
            % Young modulus and Poisson ratio
            % E = 40e9; NU = 1/3; % [Tong, Shen, Shao, Chen, 2020, EFM]
            % E = 40e9; NU = 0.2; % [Han, Li, Yu, Li, Zhang, 2022, JMPS]
            % E = 31e9; NU = 0.2; % [Galvez, Planas, Sancho, Reyes, Cendon, Casati, 2013, EFM], [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM]
            E = 30e9; NU = 0.2; % [Nguyen, Houlsby, 2007, IJNAMG], [Nguyen, 2011, IJSS], [Badnava, Msekh, Etemadi, Rabczuk, 2018, FEAD]
            % E = 24e9; NU = 0.2; % [Shi, van Dam, van Mier, Sluys, 2000, MBS], [Alfaiate, Wells, Sluys, 2002, EFM], [Nguyen, 2005, PhD thesis], [Nguyen, Houlsby, 2008, IJNAMG], [Nguyen, 2008, IJSS], [Nguyen, Korsunsky, 2008, IJSS], [Stefanou, Georgioudakis, Papadrakakis, 2014, MMUQMS], [Le, Nguyen, Bui, Sheikh, Kotousov, 2018, IJES], [Ribeiro Nogueira, Rastiello, Giry, Gatuingt, Callari, 2023, AJCE], [Liu, Chen, Yuan, 2024, AAM]
            % E = 24e9; NU = 1/3; % [Li, Lu, Huang, Yang, 2022, OE]
            % E = 20e9; NU = 0.2; % [Badnava, Mashayekhi, Kadkhodaei, 2016, IJDM]
    end
    if strcmpi(symmetry,'meanisot')
        mu = E/(1+NU)/2; % second Lame coefficient (shear modulus)
        % Elasticity matrix
        if Dim==2
            switch lower(option)
                case 'defo'
                    % Cmat = DIM3*E/(1+NU)/(1-2*NU)*...
                    %     [(1-NU),NU,0;...
                    %     NU,(1-NU),0;...
                    %     0,0,(1-2*NU)/2];
                    lambda = E*NU/(1+NU)/(1-2*NU); % first Lame coefficient
                otherwise
                    % Cmat = DIM3*E/(1-NU^2)*...
                    %     [1,NU,0;...
                    %     NU,1,0;...
                    %     0,0,(1-NU)/2];
                    lambda = E*NU/(1-NU^2); % first Lame coefficient
            end
            Cmat = DIM3*...
                [lambda+2*mu,lambda,0;...
                lambda,lambda+2*mu,0;...
                0,0,mu];
        elseif Dim==3
            % Cmat = E/(1+NU)/(1-2*NU)*...
            %     [(1-NU),NU,NU,0,0,0;...
            %     NU,(1-NU),NU,0,0,0;...
            %     NU,NU,(1-NU),0,0,0;...
            %     0,0,0,(1-2*NU)/2,0,0;...
            %     0,0,0,0,(1-2*NU)/2,0;...
            %     0,0,0,0,0,(1-2*NU)/2];
            lambda = E*NU/(1+NU)/(1-2*NU); % first Lame coefficient
            Cmat = [lambda+2*mu,lambda,lambda,0,0,0;...
                lambda,lambda+2*mu,lambda,0,0,0;...
                lambda,lambda,lambda+2*mu,0,0,0;...
                0,0,0,mu,0,0;...
                0,0,0,0,mu,0;...
                0,0,0,0,0,mu];
        end
    end
    % Energetic degradation function
    g = @(d) (1-d).^2;
    % Density
    RHO = 1;
    % RHO = 2400; % [Li, Lu, Huang, Yang, 2022, OE]
    
    % Material
    d = calc_init_dirichlet(S_phase);
    switch lower(symmetry)
        case 'isot' % almost surely isotropic material
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit,'delta',randMat.delta,'lcorr',randMat.lcorr);
        case 'meanisot' % mean isotropic material
            mat = ELAS_ANISOT('C',Cmat,'RHO',RHO,'DIM3',DIM3,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit,'delta',randMat.delta,'lcorr',randMat.lcorr);
        otherwise
            error('Wrong material symmetry class');
    end
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BU = LINE([0.0,H],[w,H]);
        BL = LINE([0.0,0.0],[w,0.0]);
        P0 = POINT([0.0,0.0]); % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Ribeiro Nogueira, Rastiello, Giry, Gatuingt, Callari, 2023, AJCE]
        % P0 = POINT([0.0,H]); % [Tong, Shen, Shao, Chen, 2020, EFM], [Li, Lu, Huang, Yang, 2022, OE]
    elseif Dim==3
        BU = PLANE([0.0,H,0.0],[w,H,0.0],[0.0,H,DIM3]);
        BL = PLANE([0.0,0.0,0.0],[w,0.0,0.0],[0.0,0.0,DIM3]);
        P0 = POINT([0.0,0.0,0.0]); % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Ribeiro Nogueira, Rastiello, Giry, Gatuingt, Callari, 2023, AJCE]
        % P0 = POINT([0.0,H,0.0]); % [Tong, Shen, Shao, Chen, 2020, EFM], [Li, Lu, Huang, Yang, 2022, OE]
    end
    
    addbc = @(S,ud) addbcDoubleEdgeCrack(S,ud,BU,BL,P0,setup);
    findddlforce = @(S) findddl(S,'UY',BU);
    
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
    switch setup
        case 1
            % [Molnar, Gravouil, 2017, FEAD], [Si, Yu, Li, Natarajan, 2023, CMAME]
            % du = 1e-4 mm during the first 400 time steps (up to u = 0.04 mm)
            % du = 1e-5 mm during the last 1000 time steps (up to u = 0.05 mm)
            % dt0 = 1e-7;
            % nt0 = 400;
            % dt1 = 1e-8;
            % nt1 = 1000;
            % t0 = linspace(dt0,nt0*dt0,nt0);
            % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
            % t = [t0,t1];
            % T = TIMEMODEL(t);
            
            % [Bhowmick, Liu, 2018, EFM]
            % du = 1e-4 mm during 500 time steps (up to u = 0.05 mm)
            % dt = 1e-7;
            % nt = 500;
            % t = linspace(dt,nt*dt,nt);
            % T = TIMEMODEL(t);
            
            % [Chen, Vasiukov, Gelebart, Park, 2019, CMAME]
            % du = 8e-6*h = 3.2e-4 mm during the first 100 time steps (up to u = 0.032 mm)
            % du = 5e-7*h = 2e-5 mm during the last 900 time steps (up to u = 0.05 mm)
            % dt0 = 3.2e-7;
            % nt0 = 100;
            % dt1 = 2e-8;
            % nt1 = 900;
            % t0 = linspace(dt0,nt0*dt0,nt0);
            % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
            % t = [t0,t1];
            % T = TIMEMODEL(t);
            
            % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            % du = 1e-4 mm during the first stage (until the phase-field reaches the threshold value)
            % du = 2e-5 mm during the last stage (as soon as the phase-field exceeds the threshold value, up to u = 0.15 mm)
            % dt0 = 1e-7;
            % dt1 = 2e-8;
            % if test
            %     dt0 = 2e-7;
            %     dt1 = 4e-8;
            % end
            % tf = 15e-5;
            % dth = 0.6;
            % T = struct('dt0',dt0,'dt1',dt1,'tf',tf,'dth',dth);
            
            % du = 1e-4 mm during the first 400 time steps (up to u = 0.04 mm)
            % du = 2e-5 mm during the last 4000 time steps (up to u = 0.12 mm)
            dt0 = 1e-7;
            nt0 = 400;
            dt1 = 2e-8;
            nt1 = 4000;
            if test
                dt0 = 1e-6;
                nt0 = 40;
                dt1 = 2e-7;
                nt1 = 400;
            end
            t0 = linspace(dt0,nt0*dt0,nt0);
            t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
            t = [t0,t1];
            T = TIMEMODEL(t);
        case 2
            % du = 1e-4 mm during the first 500 time steps (up to u = 0.05 mm)
            % du = 2e-5 mm during the last 5000 time steps (up to u = 0.15 mm)
            dt0 = 1e-7;
            nt0 = 500;
            dt1 = 2e-8;
            nt1 = 5000;
            if test
                dt0 = 1e-6;
                nt0 = 50;
                dt1 = 2e-7;
                nt1 = 500;
            end
            t0 = linspace(dt0,nt0*dt0,nt0);
            t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
            t = [t0,t1];
            T = TIMEMODEL(t);
        case 3
            % [Tong, Shen, Shao, Chen, 2020, EFM]
            % du = 1e-4 mm during 2000 time steps (up to u = 0.2 mm) for tension
            % du = -1e-4 mm during 2000 time steps (up to u = -0.2 mm) for compression
            dt = 1e-7;
            nt = 2000;
            if test
                dt = 2e-7;
                nt = 1000;
            end
            t = linspace(dt,nt*dt,nt);
            T = TIMEMODEL(t);
            
            % [Li, Lu, Huang, Yang, 2022, OE]
            % du = 1e-5 mm during 20 000 time steps (up to u = 0.2 mm)
            % dt = 1e-8;
            % nt = 2e4;
            % t = linspace(dt,nt*dt,nt);
            % T = TIMEMODEL(t);
            
            % [Liu, Chen, Yuan, 2024, AAM]
            % du = 2e-6 mm during 100 000 time steps (up to u = 0.2 mm)
            % dt = 2e-9;
            % nt = 1e5;
            % t = linspace(dt,nt*dt,nt);
            % T = TIMEMODEL(t);
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','Ca','Cb','addbc','findddlforce','findddlboundary');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','Ca','Cb','addbc','findddlforce','findddlboundary');
end

%% Solution
if solveProblem
    myparallel('start',numWorkers);
    
    %% Solution
    tTotal = tic;
    
    displayIter = false;
    displaySol  = false;
    
    nbSamples = 1;
    fun = @(S_phase,S) solvePFDetLinElas(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,...
        'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
        'displayiter',displayIter,'displaysol',displaySol);
    % fun = @(S_phase,S) solvePFDetLinElasDoubleEdgeCrack(S_phase,S,T,PFsolver,BU,BL,P0,BRight,BLeft,setup,...
    %     'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
    %     'displayiter',displayIter,'displaysol',displaySol);
    [ft,Edt,Eut,dmaxt,dt_mean,ut_mean,dt_var,ut_var,dt_sample,ut_sample] = solvePFStoLinElas(S_phase,S,T,fun,N,'nbsamples',nbSamples);
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
    
    Edt_mean = mean(Edt);
    Edt_std = std(Edt);
    Edt_ci = quantile(Edt,probs);
    
    Eut_mean = mean(Eut);
    Eut_std = std(Eut);
    Eut_ci = quantile(Eut,probs);
    
    dmaxt_mean = mean(dmaxt);
    dmaxt_std = std(dmaxt);
    dmaxt_ci = quantile(dmaxt,probs);
    
    save(fullfile(pathname,'solution.mat'),'N','ft','Edt','Eut','dmaxt',...
        'dt_mean','ut_mean','dt_var','ut_var','dt_sample','ut_sample',...
        'ft_mean','ft_std','ft_ci',...
        'Edt_mean','Edt_std','Edt_ci',...
        'Eut_mean','Eut_std','Eut_ci',...
        'dmaxt_mean','dmaxt_std','dmaxt_ci',...
        'probs','time',...
        'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
        'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw',...
        'fc','fc_mean','fc_std','fc_ci','fc_f','fc_xi','fc_bw',...
        'udc','udc_mean','udc_std','udc_ci','udc_f','udc_xi','udc_bw');
else
    load(fullfile(pathname,'solution.mat'),'N','ft','Edt','Eut','dmaxt',...
        'dt_mean','ut_mean','dt_var','ut_var','dt_sample','ut_sample',...
        'ft_mean','ft_std','ft_ci',...
        'Edt_mean','Edt_std','Edt_ci',...
        'Eut_mean','Eut_std','Eut_ci',...
        'dmaxt_mean','dmaxt_std','dmaxt_ci',...
        'probs','time',...
        'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
        'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw',...
        'fc','fc_mean','fc_std','fc_ci','fc_f','fc_xi','fc_bw',...
        'udc','udc_mean','udc_std','udc_ci','udc_f','udc_xi','udc_bw');
end

%% Outputs
filenameResults = fullfile(pathname,'results.txt');
if solveProblem
    fid = fopen(filenameResults,'w');
    fprintf(fid,'Double edge crack');
    if e==0
        fprintf(fid,' (two symmetric edge cracks)\n');
    else
        fprintf(fid,' (two asymmetric edge cracks)\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'dim      = %d\n',Dim);
    fprintf(fid,'setup    = %d\n',setup);
    fprintf(fid,'mat sym  = %s\n',symmetry);
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
    switch setup
        case {1,3}
            fprintf(fid,'mean(fmax)   = %g kN\n',fmax_mean*1e-3);
            fprintf(fid,'std(fmax)    = %g kN\n',fmax_std*1e-3);
        case 2
            fprintf(fid,'mean(fmax)   = %g N\n',fmax_mean);
            fprintf(fid,'std(fmax)    = %g N\n',fmax_std);
    end
    fprintf(fid,'disp(fmax)   = %g\n',fmax_std/fmax_mean);
    switch setup
        case {1,3}
            fprintf(fid,'%d%% ci(fmax) = [%g,%g] kN\n',(probs(2)-probs(1))*100,fmax_ci(1)*1e-3,fmax_ci(2)*1e-3);
        case 2
            fprintf(fid,'%d%% ci(fmax) = [%g,%g] N\n',(probs(2)-probs(1))*100,fmax_ci(1),fmax_ci(2));
    end
    
    fprintf(fid,'\n');
    switch setup
        case {1,3}
            fprintf(fid,'mean(fc)   = %g kN\n',fc_mean*1e-6);
            fprintf(fid,'std(fc)    = %g kN\n',fc_std*1e-6);
        case 2
            fprintf(fid,'mean(fc)   = %g N\n',fc_mean);
            fprintf(fid,'std(fc)    = %g N\n',fc_std);
    end
    fprintf(fid,'disp(fc)   = %g\n',fc_std/fc_mean);
    switch setup
        case {1,3}
            fprintf(fid,'%d%% ci(fc) = [%g,%g] kN\n',(probs(2)-probs(1))*100,fc_ci(1)*1e-3,fc_ci(2)*1e-3);
        case 2
            fprintf(fid,'%d%% ci(fc) = [%g,%g] N\n',(probs(2)-probs(1))*100,fc_ci(1),fc_ci(2));
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
end
type(filenameResults) % fprintf('%s', fileread(filenameResults))

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
    % plotDomain({D,Ca,Cb},'legend',false);
    % mysaveas(pathname,'domain',formats,renderer);
    % mymatlab2tikz(pathname,'domain.tex');
    
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
    
    % u = ut_sample(:,:,end);
    % for k=1:size(u,1)
    %     ampl = getsize(S)/max(abs(u(k,:)))/20;
    %     plotModelDeflection(S,u(k,:)','ampl',ampl,'Color','b','FaceColor',facecolordef,'FaceAlpha',facealpha,'legend',false);
    %     mysaveas(pathname,['mesh_deflected_sample_' num2str(k)],formats,renderer);
    %     
    %     figure('Name','Meshes')
    %     clf
    %     plot(S,'Color','k','FaceColor',facecolor,'FaceAlpha',facealpha);
    %     plot(S+ampl*unfreevector(S,u(k,:)'),'Color','b','FaceColor',facecolordef,'FaceAlpha',facealpha);
    %     mysaveas(pathname,['meshes_deflected_' num2str(k)],formats,renderer);
    % end
end

%% Display solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display force-displacement curve
    figure('Name','Force vs displacement')
    clf
    switch setup
        case {1,3}
            plot([0,t]*1e3,[0,ft_mean]*1e-3,'-b','LineWidth',linewidth)
            hold on
            ciplot([0,ft_ci(1,:)]*1e-3,[0,ft_ci(2,:)]*1e-3,[0,t]*1e3,'b');
        case 2
            plot([0,t]*1e3,[0,ft_mean],'-b','LineWidth',linewidth)
            hold on
            ciplot([0,ft_ci(1,:)],[0,ft_ci(2,:)],[0,t]*1e3,'b');
    end
    alpha(0.2)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    switch setup
        case {1,3}
            ylabel('Force [kN]','Interpreter',interpreter)
        case 2
            ylabel('Force [N]','Interpreter',interpreter)
    end
    legend('mean function',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'Location','NorthWest','Interpreter',interpreter)
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    colors = distinguishable_colors(N);
    figure('Name','Forces vs displacement')
    clf
    for i=1:N
        switch setup
            case {1,3}
                plot([0,t]*1e3,[0,ft(i,:)]*1e-3,'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
            case 2
                plot([0,t]*1e3,[0,ft(i,:)],'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
        end
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]'...,'Interpreter',interpreter...
        )
    switch setup
        case {1,3}
            ylabel('Force [kN]'...,'Interpreter',interpreter...
                )
        case 2
            ylabel('Force [N]'...,'Interpreter',interpreter...
                )
    end
    mysaveas(pathname,'forces_displacement',formats);
    mymatlab2tikz(pathname,'forces_displacement.tex');
    
    %% Display maximum damage-displacement curve
    figure('Name','Maximum damage vs displacement')
    clf
    plot([0,t]*1e3,[0,dmaxt_mean],'-b','LineWidth',linewidth)
    hold on
    ciplot([0,dmaxt_ci(1,:)],[0,dmaxt_ci(2,:)],[0,t]*1e3,'b');
    alpha(0.2)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Maximum damage','Interpreter',interpreter)
    legend('mean function',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'Location','NorthWest','Interpreter',interpreter)
    mysaveas(pathname,'max_damage_displacement',formats);
    mymatlab2tikz(pathname,'max_damage_displacement.tex');
    
    colors = distinguishable_colors(N);
    figure('Name','Maximum damages vs displacement')
    clf
    for i=1:N
        plot([0,t]*1e3,[0,dmaxt(i,:)],'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]'...,'Interpreter',interpreter...
        )
    ylabel('Maximum damage'...,'Interpreter',interpreter...
        )
    mysaveas(pathname,'max_damages_displacement',formats);
    mymatlab2tikz(pathname,'max_damages_displacement.tex');
    
    %% Display energy-displacement curves
    Et_ci = quantile(Eut+Edt,probs);
    figure('Name','Energies vs displacement')
    clf
    switch setup
        case 1
            plot([0,t]*1e3,[0,Eut_mean],'-b','LineWidth',linewidth)
            hold on
            plot([0,t]*1e3,[0,Edt_mean],'-r','LineWidth',linewidth)
            plot([0,t]*1e3,[0,Eut_mean+Edt_mean],'-k','LineWidth',linewidth)
            ciplot([0,Eut_ci(1,:)],[0,Eut_ci(2,:)],[0,t]*1e3,'b');
            ciplot([0,Edt_ci(1,:)],[0,Edt_ci(2,:)],[0,t]*1e3,'r');
            ciplot([0,Et_ci(1,:)],[0,Et_ci(2,:)],[0,t]*1e3,'k');
        case {2,3}
            plot([0,t]*1e3,[0,Eut_mean]*1e3,'-b','LineWidth',linewidth)
            hold on
            plot([0,t]*1e3,[0,Edt_mean]*1e3,'-r','LineWidth',linewidth)
            plot([0,t]*1e3,[0,Eut_mean+Edt_mean]*1e3,'-k','LineWidth',linewidth)
            ciplot([0,Eut_ci(1,:)]*1e3,[0,Eut_ci(2,:)]*1e3,[0,t]*1e3,'b');
            ciplot([0,Edt_ci(1,:)]*1e3,[0,Edt_ci(2,:)]*1e3,[0,t]*1e3,'r');
            ciplot([0,Et_ci(1,:)]*1e3,[0,Et_ci(2,:)]*1e3,[0,t]*1e3,'k');
    end
    alpha(0.2)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    switch setup
        case 1
            ylabel('Energy [J]','Interpreter',interpreter)
        case {2,3}
            ylabel('Energy [mJ]','Interpreter',interpreter)
    end
    legend('mean function - elastic','mean function - fracture','mean function - total',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval - elastic'],...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval - fracture'],...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval - total'],...
        'Location','NorthWest','Interpreter',interpreter)
    mysaveas(pathname,'energies_displacement',formats);
    mymatlab2tikz(pathname,'energies_displacement.tex');
    
    % Elastic energy-displacement curves
    colors = distinguishable_colors(N);
    figure('Name','Elastic energies vs displacement')
    clf
    for i=1:N
        switch setup
            case 1
                plot([0,t]*1e3,[0,Eut(i,:)],'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
            case {2,3}
                plot([0,t]*1e3,[0,Eut(i,:)]*1e3,'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
        end
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]'...,'Interpreter',interpreter...
        )
    switch setup
        case 1
            ylabel('Elastic strain energy [J]'...,'Interpreter',interpreter...
                )
        case {2,3}
            ylabel('Elastic strain energy [mJ]'...,'Interpreter',interpreter...
                )
    end
    mysaveas(pathname,'energies_elastic_displacement',formats);
    mymatlab2tikz(pathname,'energies_elastic_displacement.tex');
    
    % Fracture energy-displacement curves
    figure('Name','Fracture energies vs displacement')
    clf
    for i=1:N
        switch setup
            case 1
                plot([0,t]*1e3,[0,Edt(i,:)],'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
            case {2,3}
                plot([0,t]*1e3,[0,Edt(i,:)]*1e3,'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
        end
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]'...,'Interpreter',interpreter...
        )
    switch setup
        case 1
            ylabel('Fracture energy [J]'...,'Interpreter',interpreter...
                )
        case {2,3}
            ylabel('Fracture energy [mJ]'...,'Interpreter',interpreter...
                )
    end
    mysaveas(pathname,'energies_fracture_displacement',formats);
    mymatlab2tikz(pathname,'energies_fracture_displacement.tex');
    
    % Total energy-displacement curves
    figure('Name','Total energies vs displacement')
    clf
    for i=1:N
        switch setup
            case 1
                plot([0,t]*1e3,[0,Eut(i,:)+Edt(i,:)],'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
            case {2,3}
                plot([0,t]*1e3,[0,Eut(i,:)+Edt(i,:)]*1e3,'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
        end
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]'...,'Interpreter',interpreter...
        )
    switch setup
        case 1
            ylabel('Total energy [J]'...,'Interpreter',interpreter...
                )
        case {2,3}
            ylabel('Total energy [mJ]'...,'Interpreter',interpreter...
                )
    end
    mysaveas(pathname,'energies_total_displacement',formats);
    mymatlab2tikz(pathname,'energies_total_displacement.tex');
    
    %% Display pdf of maximum force
    figure('Name','Probability Density Estimate: Maximum force')
    clf
    switch setup
        case {1,3}
            plot(fmax_xi*1e-3,fmax_f*1e3,'-r','LineWidth',linewidth)
            hold on
            ind_fmax = find(fmax_xi>=fmax_ci(1) & fmax_xi<fmax_ci(2));
            area(fmax_xi(ind_fmax)*1e-3,fmax_f(ind_fmax)*1e3,'FaceColor','r','EdgeColor','none','FaceAlpha',0.2)
            scatter(fmax_mean*1e-3,0,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r')
        case 2
            plot(fmax_xi,fmax_f,'-r','LineWidth',linewidth)
            hold on
            ind_fmax = find(fmax_xi>=fmax_ci(1) & fmax_xi<fmax_ci(2));
            area(fmax_xi(ind_fmax),fmax_f(ind_fmax),'FaceColor','r','EdgeColor','none','FaceAlpha',0.2)
            scatter(fmax_mean,0,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r')
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    switch setup
        case {1,3}
            xlabel('$f$ [kN]','Interpreter',interpreter)
        case 2
            xlabel('$f$ [N]','Interpreter',interpreter)
    end
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
    switch setup
        case {1,3}
            plot(fc_xi*1e-3,fc_f*1e3,'-b','LineWidth',linewidth)
            hold on
            ind_fc = find(fc_xi>=fc_ci(1) & fc_xi<fc_ci(2));
            area(fc_xi(ind_fc)*1e-3,fc_f(ind_fc)*1e3,'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
            scatter(fc_mean*1e-3,0,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b')
        case 2
            plot(fc_xi,fc_f,'-b','LineWidth',linewidth)
            hold on
            ind_fc = find(fc_xi>=fc_ci(1) & fc_xi<fc_ci(2));
            area(fc_xi(ind_fc),fc_f(ind_fc),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
            scatter(fc_mean,0,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b')
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    switch setup
        case {1,3}
            xlabel('$f$ [kN]','Interpreter',interpreter)
        case 2
            xlabel('$f$ [N]','Interpreter',interpreter)
    end
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
        case 1
            tSnapshots = [0.04 0.0425 0.045 0.05 0.06 0.08]*1e-3;
        case 2
            tSnapshots = [0.055 0.06 0.065 0.07 0.08 0.10 0.12]*1e-3;
        case 3
            tSnapshots = [0.03 0.035 0.04 0.05 0.10 0.15]*1e-3;
    end
    rep = arrayfun(@(x) find(t>x-eps,1),tSnapshots);
    rep = [rep,length(T)];
    % tSnapshots = [tSnapshots,gett1(T)];
    % rep = arrayfun(@(x) find(t>x-eps,1),tSnapshots);
    
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
        % plotSolution(S,uj,'energyint','local','ampl',ampl);
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
        % plotSolution(S,uj,'energyint','local','ampl',ampl);
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
    duration = 10; % [s]
    framecount = getnbtimedof(T);
    framerate = framecount/duration;
    
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
    % evolSolution(S,uk,'energyint','local','ampl',ampl,'FrameRate',framerate,'filename','internal_energy_mean','pathname',pathname,options{:});
    
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
        % evolSolution(S,uk,'energyint','local','ampl',ampl,'FrameRate',framerate,'filename',['internal_energy_sample_' num2str(k)],'pathname',pathname,options{:});
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
