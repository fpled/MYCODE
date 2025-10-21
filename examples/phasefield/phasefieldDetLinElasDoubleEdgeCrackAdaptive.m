%% Phase-field fracture model - deterministic linear elasticity problem with double edge crack (two symmetric or asymmetric edge cracks) %%
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

Dim = 2; % space dimension Dim = 2, 3
setup = 1; % notch geometry setup = 1, 2, 3
PFmodel = 'Miehe'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 1; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
critConv = 'Energy'; % 'Solution', 'Residual', 'Energy'
initialCrack = 'GeometricNotch'; % 'GeometricCrack', 'GeometricNotch', 'InitialPhaseField'

suffix = '';

foldername = ['doubleEdgeCracksSetup' num2str(setup) '_' num2str(Dim) 'D'];
filename = ['linElas' PFmodel PFsplit PFregularization PFsolver initialCrack...
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
switch setup
    case 1
        hausd = 0.00004;
    case 2
        hausd = 0.00000075;
    case 3
        hausd = 0.01;
end
if Dim==2, hgrad = 1.1; elseif Dim==3, hgrad = 1.2; end
mmgoptions = ['-nomove -hausd ' num2str(hausd) ' -hgrad ' num2str(hgrad) ' -v -1'];
% gmshoptions = '-v 5';
% mmgoptions = '-nomove -hausd 0.01 -hgrad 1.3 -v 1';

%% Problem
if setProblem
    %% Domains and meshes
    switch setup
        case 1
            h = 40e-3; % domain height
            w = h; % w = 40e-3; % domain width
            a = 0.2*h; % a = 8e-3; % crack length
            e = 0.025*h; % e = 1e-3; % eccentricity from the middle line
            DIM3 = 1e-3; % thickness
        case 2
            h = 1e-3; % domain height
            w = h; % w = 1e-3; % domain width
            a = 0.2*h; % a = 0.2e-3; % crack length
            e = -0.05*h; % e = -0.05e-3; % eccentricity from the middle line
            DIM3 = 0.1e-3; % thickness
        case 3
            h = 120e-3; % domain height
            w = 60e-3; % domain width
            a = h/12; % a = 10e-3; % crack length
            % e = -h/48; % e = -2.5e-3; % eccentricity from the middle line
            % e = -h/24; % e = -5e-3; % eccentricity from the middle line
            e = -h/16; % e = -7.5e-3; % eccentricity from the middle line
            DIM3 = 10e-3; % thickness
    end
    
    if Dim==2
        D = DOMAIN(2,[0.0,0.0],[w,h]);
        Ca = LINE([0.0,h/2+e],[a,h/2+e]);
        Cb = LINE([w-a,h/2-e],[w,h/2-e]);
    elseif Dim==3
        D = DOMAIN(3,[0.0,0.0,0.0],[w,h,DIM3]);
        Ca = QUADRANGLE([0.0,h/2+e,0.0],[  a,h/2+e,0.0],[  a,h/2+e,DIM3],[0.0,h/2+e,DIM3]);
        Cb = QUADRANGLE([  w,h/2-e,0.0],[w-a,h/2-e,0.0],[w-a,h/2-e,DIM3],[  w,h/2-e,DIM3]);
    end
    
    switch setup
        case 1
            % clC = 1e-4; % [Molnar, Gravouil, 2017, FEAD], [Chen, Vasiukov, Gelebart, Park, 2019, CMAME], [Si, Yu, Li, Natarajan, 2023, CMAME]
            % clC = 5e-5; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            clD = 5e-4;
            clC = 1e-4;
            if test
                clD = 2e-3;
                if Dim==2
                    clC = 1.75e-4;
                elseif Dim==3
                    clC = 4e-4;
                end
            end
        case 2
            % clC = 3.3e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
            % clC = 2.5e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
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
        case 3
            % clD = 4.3e-3; % [Badnava, Msekh, Etemadi, Rabczuk, 2018, FEAD]
            % clC = 3e-3; % [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM]
            % clC = 1.5e-3; % [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM]
            % clC = 1.08e-3; % [Badnava, Msekh, Etemadi, Rabczuk, 2018, FEAD]
            % clC = 1e-3; % [Tong, Shen, Shao, Chen, 2020, EFM]
            % clC = 0.6e-3; % [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM]
            % clC = 0.5e-3; % [Li, Lu, Huang, Yang, 2022, OE]
            % clC = 0.1e-3; % [Han, Li, Yu, Li, Zhang, 2022, JMPS]
            clD = 2e-3;
            clC = 0.5e-3;
            if test
                if Dim==2
                    clD = 3e-3;
                    clC = 0.75e-3;
                elseif Dim==3
                    clD = 6e-3;
                    clC = 1.5e-3;
                    % clC = 1.75e-3;
                end
            end
    end
    switch lower(initialCrack)
        case 'geometriccrack'
            S_phase = gmshDomainWithDoubleEdgeCrack(D,Ca,Cb,clD,clC,fullfile(pathname,'gmsh_domain_double_edge_crack'));
        case 'geometricnotch'
            switch setup
                case 1
                    c = 4e-4; % crack width
                    clCtip = min(clC,c/2);
                    S_phase = gmshDomainWithDoubleEdgeNotch(D,Ca,Cb,c,clD,clCtip,fullfile(pathname,'gmsh_domain_double_edge_crack'));
                case 2
                    % c = 6.6e-6; % crack width
                    c = 1e-5; % crack width
                    clCtip = min(clC,c/2);
                    S_phase = gmshDomainWithDoubleEdgeNotch(D,Ca,Cb,c,clD,clCtip,fullfile(pathname,'gmsh_domain_double_edge_crack'));
                case 3
                    c = 2e-3; % notch width
                    S_phase = gmshDomainWithDoubleEdgeNotch(D,Ca,Cb,c,clD,clC,fullfile(pathname,'gmsh_domain_double_edge_crack'),Dim,'rectangular');
            end
        case 'initialphasefield'
            S_phase = gmshDomainWithDoubleEdgeCrack(D,Ca,Cb,clD,clC,fullfile(pathname,'gmsh_domain_double_edge_crack'),Dim,'noduplicate','refinecrack');
        otherwise
            error('Wrong model for initial crack');
    end
    
    sizemap = @(d) (clC-clD)*d+clD;
    % sizemap = @(d) clD*clC./((clD-clC)*d+clC);
    
    %% Phase-field problem
    %% Material
    switch setup
        case 1
            % Critical energy release rate (or fracture toughness)
            % gc = 3.4e3; % [Bhowmick, Liu, 2018, EFM]
            gc = 2.7e3; % [Molnar, Gravouil, 2017, FEAD], [Bhowmick, Liu, 2018, EFM], [Chen, Vasiukov, Gelebart, Park, 2019, CMAME], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Si, Yu, Li, Natarajan, 2023, CMAME]
            % Regularization parameter (width of the smeared crack)
            % l = 2e-4; % [Molnar, Gravouil, 2017, FEAD], [Chen, Vasiukov, Gelebart, Park, 2019, CMAME], [Si, Yu, Li, Natarajan, 2023, CMAME]
            l = 2e-4;
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
            % gc = 59; % [Nguyen, 2005, PhD thesis], [Nguyen, Houlsby, 2007, IJNAMG], [Nguyen, Houlsby, 2008, IJNAMG], [Nguyen, 2008, IJSS], [Nguyen, Korsunsky, 2008, IJSS], [Nguyen, 2011, IJSS], [Stefanou, Georgioudakis, Papadrakakis, 2014, MMUQMS], [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM], [Alfaiate, Wells, Sluys, 2002, EFM]
            % gc = 50; % [Galvez, Planas, Sancho, Reyes, Cendon, Casati, 2013, EFM]
            % gc = 25; % [Han, Li, Yu, Li, Zhang, 2022, JMPS]
            gc = 4.5; % [Badnava, Msekh, Etemadi, Rabczuk, 2018, FEAD]
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
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'DIM3',DIM3,'PFregularization',PFregularization);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    switch lower(initialCrack)
        case 'geometriccrack'
            if Dim==2
                Ca = POINT([  a,h/2+e]);   % crack tip
                Cb = POINT([w-a,h/2-e]); % crack tip
            elseif Dim==3
                Ca = LINE([  a,h/2+e,0.0],[  a,h/2+e,DIM3]); % crack tip
                Cb = LINE([w-a,h/2-e,0.0],[w-a,h/2-e,DIM3]); % crack tip
            end
        case 'geometricnotch'
            if Dim==2
                switch setup
                    case {1,2}
                        Ca = CIRCLE(  a-c/2,h/2+e,c/2); % circular notch
                        Cb = CIRCLE(w-a+c/2,h/2-e,c/2); % circular notch
                    case 3
                        Ca = LINE([  a,h/2+e-c/2],[  a,h/2+e+c/2]); % rectangular notch
                        Cb = LINE([w-a,h/2-e-c/2],[w-a,h/2-e+c/2]); % rectangular notch
                end
                % Ca = POINT([  a,h/2+e]);   % V notch
                % Cb = POINT([w-a,h/2-e]); % V notch
            elseif Dim==3
                switch setup
                    case {1,2}
                        Ca = CYLINDER(  a-c/2,h/2+e,0.0,c/2,DIM3); % circular notch
                        Cb = CYLINDER(w-a+c/2,h/2-e,0.0,c/2,DIM3); % circular notch
                    case 3
                        Ca = QUADRANGLE([  a,h/2+e-c/2,0.0],[  a,h/2+e+c/2,0.0],[  a,h/2+e+c/2,DIM3],[  a,h/2+e-c/2,DIM3]); % rectangular notch
                        Cb = QUADRANGLE([w-a,h/2-e-c/2,0.0],[w-a,h/2-e+c/2,0.0],[w-a,h/2-e+c/2,DIM3],[w-a,h/2-e-c/2,DIM3]); % rectangular notch
                end
                % Ca = LINE([  a,h/2+e,0.0],[  a,h/2+e,DIM3]); % V notch
                % Cb = LINE([w-a,h/2-e,0.0],[w-a,h/2-e,DIM3]); % V notch
            end
    end
    if Dim==2
        BRight = LINE([w,0.0],[w,h]);
        BLeft = LINE([0.0,0.0],[0.0,h]);
    elseif Dim==3
        BRight = PLANE([w,0.0,0.0],[w,h,0.0],[w,0.0,DIM3]);
        BLeft = PLANE([0.0,0.0,0.0],[0.0,h,0.0],[0.0,0.0,DIM3]);
    end
    
    addbcdamage = @(S_phase) addbcdamageDoubleEdgeCrack(S_phase,Ca,Cb,initialCrack);
    addbcdamageadapt = @(S_phase) addbcdamageDoubleEdgeCrackAdaptive(S_phase,Ca,Cb);
    findddlboundary = @(S_phase) union(findddl(S_phase,'T',BRight),findddl(S_phase,'T',BLeft));
    
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
    S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_domain_double_edge_crack'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
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
            % option = 'DEFO'; % [Ribeiro Nogueira, Rastiello, Giry, Gatuingt, Callari, 2023, AJCE]
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
    % Energetic degradation function
    g = @(d) (1-d).^2;
    % Density
    RHO = 1;
    % RHO = 2400; % [Li, Lu, Huang, Yang, 2022, OE]
    
    % Material
    d = calc_init_dirichlet(S_phase);
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BU = LINE([0.0,h],[w,h]);
        BL = LINE([0.0,0.0],[w,0.0]);
        P0 = POINT([0.0,0.0]); % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Ribeiro Nogueira, Rastiello, Giry, Gatuingt, Callari, 2023, AJCE]
        % P0 = POINT([0.0,h]); % [Tong, Shen, Shao, Chen, 2020, EFM], [Li, Lu, Huang, Yang, 2022, OE]
    elseif Dim==3
        BU = PLANE([0.0,h,0.0],[w,h,0.0],[0.0,h,DIM3]);
        BL = PLANE([0.0,0.0,0.0],[w,0.0,0.0],[0.0,0.0,DIM3]);
        P0 = POINT([0.0,0.0,0.0]); % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Ribeiro Nogueira, Rastiello, Giry, Gatuingt, Callari, 2023, AJCE]
        % P0 = POINT([0.0,h,0.0]); % [Tong, Shen, Shao, Chen, 2020, EFM], [Li, Lu, Huang, Yang, 2022, OE]
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
            % dthreshold = 0.6;
            % T = struct('dt0',dt0,'dt1',dt1,'tf',tf,'dthreshold',dthreshold);
            
            % du = 1e-4 mm during the first 400 time steps (up to u = 0.04 mm)
            % du = 1e-5 mm during the last 4000 time steps (up to u = 0.08 mm)
            dt0 = 1e-7;
            nt0 = 400;
            dt1 = 1e-8;
            nt1 = 4000;
            if test
                dt0 = 1e-6;
                nt0 = 40;
                dt1 = 1e-7;
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
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','D','Ca','Cb','addbc','addbcdamage','addbcdamageadapt','findddlforce','findddlboundary');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','D','Ca','Cb','addbc','addbcdamage','addbcdamageadapt','findddlforce','findddlboundary');
end

%% Solution
if solveProblem
    tTotal = tic;
    
    switch lower(PFsolver)
        case {'historyfieldelem','historyfieldnode'}
            [dt,ut,ft,St_phase,St,Ht,Edt,Eut,output] = solvePFDetLinElasAdaptive(S_phase,S,T,PFsolver,addbc,addbcdamage,addbcdamageadapt,findddlforce,findddlboundary,sizemap,...
                'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true,'filename','gmsh_domain_double_edge_crack','pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
        otherwise
            [dt,ut,ft,St_phase,St,~,Edt,Eut,output] = solvePFDetLinElasAdaptive(S_phase,S,T,PFsolver,addbc,addbcdamage,addbcdamageadapt,findddlforce,findddlboundary,sizemap,...
                'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true,'filename','gmsh_domain_double_edge_crack','pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    end
    % switch lower(PFsolver)
    %     case {'historyfieldelem','historyfieldnode'}
    %         [dt,ut,ft,St_phase,St,Ht,Edt,Eut,output] = solvePFDetLinElasDoubleEdgeCrackAdaptive(S_phase,S,T,PFsolver,Ca,Cb,BU,BL,P0,BRight,BLeft,setup,initialCrack,sizemap,...
    %             'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true,'filename','gmsh_domain_double_edge_crack','pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    %     otherwise
    %         [dt,ut,ft,St_phase,St,Ht,Edt,Eut,output] = solvePFDetLinElasDoubleEdgeCrackAdaptive(S_phase,S,T,PFsolver,Ca,Cb,BU,BL,P0,BRight,BLeft,setup,initialCrack,sizemap,...
    %             'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true,'filename','gmsh_domain_double_edge_crack','pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    % end
    
    t = gettevol(T);
    dmaxt = cellfun(@(d) max(d),dt);
    idc = find(dmaxt>=min(0.75,max(dmaxt)),1);
    fc = ft(idc);
    udc = t(idc);
    [fmax,idmax] = max(ft,[],2);
    udmax = t(idmax);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'dt','ut','ft','St_phase','St','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        save(fullfile(pathname,'solution.mat'),'Ht','-append');
    end
else
    load(fullfile(pathname,'solution.mat'),'dt','ut','ft','St_phase','St','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        load(fullfile(pathname,'solution.mat'),'Ht');
    end
end

%% Outputs
if solveProblem
    fid = fopen(fullfile(pathname,'results.txt'),'w');
    fprintf(fid,'Double edge crack');
    if e==0
        fprintf(fid,' (two symmetric edge cracks)\n');
    else
        fprintf(fid,' (two asymmetric edge cracks)\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'setup    = %d\n',setup);
    fprintf(fid,'dim      = %d\n',Dim);
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
    % mysaveas(pathname,'mesh_init',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor',facecolor,'FaceAlpha',facealpha,'legend',false);
    mysaveas(pathname,'mesh_init',formats,renderer);
    
    % u = ut{rep(end)};
    u = ut{end};
    S_final = St{end};
    
    plotModel(S_final,'Color','k','FaceColor',facecolor,'FaceAlpha',facealpha,'legend',false);
    mysaveas(pathname,'mesh_final',formats,renderer);
    
    ampl = getsize(S_final)/max(abs(u))/20;
    plotModelDeflection(S_final,u,'ampl',ampl,'Color','b','FaceColor',facecolordef,'FaceAlpha',facealpha,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor',facecolor,'FaceAlpha',facealpha);
    plot(S_final+ampl*unfreevector(S_final,u),'Color','b','FaceColor',facecolordef,'FaceAlpha',facealpha);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
end

%% Display solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display force-displacement curve
    figure('Name','Force vs displacement')
    clf
    switch setup
        case {1,3}
            plot(t*1e3,ft*1e-3,'-b','LineWidth',linewidth)
            % hold on
            % scatter(udmax*1e3,fmax*1e-3,'Marker','+','MarkerEdgeColor','b','LineWidth',linewidth)
            % scatter(udc*1e3,fc*1e-3,'Marker','+','MarkerEdgeColor','r','LineWidth',linewidth)
            % hold off
        case 2
            plot(t*1e3,ft,'-b','LineWidth',linewidth)
            % hold on
            % scatter(udmax*1e3,fmax,'Marker','+','MarkerEdgeColor','b','LineWidth',linewidth)
            % scatter(udc*1e3,fc,'Marker','+','MarkerEdgeColor','r','LineWidth',linewidth)
            % hold off
    end
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
        case 1
            tSnapshots = [0.0450 0.0475 0.050 0.0525 0.055 0.06]*1e-3;
        case 2
            tSnapshots = [0.055 0.06 0.065 0.07 0.08 0.10 0.12]*1e-3;
        case 3
            tSnapshots = [0.005 0.01 0.015 0.02 0.05 0.10 0.15]*1e-3;
    end
    rep = arrayfun(@(x) find(t>x-eps,1),tSnapshots);
    rep = [rep,length(T)];
    % tSnapshots = [tSnapshots,gett1(T)];
    % rep = arrayfun(@(x) find(t>x-eps,1),tSnapshots);
    
    for j=1:length(rep)
        dj = dt{rep(j)};
        uj = ut{rep(j)};
        Sj = St{rep(j)};
        Sj_phase = St_phase{rep(j)};
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            Hj = Ht{rep(j)};
        end
        
        plotModel(Sj,'Color','k','FaceColor',facecolor,'FaceAlpha',facealpha,'legend',false);
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
        % plotSolution(Sj,uj,'energyint','local','ampl',ampl);
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
    % evolSolutionCell(T,St,ut,'energyint','local','ampl',ampl,'FrameRate',framerate,'filename','internal_energy_density','pathname',pathname,options{:});
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
