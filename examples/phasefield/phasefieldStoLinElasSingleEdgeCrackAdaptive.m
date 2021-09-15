%% Phase field fracture model - stochastic linear elasticity problem with single edge crack %%
%%------------------------------------------------------------------------------------------%%
% [Bourdin, Francfort, Marigo, 2000, JMPS] (isotropic phase field model with no split of Bourdin et al.)
% [Miehe, Welschinger, Hofacker, 2010 IJNME] (anisotropic phase field model of Miehe et al.)
% [Miehe, Hofacker, Welschinger, 2010, CMAME] (anisotropic phase field model of Miehe et al.)
% [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME] (anisotropic phase field model of Miehe et al.)
% [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM] (anisotropic phase field model of Miehe et al.)
% [Ambati, Gerasimov, De Lorenzis, 2015, CM] (hybrid isotropic-anisotropic phase field model of Ambati et al. compared with the isotropic one of Bourdin et al. and the anisotropic ones of Amor et al. and Miehe et al.)
% [Liu, Li, Msekh, Zuo, 2016, CMS] (anisotropic phase field model of Miehe et al.)
% [Wu, Nguyen, 2018, JMPS] (hybrid isotropic-anisotropic phase field model of Wu et al.)
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM] (anisotropic phase field model of Wu et al.)
% [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US] (anisotropic phase field model of Amor et al.)
% [Wu, Nguyen, Zhou, Huang, 2020, CMAME] (anisotropic phase field model of Wu et al.)
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME] (anisotropic phase field model of Nguyen et al.)

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
numWorkers = 20;

% Deterministic model parameters
Dim = 2; % space dimension Dim = 2, 3
symmetry = 'Isotropic'; % 'Isotropic' or 'Anisotropic'. Material symmetry
ang = 30; % clockwise material orientation angle around z-axis [deg]
isotropicTest = false; % for test purposes (configuration of isotropic material with the anisotropic class). Work only for "Dim = 2" and "symmetry = 'Anisotropic'".
loading = 'Tension'; % 'Tension' or 'Shear'
PFmodel = 'Isotropic'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'

% Random model parameters
N = 5e2; % number of samples
randMat = true; % random material parameters (true or false)
randPF = true; % random phase field parameters (true or false)
correlationStructure = true;
rhoBP = 0.1; % Bravais-Pearson correlation coefficient

switch lower(symmetry)
    case 'isotropic' % isotropic material
        filename = ['phasefieldStoLinElas' symmetry 'SingleEdgeCrack' loading PFmodel];
    case 'anisotropic' % anisotropic material
        filename = ['phasefieldStoLinElas' symmetry num2str(ang) 'deg' 'SingleEdgeCrack' loading PFmode];
    otherwise
        error('Wrong material symmetry class');
end
if isotropicTest
    filename = ['phasefieldStoLinElas' 'IsotTest' 'SingleEdgeCrack' loading PFmodel];
end
if randMat
    filename = [filename 'RandMat'];
end
if randPF
    filename = [filename 'RandPF'];
end
filename = [filename 'Adaptive_' num2str(Dim) 'D_' num2str(N) 'samples'];

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','phasefield',filename);
if test
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','phasefield_test',filename);
end
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

gmshoptions = '-v 0';
mmgoptions = '-nomove -hausd 0.01 -hgrad 1.1 -v -1';
% gmshoptions = '-v 5';
% mmgoptions = '-nomove -hausd 0.01 -hgrad 1.3 -v 1';

%% Problem
if setProblem
    %% Domains and meshes
    L = 1e-3;
    a = L/2;
    if Dim==2
        e = 1;
        D = DOMAIN(2,[0.0,0.0],[L,L]);
        C = LIGNE([0.0,L/2],[a,L/2]);
    elseif Dim==3
        e = 0.1e-3;
        D = DOMAIN(3,[0.0,0.0,0.0],[L,L,e]);
        C = QUADRANGLE([0.0,L/2,0.0],[a,L/2,0.0],[a,L/2,e],[0.0,L/2,e]);
    end
    
    if Dim==2
        % clD = 6.25e-5; % [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME]
        % clD = 3e-5; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        clD = 2e-5; % [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
        % clD = 3.9e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
        % clD = 2e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        % clD = 1e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        
        % clC = 3.906e-6; % [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME]
        % clC = 2.5e-6; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        clC = 2e-6; % (shear test) [Miehe, Hofacker, Welschinger, 2010, CMAME]
        % clC = 1e-6; % (tension test) [Miehe, Welschinger, Hofacker, 2010 IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME]
        % clC = 6e-7; % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
        % clC = 3.9e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
        % clC = 2e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        % clC = 1e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        if test
            clD = 4e-5;
            clC = 1e-5;
            % clD = 1e-5;
            % clC = 1e-5;
        end
    elseif Dim==3
        clD = 4e-5;
        clC = 4e-6;
        % clD = 7.5e-6;
        % clC = 7.5e-6;
        if test
            clD = 4e-5;
            clC = 1e-5;
            % clD = 2e-5;
            % clC = 2e-5;
        end
    end
    c = clC; % crack width
    S_phase = gmshdomainwithedgesmearedcrack(D,C,c,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,'gmshoptions',gmshoptions);
    
    sizemap = @(d) (clC-clD)*d+clD;
    % sizemap = @(d) clD*clC./((clD-clC)*d+clC);
    
    %% Phase field problem
    %% Material
    switch lower(symmetry)
        case 'isotropic' % isotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 2.7e3; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
            % Regularization parameter (width of the smeared crack)
            % l = 3.75e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME]
            % l = 3e-5; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
            % l = 1.5e-5; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
            % l = 1e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
            l = 7.5e-6; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            % l = 5e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
            % l = 4e-6; % [Ambati, Gerasimov, De Lorenzis, 2015, CM]
            % eta = 0.052; w0 = 75.94; l = eta/sqrt(w0)*1e-3; % l = 6e-7; % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
        case 'anisotropic' % anisotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 10e3; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            % Regularization parameter (width of the smeared crack)
            l = 8.5e-6; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        otherwise
            error('Wrong material symmetry class');
    end
    if isotropicTest
        % Critical energy release rate (or fracture toughness)
        gc = 2.7e3; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
        % Regularization parameter (width of the smeared crack)
        l = 7.5e-6; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
    end
    % Small artificial residual stiffness
    k = 1e-10;
    % Internal energy
    H = 0;
    
    % Material
    mat_phase = FOUR_ISOT('k',gc*l,'r',gc/l+2*H);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    if Dim==2
        C = DOMAIN(2,[0.0,L/2-c/2]-[eps,eps],[a,L/2+c/2]+[eps,eps]);
    elseif Dim==3
        C = DOMAIN(3,[0.0,L/2-c/2,0.0]-[eps,eps,eps],[a,L/2+c/2,e]+[eps,eps,eps]);
    end
    S_phase = final(S_phase,'duplicate');
    S_phase = addcl(S_phase,C,'T',1);
    
    d = calc_init_dirichlet(S_phase);
    cl = sizemap(d);
    S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_domain_single_edge_crack'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    S = S_phase;
    
    S_phase = setmaterial(S_phase,mat_phase);
    S_phase = final(S_phase,'duplicate');
    S_phase = addcl(S_phase,C,'T',1);
    
    %% Stiffness matrices and sollicitation vectors
    % a_phase = BILINFORM(1,1,gc*l); % uniform values
    % % a_phase = DIFFUSIONFORM(gc*l);
    % a_phase = setfree(a_phase,0);
    % K_phase = calc_matrix(a_phase,S_phase);
    % b_phase = calc_nonhomogeneous_vector(S_phase,K_phase);
    % b_phase = -b_phase;
    % K_phase = freematrix(S_phase,K_phase);
    
    % r_phase = BILINFORM(0,0,gc/l+2*H,0); % nodal values
    % R_phase = calc_matrix(r_phase,S_phase);
    % A_phase = K_phase + M_phase;
    
    % l_phase = LINFORM(0,2*H,0); % nodal values
    % l_phase = setfree(l_phase,1);
    % b_phase = b_phase + calc_vector(l_phase,S_phase);
    
    % [A_phase,b_phase] = calc_rigi(S_phase);
    % b_phase = -b_phase + bodyload(S_phase,[],'QN',2*H);
    
    %% Linear elastic displacement field problem
    %% Materials
    % Option
    option = 'DEFO'; % plane strain [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
    % option = 'CONT'; % plane stress [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
    switch lower(symmetry)
        case 'isotropic' % isotropic material
            % Lame coefficients
            % lambda = 121.1538e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
            % mu = 80.7692e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
            lambda = 121.15e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
            mu = 80.77e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
            % Young modulus and Poisson ratio
            if Dim==2
                switch lower(option)
                    case 'defo'
                        E = mu*(3*lambda+2*mu)/(lambda+mu); %  E = 210e9;
                        NU = lambda/(lambda+mu)/2; % NU = 0.3;
                    case 'cont'
                        E = 4*mu*(lambda+mu)/(lambda+2*mu);
                        NU = lambda/(lambda+2*mu);
                end
                % E = 210e9; NU = 0.2; % [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
                % E = 210e9; NU = 0.3; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
                % kappa = 121030e6; NU=0.227; lambda=3*kappa*NU/(1+NU); mu = 3*kappa*(1-2*NU)/(2*(1+NU)); E = 3*kappa*(1-2*NU); % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
            elseif Dim==3
                E = mu*(3*lambda+2*mu)/(lambda+mu);
                NU = lambda/(lambda+mu)/2;
            end
            
        case 'anisotropic' % anisotropic material
            if Dim==2
                switch lower(option)
                    case 'defo'
                        % [Nguyen, Yvonnet, Waldmann, He, 2020,IJNME]
                        % Elasticity matrix in reference material coordinate system [Pa]
                        matElas = 1e9*[65 20 0;
                            20 260 0;
                            0 0 30];
                        theta = deg2rad(ang); % clockwise material orientation angle around z-axis [rad]
                        c = cos(theta);
                        s = sin(theta);
                        % Transition matrix for elasticity matrix from material coordinate system to global coordinate system
                        P = [c^2 s^2 -c*s;
                            s^2 c^2 c*s;
                            2*c*s -2*c*s c^2-s^2];
                        % Elasticity matrix in global coordinate system [Pa]
                        matElas = P'*matElas*P;
                    case 'cont'
                        error('Not implemented yet')
                end
                
                if isotropicTest
                    lambda = 121.15e9;
                    mu = 80.77e9;
                    if strcmpi(option,'cont')
                        E = mu*(3*lambda+2*mu)/(lambda+mu);
                        NU = lambda/(lambda+mu)/2;
                        lambda = E*NU/(1-NU^2); % first Lamé coefficient
                    end
                    matElas = [lambda+2*mu,lambda,0;...
                        lambda,lambda+2*mu,0;...
                        0,0,mu];
                end
                
            elseif Dim==3
                error('Not implemented yet')
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
        case 'isotropic' % isotropic material model for isotropic material only
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel);
        case 'anisotropic' % anisotropic material model for all symmetry classes
            mat = ELAS_ANISOT('matElas',matElas,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel);
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
    
    S = final(S,'duplicate');
    
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
            error('Wrong loading case')
    end
    
    %% Stiffness matrices and sollicitation vectors
    % [A,b] = calc_rigi(S);
    % b = -b;
    
    %% Time scheme
    switch lower(symmetry)
        case 'isotropic' % isotropic material
            if Dim==2
                switch lower(loading)
                    case 'tension'
                        % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM]
                        % du = 1e-5 mm during the first 500 time steps (up to u = 5e-3 mm)
                        % du = 1e-6 mm during the last 1300 time steps (up to u = 6.3e-3 mm)
                        dt0 = 1e-8;
                        nt0 = 500;
                        dt1 = 1e-9;
                        nt1 = 1300;
                        if test
                            dt0 = 1e-7;
                            nt0 = 50;
                            dt1 = 1e-8;
                            nt1 = 400;
                        end
                        t0 = linspace(dt0,nt0*dt0,nt0);
                        t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        t = [t0,t1];
                        
                        % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
                        % du = 1e-5 mm during 630 time steps (up to u = 6.3e-3 mm)
                        % dt = 1e-8;
                        % nt = 630;
                        % t = linspace(dt,nt*dt,nt);
                        
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
                    case 'shear'
                        % [Miehe, Welschinger, Hofacker, 2010 IJNME]
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
                        
                        % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
                        % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM],
                        % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
                        % du = 1e-5 mm during 1500 time steps (up to u = 15e-3 mm)
                        dt = 1e-8;
                        nt = 1500;
                        % nt = 2000;
                        if test
                            dt = 5e-8;
                            % nt = 300;
                            nt = 400;
                        end
                        t = linspace(dt,nt*dt,nt);
                end
            elseif Dim==3
                dt = 1e-8;
                nt = 2500;
                if test
                    dt = 1e-7;
                    nt = 250;
                end
                t = linspace(dt,nt*dt,nt);
            end
            
        case 'anisotropic' % anisotropic material
            if Dim==2
                switch lower(loading)
                    case 'tension'
                        % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
                        % du = 6e-5 mm during the first 200 time steps (up to u = 12e-3 mm)
                        % du = 2e-5 mm during the last 600 time steps (up to u = 24e-3 mm)
                        dt0 = 6e-8;
                        nt0 = 200;
                        dt1 = 2e-8;
                        nt1 = 600;
                        if test
                            dt0 = 6e-7;
                            nt0 = 20;
                            dt1 = 2e-7;
                            nt1 = 50;
                        end
                        
                    case 'shear'
                        % du = 1e-4 mm during the first 200 time steps (up to u = 20e-3 mm)
                        % du = 2e-5 mm during the last 2000 time steps (up to u = 60e-3 mm)
                        dt0 = 1e-7;
                        nt0 = 200;
                        dt1 = 2e-8;
                        nt1 = 2000;
                        if test
                            dt0 = 1e-6;
                            nt0 = 20;
                            dt1 = 2e-7;
                            nt1 = 200;
                        end
                end
                t0 = linspace(dt0,nt0*dt0,nt0);
                t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                t = [t0,t1];
                
                if isotropicTest
                    switch lower(loading)
                        case 'tension'
                            dt0 = 1e-8;
                            nt0 = 500;
                            dt1 = 1e-9;
                            nt1 = 1300;
                            if test
                                dt0 = 1e-7;
                                nt0 = 50;
                                dt1 = 1e-8;
                                nt1 = 400;
                            end
                            t0 = linspace(dt0,nt0*dt0,nt0);
                            t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                            t = [t0,t1];
                        case 'shear'
                            dt = 1e-8;
                            nt = 1500;
                            if test
                                dt = 5e-8;
                                nt = 400;
                            end
                            t = linspace(dt,nt*dt,nt);
                    end
                end
                
            elseif Dim==3
                dt = 1e-8;
                nt = 2500;
                if test
                    dt = 1e-7;
                    nt = 250;
                end
                t = linspace(dt,nt*dt,nt);
            end
        otherwise
            error('Wrong material symmetry class');
    end
    T = TIMEMODEL(t);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','D','C','BU','BL','BRight','BLeft','BFront','BBack','loading');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','D','C','BU','BL','BRight','BLeft','BFront','BBack','loading');
end

%% Solution
if solveProblem
    myparallel('start',numWorkers);
    %% Random variables
    % Material properties
    if randMat % random material parameters
        % la = -24; % la < 1/5. Parameter controlling the level of statistical fluctuation
        % deltaC1 = 1/sqrt(1-la); % coefficient of variation for bulk modulus
        % deltaC2 = 1/sqrt(1-5*la); % coefficient of variation for shear modulus
        deltaC1 = 0.1; % coefficient of variation for bulk modulus
        la = 1 - 1/deltaC1^2; % la < 1/5. Parameter controlling the level of statistical fluctuation
        deltaC2 = 1/sqrt(5/deltaC1^2 - 4); % coefficient of variation for shear modulus
        
        mC1 = E/3/(1-2*NU); % mean bulk modulus
        mC2 = mu; % mean shear modulus
        laC1 = (1-la)/mC1; % la1 > 0
        laC2 = (1-5*la)/mC2; % la2 > 0
        
        aC1 = 1-la; % a1 > 0
        bC1 = 1/laC1; % b1 > 0
        aC2 = 1-5*la; % a2 > 0
        bC2 = 1/laC2; % b2 > 0
        
        % Sample set
        C_sample(:,1) = gamrnd(aC1,bC1,N,1); % samples for bulk modulus [Pa]
        C_sample(:,2) = gamrnd(aC2,bC2,N,1); % samples for shear modulus [Pa]
        % lambda_sample = C_sample(:,1) - 2/3*C_sample(:,2); % [Pa]
        E_sample = (9*C_sample(:,1).*C_sample(:,2))./(3*C_sample(:,1)+C_sample(:,2)); % [Pa]
        NU_sample = (3*C_sample(:,1)-2*C_sample(:,2))./(6*C_sample(:,1)+2*C_sample(:,2));
    else
        E_sample = E*ones(N,1);
        NU_sample = NU*ones(N,1);
    end
    
    % Phase field properties
    if randPF % random phase field parameters
        deltaP1 = 0.1; % coefficient of variation of fracture toughness
        deltaP2 = 0.1; % coefficient of variation of regularization parameter
        aP1 = 1/deltaP1^2;
        bP1 = gc/aP1;
        aP2 = 1/deltaP2^2;
        bP2 = l/aP2;
        switch correlationStructure
            case false
                gc_sample = gamrnd(aP1,bP1,N,1); % samples for fracture toughness [N/m^2]
                l_sample = gamrnd(aP2,bP2,N,1); % samples regularization parameter [m]
            case true
                Xi = randn(N,2); % random matrix with statistically independent normalized Gaussian components
                gc_sample = gaminv(normcdf(Xi(:,1)),aP1,bP1);
                l_sample = gaminv(normcdf(rhoBP*Xi(:,1) + sqrt(1-rhoBP^2)*Xi(:,2)),aP2,bP2);
        end
    else
        gc_sample = gc*ones(N,1);
        l_sample = l*ones(N,1);
    end
    
    samples = [E_sample,NU_sample,gc_sample,l_sample];

    %% Solution
    tTotal = tic;
    
    nbSamples = 3;
    fun = @(S_phase,S,filename) solvePFDetLinElasSingleEdgeCrackAdaptive(S_phase,S,T,C,BU,BL,BRight,BLeft,BFront,BBack,loading,sizemap,'filename',filename,'pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    [ft,dt,ut,St_phase,St] = solvePFStoLinElasAdaptive(S_phase,S,T,fun,samples,'filename','gmsh_domain_single_edge_crack','pathname',pathname,'nbsamples',nbSamples);
    fmax = max(ft,[],2);
    
    time = toc(tTotal);
    
    %% Statistical outputs of solution
    probs = [0.025 0.975];
    
    mean_ft = mean(ft);
    std_ft = std(ft);
    ci_ft = quantile(ft,probs);
    
    mean_fmax = mean(fmax);
    std_fmax = std(fmax);
    ci_fmax = quantile(fmax,probs);
    
    npts = 100;
    [f_fmax,xi_fmax,bw_fmax] = ksdensity(fmax,'npoints',npts);
    
    save(fullfile(pathname,'solution.mat'),'N','dt','ut','St_phase','St',...
        'mean_ft','std_ft','ci_ft','fmax',...
        'mean_fmax','std_fmax','ci_fmax','probs','f_fmax','xi_fmax','bw_fmax','time');
else
    load(fullfile(pathname,'solution.mat'),'N','dt','ut','St_phase','St',...
        'mean_ft','std_ft','ci_ft','fmax',...
        'mean_fmax','std_fmax','ci_fmax','probs','f_fmax','xi_fmax','bw_fmax','time');
end

%% Outputs
fprintf('\n');
fprintf('dim      = %d\n',Dim);
fprintf('loading  = %s\n',loading);
if isotropicTest
    fprintf('mat sym  = isotropic test\n');
else
    fprintf('mat sym  = %s\n',symmetry);
    if strcmpi(symmetry,'anisotropic')
        fprintf('angle    = %g deg\n',ang);
    end
end
fprintf('PF model = %s\n',PFmodel);
fprintf('nb elements = %g (initial) - %g (final)\n',getnbelem(S),getnbelem(St{end}));
fprintf('nb nodes    = %g (initial) - %g (final)\n',getnbnode(S),getnbnode(St{end}));
fprintf('nb dofs     = %g (initial) - %g (final)\n',getnbddl(S),getnbddl(St{end}));
fprintf('nb time dofs = %g\n',getnbtimedof(T));
fprintf('nb samples = %g\n',N);
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

if Dim==2
    fprintf('mean(fmax)    = %g kN/mm\n',mean_fmax*1e-6);
    fprintf('std(fmax)     = %g kN/mm\n',std_fmax*1e-6);
elseif Dim==3
    fprintf('mean(fmax)    = %g kN\n',mean_fmax*1e-3);
    fprintf('std(fmax)     = %g kN\n',std_fmax*1e-3);
end
fprintf('disp(fmax)    = %g\n',std_fmax/mean_fmax);
if Dim==2
    fprintf('%d%% ci(fmax)  = [%g,%g] kN/mm\n',(probs(2)-probs(1))*100,ci_fmax(1)*1e-6,ci_fmax(2)*1e-6);
elseif Dim==3
    fprintf('%d%% ci(fmax)  = [%g,%g] kN\n',(probs(2)-probs(1))*100,ci_fmax(1)*1e-3,ci_fmax(2)*1e-3);
end

%% Display
if displayModel
    [t,rep] = gettevol(T);
    
    %% Display domains, boundary conditions and meshes
    figure('Name','Domain')
    clf
    plot(D,'FaceColor',getfacecolor(1));
    plot(C,'FaceColor','w');
    axis image
    axis off
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
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

%% Display samples of solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display force-displacement curve
    figure('Name','Force-displacement')
    clf
    ciplot(ci_ft(1,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),ci_ft(2,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),t*1e3,'b');
    alpha(0.2)
    hold on
    plot(t*1e3,mean_ft*((Dim==2)*1e-6+(Dim==3)*1e-3),'-b','Linewidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    l = legend({['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value'},'Location','NorthWest');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    %% Display pdf of critical force
    figure('Name','Probability Density Estimate: Critical force')
    clf
    plot(xi_fmax*((Dim==2)*1e-6+(Dim==3)*1e-3),f_fmax,'-b','LineWidth',linewidth)
    hold on
    ind_fmax = find(xi_fmax>=ci_fmax(1) & xi_fmax<ci_fmax(2));
    area(xi_fmax(ind_fmax)*((Dim==2)*1e-6+(Dim==3)*1e-3),f_fmax(ind_fmax),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(mean_fmax*((Dim==2)*1e-6+(Dim==3)*1e-3),0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b')
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
    mysaveas(pathname,'pdf_fmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_fmax.tex');
    
    %% Display samples of solutions at different instants
    ampl = 0;
    switch lower(symmetry)
        case 'isotropic'
            switch lower(loading)
                case 'tension'
                    rep = find(abs(t-5.5e-6)<eps | abs(t-5.75e-5)<eps | abs(t-6e-6)<eps | abs(t-6.25e-6)<eps);
                case 'shear'
                    rep = find(abs(t-1e-5)<eps | abs(t-1.25e-5)<eps | abs(t-1.35e-5)<eps | abs(t-1.5e-5)<eps);
                otherwise
                    error('Wrong loading case')
            end
        case 'anisotropic'
            switch lower(loading)
                case 'tension'
                    rep = find(abs(t-9e-6)<eps | abs(t-12e-6)<eps | abs(t-13.5e-6)<eps | abs(t-15e-6)<eps | abs(t-20e-6)<eps);
                case 'shear'
                    rep = find(abs(t-20e-6)<eps | abs(t-30e-6)<eps | abs(t-40e-6)<eps | abs(t-50e-6)<eps);
                otherwise
                    error('Wrong loading case')
            end
        otherwise
            error('Wrong material symmetry class');
    end
    rep = [rep,length(T)];
    
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
        for i=1:Dim
            evolSolutionCell(T,St(k,:),ut(k,:),'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i) '_sample_' num2str(k)],'pathname',pathname,options{:});
        end
        
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

myparallel('stop');
