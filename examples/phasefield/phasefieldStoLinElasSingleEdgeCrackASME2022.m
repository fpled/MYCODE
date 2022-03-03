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
saveParaview = true;

% test = true; % coarse mesh
test = false; % fine mesh
numWorkers = 100;
% numWorkers = 1; maxNumCompThreads(1); % mono-thread computation

% Deterministic model parameters
Dim = 2; % space dimension Dim = 2, 3
symmetry = 'Isotropic'; % 'Isotropic', 'MeanIsotropic', 'Anisotropic'. Material symmetry
ang = 30; % clockwise material orientation angle around z-axis for anisotopic material [deg]
loading = 'Shear'; % 'Tension' or 'Shear'
PFmodel = 'AnisotropicMiehe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'
PFsolver = 'HistoryFieldElem'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'

% Random model parameters
N = 100; % number of samples
% N = numWorkers;
coeff_gc = 1.0;
randMat = struct('delta',0.1,'lcorr',1e-5,'rcorr',0); % random material parameters model
randPF = struct('delta',0,'lcorr',Inf,'rcorr',0); % random phase field parameters model

switch lower(symmetry)
    case {'isotropic','meanisotropic'} % almost surely or mean isotropic material
        filename = ['phasefieldStoLinElas' symmetry 'SingleEdgeCrack' loading PFmodel PFsolver];
    case 'anisotropic' % anisotropic material
        filename = ['phasefieldStoLinElas' symmetry num2str(ang) 'deg' 'SingleEdgeCrack' loading PFmodel PFsolver];
    otherwise
        error('Wrong material symmetry class');
end
filename = [filename '_' num2str(Dim) 'D_' num2str(N) 'samples'];
if any(randMat.delta)
    filename = [filename '_RandMat_Delta' num2str(randMat.delta,'_%g') '_Lcorr' num2str(randMat.lcorr,'_%g') '_Rcorr' num2str(randMat.rcorr,'_%g')];
end
if randPF.delta
    filename = [filename '_RandPF_Delta' num2str(randPF.delta,'_%g') '_Lcorr' num2str(randPF.lcorr,'_%g') '_Rcorr' num2str(randPF.rcorr,'_%g')];
end
filename = [filename '_coeffgc' num2str(coeff_gc,'_%g')];

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

%% Problem
if setProblem
    %% Domains and meshes
    L = 1e-3;
    a = L/2;
    if Dim==2
        e = 1;
        % P1 is the crack tip position
        % P2 = [0.0,L/2];
        % P3 = [0.0,0.0];
        % P4 = [0.99*a,0.0];
        % P5 = [L,0.0];
        % P6 = [L,L];
        % P7 = [0.99*a,L];
        % P8 = [0.0,L];
        % D = POLYGON(P2,P3,P4,P5,P6,P7,P8);
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
        % clD = 2e-5; % [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
        % clD = 3.9e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
        % clD = 2e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        % clD = 1e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        
        % clC = 3.906e-6; % [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME]
        % clC = 2.5e-6; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        % clC = 2e-6; % (shear test) [Miehe, Hofacker, Welschinger, 2010, CMAME]
        % clC = 1e-6; % (tension test) [Miehe, Welschinger, Hofacker, 2010 IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME]
        % clC = 6e-7; % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
        % clC = 3.9e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
        % clC = 2e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        % clC = 1e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        clD = 2.5e-6;
        clC = 2.5e-6;
        if test
            % clD = 4e-5;
            % clC = 1e-5;
            clD = 1e-5;
            clC = 1e-5;
        end
    elseif Dim==3
        % clD = 4e-5;
        % clC = 4e-6;
        clD = 7.5e-6;
        clC = 7.5e-6;
        if test
            % clD = 4e-5;
            % clC = 1e-5;
            clD = 2e-5;
            clC = 2e-5;
        end
    end
    S_phase = gmshdomainwithedgecrack(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'));
    % S_phase = gmshdomainwithedgecrackoptim(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'));
    S = S_phase;
    
    %% Phase field problem
    %% Material
    switch lower(symmetry)
        case {'isotropic','meanisotropic'} % almost surely or mean isotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 2.7e3; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
            % Regularization parameter (width of the smeared crack)
            % l = 3.75e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME]
            % l = 3e-5; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
            % l = 1.5e-5; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
            % l = 1e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
            % l = 7.5e-6; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            % l = 5e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
            % l = 4e-6; % [Ambati, Gerasimov, De Lorenzis, 2015, CM]
            % eta = 0.052; w0 = 75.94; l = eta/sqrt(w0)*1e-3; % l = 6e-7; % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
            l = 1e-5;
        case 'anisotropic' % anisotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 10e3; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            % Regularization parameter (width of the smeared crack)
            l = 8.5e-6; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        otherwise
            error('Wrong material symmetry class');
    end
    gc = gc*coeff_gc;
    % Small artificial residual stiffness
    k = 1e-10;
    % Internal energy
    H = 0;
    
    % Material
    mat_phase = FOUR_ISOT('k',gc*l,'r',gc/l+2*H,'delta',randPF.delta,'lcorr',randPF.lcorr,'rcorr',randPF.rcorr);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
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
    % A_phase = K_phase + R_phase;
    
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
        case {'isotropic','meanisotropic'} % almost surely or mean isotropic material
            % Lame coefficients
            % lambda = 121.1538e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
            % mu = 80.7692e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
            lambda = 121.15e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
            mu = 80.77e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
            if strcmpi(symmetry,'isotropic')
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
                    % E = 210e9; NU = 0.3; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
                    % kappa = 121030e6; NU=0.227; lambda=3*kappa*NU/(1+NU); mu = 3*kappa*(1-2*NU)/(2*(1+NU)); E = 3*kappa*(1-2*NU); % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
                elseif Dim==3
                    E = mu*(3*lambda+2*mu)/(lambda+mu);
                    NU = lambda/(lambda+mu)/2;
                end
            elseif strcmpi(symmetry,'meanisotropic')
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
            
        case 'anisotropic' % anisotropic material
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
                        error('Not implemented yet')
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
        case 'isotropic' % almost surely isotropic material
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'delta',randMat.delta,'lcorr',randMat.lcorr,'rcorr',randMat.rcorr);
        case {'meanisotropic','anisotropic'} % mean isotropic or anisotropic material
            mat = ELAS_ANISOT('C',Cmat,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'delta',randMat.delta,'lcorr',randMat.lcorr,'rcorr',randMat.rcorr);
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
        case {'isotropic','meanisotropic'} % almost surely or mean isotropic material
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
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','BU','BL','BRight','BLeft','BFront','BBack','loading','symmetry','ang');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','BU','BL','BRight','BLeft','BFront','BBack','loading','symmetry','ang');
end

%% Solution
if solveProblem
    myparallel('start',numWorkers);
    
    %% Solution
    tTotal = tic;
    
    nbSamples = 1;
    fun = @(S_phase,S) solvePFDetLinElasSingleEdgeCrack(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,'display');
    [ft,dt_mean,ut_mean,dt_var,ut_var,dt_sample,ut_sample] = solvePFStoLinElas(S_phase,S,T,fun,N,'nbsamples',nbSamples);
    [fmax,idmax] = max(ft,[],2);
    t = gettevol(T);
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
    
    npts = 100;
    [fmax_f,fmax_xi,fmax_bw] = ksdensity(fmax,'npoints',npts);
    [udmax_f,udmax_xi,udmax_bw] = ksdensity(udmax,'npoints',npts);
    
    save(fullfile(pathname,'solution.mat'),'N','dt_mean','ut_mean',...
        'dt_var','ut_var','dt_sample','ut_sample','ft_mean','ft_std','ft_ci','probs',...
        'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
        'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw','time');
else
    load(fullfile(pathname,'solution.mat'),'N','dt_mean','ut_mean',...
        'dt_var','ut_var','dt_sample','ut_sample','ft_mean','ft_std','ft_ci','probs',...
        'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
        'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw','time');
end

%% Outputs
fprintf('\n');
fprintf('dim      = %d\n',Dim);
fprintf('loading  = %s\n',loading);
fprintf('mat sym  = %s\n',symmetry);
if strcmpi(symmetry,'anisotropic')
    fprintf('angle    = %g deg\n',ang);
end
fprintf('PF model = %s\n',PFmodel);
fprintf('PF solver = %s\n',PFsolver);
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('nb time dofs = %g\n',getnbtimedof(T));
fprintf('nb samples = %g\n',N);
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

if Dim==2
    fprintf('mean(fmax)    = %g kN/mm\n',fmax_mean*1e-6);
    fprintf('std(fmax)     = %g kN/mm\n',fmax_std*1e-6);
elseif Dim==3
    fprintf('mean(fmax)    = %g kN\n',fmax_mean*1e-3);
    fprintf('std(fmax)     = %g kN\n',fmax_std*1e-3);
end
fprintf('disp(fmax)    = %g\n',fmax_std/fmax_mean);
if Dim==2
    fprintf('%d%% ci(fmax)  = [%g,%g] kN/mm\n',(probs(2)-probs(1))*100,fmax_ci(1)*1e-6,fmax_ci(2)*1e-6);
elseif Dim==3
    fprintf('%d%% ci(fmax)  = [%g,%g] kN\n',(probs(2)-probs(1))*100,fmax_ci(1)*1e-3,fmax_ci(2)*1e-3);
end
fprintf('\n');

fprintf('mean(udmax)   = %g mm\n',udmax_mean*1e3);
fprintf('std(udmax)    = %g mm\n',udmax_std*1e3);
fprintf('disp(udmax)   = %g\n',udmax_std/udmax_mean);
fprintf('%d%% ci(udmax) = [%g,%g] mm\n',(probs(2)-probs(1))*100,udmax_ci(1)*1e3,udmax_ci(2)*1e3);

%% Display
if displayModel
    [t,rep] = gettevol(T);
    
    %% Display domains, boundary conditions and meshes
    plotDomain({D,C},'legend',false);
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
    figure('Name','Force-displacement')
    clf
    ciplot(ft_ci(1,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),ft_ci(2,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),t*1e3,'b');
    alpha(0.2)
    hold on
    plot(t*1e3,ft_mean*((Dim==2)*1e-6+(Dim==3)*1e-3),'-b','Linewidth',linewidth)
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
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_fmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_fmax.tex');
    
    %% Display pdf of critical displacement
    figure('Name','Probability Density Estimate: Critical displacement')
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
    ylabel('$p_{U_c}(u)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_udmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_udmax.tex');
    
    %% Display means, variances and samples of solutions at different instants
    ampl = 0;
    switch lower(symmetry)
        case {'isotropic','meanisotropic'} % almost surely or mean isotropic material
            switch lower(loading)
                case 'tension'
                    rep = find(abs(t-5.5e-6)<eps | abs(t-5.75e-5)<eps | abs(t-6e-6)<eps | abs(t-6.25e-6)<eps);
                case 'shear'
                    rep = find(abs(t-1e-5)<eps | abs(t-1.25e-5)<eps | abs(t-1.35e-5)<eps | abs(t-1.5e-5)<eps);
                otherwise
                    error('Wrong loading case')
            end
        case 'anisotropic' % anisotropic material
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
            plotSolution(S,uj_var,'displ_var',i,'ampl',ampl);
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
    uk = TIMEMATRIX(reshape(ut_mean(:,:),sz_u),T);
    uk_var = TIMEMATRIX(reshape(ut_var(:,:),sz_u),T);
    
    evolSolution(S_phase,dk,'FrameRate',framerate,'filename','damage_mean','pathname',pathname,options{:});
    evolSolution(S_phase,dk_var,'FrameRate',framerate,'filename','damage_var','pathname',pathname,options{:});
    for i=1:Dim
        evolSolution(S,uk,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i) '_mean'],'pathname',pathname,options{:});
        evolSolution(S,uk_var,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i) '_var'],'pathname',pathname,options{:});
    end
    
    % for i=1:(Dim*(Dim+1)/2)
    %     evolSolution(S,uk,'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['epsilon_mean_' num2str(i)],'pathname',pathname,options{:});
    %     evolSolution(S,uk,'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['sigma_mean_' num2str(i)],'pathname',pathname,options{:});
    % end
    %
    % evolSolution(S,uk,'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename','epsilon_von_mises_mean','pathname',pathname,options{:});
    % evolSolution(S,uk,'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename','sigma_von_mises_mean','pathname',pathname,options{:});
    % evolSolution(S,uk,'energyint','','ampl',ampl,'FrameRate',framerate,'filename','internal_energy_mean','pathname',pathname,options{:});

    for k=1:size(dt_sample,1)
        dk = TIMEMATRIX(reshape(dt_sample(k,:,:),sz_d),T);
        uk = TIMEMATRIX(reshape(ut_sample(k,:,:),sz_u),T);
        
        evolSolution(S_phase,dk,'FrameRate',framerate,'filename',['damage_sample_' num2str(k)],'pathname',pathname,options{:});
        for i=1:Dim
            evolSolution(S,uk,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i) '_sample_' num2str(k)],'pathname',pathname,options{:});
        end
        
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
