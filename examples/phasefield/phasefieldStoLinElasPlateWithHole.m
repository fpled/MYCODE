%% Phase-field fracture model - stochastic linear elasticity problem %%
%  Plate with a central circular hole under compression              %%
%%-------------------------------------------------------------------%%
% [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS] (experimental tests)
% [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF] (anisotropic phase-field model of Miehe et al.)
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME] (anisotropic phase-field model of He et al.)
% [Luo, Chen, Wang, Li, 2022, CM] (anisotropic phase-field model of He et al. with anisotropic fracture surface energy)
% [Vu, Le Quang, He, 2024, AES] (anisotropic phase-field model of He et al. with improved degradation function of Sargado et al.)

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
symmetry = 'Isot'; % 'Isot', 'MeanIsot', 'Anisot'. Material symmetry
ang = 45; % clockwise material orientation angle around z-axis for anisotopic material [deg]
PFmodel = 'Miehe'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 1; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
critConv = 'Energy'; % 'Solution', 'Residual', 'Energy'
FEmesh = 'Optim'; % 'Unif' or 'Optim'

% Random model parameters
% N = 500; % number of samples
N = numWorkers;
L = 15e-3; H = 2*L; lcorr = min(L,H)/10; % lcorr = 15e-4;
randMat = struct('delta',0.2,'lcorr',lcorr); % random material parameters model
randPF = struct('aGc',0,'bGc',0,'lcorr',Inf); % random phase-field parameters model

suffix = '';

foldername = ['plateWithHole_' num2str(Dim) 'D'];
filename = ['linElas' symmetry];
if strcmpi(symmetry,'anisot') % anisotropic material
    filename = [filename num2str(ang) 'deg'];
end
filename = [filename PFmodel PFsplit PFregularization PFsolver...
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
    % 2D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Vu, Le Quang, He, 2024, AES]
    % L = 20e-3; % length
    % H = 30e-3; % height
    % r = 4e-3; % radius of the hole
    
    % 2D and 3D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], 2D [Vu, Le Quang, He, 2024, AES]
    % L = 100e-3; % length
    % H = 65e-3; % height
    % e = 40e-3; % thickness [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
    % % e = 20e-3; % thickness [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS]
    % % r = 1.5e-3; % radius of the hole [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS]
    % % r = 2e-3; % radius of the hole [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS]
    % r = 2.5e-3; % radius of the hole [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS], [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
    % % r = 3e-3; % radius of the hole [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS]
    
    % 2D [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM]
    L = 15e-3; % length
    H = 2*L; % height
    r = 3e-3; % radius of the hole
    
    if Dim==2
        e = 1; % thickness
        D = DOMAIN(2,[0.0,0.0],[L,H]);
        C = CIRCLE(L/2,H/2,r);
    elseif Dim==3
        e = 1e-3; % thickness
        D = DOMAIN(3,[0.0,0.0,0.0],[L,H,e]);
        C = CIRCLE(L/2,H/2,0.0,r);
        C = CYLINDER(C,e); % CYLINDER(L/2,H/2,0.0,r,e);
    end
    
    switch lower(FEmesh)
        case 'unif' % uniform mesh
            % 2D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
            % cl = 2e-3;
            % cl = 1e-3;
            % cl = 1/15*1e-3;
            % cl = 0.05e-3;
            % cl = 1/30*1e-3;
            % cl = 0.025e-3;
            % cl = 1/60*1e-3;
            % cl = 0.0125e-3;
            % cl = 0.01e-3;
            
            % 2D [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            cl = 0.06e-3;
            if test
                % cl = 2e-3; % 2D and 3D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
                cl = 0.24e-3; % 2D [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            end
            clD = cl; % characteristic length for domain
            clC = cl; % characteristic length for circular hole
            B = [];
        case 'optim' % optimized mesh
            % 2D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], 2D [Vu, Le Quang, He, 2024, AES]
            % clD = 1e-3; % characteristic length for domain
            % clC = 0.01e-3; % characteristic length for circular hole
            
            % 2D and 3D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], 2D [Vu, Le Quang, He, 2024, AES]
            % clD = 0.25e-3; % characteristic length for domain
            % clC = 0.05e-3; % characteristic length for circular hole
            % if test
            %     clD = 10e-3; % 2D and 3D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
            %     clC = 1e-3; % 2D and 3D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
            % end
            
            % 2D [Luo, Chen, Wang, Li, 2022, CM]
            % clD = 0.12e-3; % characteristic length for domain
            % clC = 0.024e-3; % characteristic length for circular hole
            
            % 2D [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            clD = 0.24e-3; % characteristic length for domain
            clC = 0.06e-3; % characteristic length for circular hole
            if test
                if Dim==2
                    clD = 0.48e-3;
                    clC = 0.24e-3;
                elseif Dim==3
                    clD = 0.96e-3;
                    clC = 0.36e-3;
                end
            end
            VIn = clC; VOut = clD;
            switch lower(PFmodel)
                case {'bourdin','amor','heamor'}
                    XMin = 0; XMax = L;
                    YMin = H/2-3/2*r; YMax = H/2+3/2*r;
                    Thickness = H/2-3/2*r;
                otherwise
                    XMin = L/2-3/2*r; XMax = L/2+3/2*r;
                    YMin = 0; YMax = H;
                    Thickness = L/2-3/2*r;
            end
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
    if Dim==2
        S_phase = gmshDomainWithHole(D,C,clD,clC,fullfile(pathname,'gmsh_plate_with_hole'),Dim,'Box',B);
    elseif Dim==3
        S_phase = gmshDomainWithHole(D,C,clD,clC,fullfile(pathname,'gmsh_plate_with_hole'),Dim,'Box',B,'extrude');
    end
    S = S_phase;
    
    %% Phase-field problem
    %% Material
    switch lower(symmetry)
        case {'isot','meanisot'} % almost surely or mean isotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 1.4; % [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS], [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM], [Vu, Le Quang, He, 2024, AES]
            % Regularization parameter (width of the smeared crack)
            % l = 0.5e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
            % l = 0.45e-3; % [Vu, Le Quang, He, 2024, AES]
            % l = 0.3e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Vu, Le Quang, He, 2024, AES]
            % l = 0.2e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Vu, Le Quang, He, 2024, AES]
            l = 0.12e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM]
            % l = 0.1e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Vu, Le Quang, He, 2024, AES]
            % l = 0.06e-3; % [Vu, Le Quang, He, 2024, AES]
            % l = 0.05e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
            % l = 0.025e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
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
    [K,R,Qn] = setphasefieldparam(gc,l,PFregularization);
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'DIM3',e,'PFregularization',PFregularization,'aGc',randPF.aGc,'bGc',randPF.bGc,'lcorr',randPF.lcorr);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BRight = LINE([L,0.0],[L,H]);
        BLeft = LINE([0.0,0.0],[0.0,H]);
    elseif Dim==3
        BRight = PLANE([L,0.0,0.0],[L,H,0.0],[L,0.0,e]);
        BLeft = PLANE([0.0,0.0,0.0],[0.0,H,0.0],[0.0,0.0,e]);
    end
    
    findddlboundary = @(S_phase) union(findddl(S_phase,'T',BRight),findddl(S_phase,'T',BLeft));
    
    S_phase = final(S_phase);
    
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
    option = 'DEFO'; % plane strain [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS], [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM], [Vu, Le Quang, He, 2024, AES]
    % option = 'CONT'; % plane stress [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
    switch lower(symmetry)
        case {'isot','meanisot'} % almost surely or mean isotropic material
            % Young modulus and Poisson ratio
            E = 12e9; % [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS], [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM], [Vu, Le Quang, He, 2024, AES]
            NU = 0.3; % [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS], [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM], [Vu, Le Quang, He, 2024, AES]
            if strcmpi(symmetry,'meanisot')
                mu = E/(1+nu)/2; % second Lame coefficient (shear modulus)
                % Elasticity matrix
                if Dim==2
                    switch lower(option)
                        case 'defo'
                            % Cmat = e*E/(1+nu)/(1-2*nu)*...
                            %     [(1-nu),nu,0;...
                            %     nu,(1-nu),0;...
                            %     0,0,(1-2*nu)/2];
                            lambda = E*nu/(1+nu)/(1-2*nu); % first Lame coefficient
                        otherwise
                            % Cmat = e*E/(1-nu^2)*...
                            %     [1,nu,0;...
                            %     nu,1,0;...
                            %     0,0,(1-nu)/2];
                            lambda = E*nu/(1-nu^2); % first Lame coefficient
                    end
                    Cmat = e*...
                        [lambda+2*mu,lambda,0;...
                        lambda,lambda+2*mu,0;...
                        0,0,mu];
                elseif Dim==3
                    % Cmat = E/(1+nu)/(1-2*nu)*...
                    %     [(1-nu),nu,nu,0,0,0;...
                    %     nu,(1-nu),nu,0,0,0;...
                    %     nu,nu,(1-nu),0,0,0;...
                    %     0,0,0,(1-2*nu)/2,0,0;...
                    %     0,0,0,0,(1-2*nu)/2,0;...
                    %     0,0,0,0,0,(1-2*nu)/2];
                    lambda = E*nu/(1+nu)/(1-2*nu); % first Lame coefficient
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
        BU = LINE([0.0,H],[L,H]);
        BL = LINE([0.0,0.0],[L,0.0]);
    elseif Dim==3
        BU = PLANE([0.0,H,0.0],[L,H,0.0],[0.0,H,e]);
        BL = PLANE([0.0,0.0,0.0],[L,0.0,0.0],[0.0,0.0,e]);
    end
    PD = getvertices(D);
    P1 = POINT(PD{1});
    P2 = POINT(PD{2});
    
    addbc = @(S,ud) addbcPlateWithHole(S,ud,BU,BL,P1,P2);
    findddlforce = @(S) findddl(S,'UY',BU);
    
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
    % 2D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], 2D [Vu, Le Quang, He, 2024, AES]
    % % du = 5e-3 mm during 11 time steps (up to u = 0.055 mm)
    % % du = 1.5e-3 mm during 37 time steps (up to u = 0.0555 mm)
    % % du = 1e-3 mm during 55 time steps (up to u = 0.055 mm)
    % % du = 2e-4 mm during 275 time steps (up to u = 0.055 mm)
    % % du = 3e-4 mm during 183 time steps (up to u = 0.0549 mm)
    % du = 1e-4 mm during 550 time steps (up to u = 0.055 mm)
    % % du = 5e-5 mm during 1100 time steps (up to u = 0.055 mm)
    % % du = 3e-5 mm during 1833 time steps (up to u = 0.05499 mm)
    % % du = 2e-5 mm during 2750 time steps (up to u = 0.055 mm)
    % dt = 1e-7;
    % nt = 550;
    % t = linspace(dt,nt*dt,nt);
    % T = TIMEMODEL(t);
    
    % 2D and 3D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
    % du = 1e-3 mm during the first stage (until the phase-field reaches the threshold value)
    % du = 1e-4 mm during the last stage (as soon as the phase-field exceeds the threshold value, up to u = 0.2 mm)
    % dt0 = 1e-6;
    % dt1 = 1e-7;
    % if test
    %     dt0 = 2e-6;
    %     dt1 = 2e-7;
    % end
    % tf = 2e-4;
    % % dth = 0.9;
    % dth = 0.6;
    
    % 2D [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
    % du = 8e-5 mm during the first stage (until the phase-field reaches the threshold value)
    % du = 2e-5 mm during the last stage (as soon as the phase-field exceeds the threshold value, up to u = 25e-3 mm)
    % dt0 = 8e-8;
    % dt1 = 2e-8;
    % if test
    %     dt0 = 16e-8;
    %     dt1 = 4e-8;
    % end
    % tf = 25e-6;
    % dth = 0.6;
    % T = struct('dt0',dt0,'dt1',dt1,'tf',tf,'dth',dth);
    
    % du = 2e-5 mm during 1250 time steps (u = 25e-3 mm)
    dt = 2e-8;
    nt = 1250;
    if test
        dt = 5e-8;
        nt = 500;
    end
    t = linspace(dt,nt*dt,nt);
    T = TIMEMODEL(t);
    
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
    
    displayIter = false;
    displaySol  = false;
    
    nbSamples = 1;
    fun = @(S_phase,S) solvePFDetLinElas(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,...
        'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
        'displayiter',displayIter,'displaysol',displaySol);
    % fun = @(S_phase,S) solvePFDetLinElasPlateWithHole(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P1,P2,...
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
if solveProblem
    filenameResults = fullfile(pathname,'results.txt');
    fid = fopen(filenameResults,'w');
    fprintf(fid,'Plate with hole\n');
    fprintf(fid,'\n');
    fprintf(fid,'dim      = %d\n',Dim);
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
    type(filenameResults) % fprintf('%s', fileread(filenameResults))
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

%% Display statistics of solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display force-displacement curve
    figure('Name','Force vs displacement')
    clf
    plot(t*1e3,ft_mean*((Dim==2)*1e-6+(Dim==3)*1e-3),'-b','LineWidth',linewidth)
    hold on
    ciplot(ft_ci(1,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),ft_ci(2,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),t*1e3,'b');
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
        plot(t*1e3,ft(i,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
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
    
    %% Display maximum damage-displacement curve
    figure('Name','Maximum damage vs displacement')
    clf
    plot(t*1e3,dmaxt_mean,'-b','LineWidth',linewidth)
    hold on
    ciplot(dmaxt_ci(1,:),dmaxt_ci(2,:),t*1e3,'b');
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
        plot(t*1e3,dmaxt(i,:),'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
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
    
    %% Display energy-displacement curve
    Et_ci = quantile(Eut+Edt,probs);
    figure('Name','Energies vs displacement')
    clf
    plot(t*1e3,Eut_mean,'-b','LineWidth',linewidth)
    hold on
    plot(t*1e3,Edt_mean,'-r','LineWidth',linewidth)
    plot(t*1e3,Eut_mean+Edt_mean,'-k','LineWidth',linewidth)
    ciplot(Eut_ci(1,:),Eut_ci(2,:),t*1e3,'b');
    ciplot(Edt_ci(1,:),Edt_ci(2,:),t*1e3,'r');
    ciplot(Et_ci(1,:),Et_ci(2,:),t*1e3,'k');
    alpha(0.2)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Energy [J]','Interpreter',interpreter)
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
        plot(t*1e3,Eut(i,:),'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]'...,'Interpreter',interpreter...
        )
    ylabel('Elastic strain energy [J]'...,'Interpreter',interpreter...
        )
    mysaveas(pathname,'energies_elastic_displacement',formats);
    mymatlab2tikz(pathname,'energies_elastic_displacement.tex');
    
    % Fracture energy-displacement curves
    figure('Name','Fracture energies vs displacement')
    clf
    for i=1:N
        plot(t*1e3,Edt(i,:),'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]'...,'Interpreter',interpreter...
        )
    ylabel('Fracture energy [J]'...,'Interpreter',interpreter...
        )
    mysaveas(pathname,'energies_fracture_displacement',formats);
    mymatlab2tikz(pathname,'energies_fracture_displacement.tex');
    
    % Total energy-displacement curves
    figure('Name','Total energies vs displacement')
    clf
    for i=1:N
        plot(t*1e3,Eut(i,:)+Edt(i,:),'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]'...,'Interpreter',interpreter...
        )
    ylabel('Total energy [J]'...,'Interpreter',interpreter...
        )
    mysaveas(pathname,'energies_total_displacement',formats);
    mymatlab2tikz(pathname,'energies_total_displacement.tex');
    
    %% Display pdf of maximum force
    figure('Name','Probability Density Estimate: Maximum force')
    clf
    plot(fmax_xi*((Dim==2)*1e-6+(Dim==3)*1e-3),fmax_f*((Dim==2)*1e6+(Dim==3)*1e3),'-r','LineWidth',linewidth)
    hold on
    ind_fmax = find(fmax_xi>=fmax_ci(1) & fmax_xi<fmax_ci(2));
    area(fmax_xi(ind_fmax)*((Dim==2)*1e-6+(Dim==3)*1e-3),fmax_f(ind_fmax)*((Dim==2)*1e6+(Dim==3)*1e3),'FaceColor','r','EdgeColor','none','FaceAlpha',0.2)
    scatter(fmax_mean*((Dim==2)*1e-6+(Dim==3)*1e-3),0,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r')
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
    plot(fc_xi*((Dim==2)*1e-6+(Dim==3)*1e-3),fc_f*((Dim==2)*1e6+(Dim==3)*1e3),'-b','LineWidth',linewidth)
    hold on
    ind_fc = find(fc_xi>=fc_ci(1) & fc_xi<fc_ci(2));
    area(fc_xi(ind_fc)*((Dim==2)*1e-6+(Dim==3)*1e-3),fc_f(ind_fc)*((Dim==2)*1e6+(Dim==3)*1e3),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(fc_mean*((Dim==2)*1e-6+(Dim==3)*1e-3),0,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b')
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
    % 2D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
    % tSnapshots = [0.04 0.06 0.08 0.09]*1e-3;
    
    % 2D [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
    % tSnapshots = [18.5 24.6]*1e-6;
    
    % 2D [Luo, Chen, Wang, Li, 2022, CM]
    % tSnapshots = [17.56 19.00 22.00 25.00]*1e-6;
    
    tSnapshots = [18 18.5 19 20 21 22]*1e-6;
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
