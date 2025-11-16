%% Phase-field fracture model - deterministic linear elasticity problem %%
%  Plate with a central circular hole under compression                 %%
%%----------------------------------------------------------------------%%
% [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS] (experimental tests)
% [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF] (anisotropic phase-field model of Miehe et al.)
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME] (anisotropic phase-field model of He et al.)
% [Luo, Chen, Wang, Li, 2022, CM] (anisotropic phase-field model of He et al. with anisotropic fracture surface energy)

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
symmetry = 'Isot'; % 'Isot' or 'Anisot'. Material symmetry
ang = 45; % clockwise material orientation angle around z-axis for anisotopic material [deg]
PFmodel = 'Miehe'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 1; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
critConv = 'Energy'; % 'Solution', 'Residual', 'Energy'
FEmesh = 'Optim'; % 'Unif' or 'Optim'
heff = 1; % self-healing efficiency
dact = 0.5; % damage activation threshold
ratiohcgc = 1; % ratio Hc/Gc
healing = (heff~=0);

% PFmodels = {'Bourdin','Amor','Miehe','HeAmor','HeFreddi','Zhang'};
% PFsplits = {'Strain','Stress'};
% PFregularizations = {'AT1','AT2'};
% PFsolvers = {'HistoryFieldElem','BoundConstrainedOptim'};
% maxIters = [1,Inf];

% for iPFmodel=1:length(PFmodels)
% PFmodel = PFmodels{iPFmodel};
% for iPFsplit=1:length(PFsplits)
% PFsplit = PFsplits{iPFsplit};
% for iPFRegularization=1:length(PFregularizations)
% PFregularization = PFregularizations{iPFRegularization};
% for iPFsolver=1:length(PFsolvers)
% PFsolver = PFsolvers{iPFsolver};
% for imaxIter=1:length(maxIters)
% maxIter = maxIters(imaxIter);
% close all

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
filename = [filename 'Mesh' FEmesh];
if healing
    filename = [filename 'Heff' num2str(heff) 'Dact' num2str(dact) 'RatioHcGc' num2str(ratiohcgc)];
end
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
    % 2D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
    % L = 20e-3; % length
    % H = 30e-3; % height
    % r = 4e-3; % radius of the hole
    
    % 2D and 3D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
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
            % 2D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
            % clD = 1e-3; % characteristic length for domain
            % clC = 0.01e-3; % characteristic length for circular hole
            
            % 2D and 3D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
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
    S_healing = S_phase;
    
    %% Phase-field problem
    %% Material
    switch lower(symmetry)
        case 'isot' % isotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 1.4; % [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS], [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM]
            % Regularization parameter (width of the smeared crack)
            % l = 0.5e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
            % l = 0.3e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
            % l = 0.2e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
            l = 0.12e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM]
            % l = 0.1e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
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
    if healing
        % Healing toughness
        hc = ratiohcgc*gc;
        % Regularization parameter (width of the healed zone)
        lh = l;
    end
    % Small artificial residual stiffness
    % k = 1e-12;
    k = 0;
    
    % Material
    [K,R,Qn] = setphasefieldparam(gc,l,PFregularization);
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'DIM3',e,'PFregularization',PFregularization);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    if healing
        [Kh,Rh,Qnh] = sethealingfieldparam(hc,lh,PFregularization);
        mat_healing = FOUR_ISOT('k',Kh,'r',Rh,'qn',Qnh,'DIM3',e,'PFregularization',PFregularization);
        mat_healing = setnumber(mat_healing,1);
        S_healing = setmaterial(S_healing,mat_healing);
    end
    
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
    if healing
        S_healing = final(S_healing);
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
    option = 'DEFO'; % plane strain [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS], [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM]
    % option = 'CONT'; % plane stress [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
    switch lower(symmetry)
        case 'isot' % isotropic material
            % Lame coefficients
            % Young modulus and Poisson ratio
            E = 12e9; % [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS], [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM]
            NU = 0.3; % [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS], [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM]
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
    if healing
        % deff = @(d,h) max(min(d-heff*h,1),0);
        deff = @(d,h) d-heff*h;
        fundact = @(d,h) deff(d,h)>dact;
        g = @(d,h) (1-deff(d,h)).^2;
    else
        g = @(d) (1-d).^2;
    end
    % Density
    RHO = 1;
    
    % Material
    d = calc_init_dirichlet(S_phase);
    if healing
        h = calc_init_dirichlet(S_healing);
    end
    switch lower(symmetry)
        case 'isot' % isotropic material
            if healing
                mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'h',h,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
            else
                mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
            end
        case 'anisot' % anisotropic material
            if healing
                mat = ELAS_ANISOT('C',Cmat,'RHO',RHO,'DIM3',e,'d',d,'h',h,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
            else
                mat = ELAS_ANISOT('C',Cmat,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
            end
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
    % 2D [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
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
    dt0 = 8e-8;
    dt1 = 2e-8;
    if test
        dt0 = 16e-8;
        dt1 = 4e-8;
    end
    tf = 25e-6;
    dth = 0.6;
    T = struct('dt0',dt0,'dt1',dt1,'tf',tf,'dth',dth);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','addbc','findddlforce','findddlboundary');
    if healing
        save(fullfile(pathname,'problem.mat'),'S_healing','deff','heff','dact','fundact','-append');
    end
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','addbc','findddlforce','findddlboundary');
    if healing
        load(fullfile(pathname,'problem.mat'),'S_healing','deff','heff','dact','fundact');
    end
end

%% Solution
if solveProblem
    tTotal = tic;
    
    displayIter = true;
    displaySol  = false;
    
    if healing
        switch lower(PFsolver)
            case {'historyfieldelem','historyfieldnode'}
                [dt,ht,ut,ft,Ht,Edt,Eht,Eut,output] = solvePFSHDetLinElasThreshold(S_phase,S_healing,S,T,PFsolver,addbc,findddlforce,findddlboundary,...
                    'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
                    'heff',heff,'deff',deff,'dact',dact,'fundact',fundact,...
                    'displayiter',displayIter,'displaysol',displaySol);
            otherwise
                [dt,ht,ut,ft,~,Edt,Eht,Eut,output] = solvePFSHDetLinElasThreshold(S_phase,S_healing,S,T,PFsolver,addbc,findddlforce,findddlboundary,...
                    'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
                    'heff',heff,'deff',deff,'dact',dact,'fundact',fundact,...
                    'displayiter',displayIter,'displaysol',displaySol);
        end
    else
        switch lower(PFsolver)
            case {'historyfieldelem','historyfieldnode'}
                [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElasThreshold(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,...
                    'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
                    'displayiter',displayIter,'displaysol',displaySol);
            otherwise
                [dt,ut,ft,~,Edt,Eut,output] = solvePFDetLinElasThreshold(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,...
                    'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
                    'displayiter',displayIter,'displaysol',displaySol);
        end
        % switch lower(PFsolver)
        %     case {'historyfieldelem','historyfieldnode'}
        %         [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElasPlateWithHoleThreshold(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P1,P2,...
        %             'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
        %             'displayiter',displayIter,'displaysol',displaySol);
        %     otherwise
        %         [dt,ut,ft,~,Edt,Eut,output] = solvePFDetLinElasPlateWithHoleThreshold(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P1,P2,...
        %             'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
        %             'displayiter',displayIter,'displaysol',displaySol);
        % end
    end
    
    T = gettimemodel(dt);
    t = gettevol(T);
    dt_val = getvalue(dt);
    dmaxt = max(dt_val);
    idc = find(dmaxt>=min(0.75,max(dmaxt)),1);
    if healing
        ht_val = getvalue(ht);
        Dt_val = deff(dt_val,ht_val);
        hmaxt = max(ht_val);
        Dmaxt = max(Dt_val);
        idc = find(Dmaxt>=min(0.75,max(Dmaxt)),1);
    end
    fc = ft(idc);
    udc = t(idc);
    [fmax,idmax] = max(ft,[],2);
    udmax = t(idmax);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'dt','ut','ft','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','T','time');
    if healing
        save(fullfile(pathname,'solution.mat'),'ht','Eht','hmaxt','Dmaxt','-append');
    end
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        save(fullfile(pathname,'solution.mat'),'Ht','-append');
    end
else
    load(fullfile(pathname,'solution.mat'),'dt','ut','ft','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','T','time');
    if healing
        load(fullfile(pathname,'solution.mat'),'ht','Eht','hmaxt','Dmaxt');
    end
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        load(fullfile(pathname,'solution.mat'),'Ht');
    end
end

%% Outputs
if solveProblem
    fid = fopen(fullfile(pathname,'results.txt'),'w');
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
    plot(t*1e3,ft*((Dim==2)*1e-6+(Dim==3)*1e-3),'-b','LineWidth',linewidth)
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
    
    if healing
        %% Display maximum healing-displacement curve
        figure('Name','Maximum healing vs displacement')
        clf
        plot(t*1e3,hmaxt,'-b','LineWidth',linewidth)
        grid on
        box on
        set(gca,'FontSize',fontsize)
        xlabel('Displacement [mm]','Interpreter',interpreter)
        ylabel('Maximum healing','Interpreter',interpreter)
        mysaveas(pathname,'max_healing_displacement',formats);
        mymatlab2tikz(pathname,'max_healing_displacement.tex');
        
         %% Display maximum effective damage-displacement curve
        figure('Name','Maximum effective damage vs displacement')
        clf
        plot(t*1e3,Dmaxt,'-b','LineWidth',linewidth)
        grid on
        box on
        set(gca,'FontSize',fontsize)
        xlabel('Displacement [mm]','Interpreter',interpreter)
        ylabel('Maximum effective damage','Interpreter',interpreter)
        mysaveas(pathname,'max_damage_eff_displacement',formats);
        mymatlab2tikz(pathname,'max_damage_eff_displacement.tex');
    end
    
    %% Display energy-displacement curves
    figure('Name','Energies vs displacement')
    clf
    plot(t*1e3,Eut,'-b','LineWidth',linewidth)
    hold on
    plot(t*1e3,Edt,'-r','LineWidth',linewidth)
    if healing
        plot(t*1e3,Eht,'-g','LineWidth',linewidth)
        plot(t*1e3,Eut+Edt+Eht,'-k','LineWidth',linewidth)
    else
        plot(t*1e3,Eut+Edt,'-k','LineWidth',linewidth)
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Energy [J]','Interpreter',interpreter)
    if healing
        leg = {'$\Psi_u$','$\Psi_c$','$\Psi_h$','$\Psi_{\mathrm{tot}}$'};
    else
        leg = {'$\Psi_u$','$\Psi_c$','$\Psi_{\mathrm{tot}}$'};
    end
    legend(leg{:},'Location','NorthWest','Interpreter',interpreter)
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
        dj = getmatrixatstep(dt,rep(j));
        if healing
            hj = getmatrixatstep(ht,rep(j));
            Dj = deff(dj,hj);
        end
        uj = getmatrixatstep(ut,rep(j));
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            Hj = getmatrixatstep(Ht,rep(j));
        end
        
        plotSolution(S_phase,dj);
        mysaveas(pathname,['damage_t' num2str(rep(j))],formats,renderer);
        
        if healing
            plotSolution(S_healing,hj);
            mysaveas(pathname,['healing_t' num2str(rep(j))],formats,renderer);
            plotSolution(S_phase,Dj);
            mysaveas(pathname,['damage_eff_t' num2str(rep(j))],formats,renderer);
        end
        
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
    if healing
        evolSolution(S_healing,ht,'FrameRate',framerate,'filename','healing','pathname',pathname,options{:});
        Dt = deff(dt,ht);
        evolSolution(S_phase,Dt,'FrameRate',framerate,'filename','damage_eff','pathname',pathname,options{:});
    end
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
        if healing
            hi = getmatrixatstep(ht,rep(i));
            Di = deff(di,hi);
        end
        ui = getmatrixatstep(ut,rep(i));
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            Hi = getmatrixatstep(Ht,rep(i));
        end
        
        if healing
            switch lower(PFsolver)
                case 'historyfieldelem'
                    write_vtk_mesh(S,{di,hi,Di,ui},{Hi},...
                        {'damage','healing','effective damage','displacement'},{'internal energy density history'},...
                        pathname,'solution',1,i-1);
                case 'historyfieldnode'
                    write_vtk_mesh(S,{di,hi,Di,ui,Hi},[],...
                        {'damage','healing','effective damage','displacement','internal energy density history'},[],...
                        pathname,'solution',1,i-1);
                otherwise
                    write_vtk_mesh(S,{di,hi,Di,ui},[],...
                        {'damage','healing','effective damage','displacement'},[],...
                        pathname,'solution',1,i-1);
            end
        else
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
    end
    make_pvd_file(pathname,'solution',1,length(T));
end

% myparallel('stop');

% end
% end
% end
% end
% end
