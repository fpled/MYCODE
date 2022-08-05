%% Phase field fracture model - stochastic linear elasticity problem   %%
%  Asymmetric notched plate with three holes under three-point bending %%
%%---------------------------------------------------------------------%%
% [Ingraffea, Grigoriu, 1990] (experimental tests)
% [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM] (SIF-based method with local remeshing and special FE)
% [Ventura, Xu, Belytschko, 2002, IJNME] (vector level set method with discontinuous enrichment in meshless method)
% [Guidault, Allix, Champaney, Cornuault, 2008, CMAME] (MsXFEM)
% [Miehe, Welschinger, Hofacker, 2010, IJNME] (anisotropic phase field model of Miehe et al.)
% [Miehe, Hofacker, Welschinger, 2010, CMAME] (anisotropic phase field model of Miehe et al.)
% [Häusler, Lindhorst, Horst, 2011, IJNME] (XFEM)
% [Geniaut, Galenne, 2012, IJSS] (XFEM)
% [Passieux, Rethore, Gravouil, Baietto, 2013, CM] (XFEM)
% [Ambati, Gerasimov, De Lorenzis, 2015, CM] (hybrid isotropic-anisotropic phase field model of Ambati et al. compared with the anisotropic one of Miehe et al.)
% [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME] (isotropic phase field model with no split of Bourdin et al. compared to experimental data of [Winkler PhD thesis, 2001])
% [Wu, Nguyen, 2018, JMPS] (hybrid isotropic-anisotropic phase field model of Wu et al.)
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM] (anisotropic phase field model of Wu et al.)

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

test = true; % coarse mesh
% test = false; % fine mesh

numWorkers = 4;
% numWorkers = 1; maxNumCompThreads(1); % mono-thread computation

% Deterministic model parameters
Dim = 2; % space dimension Dim = 2
setup = 2; % notch geometry setup = 1, 2, 3, 4, 5
PFmodel = 'AnisotropicMiehe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicSpectral', 'AnisotropicHe'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 100; % maximum number of iterations at each loading increment
initialCrack = 'GeometricNotch'; % 'GeometricCrack', 'GeometricNotch', 'InitialPhaseField'

% Random model parameters
% N = 500; % number of samples
N = numWorkers;
randMat = struct('delta',0.1,'lcorr',1e-4); % random material parameters model
randPF = struct('aGc',0,'bGc',0,'lcorr',Inf); % random phase field parameters model

filename = ['phasefieldStoLinElasAsymmetricNotchedPlateSetup' num2str(setup) PFmodel PFsplit PFregularization PFsolver initialCrack 'MaxIter' num2str(maxIter) 'Adaptive'];
filename = [filename '_' num2str(N) 'samples'];
if any(randMat.delta)
    filename = [filename '_RandMat_Delta' num2str(randMat.delta,'_%g') '_Lcorr' num2str(randMat.lcorr,'_%g')];
end
if any(randPF.aGc) && any(randPF.bGc)
    gcbounds = [randPF.aGc(:),randPF.bGc(:)]';
    filename = [filename '_RandPF_Gc' num2str(gcbounds(:)','_%g') '_Lcorr' num2str(randPF.lcorr,'_%g')];
end

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
formats = {'epsc'};
renderer = 'OpenGL';

gmshoptions = '-v 0';
mmgoptions = '-nomove -hausd 0.01 -hgrad 1.1 -v -1';
% gmshoptions = '-v 5';
% mmgoptions = '-nomove -hausd 0.01 -hgrad 1.3 -v 1';

%% Problem
if setProblem
    %% Domains and meshes
    unit = 1e-3; % for mm % [Guidault, Allix, Champaney, Cornuault, 2008, CMAME], [Miehe, Welschinger, Hofacker, 2010, IJNME],
    % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Passieux, Rethore, Gravouil, Baietto, 2013, CM]
    % unit = 25.4e-3; % for inch % [Ingraffea, Grigoriu, 1990], [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM],
    % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    switch setup
        case 1 % [Ingraffea, Grigoriu, 1990], [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM], [Ventura, Xu, Belytschko, 2002, IJNME],
            % [Guidault, Allix, Champaney, Cornuault, 2008, CMAME], [Häusler, Lindhorst, Horst, 2011, IJNME], [Geniaut, Galenne,2012, IJSS],
            % [Passieux, Rethore, Gravouil, Baietto, 2013, CM]
            a = 1.5*unit; % crack length
            b = 5*unit; % crack offset from the centerline
        case 2 % [Ingraffea, Grigoriu, 1990], [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM], [Ventura, Xu, Belytschko, 2002, IJNME],
            % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Häusler, Lindhorst, Horst, 2011, IJNME],
            % [Geniaut, Galenne, 2012, IJSS], [Passieux, Rethore, Gravouil, Baietto, 2013, CM], [Ambati, Gerasimov, De Lorenzis, 2015, CM],
            % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME],
            % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
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
    
    clD = 0.1*unit; % characteristic length for domain
    % cl = 0.01*unit; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    cl = 0.025*unit/2; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME]
    % cl = 0.01*unit/2; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    if test
        clD = 0.2*unit;
        cl = 0.025*unit;
%     else
%         clD = min(min(min(randMat.lcorr),min(randPF.lcorr))/4,clD);
%         cl = min(min(min(randMat.lcorr),min(randPF.lcorr))/4,cl);
    end
    clC = cl; % characteristic length for edge crack/notch
    clH = cl; % characteristic length for circular holes
    switch lower(initialCrack)
        case 'geometriccrack'
            S_phase = gmshasymmetricnotchedplatewithedgecrack(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_domain_asymmetric_notched_plate'),Dim,'duplicate');
        case 'geometricnotch'
            c = 0.025*unit; % crack width
            S_phase = gmshasymmetricnotchedplatewithedgesmearedcrack(a,b,c,clD,clC,clH,unit,fullfile(pathname,'gmsh_domain_asymmetric_notched_plate'));
        case 'initialphasefield'
            S_phase = gmshasymmetricnotchedplatewithedgecrack(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_domain_asymmetric_notched_plate'),Dim,lower(initialCrack));
        otherwise
            error('Wrong model for initial crack');
    end
    
    sizemap = @(d) (clC-clD)*d+clD;
    % sizemap = @(d) clD*clC./((clD-clC)*d+clC);
    
    %% Phase field problem
    %% Material
    % Critical energy release rate (or fracture toughness)
    gc = 1e3;
    % gc = 304.321; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % gc = 315; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % Regularization parameter (width of the smeared crack)
    % l = 0.05*unit; % [Wu, Nguyen, 2018, JMPS]
    l = 0.025*unit; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % l = 0.01*unit; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % Small artificial residual stiffness
    k = 1e-12;
    % k = 0;
    
    % Material
    switch lower(PFregularization)
        case 'at1'
            % c0 = 8/3;
            K = 3/4*gc*l; % K = 2*(gc*l)/c0;
            R = 0;
            Qn = -3/8*gc/l; % Qn = -(gc/l)/c0;
        case 'at2'
            % c0 = 2;
            K = gc*l; % K = 2*(gc*l)/c0;
            R = gc/l; % R = 2*(gc/l)/c0;
            Qn = 0;
        otherwise
            error('Wrong regularization model');
    end
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'PFregularization',PFregularization,'aGc',randPF.aGc,'bGc',randPF.bGc,'lcorr',randPF.lcorr);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    R = 2*unit;
    BU = CIRCLE(0.0,h,R);
    BL = CIRCLE(-ls,-h,R);
    BR = CIRCLE(ls,-h,R);
    H1 = CIRCLE(-lh,h-ph-2*dh,r);
    H2 = CIRCLE(-lh,h-ph-dh,r);
    H3 = CIRCLE(-lh,h-ph,r);
    
    switch lower(initialCrack)
        case 'geometriccrack'
            C = POINT([-b,-h+a]);
            S_phase = final(S_phase,'duplicate');
        case 'geometricnotch'
            C = LIGNE([-b-c/2,-h+a],[-b+c/2,-h+a]);
            S_phase = final(S_phase);
        case 'initialphasefield'
            S_phase = final(S_phase);
        otherwise
            error('Wrong model for initial crack');
    end
    S_phase = addcl(S_phase,C,'T',1);
    S_phase = addcl(S_phase,BU,'T');
    S_phase = addcl(S_phase,BL,'T');
    S_phase = addcl(S_phase,BR,'T');
    S_phase = addcl(S_phase,H1,'T',1);
    S_phase = addcl(S_phase,H2,'T',1);
    S_phase = addcl(S_phase,H3,'T',1);
    
    d = calc_init_dirichlet(S_phase);
    cl = sizemap(d);
    S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_domain_asymmetric_notched_plate'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    S = S_phase;
    
    S_phase = setmaterial(S_phase,mat_phase);
    switch lower(initialCrack)
        case 'geometriccrack'
            S_phase = final(S_phase,'duplicate');
        case 'geometricnotch'
            S_phase = final(S_phase);
        case 'initialphasefield'
            S_phase = final(S_phase);
            S_phase = addcl(S_phase,C,'T',1);
        otherwise
            error('Wrong model for initial crack');
    end
    S_phase = addcl(S_phase,BU,'T');
    S_phase = addcl(S_phase,BL,'T');
    S_phase = addcl(S_phase,BR,'T');
    
    %% Stiffness matrices and sollicitation vectors
    % a_phase = BILINFORM(1,1,K); % uniform values
    % % a_phase = DIFFUSIONFORM(K);
    % a_phase = setfree(a_phase,0);
    % K_phase = calc_matrix(a_phase,S_phase);
    % b_phase = calc_nonhomogeneous_vector(S_phase,K_phase);
    % b_phase = -b_phase;
    % K_phase = freematrix(S_phase,K_phase);
    
    % r_phase = BILINFORM(0,0,R,0); % nodal values
    % R_phase = calc_matrix(r_phase,S_phase);
    % A_phase = K_phase + R_phase;
    
    % l_phase = LINFORM(0,Qn,0); % nodal values
    % l_phase = setfree(l_phase,1);
    % b_phase = b_phase + calc_vector(l_phase,S_phase);
    
    % [A_phase,b_phase] = calc_rigi(S_phase);
    % b_phase = -b_phase + bodyload(S_phase,[],'QN',Qn);
    
    %% Linear elastic displacement field problem
    %% Materials
    % Option
    option = 'DEFO'; % plane strain [Guidault, Allix, Champaney, Cornuault, 2008, CMAME], [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM]
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
    % E = 3.102e9; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % E = 3.275e9; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % NU = 0.35; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % E = 200e9; % [Guidault, Allix, Champaney, Cornuault, 2008, CMAME]
    % NU = 0.3; % [Guidault, Allix, Champaney, Cornuault, 2008, CMAME]
    % Energetic degradation function
    g = @(d) (1-d).^2;
    % Thickness
    DIM3 = 1;
    % DIM3 = 0.5*unit; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % Density
    RHO = 1;
    
    % Material
    d = calc_init_dirichlet(S_phase);
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit,'delta',randMat.delta,'lcorr',randMat.lcorr);
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    PU = POINT([0.0,h]);
    PL = POINT([-ls,-h]);
    PR = POINT([ls,-h]);
    
    switch lower(initialCrack)
        case 'geometriccrack'
            S = final(S,'duplicate');
        otherwise
            S = final(S);
    end
    
    ud = 0;
    S = addcl(S,PU,'UY',ud);
    S = addcl(S,PL,{'UX','UY'});
    S = addcl(S,PR,'UY');
    
    %% Stiffness matrices and sollicitation vectors
    % [A,b] = calc_rigi(S);
    % b = -b;
    
    %% Time scheme
    % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % dt = 1e-4*unit;
    % nt = 2500;
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
    save(fullfile(pathname,'problem.mat'),'unit','T','S_phase','S','sizemap','C','BU','BL','BR','H1','H2','H3','PU','PL','PR');
else
    load(fullfile(pathname,'problem.mat'),'unit','T','S_phase','S','sizemap','C','BU','BL','BR','H1','H2','H3','PU','PL','PR');
end

%% Solution
if solveProblem
    myparallel('start',numWorkers);
    
    %% Solution
    tTotal = tic;
    
    nbSamples = 1;
    fun = @(S_phase,S,filename) solvePFDetLinElasAsymmetricNotchedPlateAdaptive(S_phase,S,T,PFsolver,C,BU,BL,BR,H1,H2,H3,PU,PL,PR,sizemap,...
        'filename',filename,'pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    [ft,dt,ut,St_phase,St] = solvePFStoLinElasAdaptive(S_phase,S,T,fun,N,'filename','gmsh_domain_asymmetric_notched_plate','pathname',pathname,'nbsamples',nbSamples);
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
    
    save(fullfile(pathname,'solution.mat'),'N','ft','dt','ut','St_phase','St',...
        'ft_mean','ft_std','ft_ci','probs','time',...
        'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
        'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw');
else
    load(fullfile(pathname,'solution.mat'),'N','ft','dt','ut','St_phase','St',...
        'ft_mean','ft_std','ft_ci','probs','time',...
        'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
        'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw');
end

%% Outputs
fprintf('\n');
fprintf('setup    = %d\n',setup);
fprintf('PF model = %s\n',PFmodel);
fprintf('PF split = %s\n',PFsplit);
fprintf('PF regularization = %s\n',PFregularization);
fprintf('PF solver = %s\n',PFsolver);
fprintf('nb elements = %g (initial) - %g (final)\n',getnbelem(S),getnbelem(St{end}));
fprintf('nb nodes    = %g (initial) - %g (final)\n',getnbnode(S),getnbnode(St{end}));
fprintf('nb dofs     = %g (initial) - %g (final)\n',getnbddl(S),getnbddl(St{end}));
fprintf('nb time dofs = %g\n',getnbtimedof(T));
fprintf('nb samples = %g\n',N);
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

fprintf('mean(fmax)   = %g kN/mm\n',fmax_mean*1e-6);
fprintf('std(fmax)    = %g kN/mm\n',fmax_std*1e-6);
fprintf('disp(fmax)   = %g\n',fmax_std/fmax_mean);
fprintf('%d%% ci(fmax) = [%g,%g] kN/mm\n',(probs(2)-probs(1))*100,fmax_ci(1)*1e-6,fmax_ci(2)*1e-6);
fprintf('\n');

fprintf('mean(udmax)   = %g mm\n',udmax_mean*1e3);
fprintf('std(udmax)    = %g mm\n',udmax_std*1e3);
fprintf('disp(udmax)   = %g\n',udmax_std/udmax_mean);
fprintf('%d%% ci(udmax) = [%g,%g] mm\n',(probs(2)-probs(1))*100,udmax_ci(1)*1e3,udmax_ci(2)*1e3);

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

%% Display samples of solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display force-displacement curve
    figure('Name','Force-displacement')
    clf
    plot(t*1e3,ft_mean*1e-6,'-b','Linewidth',linewidth)
    hold on
    ciplot(ft_ci(1,:)*1e-6,ft_ci(2,:)*1e-6,t*1e3,'b');
    alpha(0.2)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    l = legend({'mean function',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval']},...
        'Location','NorthWest');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    figure('Name','Force-displacement')
    clf
    color = distinguishable_colors(N);
    for i=1:N
        plot(t*1e3,ft(i,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),'LineStyle','-','Color',color(i,:),'Linewidth',linewidth)
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    mysaveas(pathname,'forces_displacement',formats);
    mymatlab2tikz(pathname,'forces_displacement.tex');
    
    %% Display pdf of maximum force
    figure('Name','Probability Density Estimate: Maximum force')
    clf
    plot(fmax_xi*1e-6,fmax_f,'-b','LineWidth',linewidth)
    hold on
    ind_fmax = find(fmax_xi>=fmax_ci(1) & fmax_xi<fmax_ci(2));
    area(fmax_xi(ind_fmax)*1e-6,fmax_f(ind_fmax),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(fmax_mean*1e-6,0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$f$ [kN]','Interpreter',interpreter)
    ylabel('$p_{F_{\mathrm{max}}}(f)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_fmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_fmax.tex');
    
    %% Display pdf of maximum displacement
    figure('Name','Probability Density Estimate: Maximum displacement')
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
    ylabel('$p_{U_{D,\mathrm{max}}}(u)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_udmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_udmax.tex');
    
    %% Display samples of solutions at different instants
    ampl = 0;
    switch setup
        case {1,4,5}
            rep = find(abs(t-0.190*unit)<eps | abs(t-0.201*unit)<eps | abs(t-0.203*unit)<eps | abs(t-0.204*unit)<eps | abs(t-0.205*unit)<eps | abs(t-0.207*unit)<eps);
        case {2,3}
            rep = find(abs(t-0.210*unit)<eps | abs(t-0.215*unit)<eps | abs(t-0.218*unit)<eps | abs(t-0.219*unit)<eps | abs(t-0.220*unit)<eps | abs(t-0.222*unit)<eps);
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
