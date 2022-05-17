%% Phase field fracture model - deterministic linear elasticity problem %%
%  Asymmetric notched plate with three holes under three-point bending  %%
%%----------------------------------------------------------------------%%
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
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displayModel = false;
displaySolution = false;
makeMovie = false;
saveParaview = true;

test = true; % coarse mesh
% test = false; % fine mesh

Dim = 2; % space dimension Dim = 2
setup = 2; % notch geometry setup = 1, 2, 3, 4, 5
PFmodel = 'AnisotropicMiehe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
initialCrack = 'GeometricCrack'; % 'GeometricCrack', 'GeometricNotch', 'InitialPhaseField'

filename = ['phasefieldDetLinElasAsymmetricNotchedPlateSetup' num2str(setup) PFmodel PFsolver initialCrack 'Adaptive'];

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
    % Internal energy
    H = 0;
    
    % Material
    mat_phase = FOUR_ISOT('k',gc*l,'r',gc/l+2*H);
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
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel);
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
    tTotal = tic;
    
    switch lower(PFsolver)
        case {'historyfieldelem','historyfieldnode'}
            [dt,ut,ft,St_phase,St,Ht] = solvePFDetLinElasAsymmetricNotchedPlateAdaptive(S_phase,S,T,PFsolver,C,BU,BL,BR,H1,H2,H3,PU,PL,PR,sizemap,...
                'pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions,'display');
        otherwise
            [dt,ut,ft,St_phase,St] = solvePFDetLinElasAsymmetricNotchedPlateAdaptive(S_phase,S,T,PFsolver,C,BU,BL,BR,H1,H2,H3,PU,PL,PR,sizemap,...
                'pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions,'display');
    end
    [fmax,idmax] = max(ft,[],2);
    t = gettevol(T);
    udmax = t(idmax);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'dt','ut','ft','St_phase','St','fmax','udmax','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        save(fullfile(pathname,'solution.mat'),'Ht','-append');
    end
else
    load(fullfile(pathname,'solution.mat'),'dt','ut','ft','St_phase','St','fmax','udmax','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        load(fullfile(pathname,'solution.mat'),'Ht');
    end
end

%% Outputs
fprintf('\n');
fprintf('setup    = %d\n',setup);
fprintf('PF model = %s\n',PFmodel);
fprintf('PF solver = %s\n',PFsolver);
fprintf('nb elements = %g (initial) - %g (final)\n',getnbelem(S),getnbelem(St{end}));
fprintf('nb nodes    = %g (initial) - %g (final)\n',getnbnode(S),getnbnode(St{end}));
fprintf('nb dofs     = %g (initial) - %g (final)\n',getnbddl(S),getnbddl(St{end}));
fprintf('nb time dofs = %g\n',getnbtimedof(T));
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

fprintf('fmax  = %g kN/mm\n',fmax*1e-6);
fprintf('udmax = %g mm\n',udmax*1e3);

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
    
    %% Display force-displacement curve
    figure('Name','Force-displacement')
    clf
    plot(t*1e3,ft*1e-6,'-b','Linewidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    %% Display solutions at different instants
    ampl = 0;
    switch setup
        case {1,4,5}
            rep = find(abs(t-0.190*unit)<eps | abs(t-0.201*unit)<eps | abs(t-0.203*unit)<eps | abs(t-0.204*unit)<eps | abs(t-0.205*unit)<eps | abs(t-0.207*unit)<eps);
        case {2,3}
            rep = find(abs(t-0.210*unit)<eps | abs(t-0.215*unit)<eps | abs(t-0.218*unit)<eps | abs(t-0.219*unit)<eps | abs(t-0.220*unit)<eps | abs(t-0.222*unit)<eps);
    end
    rep = [rep,length(T)];
    
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
