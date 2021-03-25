%% Phase field fracture model - stochastic linear elasticity problem with single edge crack %%
%%------------------------------------------------------------------------------------------%%
% [Bourdin, Francfort, Marigo, 2000, JMPS] (isotropic phase field model with no split of Bourdin et al.)
% [Miehe, Welschinger, Hofacker, 2010 IJNME] (anisotropic phase field model of Miehe et al.)
% [Miehe, Hofacker, Welschinger, 2010, CMAME] (anisotropic phase field model of Miehe et al.)
% [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME] (anisotropic phase field model of Miehe et al.)
% [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM] (anisotropic phase field model of Miehe et al.)
% [Ambati, Gerasimov, De Lorenzis, 2015, CM] (hybrid isotropic-anisotropic phase field model of Ambati et al. compared with the isotropic one of Bourdin et al. and the anisotropic ones of Amor et al. and Miehe et al.)
% [Liu, Li, Msekh, Zuo, 2016, CMS] (anisotropic phase field model of Miehe et al.)
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM] (anisotropic phase field model of Wu et al.)
% [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US] (anisotropic phase field model of Amor et al.)

% clc
clearvars
close all
% rng('default');
myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = false;

test = true; % coarse mesh
% test = false; % fine mesh

Dim = 2; % space dimension Dim = 2, 3
loading = 'Shear'; % 'Tension' or 'Shear'
PFmodel = 'Isotropic'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe'
randMat = true; % random material parameters (true or false)
randPF = false; % random phase field parameters (true or false)

filename = ['phasefieldStoLinElasSingleEdgeCrack' loading PFmodel '_' num2str(Dim) 'D'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','phasefield',filename);
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
        D = DOMAIN(2,[0.0,0.0],[L,L]);
        C = LIGNE([0.0,L/2],[a,L/2]);
    elseif Dim==3
        e = 0.1e-3;
        D = DOMAIN(3,[0.0,0.0,0.0],[L,L,e]);
        C = QUADRANGLE([0.0,L/2,0.0],[a,L/2,0.0],[a,L/2,e],[0.0,L/2,e]);
    end
    
    if Dim==2
        % clD = 2e-5; % [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
        % clC = 2e-6; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
        % clC = 1e-6; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
        % clC = 6e-7; % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
        clD = 3.9e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM]
        clC = 3.9e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM]
        if test
            clD = 1.5e-5;
            clC = 1.5e-5;
            % clD = 4e-5;
            % clC = 1e-5;
        end
    elseif Dim==3
        clD = 4e-5;
        clC = 4e-6;
        if test
            clD = 4e-5;
            clC = 1e-5;
        end
    end
    S_phase = gmshdomainwithedgecrack(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'));
    S = S_phase;
    
    %% Phase field problem
    %% Material
    % Critical energy release rate (or fracture toughness)
    gc = 2.7e3;
    % Regularization parameter (width of the smeared crack)
    % l = 1e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME]
    % l = 1.5e-5; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM]
    % l = 3.75e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME]
    l = 7.5e-6; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS]
    % l = 4e-6; % [Ambati, Gerasimov, De Lorenzis, 2015, CM]
    % eta = 0.052; w0 = 75.94; l = eta/sqrt(w0)*1e-3; % l = 6e-7; % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
    % Small artificial residual stiffness
    k = 1e-10;
    % Internal energy
    H = 0;
    
    % Material
    mat_phase = FOUR_ISOT('k',gc*l,'r',gc/l+2*H);
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
    % A_phase = K_phase + M_phase;
    
    % l_phase = LINFORM(0,2*H,0); % nodal values
    % l_phase = setfree(l_phase,1);
    % b_phase = b_phase + calc_vector(l_phase,S_phase);
    
    % [A_phase,b_phase] = calc_rigi(S_phase);
    % b_phase = -b_phase + bodyload(S_phase,[],'QN',2*H);
    
    %% Linear elastic displacement field problem
    %% Materials
    % Option
    option = 'DEFO'; % plane strain
    % option = 'CONT'; % plane stress [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS]
    % Lame coefficients
    % lambda = 121.1538e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
    % mu = 80.7692e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
    lambda = 121.15e9;
    mu = 80.77e9;
    % Young modulus and Poisson ratio
    if Dim==2
        switch lower(option)
            case 'defo'
                E = mu*(3*lambda+2*mu)/(lambda+mu);
                NU = lambda/(lambda+mu)/2;
            case 'cont'
                E = 4*mu*(lambda+mu)/(lambda+2*mu);
                NU = lambda/(lambda+2*mu);
        end
        % E = 210e9; NU = 0.2; % [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM]
        % kappa = 121030e6; NU=0.227; lambda=3*kappa*NU/(1+NU); mu = 3*kappa*(1-2*NU)/(2*(1+NU)); E = 3*kappa*(1-2*NU); % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
    elseif Dim==3
        E = mu*(3*lambda+2*mu)/(lambda+mu);
        NU = lambda/(lambda+mu)/2;
    end
    % Energetic degradation function
    g = @(d) (1-d).^2;
    % Density
    RHO = 1;
    
    % Material
    d = calc_init_dirichlet(S_phase);
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel);
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
            S = addcl(S,BU,'UY',ud);
            S = addcl(S,BL,'UY');
            if Dim==2
                S = addcl(S,POINT([0.0,0.0]),'UX');
            elseif Dim==3
                S = addcl(S,POINT([0.0,0.0,0.0]),{'UX','UZ'});
            end
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
    if Dim==2
        switch lower(loading)
            case 'tension'
                % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM]
                % du = 1e-5 mm during the first 500 time steps (up to u = 5e-3 mm)
                % du = 1e-6 mm during the last 1300 time steps (up to u = 6.3e-3 mm)
                dt0 = 1e-8;
                nt0 = 500;
                if test
                    dt0 = 1e-7;
                    nt0 = 50;
                end
                t0 = linspace(dt0,nt0*dt0,nt0);
                dt1 = 1e-9;
                nt1 = 1300;
                if test
                    dt1 = 1e-8;
                    nt1 = 130;
                end
                %
                t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                t = [t0,t1];
                
                % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
                % dt = 1e-8;
                % nt = 630;
                % t = linspace(dt,nt*dt,nt);
                
                % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
                % dt = 1e-7;
                % nt = 63;
                % t = linspace(dt,nt*dt,nt);
                
                % [Liu, Li, Msekh, Zuo, 2016, CMS]
                % du = 1e-4 mm during the first 50 time steps (up to u = 5e-3 mm)
                % du = 1e-6 mm during the last 1300 time steps (up to u = 6.3e-3 mm)
                % dt0 = 1e-7;
                % nt0 = 50;
                % t0 = linspace(dt0,nt0*dt0,nt0);
                % dt1 = 1e-9;
                % nt1 = 1300;
                % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                % t = [t0,t1];
            case 'shear'
                % [Miehe, Welschinger, Hofacker, 2010 IJNME]
                % du = 1e-4 mm during the first 100 time steps (up to u = 10e-3 mm)
                % du = 1e-6 mm during the last 10 000 time steps (up to u = 20e-3 mm)
                % dt0 = 1e-7;
                % nt0 = 100;
                % t0 = linspace(dt0,nt0*dt0,nt0);
                % dt1 = 1e-9;
                % nt1 = 10000;
                % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                % t = [t0,t1];
                
                % [Liu, Li, Msekh, Zuo, 2016, CMS]
                % du = 1e-4 mm during the first 50 time steps (up to u = 5e-3 mm)
                % du = 1e-5 mm during the last 1500 time steps (up to u = 20e-3 mm)
                % dt0 = 1e-7;
                % nt0 = 50;
                % t0 = linspace(dt0,nt0*dt0,nt0);
                % dt1 = 1e-8;
                % nt1 = 1500;
                % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                % t = [t0,t1];
                
                % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
                % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM], [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
                dt = 1e-8;
                nt = 2000;
                if test
                    dt = 1e-7;
                    nt = 200;
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
    T = TIMEMODEL(t);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','BU','BL','BRight','BLeft','BFront','BBack','loading');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','BU','BL','BRight','BLeft','BFront','BBack','loading');
end

%% Solution
if solveProblem
    
    % Number of samples
    if test
        N = 10;
    else
        N = 5e2;
    end
    
    %% Random variables
    % Material properties
    if randMat % random material parameters
        % la = -24; % la < 1/5. Parameter controlling the level of statistical fluctuation
        % delta1 = 1/sqrt(1-la); % coefficient of variation for bulk modulus
        % delta2 = 1/sqrt(1-5*la); % coefficient of variation for shear modulus
        delta1 = 0.1; % coefficient of variation for bulk modulus
        la = 1 - 1/delta1^2; % la < 1/5. Parameter controlling the level of statistical fluctuation
        delta2 = 1/sqrt(5/delta1^2 - 4); % coefficient of variation for shear modulus
        
        mC1 = E/3/(1-2*NU); % mean bulk modulus
        mC2 = mu; % mean shear modulus
        la1 = (1-la)/mC1; % la1 > 0
        la2 = (1-5*la)/mC2; % la2 > 0
        
        a1 = 1-la; % a1 > 0
        b1 = 1/la1; % b1 > 0
        a2 = 1-5*la; % a2 > 0
        b2 = 1/la2; % b2 > 0
        
        % Sample set
        C_sample(:,1) = gamrnd(a1,b1,N,1); % [Pa]
        C_sample(:,2) = gamrnd(a2,b2,N,1); % [Pa]
        % lambda_sample = C_sample(:,1) - 2/3*C_sample(:,2); % [Pa]
        E_sample = (9*C_sample(:,1).*C_sample(:,2))./(3*C_sample(:,1)+C_sample(:,2)); % [Pa]
        NU_sample = (3*C_sample(:,1)-2*C_sample(:,2))./(6*C_sample(:,1)+2*C_sample(:,2));
    else
        E_sample = E*ones(N,1);
        NU_sample = NU*ones(N,1);
    end
    
    % Phase field properties
    if randPF % random phase field parameters
        delta3 = 0.1; % coefficient of variation of fracture toughness
        delta4 = 0.1; % coefficient of variation of regularization parameter
        a3 = 1/delta3^2;
        b3 = gc/a3;
        a4 = 1/delta4^2;
        b4 = l/a4;
        gc_sample = gamrnd(a3,b3,N,1);
        l_sample = gamrnd(a4,b4,N,1);
    else
        gc_sample = gc*ones(N,1);
        l_sample = l*ones(N,1);
    end
    
    samples = [E_sample,NU_sample,gc_sample,l_sample];

    %% Solution
    tTotal = tic;
    
    fun = @(S,S_phase) solvePFDetLinElasSingleEdgeCrack(S,S_phase,T,BU,BL,BRight,BLeft,BFront,BBack,loading);
    [Ht,dt,ut,ft] = solvePFStoLinElas(S,S_phase,T,fun,samples,'display');
    fmax = max(ft,[],2);
    
    time = toc(tTotal);
    
    %% Statistical outputs of solution
    sz_phase = [getnbddl(S_phase),getnbtimedof(T)];
    sz = [getnbddl(S),getnbtimedof(T)];

    mean_Ht = mean(Ht);
    mean_Ht = reshape(mean_Ht,sz_phase);
    mean_Ht = TIMEMATRIX(mean_Ht,T);
    
    mean_dt = mean(dt);
    mean_dt = reshape(mean_dt,sz_phase);
    mean_dt = TIMEMATRIX(mean_dt,T);
    
    mean_ut = mean(ut);
    mean_ut = reshape(mean_ut,sz);
    mean_ut = TIMEMATRIX(mean_ut,T);
    
    probs = [0.025 0.975];
    
    mean_ft = mean(ft);
    std_ft = std(ft);
    ci_ft = quantile(ft,probs);
    
    mean_fmax = mean(fmax);
    std_fmax = std(fmax);
    ci_fmax = quantile(fmax,probs);
    
    npts = 100;
    [f_fmax,xi_fmax,bw_fmax] = ksdensity(fmax,'npoints',npts);

    save(fullfile(pathname,'solution.mat'),'N','mean_Ht','mean_dt','mean_ut',...
        'mean_ft','std_ft','ci_ft','fmax',...
        'mean_fmax','std_fmax','ci_fmax','probs','f_fmax','xi_fmax','bw_fmax','time');
else
    load(fullfile(pathname,'solution.mat'),'N','mean_Ht','mean_dt','mean_ut',...
        'mean_ft','std_ft','ci_ft','fmax',...
        'mean_fmax','std_fmax','ci_fmax','probs','f_fmax','xi_fmax','bw_fmax','time');
end

%% Outputs
fprintf('\n');
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('nb time dofs = %g\n',getnbtimedof(T));
fprintf('nb samples = %g\n',N);
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

fprintf('mean(fmax)    = %e\n',mean_fmax);
fprintf('std(fmax)     = %e\n',std_fmax);
fprintf('disp(fmax)    = %e\n',std_fmax/mean_fmax);
fprintf('%d%% ci(fmax)  = [%e,%e]\n',(probs(2)-probs(1))*100,ci_fmax(1),ci_fmax(2));

%% Display
if displaySolution
    [t,rep] = gettevol(T);
    mean_u = getmatrixatstep(mean_ut,rep(end));
    
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
    
    ampl = getsize(S)/max(abs(mean_u))/20;
    plotModelDeflection(S,mean_u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*unfreevector(S,mean_u),'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
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
    if Dim==2
        ylabel('Force [kN/mm]','Interpreter',interpreter)
    elseif Dim==3
        ylabel('Force [kN]','Interpreter',interpreter)
    end
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
    if Dim==2
        xlabel('$f$ [kN/mm]','Interpreter',interpreter)
    elseif Dim==3
        xlabel('$f$ [kN]','Interpreter',interpreter)
    end
    ylabel('$p_{F_c}(f)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_fmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_fmax.tex');
    
    %% Display evolution of mean solutions
    ampl = 0;
    % ampl = getsize(S)/max(max(abs(getvalue(mean_ut))))/20;
    
    options = {'plotiter',true,'plottime',false};
    framerate = 80;
    
%     evolSolution(S_phase,mean_Ht,'FrameRate',framerate,'filename','mean_internal_energy','pathname',pathname,options{:});
    
%     evolSolution(S_phase,mean_dt,'FrameRate',framerate,'filename','mean_damage','pathname',pathname,options{:});
%     for i=1:Dim
%         evolSolution(S,mean_ut,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['mean_displacement_' num2str(i)],'pathname',pathname,options{:});
%     end
    
%     for i=1:(Dim*(Dim+1)/2)
%         evolSolution(S,mean_ut,'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['mean_epsilon_' num2str(i)],'pathname',pathname,options{:});
%         evolSolution(S,mean_ut,'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['mean_sigma_' num2str(i)],'pathname',pathname,options{:});
%     end
%     
%     evolSolution(S,mean_ut,'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename','mean_epsilon_von_mises','pathname',pathname,options{:});
%     evolSolution(S,mean_ut,'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename','mean_sigma_von_mises','pathname',pathname,options{:});
    
    %% Display mean solutions at differents instants
    switch lower(loading)
        case 'tension'
            rep = find(abs(t-5.5e-6)<eps | abs(t-5.75e-5)<eps | abs(t-6e-6)<eps | abs(t-6.25e-6)<eps);
        case 'shear'
            rep = find(abs(t-1e-5)<eps | abs(t-1.25e-5)<eps | abs(t-1.35e-5)<eps | abs(t-1.5e-5)<eps);
        otherwise
            error('Wrong loading case')
    end
    for j=1:length(rep)
        close all
        mean_Hj = getmatrixatstep(mean_Ht,rep(j));
        mean_dj = getmatrixatstep(mean_dt,rep(j));
        mean_uj = getmatrixatstep(mean_ut,rep(j));
        
%         plotSolution(S_phase,mean_Hj);
%         mysaveas(pathname,['mean_internal_energy_t' num2str(rep(j))],formats,renderer);
        
        plotSolution(S_phase,mean_dj);
        mysaveas(pathname,['mean_damage_t' num2str(rep(j))],formats,renderer);
        
        for i=1:Dim
            plotSolution(S,mean_uj,'displ',i,'ampl',ampl);
            mysaveas(pathname,['mean_displacement_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        end
        
%         for i=1:(Dim*(Dim+1)/2)
%             plotSolution(S,mean_uj,'epsilon',i,'ampl',ampl);
%             mysaveas(pathname,['mean_epsilon_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
%             
%             plotSolution(S,mean_uj,'sigma',i,'ampl',ampl);
%             mysaveas(pathname,['mean_sigma_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
%         end
%         
%         plotSolution(S,mean_uj,'epsilon','mises','ampl',ampl);
%         mysaveas(pathname,['mean_epsilon_von_mises_t' num2str(rep(j))],formats,renderer);
%         
%         plotSolution(S,mean_uj,'sigma','mises','ampl',ampl);
%         mysaveas(pathname,['mean_sigma_von_mises_t' num2str(rep(j))],formats,renderer);
    end
    
end

%% Save mean solutions
[t,rep] = gettevol(T);
for i=1:length(T)
    mean_Hi = getmatrixatstep(mean_Ht,rep(i));
    mean_di = getmatrixatstep(mean_dt,rep(i));
    mean_ui = getmatrixatstep(mean_ut,rep(i));
    
    write_vtk_mesh(S,{mean_Hi,mean_di,mean_ui},[],...
        {'internal energy','damage','displacement'},[],...
        pathname,'mean_solution',1,i-1);
end
make_pvd_file(pathname,'mean_solution',1,length(T));

myparallel('stop');
