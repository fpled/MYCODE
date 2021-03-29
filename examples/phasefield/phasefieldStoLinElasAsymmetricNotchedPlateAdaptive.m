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
% [Ambati, Gerasimov, De Lorenzis, 2015, CM] (hybrid isotropic-anisotropic phase field model of Ambati et al. compared with the isotropic one of Bourdin et al. and the anisotropic ones of Amor et al. and Miehe et al.)
% [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME] (isotropic phase field model with no split of Bourdin et al. compared to experimental data of [Winkler PhD thesis, 2001])
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM] (anisotropic phase field model of Wu et al.)

% clc
clearvars
close all
% rng('default');
myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = false;

test = true; % coarse mesh and small number of samples
% test = false; % fine mesh and high number of samples

setup = 1; % notch geometry setup = 1, 2
PFmodel = 'AnisotropicMiehe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe'
randMat = true; % random material parameters (true or false)
randPF = true; % random phase field parameters (true or false)

filename = ['phasefieldStoLinElasAsymmetricNotchedPlateSetup' num2str(setup) PFmodel];
if randMat
    filename = [filename 'RandMat'];
end
if randPF
    filename = [filename 'RandPF'];
end
filename = [filename 'Adaptive'];
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

gmshoptions = '-v 0';
mmgoptions = '-nomove -v -1';
% gmshoptions = '-v 5';
% mmgoptions = '-nomove -v 1';

%% Problem
if setProblem
    %% Domains and meshes
    unit = 1e-3; % for mm
    % unit = 25.4e-3; % for inch % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    switch setup
        case 1
            a = 1.5*unit; % crack length
            b = 5*unit; % crack offset from the centerline
            % a = 2.5*unit; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
            % b = 6*unit; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
        case 2
            a = 1*unit; % crack length
            b = 6*unit; % crack offset from the centerline
    end
    h = 4*unit;
    C = LIGNE([-b,-h],[-b,-h+a]);
    clD = 0.1*unit; % characteristic length for domain
    % cl = 0.01*unit; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM]
    cl = 0.025*unit/2; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME]
    % cl = 0.01*unit/2; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME]
    % cl = 0.01*unit/5; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    if test
        clD = 0.2*unit;
        cl = 0.1*unit;
    end
    clC = cl; % characteristic length for edge crack/notch
    clH = cl; % characteristic length for circular holes
    % S_phase = gmshasymmetricnotchedplate(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_domain_asymmetric_notched_plate'),2,'gmshoptions',gmshoptions);
    S_phase = gmshasymmetricnotchedplatewithedgesmearedcrack(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_domain_asymmetric_notched_plate'),2,'gmshoptions',gmshoptions);
    
    c = a/1e3;
    CL = LIGNE([-b-c/2,-h],[-b,-h+a]);
    CR = LIGNE([-b+c/2,-h],[-b,-h+a]);
    
    % sizemap = @(d) (clC-clD)*d+clD;
    sizemap = @(d) clD*clC./((clD-clC)*d+clC);
    
    %% Phase field problem
    %% Material
    % Critical energy release rate (or fracture toughness)
    gc = 1e3;
    % gc = 304.321; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % Regularization parameter (width of the smeared crack)
    l = 0.025*unit; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM]
    % l = 0.01*unit; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % Small artificial residual stiffness
    k = 1e-10;
    % Internal energy
    H = 0;
    
    % Material
    mat_phase = FOUR_ISOT('k',gc*l,'r',gc/l+2*H);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    BU = CIRCLE(0.0,h,2.5*unit);
    BL = CIRCLE(-9*unit,-h,2.5*unit);
    BR = CIRCLE(9*unit,-h,2.5*unit);
    
    S_phase = final(S_phase,'duplicate');
    S_phase = addcl(S_phase,CL,'T',1);
    S_phase = addcl(S_phase,CR,'T',1);
    S_phase = addcl(S_phase,BU,'T');
    S_phase = addcl(S_phase,BL,'T');
    S_phase = addcl(S_phase,BR,'T');
    
    d = calc_init_dirichlet(S_phase);
    cl = sizemap(d);
    S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_domain_asymmetric_notched_plate'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    S = S_phase;
    
    S_phase = setmaterial(S_phase,mat_phase);
    S_phase = final(S_phase,'duplicate');
    S_phase = addcl(S_phase,CL,'T',1);
    S_phase = addcl(S_phase,CR,'T',1);
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
    % Lame coefficients
    lambda = 12e9;
    mu = 8e9;
    % Young modulus and Poisson ratio
    switch lower(option)
        case 'defo'
            E = mu*(3*lambda+2*mu)/(lambda+mu);
            NU = lambda/(lambda+mu)/2;
        case 'cont'
            E = 4*mu*(lambda+mu)/(lambda+2*mu);
            NU = lambda/(lambda+2*mu);
    end
    % E = 3.102e9; % NU = 0.35; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    % Energetic degradation function
    g = @(d) (1-d).^2;
    % Thickness
    DIM3 = 1;
    % Density
    RHO = 1;
    
    % Material
    d = calc_init_dirichlet(S_phase);
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel);
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    B = LIGNE([-10*unit,h],[10*unit,h]);
    PU = POINT([0.0,h]);
    PL = POINT([-9*unit,-h]);
    PR = POINT([9*unit,-h]);
    
    S = final(S,'duplicate');
    
    ud = 0;
    S = addcl(S,PU,'UY',ud);
    S = addcl(S,PL,{'UX','UY'});
    S = addcl(S,PR,'UY');
    
    %% Stiffness matrices and sollicitation vectors
    % [A,b] = calc_rigi(S);
    % b = -b;
    
    %% Time scheme
    % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM]
    % dt = 1e-4*unit;
    % nt = 2500;
    % t = linspace(dt,nt*dt,nt);
    
    % [Ambati, Gerasimov, De Lorenzis, 2015, CM]
    % du = 1e-3 mm during the first 200 time steps (up to u = 0.2 mm)
    % du = 1e-4 mm during the last  500 time steps (up to u = 0.25 mm)
    dt0 = 1e-3*unit;
    nt0 = 200;
    t0 = linspace(dt0,nt0*dt0,nt0);
    dt1 = 1e-4*unit;
    nt1 = 500;
    t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
    t = [t0,t1];
    
    T = TIMEMODEL(t);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'unit','T','S_phase','S','sizemap','C','CL','CR','B','BU','BL','BR','PU','PL','PR');
else
    load(fullfile(pathname,'problem.mat'),'unit','T','S_phase','S','sizemap','C','CL','CR','B','BU','BL','BR','PU','PL','PR');
end

%% Solution
if solveProblem
    
    % Number of samples
    if test
        N = 8;
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
    
    nbSamples = 3;
    fun = @(S_phase,S,filename) solvePFDetLinElasAsymmetricNotchedPlateAdaptive(S_phase,S,T,CL,CR,BU,BL,BR,PU,PL,PR,sizemap,'filename',filename,'pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    [ft,dt,ut,Ht,St_phase,St] = solvePFStoLinElasAdaptive(S_phase,S,T,fun,samples,'filename','gmsh_domain_asymmetric_notched_plate','pathname',pathname,'nbsamples',nbSamples);
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

    save(fullfile(pathname,'solution.mat'),'N','dt','ut','Ht','St_phase','St',...
        'mean_ft','std_ft','ci_ft','fmax',...
        'mean_fmax','std_fmax','ci_fmax','probs','f_fmax','xi_fmax','bw_fmax','time');
else
    load(fullfile(pathname,'solution.mat'),'N','dt','ut','Ht','St_phase','St',...
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

fprintf('mean(fmax)    = %g kN/mm\n',mean_fmax*1e-6);
fprintf('std(fmax)     = %g kN/mm\n',std_fmax*1e-6);
fprintf('disp(fmax)    = %g\n',std_fmax/mean_fmax);
fprintf('%d%% ci(fmax)  = [%g,%g] kN/mm\n',(probs(2)-probs(1))*100,ci_fmax(1)*1e-6,ci_fmax(2)*1e-6);

%% Display
if displaySolution
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
    
    % DO NOT WORK WITH MESH ADAPTATION
    % u = getmatrixatstep(ut,rep(end));
%     u = ut(:,rep(end));
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
    
    %% Display force-displacement curve
    figure('Name','Force-displacement')
    clf
    ciplot(ci_ft(1,:)*1e-6,ci_ft(2,:)*1e-6,t*1e3,'b');
    alpha(0.2)
    hold on
    plot(t*1e3,mean_ft*1e-6,'-b','Linewidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN/mm]','Interpreter',interpreter)
    l = legend({['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value'},'Location','NorthWest');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    %% Display pdf of critical force
    figure('Name','Probability Density Estimate: Critical force')
    clf
    plot(xi_fmax*1e-6,f_fmax,'-b','LineWidth',linewidth)
    hold on
    ind_fmax = find(xi_fmax>=ci_fmax(1) & xi_fmax<ci_fmax(2));
    area(xi_fmax(ind_fmax)*1e-6,f_fmax(ind_fmax),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(mean_fmax*1e-6,0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$f$ [kN/mm]','Interpreter',interpreter)
    ylabel('$p_{F_c}(f)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_fmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_fmax.tex');
    
    %% Display evolution of samples of solutions
    ampl = 0;
    % DO NOT WORK WITH MESH ADAPTATION
    % ampl = getsize(S)/max(max(abs(getvalue(ut))))/20;
    % umax = cellfun(@(u) max(abs(u)),ut,'UniformOutput',false);
    % ampl = getsize(S)/max([umax{:}])/20;
    
    options = {'plotiter',true,'plottime',false};
    framerate = 80;
    
%     for k=1:size(St,1)
%         evolModel(T,St(k,:),'FrameRate',framerate,'filename',['mesh_sample_' num2str(k)],'pathname',pathname,options{:});
%         
%         evolSolutionCell(T,St_phase(k,:),dt(k,:),'FrameRate',framerate,'filename',['damage_sample_' num2str(k)],'pathname',pathname,options{:});
%         for i=1:2
%             evolSolutionCell(T,St(k,:),ut(k,:),'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i) '_sample_' num2str(k)],'pathname',pathname,options{:});
%         end
%         
% %         for i=1:3
% %             evolSolutionCell(T,St(k,:),ut(k,:),'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['epsilon_' num2str(i) '_sample_' num2str(k)],'pathname',pathname,options{:});
% %             evolSolutionCell(T,St(k,:),ut(k,:),'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['sigma_' num2str(i) '_sample_' num2str(k)],'pathname',pathname,options{:});
% %         end
% %         
% %         evolSolutionCell(T,St(k,:),ut(k,:),'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename',['epsilon_von_mises_sample_' num2str(k)],'pathname',pathname,options{:});
% %         evolSolutionCell(T,St(k,:),ut(k,:),'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename',['sigma_von_mises_sample_' num2str(k)],'pathname',pathname,options{:});
% %         
% %         evolSolutionCell(T,St_phase(k,:),Ht(k,:),'FrameRate',framerate,'filename',['internal_energy_sample_' num2str(k)],'pathname',pathname,options{:});
%     end
    
    %% Display samples of solutions at different instants
    rep = find(abs(t-0.210*unit)<eps | abs(t-0.215*unit)<eps | abs(t-0.218*unit)<eps | abs(t-0.220*unit)<eps | abs(t-0.222*unit)<eps);
    for k=1:size(St,1)
    for j=1:length(rep)
        close all
        % DO NOT WORK WITH MESH ADAPTATION
        % dj = getmatrixatstep(dt,rep(j));
        % uj = getmatrixatstep(ut,rep(j));
        % Hj = getmatrixatstep(Ht,rep(j));
        dj = dt{k,rep(j)};
        uj = ut{k,rep(j)};
        Hj = Ht{k,rep(j)};
        Sj = St{k,rep(j)};
        Sj_phase = St_phase{k,rep(j)};
        
        plotModel(Sj,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
        mysaveas(pathname,['mesh_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        
        plotSolution(Sj_phase,dj);
        mysaveas(pathname,['damage_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        
        for i=1:2
            plotSolution(Sj,uj,'displ',i,'ampl',ampl);
            mysaveas(pathname,['displacement_' num2str(i) '_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        end
        
%         for i=1:3
%             plotSolution(Sj,uj,'epsilon',i,'ampl',ampl);
%             mysaveas(pathname,['epsilon_' num2str(i) '_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
%             
%             plotSolution(Sj,uj,'sigma',i,'ampl',ampl);
%             mysaveas(pathname,['sigma_' num2str(i) '_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
%         end
%         
%         plotSolution(Sj,uj,'epsilon','mises','ampl',ampl);
%         mysaveas(pathname,['epsilon_von_mises_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
%         
%         plotSolution(Sj,uj,'sigma','mises','ampl',ampl);
%         mysaveas(pathname,['sigma_von_mises_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
        
%         plotSolution(Sj_phase,Hj);
%         mysaveas(pathname,['internal_energy_sample_' num2str(k) '_t' num2str(rep(j))],formats,renderer);
    end
    end
    
end

%% Save samples of solutions
[t,rep] = gettevol(T);
for k=1:size(St,1)
for i=1:length(T)
    % DO NOT WORK WITH MESH ADAPTATION
    % di = getmatrixatstep(dt,rep(i));
    % ui = getmatrixatstep(ut,rep(i));
    % Hi = getmatrixatstep(Ht,rep(i));
    di = dt{k,rep(i)};
    ui = ut{k,rep(i)};
    Hi = Ht{k,rep(i)};
    Si = St{k,rep(i)};
    % Si_phase = St_phase{k,rep(i)};
    
    write_vtk_mesh(Si,{di,ui,Hi},[],...
        {'damage','displacement','internal energy'},[],...
        pathname,['solution_sample_' num2str(k)],1,i-1);
end
make_pvd_file(pathname,['solution_sample_' num2str(k)],1,length(T));
end

myparallel('stop');
