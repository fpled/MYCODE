%% Phase field fracture model - deterministic linear elasticity problem %%
%  Plate with a central circular hole under compression                 %%
%%----------------------------------------------------------------------%%
% [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS] (experimental tests)
% [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF] (anisotropic phase field model of Miehe et al.)
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME] (anisotropic phase field model of Nguyen et al.)

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

% test = true; % coarse mesh
test = false; % fine mesh

Dim = 2; % space dimension Dim = 2
symmetry = 'Isotropic'; % 'Isotropic' or 'Anisotropic'. Material symmetry
ang = 30; % clockwise material orientation angle around z-axis for anisotopic material [deg]
PFmodel = 'AnisotropicMiehe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicSpectral', 'AnisotropicHe'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'HistoryFieldElem'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = Inf; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment

% PFmodels = {'Isotropic','AnisotropicAmor','AnisotropicMiehe', 'AnisotropicSpectral','AnisotropicHe'};
% PFsplits = {'Strain','Stress'};
% PFregularizations = {'AT1','AT2'};
% PFsolvers = {'HistoryFieldElem','BoundConstrainedOptim'};
% maxIters = [1,100];

% for iPFmodel=1:length(PFmodels)
% PFmodel = PFmodels{iPFmodel};
% for iPFsplit=1:length(PFsplits)
% PFsplit = PFsplits{iPFsplit};
% for iPFRegularization=1:length(PFregularizations)
% PFregularization = PFregularizations{iPFRegularization};
% for iPFsolver=1:length(PFsolvers) 
% PFsolver = PFsolvers{iPFsolver};
% close all
% for imaxIter=1:length(maxIters) 
% maxIter = maxIters(imaxIter);

filename = ['phasefieldDetLinElas' symmetry];
if strcmpi(symmetry,'anisotropic') % anisotropic material
    filename = [filename num2str(ang) 'deg'];
end
filename = [filename 'PlatewithHole' PFmodel PFsplit PFregularization PFsolver...
    'MaxIter' num2str(maxIter) 'Tol' num2str(tolConv) '_' num2str(Dim) 'D'];

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
formats = {'epsc','fig'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    if Dim==2
        % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        L = 15e-3; % length
        h = 2*L; % height
        e = 1; % thickness
        r = 3e-3; % radius of the hole
        % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
        % L = 100e-3; % length
        % h = 65e-3; % height
        % e = 40e-3; % thickness
        % r = 2.5e-3; % radius of the hole
        D = DOMAIN(2,[0.0,0.0],[L,h]);
        C = CIRCLE(L/2,h/2,r);
    elseif Dim==3 % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
        L = 100e-3; % length
        h = 65e-3; % height
        e = 40e-3; % thickness
        r = 2.5e-3; % radius of the hole
        D = DOMAIN(3,[0.0,0.0,0.0],[L,h,e]);
        % C = CYLINDER
    end
    
    if Dim==2
        % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        clD = 0.06e-3; % characteristic length for domain
        clC = 0.06e-3; % characteristic length for circular hole
        % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
        % clD = 0.25e-3; % characteristic length for domain
        % clC = 0.05e-3; % characteristic length for circular hole
        if test
            clD = 0.25e-3;
            clC = 0.12e-3;
        end
    elseif Dim==3
        % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
        clD = 0.25e-3; % characteristic length for domain
        clC = 0.05e-3; % characteristic length for circular hole
        if test
            clD = 0.25e-3;
            clC = 0.12e-3;
        end
    end
    S_phase = gmshdomainwithhole(D,C,clD,clC,fullfile(pathname,'gmsh_domain_with_hole'));
    % S_phase = gmsh2femobject(2,fullfile(getfemobjectoptions('path'),'MYCODE','examples','phasefield','gmsh_domain_with_hole_Matthieu_Noel.msh'),2);
    S = S_phase;
    
    %% Phase field problem
    %% Material
    switch lower(symmetry)
        case 'isotropic' % isotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 1.4; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            % Regularization parameter (width of the smeared crack)
            l = 0.12e-3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        case 'anisotropic' % anisotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 10e3; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            % Regularization parameter (width of the smeared crack)
            l = 8.5e-6; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        otherwise
            error('Wrong material symmetry class');
    end
    % Small artificial residual stiffness
    % k = 1e-12;
    k = 0;
    
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
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'PFregularization',PFregularization);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
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
    option = 'DEFO'; % plane strain [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
    % option = 'CONT'; % plane stress [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
    switch lower(symmetry)
        case 'isotropic' % isotropic material
            % Lame coefficients
            % Young modulus and Poisson ratio
            E = 12e9; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            NU = 0.3; % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
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
        case 'isotropic' % isotropic material
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
        case 'anisotropic' % anisotropic material
            mat = ELAS_ANISOT('C',Cmat,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
        otherwise
            error('Wrong material symmetry class');
    end
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BU = LIGNE([0.0,h],[L,h]);
        BL = LIGNE([0.0,0.0],[L,0.0]);
        BRight = LIGNE([L,0.0],[L,h]);
        BLeft = LIGNE([0.0,0.0],[0.0,h]);
    elseif Dim==3
        BU = PLAN([0.0,h,0.0],[L,h,0.0],[0.0,h,e]);
        BL = PLAN([0.0,0.0,0.0],[L,0.0,0.0],[0.0,0.0,e]);
        BRight = PLAN([L,0.0,0.0],[L,h,0.0],[L,0.0,e]);
        BLeft = PLAN([0.0,0.0,0.0],[0.0,h,0.0],[0.0,0.0,e]);
    end
    P0 = getvertices(D);
    P0 = POINT(P0{1});
    
    S = final(S);
    
    ud = 0;
    if Dim==2
        S = addcl(S,BU,'UY',ud);
    elseif Dim==3
        S = addcl(S,BU,'UY',ud);
    end
    S = addcl(S,BL,'UY');
    S = addcl(S,P0,'UX');
    
    %% Stiffness matrices and sollicitation vectors
    % [A,b] = calc_rigi(S);
    % b = -b;
    
    %% Time scheme
    if Dim==2
        % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
        % du = 1e-3 mm during the first stage (until the phase field reaches the threshold value)
        % du = 1e-4 mm during the last stage (as soon as the phase field exceeds the threshold value)
        % dt0 = 1e-6;
        % dt1 = 1e-7;
        % tf = 25e-6;
        % dthreshold = 0.9;
        
        % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        % du = 8e-5 mm during the first stage (until the phase field reaches the threshold value)
        % du = 2e-5 mm during the last stage (as soon as the phase field exceeds the threshold value)
        dt0 = 8e-8;
        dt1 = 2e-8;
        if test
            dt0 = 16e-8;
            dt1 = 4e-8;
        end
        tf = 25e-6;
        dthreshold = 0.6;
    elseif Dim==3
        % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
        % du = 1e-3 mm during the first stage (until the phase field reaches the threshold value)
        % du = 1e-4 mm during the last stage (as soon as the phase field exceeds the threshold value)
        dt0 = 1e-6;
        dt1 = 1e-7;
        if test
            dt0 = 2e-6;
            dt1 = 2e-7;
        end
        tf = 25e-6;
        dthreshold = 0.9;
    end
    
%     t0 = linspace(dt0,nt0*dt0,nt0);
%     t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
%     t = [t0,t1];
%     
%     T = TIMEMODEL(t);
    T = struct('dt0',dt0,'dt1',dt1,'tf',tf,'dthreshold',dthreshold);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','BU','BL','BLeft','BRight','symmetry','ang');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','BU','BL','BLeft','BRight','symmetry','ang');
end

%% Solution
if solveProblem
    tTotal = tic;
    
    switch lower(PFsolver)
        case {'historyfieldelem','historyfieldnode'}
            [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElasPlatewithHoleThreshold(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P0,'maxiter',maxIter,'tol',tolConv);
        otherwise
            [dt,ut,ft,~,Edt,Eut,output] = solvePFDetLinElasPlatewithHoleThreshold(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P0,'maxiter',maxIter,'tol',tolConv);
    end
    [fmax,idmax] = max(ft,[],2);
    T = gettimemodel(dt);
    t = gettevol(T);
    udmax = t(idmax);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'dt','ut','ft','Edt','Eut','output','fmax','udmax','T','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        save(fullfile(pathname,'solution.mat'),'Ht','-append');
    end
else
    load(fullfile(pathname,'solution.mat'),'dt','ut','ft','Edt','Eut','output','fmax','udmax','T','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        load(fullfile(pathname,'solution.mat'),'Ht');
    end
end

%% Outputs
fprintf('\n');
fprintf('dim      = %d\n',Dim);
fprintf('mat sym  = %s\n',symmetry);
if strcmpi(symmetry,'anisotropic')
    fprintf('angle    = %g deg\n',ang);
end
fprintf('PF model = %s\n',PFmodel);
fprintf('PF split = %s\n',PFsplit);
fprintf('PF regularization = %s\n',PFregularization);
fprintf('PF solver = %s\n',PFsolver);
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('nb time dofs = %g\n',getnbtimedof(T));
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

fprintf('fmax  = %g kN/mm\n',fmax*1e-6);
fprintf('udmax = %g mm\n',udmax*1e3);

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
    
    u = getmatrixatstep(ut,rep(end));
    ampl = getsize(S)/max(abs(u))/20;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*unfreevector(S,u),'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
end

%% Display solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display force-displacement curve
    figure('Name','Force vs displacement')
    clf
    plot(t*1e3,ft*((Dim==2)*1e-6+(Dim==3)*1e-3),'-b','Linewidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    %% Display energy-displacement curves
    figure('Name','Energies vs displacement')
    clf
    plot(t*1e3,Eut,'-b','Linewidth',linewidth)
    hold on
    plot(t*1e3,Edt,'-r','Linewidth',linewidth)
    plot(t*1e3,Eut+Edt,'-k','Linewidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Energy [J]','Interpreter',interpreter)
    l = legend('$\Psi_u$','$\Psi_c$','$\Psi_{\mathrm{tot}}$',...
        'Location','NorthWest');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'energies_displacement',formats);
    mymatlab2tikz(pathname,'energies_displacement.tex');

    %% Display outputs of iterative resolution
    figure('Name','Number of iterations vs displacement')
    clf
    plot(t*1e3,output.iteration,'-b','Linewidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Number of iterations','Interpreter',interpreter)
    mysaveas(pathname,'nb_iterations_displacement',formats);
    mymatlab2tikz(pathname,'nb_iterations_displacement.tex');
    
    figure('Name','Computing time vs displacement')
    clf
    plot(t*1e3,output.time,'-r','Linewidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Computing time [s]','Interpreter',interpreter)
    mysaveas(pathname,'cpu_time_displacement',formats);
    mymatlab2tikz(pathname,'cpu_time_displacement.tex');
    
    figure('Name','Error vs displacement')
    clf
    plot(t*1e3,output.error,'-k','Linewidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Error','Interpreter',interpreter)
    mysaveas(pathname,'error_displacement',formats);
    mymatlab2tikz(pathname,'error_displacement.tex');
    
    %% Display solutions at different instants
    ampl = 0;
    rep = find(abs(t-18.5e-6)<eps | abs(t-24.6e-6)<eps);
    rep = [rep,length(T)];
    
    for j=1:length(rep)
        dj = getmatrixatstep(dt,rep(j));
        uj = getmatrixatstep(ut,rep(j));
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            Hj = getmatrixatstep(Ht,rep(j));
        end
        
        plotSolution(S_phase,dj);
        mysaveas(pathname,['damage_t' num2str(rep(j))],formats,renderer);
        
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
        % plotSolution(S,uj,'energyint','','ampl',ampl);
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
    % evolSolution(S,ut,'energyint','','ampl',ampl,'FrameRate',framerate,'filename','internal_energy_density','pathname',pathname,options{:});
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
        ui = getmatrixatstep(ut,rep(i));
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            Hi = getmatrixatstep(Ht,rep(i));
        end
        
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
    make_pvd_file(pathname,'solution',1,length(T));
end

% myparallel('stop');

% end
% end
% end
% end
% end