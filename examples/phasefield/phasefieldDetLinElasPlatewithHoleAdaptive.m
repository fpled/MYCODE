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
saveParaview = true;

test = true; % coarse mesh
% test = false; % fine mesh

Dim = 2; % space dimension Dim = 2
symmetry = 'Isotropic'; % 'Isotropic' or 'Anisotropic'. Material symmetry
ang = 30; % clockwise material orientation angle around z-axis for anisotopic material [deg]
PFmodel = 'AnisotropicMiehe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'

switch lower(symmetry)
    case 'isotropic' % isotropic material
        filename = ['phasefieldDetLinElas' symmetry 'PlatewithHole' PFmodel PFsolver 'Adaptive_' num2str(Dim) 'D'];
    case 'anisotropic' % anisotropic material
        filename = ['phasefieldDetLinElas' symmetry num2str(ang) 'deg' 'PlatewithHole' PFmodel PFsolver 'Adaptive_' num2str(Dim) 'D'];
    otherwise
        error('Wrong material symmetry class');
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
        % clD = 0.06e-3; % characteristic length for domain
        % clC = 0.06e-3; % characteristic length for circular hole
        % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
        clD = 0.25e-3; % characteristic length for domain
        clC = 0.05e-3; % characteristic length for circular hole
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
    
    sizemap = @(d) (clC-clD)*d+clD;
    % sizemap = @(d) clD*clC./((clD-clC)*d+clC);
    
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
    % k = 1e-10;
    k = 0;
    % Internal energy
    H = 0;
    
    % Material
    mat_phase = FOUR_ISOT('k',gc*l,'r',gc/l+2*H);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    S_phase = final(S_phase);
    S_phase = addcl(S_phase,C,'T',1);
    
    d = calc_init_dirichlet(S_phase);
    cl = sizemap(d);
    S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_domain_with_hole'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    S = S_phase;
    
    S_phase = setmaterial(S_phase,mat_phase);
    S_phase = final(S_phase);
    
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
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel);
        case 'anisotropic' % anisotropic material
            mat = ELAS_ANISOT('C',Cmat,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel);
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
    elseif Dim==3
        BU = PLAN([0.0,h,0.0],[L,h,0.0],[0.0,h,e]);
        BL = PLAN([0.0,0.0,0.0],[L,0.0,0.0],[0.0,0.0,e]);
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
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','D','C','BU','BL','symmetry','ang');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','D','C','BU','BL','symmetry','ang');
end

%% Solution
if solveProblem
    tTotal = tic;
    
    switch lower(PFsolver)
        case {'historyfieldelem','historyfieldnode'}
            [dt,ut,ft,T,St_phase,St,Ht] = solvePFDetLinElasPlatewithHoleAdaptiveThreshold(S_phase,S,T,PFsolver,BU,BL,P0,C,sizemap,...
                'pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions,'display');
        otherwise
            [dt,ut,ft,T,St_phase,St] = solvePFDetLinElasPlatewithHoleAdaptiveThreshold(S_phase,S,T,PFsolver,BU,BL,P0,C,sizemap,...
                'pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions,'display');
    end
    [fmax,idmax] = max(ft,[],2);
    t = gettevol(T);
    udmax = t(idmax);

    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'dt','ut','ft','T','St_phase','St','fmax','udmax','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        save(fullfile(pathname,'solution.mat'),'Ht','-append');
    end
else
    load(fullfile(pathname,'solution.mat'),'dt','ut','ft','T','St_phase','St','fmax','udmax','time');
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
    plot(t*1e3,ft*((Dim==2)*1e-6+(Dim==3)*1e-3),'-b','Linewidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    %% Display solutions at different instants
    ampl = 0;
    rep = find(abs(t-18.5e-6)<eps | abs(t-24.6e-6)<eps);
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
