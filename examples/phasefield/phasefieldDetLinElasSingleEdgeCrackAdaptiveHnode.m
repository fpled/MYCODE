%% Phase field fracture model - deterministic linear elasticity problem with single edge crack %%
%%---------------------------------------------------------------------------------------------%%
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
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = false;

test = true; % coarse mesh
% test = false; % fine mesh

Dim = 2; % space dimension Dim = 2, 3
loading = 'Shear'; % 'Tension' or 'Shear'
PFmodel = 'Isotropic'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'

filename = ['phasefieldDetLinElasSingleEdgeCrack' loading PFmodel 'AdaptiveHnode_' num2str(Dim) 'D'];
if test
    filename = [filename '_test'];
end
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
        % clC = 3.906e-6; % [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME]
        % clC = 2.5e-6; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        clC = 2e-6; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
        % clC = 1e-6; % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME]
        % clC = 6e-7; % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
        % clD = 3.9e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
        % clC = 3.9e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
        % clD = 2e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        % clC = 2e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        % clD = 1e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
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
    % Critical energy release rate (or fracture toughness)
    gc = 2.7e3;
    % Regularization parameter (width of the smeared crack)
    % l = 3.75e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME]
    % l = 3e-5; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % l = 1.5e-5; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
    % l = 1e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
    l = 7.5e-6; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
    % l = 5e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
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
    % Lame coefficients
    % lambda = 121.1538e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
    % mu = 80.7692e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
    lambda = 121.15e9;
    mu = 80.77e9;
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
                    nt1 = 300;
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
                % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM], 
                % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
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
    T = TIMEMODEL(t);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','D','C','BU','BL','BRight','BLeft','BFront','BBack','loading');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','D','C','BU','BL','BRight','BLeft','BFront','BBack','loading');
end

%% Solution
if solveProblem
    tTotal = tic;
    
    [dt,ut,ft,St_phase,St] = solvePFDetLinElasSingleEdgeCrackAdaptiveHnode(S_phase,S,T,C,BU,BL,BRight,BLeft,BFront,BBack,loading,sizemap,...
        'filename','gmsh_domain_single_edge_crack','pathname',pathname,'gmshoptions',gmshoptions,'mmgoptions',mmgoptions,'display');
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'dt','ut','ft','St_phase','St','time');
else
    load(fullfile(pathname,'solution.mat'),'dt','ut','ft','St_phase','St','time');
end

%% Outputs
fprintf('\n');
fprintf('dim      = %d\n',Dim);
fprintf('loading  = %s\n',loading);
fprintf('PF model = %s\n',PFmodel);
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('nb time dofs = %g\n',getnbtimedof(T));
fprintf('elapsed time = %f s\n',time);

%% Display
if displaySolution
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
    
    % DO NOT WORK WITH MESH ADAPTATION
    % u = getmatrixatstep(ut,rep(end));
%     u = ut{rep(end)};
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
    
    if Dim~=3
    %% Display evolution of solutions
    ampl = 0;
    % DO NOT WORK WITH MESH ADAPTATION
    % ampl = getsize(S)/max(max(abs(getvalue(ut))))/20;
    % umax = cellfun(@(u) max(abs(u)),ut,'UniformOutput',false);
    % ampl = getsize(S)/max([umax{:}])/20;
    
    options = {'plotiter',true,'plottime',false};
    framerate = 80;
    
    evolModel(T,St,'FrameRate',framerate,'filename','mesh','pathname',pathname,options{:});
    
    evolSolutionCell(T,St_phase,dt,'FrameRate',framerate,'filename','damage','pathname',pathname,options{:});
    for i=1:Dim
        evolSolutionCell(T,St,ut,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i)],'pathname',pathname,options{:});
    end
    
%     for i=1:(Dim*(Dim+1)/2)
%         evolSolutionCell(T,St,ut,'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['epsilon_' num2str(i)],'pathname',pathname,options{:});
%         evolSolutionCell(T,St,ut,'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['sigma_' num2str(i)],'pathname',pathname,options{:});
%     end
%     
%     evolSolutionCell(T,St,ut,'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename','epsilon_von_mises','pathname',pathname,options{:});
%     evolSolutionCell(T,St,ut,'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename','sigma_von_mises','pathname',pathname,options{:});
%     evolSolutionCell(T,St,ut,'energyint','','ampl',ampl,'FrameRate',framerate,'filename','internal_energy','pathname',pathname,options{:});
    
    %% Display solutions at different instants
    switch lower(loading)
        case 'tension'
            rep = find(abs(t-5.5e-6)<eps | abs(t-5.75e-5)<eps | abs(t-6e-6)<eps | abs(t-6.25e-6)<eps);
        case 'shear'
            rep = find(abs(t-1e-5)<eps | abs(t-1.25e-5)<eps | abs(t-1.35e-5)<eps | abs(t-1.5e-5)<eps);
        otherwise
            error('Wrong loading case')
    end
    rep = [rep,length(T)];
    for j=1:length(rep)
        % DO NOT WORK WITH MESH ADAPTATION
        % dj = getmatrixatstep(dt,rep(j));
        % uj = getmatrixatstep(ut,rep(j));
        dj = dt{rep(j)};
        uj = ut{rep(j)};
        Sj = St{rep(j)};
        Sj_phase = St_phase{rep(j)};
        
        plotModel(Sj,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
        mysaveas(pathname,['mesh_t' num2str(rep(j))],formats,renderer);
        
        plotSolution(Sj_phase,dj);
        mysaveas(pathname,['damage_t' num2str(rep(j))],formats,renderer);
        
        for i=1:Dim
            plotSolution(Sj,uj,'displ',i,'ampl',ampl);
            mysaveas(pathname,['displacement_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        end
        
%         for i=1:(Dim*(Dim+1)/2)
%             plotSolution(Sj,uj,'epsilon',i,'ampl',ampl);
%             mysaveas(pathname,['epsilon_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
%             
%             plotSolution(Sj,uj,'sigma',i,'ampl',ampl);
%             mysaveas(pathname,['sigma_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
%         end
%         
%         plotSolution(Sj,uj,'epsilon','mises','ampl',ampl);
%         mysaveas(pathname,['epsilon_von_mises_t' num2str(rep(j))],formats,renderer);
%         
%         plotSolution(Sj,uj,'sigma','mises','ampl',ampl);
%         mysaveas(pathname,['sigma_von_mises_t' num2str(rep(j))],formats,renderer);
%         
%         plotSolution(Sj,uj,'energyint','','ampl',ampl);
%         mysaveas(pathname,['internal_energy_t' num2str(rep(j))],formats,renderer);
    end
    end
    
end

%% Save solutions
[t,rep] = gettevol(T);
for i=1:length(T)
    % DO NOT WORK WITH MESH ADAPTATION
    % di = getmatrixatstep(dt,rep(i));
    % ui = getmatrixatstep(ut,rep(i));
    di = dt{rep(i)};
    ui = ut{rep(i)};
    Si = St{rep(i)};
    % Si_phase = St_phase{rep(i)};
    
    write_vtk_mesh(Si,{di,ui},[],...
        {'damage','displacement'},[],...
        pathname,'solution',1,i-1);
end
make_pvd_file(pathname,'solution',1,length(T));

% myparallel('stop');