%% Monoscale deterministic linear elasticity problem with curve crack %%
%%--------------------------------------------------------------------%%
% [Bourdin, Francfort, Marigo, 2000, JMPS]
% [Amor, Marigo, Maurini, 2009, JMPS]
% [Miehe, Hofacker, Welschinger, 2010, CMAME]
% [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
% [Wu, Nguyen, Nguyen, Sutula, Borad, Sinaie, 2018, AAM]


% clc
clearvars
close all
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

loading = 'Shear'; % 'Pull' or 'Shear'
filename = ['linElasCurvedCracks' loading '_bis'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','monoscaleDet',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    L = 1e-3;
    a = L/2;
    D = DOMAIN(2,[0.0,0.0],[L,L]);
    
    P = [a,L/2];
    B = LIGNE([0.0,L/2],P);
    
    option = 'DEFO'; % plane stress
    clD = 2e-5;
    % clP = 6e-7;
    clP = 5e-6;
    S_phase = gmshdomainwithcurvedcrack(D,P,clD,clP,fullfile(pathname,'gmsh_domain_curved_crack'));
    
%     load(fullfile(getfemobjectoptions('path'),'MYCODE',...
%         'examples','monoscale','micro2Dtest_coarse.mat'));
%     node = NODE([Nx,1-Ny]*1e-3,1:size(Nx,1));
%     elem = Connect;
%     elemtype = 'TRI3';
%     S_phase = MODEL('PLAN');
%     S_phase = addnode(S_phase,node);
%     S_phase = addelem(S_phase,elemtype,elem,'option',option)
    
    S_phase = concatgroupelem(S_phase);
    S = setoption(S_phase,option);
    
    %% Phase field problem
    %% Material
    % Fracture toughness
    gc = 2700;
    % Regularization parameter (width of the smeared crack)
    l = 7.5e-6;
    % Small parameter
    k = 1e-10;
    % Internal energy
    H = 0;
    
    % Material
    mat_phase = FOUR_ISOT('k',gc*l,'r',gc/l+2*H);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    S_phase = final(S_phase,'duplicate');
    S_phase = addcl(S_phase,B,'T',1);
    
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
    % Lame coefficients
    lambda = 121.15e9;
    mu = 80.77e9;
    % Young modulus and Poisson ratio
    switch lower(option)
        case 'defo'
            E = mu*(3*lambda+2*mu)/(lambda+mu);
            NU = lambda/(lambda+mu)/2;
        case 'cont'
            E = 4*mu*(lambda+mu)/(lambda+2*mu);
            NU = lambda/(lambda+2*mu);
    end
    % Degradation/Damage function
    g = @(d) (1-d).^2;
    % Thickness
    DIM3 = 1;
    % Density
    RHO = 1;
    
    % Material
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
    mat = setnumber(mat,1);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    LU = LIGNE([0.0,L],[L,L]);
    LL = LIGNE([0.0,0.0],[L,0.0]);
    LRight = LIGNE([L,0.0],[L,L]);
    LLeft = LIGNE([0.0,0.0],[0.0,L]);
    
    S = final(S,'duplicate');
    
    ud = 0;
    switch lower(loading)
        case 'pull'
            S = addcl(S,LU,'UY',ud);
        case 'shear'
            S = addcl(S,LU,{'UX','UY'},[ud;0]);
            S = addcl(S,LRight,'UY');
            S = addcl(S,LLeft,'UY');
        otherwise
            error('Wrong loading case')
    end
    S = addcl(S,LL);
    
    %% Stiffness matrices and sollicitation vectors
    % [A,b] = calc_rigi(S);
    % b = -b;
    
    %% Time scheme
    dt = 2e-8;
    nt = 1500;
    t0 = 2e-8;
    t1 = nt*dt;
    T = TIMEMODEL(t0,t1,nt-1);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','B','LU','LL','LRight','LLeft','gc','l','E','g');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','B','LU','LL','LRight','LLeft','gc','l','E','g');
end

%% Solution
if solveProblem
    
    tTotal = tic;
    
    t = gett(T);
    nt = getnt(T);
    
    Ht = cell(1,length(T));
    dt = cell(1,length(T));
    ut = cell(1,length(T));
    
    sz_H = getnbelem(S);
    sz_d = getnbddl(S_phase);
    sz_u = getnbddl(S);
    H = FEELEMFIELD(zeros(sz_H,1),S);
    d = zeros(sz_d,1);
    u = zeros(sz_u,1);
    
    fprintf('\n+------------+------------+------------+------------+\n');
    fprintf('|    Iter    |  norm(H)   |  norm(d)   |  norm(u)   |\n');
    fprintf('+------------+------------+------------+------------+\n');
    
    for i=1:length(T)
        
        % Internal energy field
        mats = MATERIALS(S);
        for m=1:length(mats)
            mats{m} = setparam(mats{m},'E',E);
        end
        S = actualisematerials(S,mats);
        
        h_old = getvalue(H);
        H = calc_energyint(S,u);
        h = getvalue(H);
        for p=1:getnbgroupelem(S)
            he = double(h{p});
            he_old = double(h_old{p});
            rep = find(he <= he_old);
            he(rep) = he_old(rep);
            h{p} = he;
        end
        H = FEELEMFIELD(h,'storage',getstorage(H),'type',gettype(H),'ddl',getddl(H));
        
        % Phase field
        mats_phase = MATERIALS(S_phase);
        for m=1:length(mats_phase)
            mats_phase{m} = setparam(mats_phase{m},'r',gc/l+2*H);
        end
        S_phase = actualisematerials(S_phase,mats_phase);
        
        [A_phase,b_phase] = calc_rigi(S_phase);
        b_phase = -b_phase + bodyload(S_phase,[],'QN',2*H);
        
        d = A_phase\b_phase;
        d = unfreevector(S_phase,d);
        
        % Displacement field
        mats = MATERIALS(S);
        for m=1:length(mats)
            mats{m} = setparam(mats{m},'E',FENODEFIELD(E.*(g(d)+k)));
        end
        S = actualisematerials(S,mats);
        
        S = removebc(S);
        ud = t(i);
        switch lower(loading)
            case 'pull'
                S = addcl(S,LU,'UY',ud);
            case 'shear'
                S = addcl(S,LU,{'UX','UY'},[ud;0]);
                S = addcl(S,LRight,'UY');
                S = addcl(S,LLeft,'UY');
            otherwise
                error('Wrong loading case')
        end
        S = addcl(S,LL);
        
        [A,b] = calc_rigi(S);
        b = -b;
        
        u = A\b;
        u = unfreevector(S,u);
        
        % Update fields
        Ht{i} = H;
        dt{i} = d;
        ut{i} = u;
        
        fprintf('| %10d | %9.4e | %9.4e | %9.4e |\n',i,norm(squeeze(double(Ht{i}))),norm(dt{i}),norm(ut{i}));
        
    end
    
    fprintf('+------------+------------+------------+------------+\n');
    
    Ht = TIMEMATRIX(cellfun(@(H) squeeze(double(H)),Ht,'UniformOutput',false),T,[sz_H,1]);
    dt = TIMEMATRIX(dt,T,[sz_d,1]);
    ut = TIMEMATRIX(ut,T,[sz_u,1]);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'Ht','dt','ut','time');
else
    load(fullfile(pathname,'solution.mat'),'Ht','dt','ut','time');
end

%% Outputs
fprintf('\n');
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('nb time dofs = %g\n',getnbtimedof(T));
fprintf('elapsed time = %f s\n',time);

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain({D,B},'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 0.5;
    v = calc_init_dirichlet(S);
    [hN,legN] = vectorplot(S,'U',v,ampl,'r','LineWidth',1);
    % legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    % plotModel(S,'legend',false);
    % mysaveas(pathname,'mesh',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    [t,rep] = gettevol(T);
    u = getmatrixatstep(ut,rep(end));
    ampl = getsize(S)/max(abs(u))/20;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*unfreevector(S,u),'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    %% Display evolution of solutions
    ampl = 0;
    % ampl = getsize(S)/max(max(abs(getvalue(ut))))/20;
    
    options = {'plotiter',true,'plottime',false};
    
%     figure('Name','Solution H')
%     clf
%     T = setevolparam(T,'colorbar',true,'FontSize',fontsize,options{:});
%     frame = evol(T,Ht,S_phase,'rescale',true);
%     saveMovie(frame,'filename','internal_energy','pathname',pathname);
    
    evolSolution(S_phase,dt,'filename',['damage_' num2str(i)],'pathname',pathname,options{:});
    for i=1:2
        evolSolution(S,ut,'displ',i,'ampl',ampl,'filename',['displacement_' num2str(i)],'pathname',pathname,options{:});
    end
    
%     for i=1:3
%         evolSolution(S,ut,'epsilon',i,'ampl',ampl,'filename',['epsilon_' num2str(i)],'pathname',pathname,options{:});
%         evolSolution(S,ut,'sigma',i,'ampl',ampl,'filename',['sigma_' num2str(i)],'pathname',pathname,options{:});
%     end
%     
%     evolSolution(S,ut,'epsilon','mises','ampl',ampl,'filename','epsilon_von_mises','pathname',pathname,options{:});
%     evolSolution(S,ut,'sigma','mises','ampl',ampl,'filename','sigma_von_mises','pathname',pathname,options{:});
    
    %% Display solutions at differents instants
    rep = [500,650];
    for j=1:length(rep)
        close all
        Hj = getmatrixatstep(Ht,rep(j));
        dj = getmatrixatstep(dt,rep(j));
        uj = getmatrixatstep(ut,rep(j));
        
%         figure('Name','Solution H')
%         clf
%         plot(S_phase,Hj);
%         colorbar
%         set(gca,'FontSize',fontsize)
%         mysaveas(pathname,['internal_energy_t' num2str(rep(j))],formats,renderer);
        
        plotSolution(S_phase,dj);
        mysaveas(pathname,['damage_t' num2str(rep(j))],formats,renderer);
        
        for i=1:2
            plotSolution(S,uj,'displ',i,'ampl',ampl);
            mysaveas(pathname,['displacement_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        end
        
%         for i=1:3
%             plotSolution(S,uj,'epsilon',i,'ampl',ampl);
%             mysaveas(pathname,['epsilon_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
%             
%             plotSolution(S,uj,'sigma',i,'ampl',ampl);
%             mysaveas(pathname,['sigma_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
%         end
%         
%         plotSolution(S,uj,'epsilon','mises','ampl',ampl);
%         mysaveas(pathname,['epsilon_von_mises_t' num2str(rep(j))],formats,renderer);
%         
%         plotSolution(S,uj,'sigma','mises','ampl',ampl);
%         mysaveas(pathname,['sigma_von_mises_t' num2str(rep(j))],formats,renderer);
    end
    
end

% myparallel('stop');
