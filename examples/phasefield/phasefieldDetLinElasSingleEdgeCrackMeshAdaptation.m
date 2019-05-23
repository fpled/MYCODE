%% Phase field fracture model - deterministic linear elasticity problem with single edge crack %%
%%---------------------------------------------------------------------------------------------%%
% [Bourdin, Francfort, Marigo, 2000, JMPS]
% [Amor, Marigo, Maurini, 2009, JMPS]
% [Miehe, Hofacker, Welschinger, 2010, CMAME]
% [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME]
% [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
% [Wu, Nguyen, Nguyen, Sutula, Borad, Sinaie, 2018, AAM]
% [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]


% clc
clearvars
close all
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

Dim = 2; % space dimension Dim = 2, 3
loading = 'Shear'; % 'Pull' or 'Shear'
filename = ['phasefieldDetLinElasSingleEdgeCrack' loading '_' num2str(Dim) 'D_MeshAdaptation'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','phasefield',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
formats = {'fig','epsc'};
renderer = 'OpenGL';

gmshoptions = '-v 0';
mmgoptions = '-v 0';
% gmshoptions = '-v 5';
% mmgoptions = '-v 1';

%% Problem
if setProblem
    %% Domains and meshes
    L = 1e-3;
    a = L/2;
    if Dim==2
        D = DOMAIN(2,[0.0,0.0],[L,L]);
        C = LIGNE([0.0,L/2],[a,L/2]);
    elseif Dim==3
        D = DOMAIN(3,[0.0,0.0,0.0],[L,L,L]);
        C = QUADRANGLE([0.0,0.0,L/2],[a,0.0,L/2],[a,L,L/2],[0.0,L,L/2]);
    end
    
    if Dim==2
        clD = 2e-5;
        % clC = 6e-7;
        clC = 2e-6;
    elseif Dim==3
        clD = 5e-5;
        clC = 5e-6;
    end
    % S_phase = gmshdomainwithedgecrack(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,'gmshoptions',gmshoptions);
    S_phase = gmshdomainwithedgesmearedcrack(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,'gmshoptions',gmshoptions);
    
    if Dim==2
        CU = LIGNE([0.0,L/2+clC/2],[a,L/2]);
        CL = LIGNE([0.0,L/2-clC/2],[a,L/2]);
    elseif Dim==3
        CU = QUADRANGLE([0.0,0.0,L/2+clC/2],[a,0.0,L/2],[a,L,L/2],[0.0,L,L/2+clC/2]);
        CL = QUADRANGLE([0.0,0.0,L/2-clC/2],[a,0.0,L/2],[a,L,L/2],[0.0,L,L/2-clC/2]);
    end
    
    % sizemap = @(d) (clC-clD)*d+clD;
    sizemap = @(d) clD*clC./((clD-clC)*d+clC);
    
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
    % S_phase = addcl(S_phase,C,'T',1);
    S_phase = addcl(S_phase,CU,'T',1);
    S_phase = addcl(S_phase,CL,'T',1);
    
    d = calc_init_dirichlet(S_phase);
    cl = sizemap(d);
    S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_domain_single_edge_crack'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    S = S_phase;
    
    S_phase = setmaterial(S_phase,mat_phase);
    S_phase = final(S_phase,'duplicate');
    % S_phase = addcl(S_phase,C,'T',1);
    S_phase = addcl(S_phase,CU,'T',1);
    S_phase = addcl(S_phase,CL,'T',1);
    
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
    % Energetic degradation function
    g = @(d) (1-d).^2;
    % Thickness
    DIM3 = 1;
    % Density
    RHO = 1;
    
    % Material
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
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
        BU = PLAN([0.0,0.0,L],[L,0.0,L],[0.0,L,L]);
        BL = PLAN([0.0,0.0,0.0],[L,0.0,0.0],[0.0,L,0.0]);
        BRight = PLAN([L,0.0,0.0],[L,L,0.0],[L,0.0,L]);
        BLeft = PLAN([0.0,0.0,0.0],[0.0,L,0.0],[0.0,0.0,L]);
        BFront = PLAN([0.0,0.0,0.0],[L,0.0,0.0],[0.0,0.0,L]);
        BBack = PLAN([0.0,L,0.0],[L,L,0.0],[0.0,L,L]);
    end
    
    S = final(S,'duplicate');
    
    ud = 0;
    switch lower(loading)
        case 'pull'
            S = addcl(S,BU,'UY',ud);
        case 'shear'
            if Dim==2
                S = addcl(S,BU,{'UX','UY'},[ud;0]);
                S = addcl(S,BRight,'UY');
                S = addcl(S,BLeft,'UY');
            elseif Dim==3
                S = addcl(S,BU,{'UX','UY','UZ'},[ud;0;0]);
                S = addcl(S,BRight,{'UY','UZ'});
                S = addcl(S,BLeft,{'UY','UZ'});
                S = addcl(S,BFront,{'UY','UZ'});
                S = addcl(S,BBack,{'UY','UZ'});
            end
        otherwise
            error('Wrong loading case')
    end
    S = addcl(S,BL);
    
    %% Stiffness matrices and sollicitation vectors
    % [A,b] = calc_rigi(S);
    % b = -b;
    
    %% Time scheme
    if Dim==2
        dt = 2e-8;
        nt = 1500;
    elseif Dim==3
        dt = 2e-8;
        nt = 2500;
    end
    t0 = dt;
    t1 = nt*dt;
    T = TIMEMODEL(t0,t1,nt-1);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','D','C','CU','CL','BU','BL','BRight','BLeft','BFront','BBack','gc','l','E','g');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','D','C','CU','CL','BU','BL','BRight','BLeft','BFront','BBack','gc','l','E','g');
end

%% Solution
if solveProblem
    
    tTotal = tic;
    
    t = gett(T);
    
    Ht = cell(1,length(T));
    dt = cell(1,length(T));
    ut = cell(1,length(T));
    St_phase = cell(1,length(T));
    St = cell(1,length(T));
    
    sz_phase = getnbddl(S_phase);
    sz = getnbddl(S);
    H = zeros(sz_phase,1);
    u = zeros(sz,1);
    
    fprintf('\n+----------+----------+----------+------------+------------+------------+\n');
    fprintf('|   Iter   | Nb nodes | Nb elems |  norm(H)   |  norm(d)   |  norm(u)   |\n');
    fprintf('+----------+----------+----------+------------+------------+------------+\n');
    fprintf('| %8d | %8d | %8d | %9.4e | %9.4e | %9.4e |\n',0,getnbnode(S),getnbelem(S),0,0,0);
    
    for i=1:length(T)
        
        % Internal energy field
        mats = MATERIALS(S);
        for m=1:length(mats)
            mats{m} = setparam(mats{m},'E',E);
        end
        S = actualisematerials(S,mats);
        
        h_old = double(H);
        H = FENODEFIELD(calc_energyint(S,u,'node'));
        h = double(H);
        rep = find(h <= h_old);
        h(rep) = h_old(rep);
        H = setvalue(H,h);
        
        % Phase field
        mats_phase = MATERIALS(S_phase);
        for m=1:length(mats_phase)
            mats_phase{m} = setparam(mats_phase{m},'r',FENODEFIELD(gc/l+2*H));
        end
        S_phase = actualisematerials(S_phase,mats_phase);
        
        [A_phase,b_phase] = calc_rigi(S_phase);
        b_phase = -b_phase + bodyload(S_phase,[],'QN',FENODEFIELD(2*H));
        
        d = A_phase\b_phase;
        d = unfreevector(S_phase,d);
        
        % Mesh adaptation
        S_phase_old = S_phase;
        cl = sizemap(d);
        S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_domain_single_edge_crack'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
        S = S_phase;
        
        for m=1:length(mats_phase)
            S_phase = setmaterial(S_phase,mats_phase{m},m);
        end
        % S_phase = actualisematerials(S_phase,mats_phase);
        S_phase = final(S_phase,'duplicate');
        S_phase = addcl(S_phase,CU,'T',1);
        S_phase = addcl(S_phase,CL,'T',1);
        
        P_phase = calcProjection(S_phase,S_phase_old,[],'free',false);
        d = P_phase'*d;
        h = P_phase'*h;
        H = setvalue(H,h);
        
        % Displacement field
        for m=1:length(mats)
            mats{m} = setparam(mats{m},'E',FENODEFIELD(E.*(g(d)+k)));
            S = setmaterial(S,mats{m},m);
        end
        % S = actualisematerials(S,mats);
        S = final(S,'duplicate');
        S = removebc(S);
        ud = t(i);
        switch lower(loading)
            case 'pull'
                S = addcl(S,BU,'UY',ud);
            case 'shear'
                if Dim==2
                    S = addcl(S,BU,{'UX','UY'},[ud;0]);
                    S = addcl(S,BRight,'UY');
                    S = addcl(S,BLeft,'UY');
                elseif Dim==3
                    S = addcl(S,BU,{'UX','UY','UZ'},[ud;0;0]);
                    S = addcl(S,BRight,{'UY','UZ'});
                    S = addcl(S,BLeft,{'UY','UZ'});
                    S = addcl(S,BRight,{'UY','UZ'});
                    S = addcl(S,BLeft,{'UY','UZ'});
                end
            otherwise
                error('Wrong loading case')
        end
        S = addcl(S,BL);
        
        [A,b] = calc_rigi(S);
        b = -b;
        
        u = A\b;
        u = unfreevector(S,u);
        
        % Update fields
        Ht{i} = double(H);
        dt{i} = d;
        ut{i} = u;
        St_phase{i} = S_phase;
        St{i} = S;
        
        fprintf('| %8d | %8d | %8d | %9.4e | %9.4e | %9.4e |\n',i,getnbnode(S),getnbelem(S),norm(Ht{i}),norm(dt{i}),norm(ut{i}));
        
    end
    
    fprintf('+----------+----------+----------+------------+------------+------------+\n');
    
    % DO NOT WORK WITH MESH ADAPTATION
    % Ht = TIMEMATRIX(Ht,T,[sz_phase,1]);
    % dt = TIMEMATRIX(dt,T,[sz_phase,1]);
    % ut = TIMEMATRIX(ut,T,[sz,1]);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'Ht','dt','ut','St','St_phase','time');
else
    load(fullfile(pathname,'solution.mat'),'Ht','dt','ut','St','St_phase','time');
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
    % legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions_damage',formats,renderer);
    
    % plotModel(S,'legend',false);
    % mysaveas(pathname,'mesh',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    [t,rep] = gettevol(T);
    % DO NOT WORK WITH MESH ADAPTATION
    % u = getmatrixatstep(ut,rep(end));
    u = ut{rep(end)};
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
    % DO NOT WORK WITH MESH ADAPTATION
    % ampl = getsize(S)/max(max(abs(getvalue(ut))))/20;
    % umax = cellfun(@(u) max(abs(u)),ut,'UniformOutput',false);
    % ampl = getsize(S)/max([umax{:}])/20;
    
    options = {'plotiter',true,'plottime',false};
    
    evolModel(T,St,'filename','mesh','pathname',pathname,options{:});
    
%     evolSolutionCell(T,St_phase,Ht,'filename','internal_energy','pathname',pathname,options{:});
    
    evolSolutionCell(T,St_phase,dt,'filename','damage','pathname',pathname,options{:});
    for i=1:Dim
        evolSolutionCell(T,St,ut,'displ',i,'ampl',ampl,'filename',['displacement_' num2str(i)],'pathname',pathname,options{:});
    end
    
%     for i=1:(Dim*(Dim+1)/2)
%         evolSolutionCell(T,St,ut,'epsilon',i,'ampl',ampl,'filename',['epsilon_' num2str(i)],'pathname',pathname,options{:});
%         evolSolutionCell(T,St,ut,'sigma',i,'ampl',ampl,'filename',['sigma_' num2str(i)],'pathname',pathname,options{:});
%     end
%     
%     evolSolutionCell(T,St,ut,'epsilon','mises','ampl',ampl,'filename','epsilon_von_mises','pathname',pathname,options{:});
%     evolSolutionCell(T,St,ut,'sigma','mises','ampl',ampl,'filename','sigma_von_mises','pathname',pathname,options{:});
    
    %% Display solutions at differents instants
    rep = [500,650];
    for j=1:length(rep)
        close all
        % DO NOT WORK WITH MESH ADAPTATION
        % Hj = getmatrixatstep(Ht,rep(j));
        % dj = getmatrixatstep(dt,rep(j));
        % uj = getmatrixatstep(ut,rep(j));
        Hj = Ht{rep(j)};
        dj = dt{rep(j)};
        uj = ut{rep(j)};
        Sj = St{rep(j)};
        Sj_phase = St_phase{rep(j)};
        
%         plotSolution(Sj_phase,Hj);
%         mysaveas(pathname,['internal_energy_t' num2str(rep(j))],formats,renderer);
        
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
    end
    
end

%% Save solutions
[t,rep] = gettevol(T);
for i=1:length(T)
    Hi = getmatrixatstep(Ht,rep(i));
    di = getmatrixatstep(dt,rep(i));
    ui = getmatrixatstep(ut,rep(i));
    
    write_vtk_mesh(S,{Hi,di,ui},[],...
        {'internal energy','damage','displacement'},[],...
        pathname,'solution',1,i-1);
end
make_pvd_file(pathname,'solution',1,length(T));

% myparallel('stop');
