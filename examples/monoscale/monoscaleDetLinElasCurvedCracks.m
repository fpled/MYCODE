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
filename = ['linElasCurvedCracks' loading];
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
    % clD = 2e-5;
    % clP = 6e-7;
    % clD = 2e-5;
    % clP = 2e-5;
    clD = 2e-5;
    clP = 2e-5;
    % S_phase = gmshdomainwithcurvedcrackmodif(D,P,clD,clP,fullfile(pathname,'gmsh_domain_curved_crack'));
    load('./micro2Dtest_coarse.mat');
    node = NODE([Nx,1-Ny]*1e-3,1:size(Nx,1));
    elem = Connect;
    elemtype = 'TRI3';
    S_phase = MODEL('PLAN');
    S_phase = addnode(S_phase,node);
    S_phase = addelem(S_phase,elemtype,elem,'option',option);
    
    S_phase = concatgroupelem(S_phase);
    S = setoption(S_phase,option);
    
    %% Initialization
    H = zeros(S_phase.nbnode,1);
    % H = 0;
    
    %% Phase field problem
    %% Material
    % Fracture toughness
    gc = 2700;
    % Regularization parameter (width of the smeared crack)
    l = 7.5e-6;
    % l = 1.5e-5;
    % Small parameter
    k = 1e-10;
    
    % Material
    mat_phase = FOUR_ISOT('k',gc*l,'r',FENODEFIELD(gc/l+2*H));
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    S_phase = final(S_phase,'duplicate');
    % S_phase = addcl(S_phase,B,'T',1);
    
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
    
    % l_phase = LINFORM(0,2*H,0);
    % l_phase = setfree(l_phase,1);
    % b_phase = b_phase + calc_vector(l_phase,S_phase);
    
    [A_phase,b_phase] = calc_rigi(S_phase);
    b_phase = -b_phase + bodyload(S_phase,[],'QN',FENODEFIELD(2*H)); 
    
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
    E = 210e9;
    NU = 0.3;
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
    [A,b] = calc_rigi(S);
    b = -b;
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S_phase','S','D','B','A','b');
else
    load(fullfile(pathname,'problem.mat'),'S_phase','S','D','B','A','b');
end

%% Solution
if solveProblem
    
    for n=1:1000
        t = tic;
        % Phase field solution
%         a_phase = BILINFORM(1,1,gc*l); % uniform values
%         % a_phase = DIFFUSIONFORM(gc*l);
%         a_phase = setfree(a_phase,0);
%         K_phase = calc_matrix(a_phase,S_phase);
%         b_phase = calc_nonhomogeneous_vector(S_phase,K_phase);
%         b_phase = -b_phase;
%         K_phase = freematrix(S_phase,K_phase);
%         
%         r_phase = BILINFORM(0,0,gc/l+2*H,0); % nodal values
%         R_phase = calc_matrix(r_phase,S_phase);
%         A_phase = K_phase + R_phase;
%         
%         l_phase = LINFORM(0,2*H,0);
%         l_phase = setfree(l_phase,1);
%         b_phase = b_phase + calc_vector(l_phase,S_phase);
        
        mat_phase = setparam(mat_phase,'r',FENODEFIELD(gc/l+2*H));
        S_phase = setmaterial(S_phase,mat_phase);
        
        [A_phase,b_phase] = calc_rigi(S_phase);
        b_phase = -b_phase + bodyload(S_phase,[],'QN',FENODEFIELD(2*H));
        
        d = A_phase\b_phase;
        
        % Displacement field solution
        d = unfreevector(S_phase,d);
        mat = setparam(mat,'E',FENODEFIELD(E.*(g(d)+k)));
        S = setmaterial(S,mat);
        
        S = removebc(S);
        ud = ud + 2e-8;
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
        
        e = calc_epsilon(S,u,'node');
        s = calc_sigma(S,u,'node');
        H_old = H;
        H = squeeze(1/2*sum(double(s.*e),1));
        rep = find(H <= H_old);
        H(rep) = H_old(rep);
        
        if ~mod(n,10)
            plotSolution(S_phase,d);
            % caxis([0,1])
%             figure('Name','Solution u')
%             clf
%             subplot(1,2,1)
%             plot_sol(S,u,'displ',1)
%             colorbar
%             subplot(1,2,2)
%             plot_sol(S,u,'displ',2)
%             colorbar
%             
%             figure('Name','Solution sigma')
%             clf
%             subplot(1,3,1)
%             plot(s,S,'compo','SMXX')
%             colorbar
%             subplot(1,3,2)
%             plot(s,S,'compo','SMYY')
%             colorbar
%             subplot(1,3,3)
%             plot(s,S,'compo','SMXY')
%             colorbar
            
%             plotSolution(S_phase,H);
            
            pause(1)
        end
        
        time = toc(t);
        
        fprintf('n = %d, max(d) = %g, norm(u) = %g\n',n,max(d),norm(u));
        
    end
    
    plotSolution(S_phase,d);
    
    save(fullfile(pathname,'solution.mat'),'u','e','s','time');
else
    load(fullfile(pathname,'solution.mat'),'u','e','s','time');
end

%% Outputs
fprintf('\n');
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
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
    
    ampl = getsize(S)/max(abs(u))/20;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*unfreevector(S,u),'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    %% Display solutions
    % ampl = 0;
    ampl = getsize(S)/max(abs(u))/20;
    
    for i=1:2
        plotSolution(S,u,'displ',i,'ampl',ampl);
        mysaveas(pathname,['u_' num2str(i)],formats,renderer);
    end
    
    for i=1:3
        plotSolution(S,u,'epsilon',i,'ampl',ampl);
        mysaveas(pathname,['eps_' num2str(i)],formats,renderer);
        
        plotSolution(S,u,'sigma',i,'ampl',ampl);
        mysaveas(pathname,['sig_' num2str(i)],formats,renderer);
    end
    
    plotSolution(S,u,'epsilon','mises','ampl',ampl);
    mysaveas(pathname,'eps_von_mises',formats,renderer);
    
    plotSolution(S,u,'sigma','mises','ampl',ampl);
    mysaveas(pathname,'sig_von_mises',formats,renderer);
    
    % u = unfreevector(S,u);
    %
    % figure('Name','Solution eps_xx')
    % clf
    % plot(e,S+ampl*u,'compo','EPXX')
    % colorbar
    % set(gca,'FontSize',fontsize)
    % mysaveas(pathname,'eps_xx',formats,renderer);
    %
    % figure('Name','Solution eps_yy')
    % clf
    % plot(e,S+ampl*u,'compo','EPYY')
    % colorbar
    % set(gca,'FontSize',fontsize)
    % mysaveas(pathname,'eps_yy',formats,renderer);
    %
    % figure('Name','Solution eps_xy')
    % clf
    % plot(e,S+ampl*u,'compo','EPXY')
    % colorbar
    % set(gca,'FontSize',fontsize)
    % mysaveas(pathname,'eps_xy',formats,renderer);
    
    % figure('Name','Solution sig_xx')
    % clf
    % plot(s,S+ampl*u,'compo','SMXX')
    % colorbar
    % set(gca,'FontSize',fontsize)
    % mysaveas(pathname,'sig_xx',formats,renderer);
    %
    % figure('Name','Solution sig_yy')
    % clf
    % plot(s,S+ampl*u,'compo','SMYY')
    % colorbar
    % set(gca,'FontSize',fontsize)
    % mysaveas(pathname,'sig_yy',formats,renderer);
    %
    % figure('Name','Solution sig_xy')
    % clf
    % plot(s,S+ampl*u,'compo','SMXY')
    % colorbar
    % set(gca,'FontSize',fontsize)
    % mysaveas(pathname,'sig_xy',formats,renderer);
end

% myparallel('stop');
