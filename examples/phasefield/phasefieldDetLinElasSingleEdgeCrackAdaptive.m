%% Phase field fracture model - deterministic linear elasticity problem with single edge crack %%
%%---------------------------------------------------------------------------------------------%%
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
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = false;

Dim = 3; % space dimension Dim = 2, 3
loading = 'Shear'; % 'Tension' or 'Shear'
PFmodel = 'AnisotropicMiehe'; % 'Isotropic', 'AnisotropicAmor' or 'AnisotropicMiehe'
filename = ['phasefieldDetLinElasSingleEdgeCrack' loading PFmodel 'Adaptive_' num2str(Dim) 'D'];
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
        clD = 2e-5; % [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
        clC = 2e-6; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
        % clC = 1e-6; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
        % clC = 6e-7; % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
        % clD = 3.9e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM]
        % clC = 3.9e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM]
        clD = 1.5e-5; % test
        clC = 1.5e-5; % test
%         clD = 4e-5; % test
%         clC = 1e-5; % test
    elseif Dim==3
        clD = 4e-5;
        clC = 4e-6;
        clD = 4e-5; % test
        clC = 1e-5; % test
    end
    % S_phase = gmshdomainwithedgecrack(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,'gmshoptions',gmshoptions);
    S_phase = gmshdomainwithedgesmearedcrack(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,'gmshoptions',gmshoptions);
    
    if Dim==2
        CU = LIGNE([0.0,L/2+clC/2],[a,L/2]);
        CL = LIGNE([0.0,L/2-clC/2],[a,L/2]);
    elseif Dim==3
        CU = QUADRANGLE([0.0,L/2+clC/2,0.0],[a,L/2,0.0],[a,L/2,e],[0.0,L/2+clC/2,e]);
        CL = QUADRANGLE([0.0,L/2-clC/2,0.0],[a,L/2,0.0],[a,L/2,e],[0.0,L/2-clC/2,e]);
    end
    
    % sizemap = @(d) (clC-clD)*d+clD;
    sizemap = @(d) clD*clC./((clD-clC)*d+clC);
    
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
    S_phase = addcl(S_phase,CU,'T',1);
    S_phase = addcl(S_phase,CL,'T',1);
    
    d = calc_init_dirichlet(S_phase);
    cl = sizemap(d);
    S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_domain_single_edge_crack'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    S = S_phase;
    
    S_phase = setmaterial(S_phase,mat_phase);
    S_phase = final(S_phase,'duplicate');
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
    % option = 'CONT'; % plane stress [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS]
    % Lame coefficients
    % lambda = 121.1538e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
    % mu = 80.7692e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
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
    % E = 210e9; NU = 0.2; % [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM]
    % kappa = 121030e6; NU=0.227; lambda=3*kappa*NU/(1+NU); mu = 3*kappa*(1-2*NU)/(2*(1+NU)); E = 3*kappa*(1-2*NU); % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
    % Energetic degradation function
    g = @(d) (1-d).^2;
    % Density
    RHO = 1;
    
    % Material
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
                dt0 = 1e-7; % test
                nt0 = 50; % test
                t0 = linspace(dt0,nt0*dt0,nt0);
                dt1 = 1e-9;
                nt1 = 1300;
                dt1 = 1e-8; % test
                nt1 = 130; % test
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
                dt = 1e-7; % test
                nt = 200; % test
                t = linspace(dt,nt*dt,nt);
        end
    elseif Dim==3
        dt = 1e-8;
        nt = 2500;
        dt = 1e-7; % test
        nt = 250; % test
        t = linspace(dt,nt*dt,nt);
    end
    T = TIMEMODEL(t);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','D','C','CU','CL','BU','BL','BRight','BLeft','BFront','BBack','gc','l');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','D','C','CU','CL','BU','BL','BRight','BLeft','BFront','BBack','gc','l');
end

%% Solution
if solveProblem
    
    tTotal = tic;
    
    t = gett(T);
    
    Ht = cell(1,length(T));
    dt = cell(1,length(T));
    ut = cell(1,length(T));
    ft = zeros(1,length(T));
    St_phase = cell(1,length(T));
    St = cell(1,length(T));
    
    sz_phase = getnbddl(S_phase);
    sz = getnbddl(S);
    H = zeros(sz_phase,1);
    u = zeros(sz,1);
    
    fprintf('\n+----------+-----------+-----------+----------+----------+------------+------------+------------+\n');
    fprintf('|   Iter   |  u [mm]   |  f [kN]   | Nb nodes | Nb elems |  norm(H)   |  norm(d)   |  norm(u)   |\n');
    fprintf('+----------+-----------+-----------+----------+----------+------------+------------+------------+\n');
    fprintf('| %8d | %6.3e | %6.3e | %8d | %8d | %9.4e | %9.4e | %9.4e |\n',0,0,0,getnbnode(S),getnbelem(S),0,0,0);
    
    for i=1:length(T)
        
        % Internal energy field
        h_old = double(H);
        H = FENODEFIELD(calc_energyint(S,u,'node','positive'));
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
        mats = MATERIALS(S);
        S_phase_old = S_phase;
        % S_old = S;
        cl = sizemap(d);
        S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_domain_single_edge_crack'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
        S = S_phase;
        
        for m=1:length(mats_phase)
            S_phase = setmaterial(S_phase,mats_phase{m},m);
        end
        S_phase = final(S_phase,'duplicate');
        S_phase = addcl(S_phase,CU,'T',1);
        S_phase = addcl(S_phase,CL,'T',1);
        
        % P_phase = calcProjection(S_phase,S_phase_old,[],'free',false);
        P_phase = calcProjection(S_phase,S_phase_old,[],'free',false,'full',true);
        d = P_phase'*d;
        h = P_phase'*h;
        H = setvalue(H,h);
        
        % Displacement field
        for m=1:length(mats)
            mats{m} = setparam(mats{m},'d',d);
            S = setmaterial(S,mats{m},m);
        end
        S = final(S,'duplicate');
        % S = removebc(S);
        ud = t(i);
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
        
        % P = calcProjection(S,S_old,[],'free',false);
        % P = calcProjection(S,S_old,[],'free',false,'full',true);
        P = kron(P_phase,eye(Dim));
        u = P'*u;
        for m=1:length(mats)
            mats{m} = setparam(mats{m},'u',u);
        end
        S = actualisematerials(S,mats);
        
        [A,b] = calc_rigi(S,'nofree');
        b = -b;
        
        u = freematrix(S,A)\b;
        u = unfreevector(S,u);
        
        switch lower(loading)
            case 'tension'
                numddl = findddl(S,'UY',BU);
            case 'shear'
                numddl = findddl(S,'UX',BU);
            otherwise
                error('Wrong loading case')
        end
        f = A(numddl,:)*u;
        f = sum(f);
        
        % Update fields
        Ht{i} = double(H);
        dt{i} = d;
        ut{i} = u;
        ft(i) = f;
        St_phase{i} = S_phase;
        St{i} = S;
        
        fprintf('| %8d | %6.3e | %6.3e | %8d | %8d | %9.4e | %9.4e | %9.4e |\n',i,t(i)*1e3,ft(i)*((Dim==2)*1e-6+(Dim==3)*1e-3),getnbnode(S),getnbelem(S),norm(Ht{i}),norm(dt{i}),norm(ut{i}));
        
    end
    
    fprintf('+----------+-----------+-----------+----------+----------+------------+------------+------------+\n');
    
    % DO NOT WORK WITH MESH ADAPTATION
    % Ht = TIMEMATRIX(Ht,T,[sz_phase,1]);
    % dt = TIMEMATRIX(dt,T,[sz_phase,1]);
    % ut = TIMEMATRIX(ut,T,[sz,1]);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'Ht','dt','ut','ft','St','St_phase','time');
else
    load(fullfile(pathname,'solution.mat'),'Ht','dt','ut','ft','St','St_phase','time');
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
    [t,rep] = gettevol(T);
    % DO NOT WORK WITH MESH ADAPTATION
    % u = getmatrixatstep(ut,rep(end));
%     u = ut{rep(end)};
    u = ut{end};
    S_final = St{end};
    
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
    mysaveas(pathname,'mesh_init',formats,renderer);
    
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
    if Dim==2
        ylabel('Force [kN/mm]','Interpreter',interpreter)
    elseif Dim==3
        ylabel('Force [kN]','Interpreter',interpreter)
    end
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    %% Display evolution of solutions
    ampl = 0;
    % DO NOT WORK WITH MESH ADAPTATION
    % ampl = getsize(S)/max(max(abs(getvalue(ut))))/20;
    % umax = cellfun(@(u) max(abs(u)),ut,'UniformOutput',false);
    % ampl = getsize(S)/max([umax{:}])/20;
    
    options = {'plotiter',true,'plottime',false};
    framerate = 80;
    
%     evolModel(T,St,'FrameRate',framerate,'filename','mesh','pathname',pathname,options{:});
    
%     evolSolutionCell(T,St_phase,Ht,'FrameRate',framerate,'filename','internal_energy','pathname',pathname,options{:});
    
%     evolSolutionCell(T,St_phase,dt,'FrameRate',framerate,'filename','damage','pathname',pathname,options{:});
%     for i=1:Dim
%         evolSolutionCell(T,St,ut,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i)],'pathname',pathname,options{:});
%     end
    
%     for i=1:(Dim*(Dim+1)/2)
%         evolSolutionCell(T,St,ut,'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['epsilon_' num2str(i)],'pathname',pathname,options{:});
%         evolSolutionCell(T,St,ut,'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['sigma_' num2str(i)],'pathname',pathname,options{:});
%     end
%     
%     evolSolutionCell(T,St,ut,'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename','epsilon_von_mises','pathname',pathname,options{:});
%     evolSolutionCell(T,St,ut,'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename','sigma_von_mises','pathname',pathname,options{:});
    
    %% Display solutions at differents instants
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
        % DO NOT WORK WITH MESH ADAPTATION
        % Hj = getmatrixatstep(Ht,rep(j));
        % dj = getmatrixatstep(dt,rep(j));
        % uj = getmatrixatstep(ut,rep(j));
        Hj = Ht{rep(j)};
        dj = dt{rep(j)};
        uj = ut{rep(j)};
        Sj = St{rep(j)};
        Sj_phase = St_phase{rep(j)};
        
        plotModel(Sj,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
        mysaveas(pathname,['mesh_t' num2str(rep(j))],formats,renderer);
        
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
    % DO NOT WORK WITH MESH ADAPTATION
    % Hi = getmatrixatstep(Ht,rep(i));
    % di = getmatrixatstep(dt,rep(i));
    % ui = getmatrixatstep(ut,rep(i));
    Hi = Ht{rep(i)};
    di = dt{rep(i)};
    ui = ut{rep(i)};
    Si = St{rep(i)};
    % Si_phase = St_phase{rep(i)};
    
    write_vtk_mesh(Si,{Hi,di,ui},[],...
        {'internal energy','damage','displacement'},[],...
        pathname,'solution',1,i-1);
end
make_pvd_file(pathname,'solution',1,length(T));

% myparallel('stop');
