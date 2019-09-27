%% Phase field fracture model - deterministic linear elasticity problem  %%
%  Asymmetric notched plate with three holes under three-point bending   %%
%%-----------------------------------------------------------------------%%
% [Ingraffea, Grigoriu, 1990]
% [Bittencourt, Wawrzynek, Ingraffea, Sousa, 1996, EFM]
% [Ventura, Xu, Belytschko, 2002, IJNME]
% [Guidault, Allix, Champaney, Cornuault, 2008, CMAME]
% [Miehe, Welschinger, Hofacker, 2010, IJNME]
% [Miehe, Hofacker, Welschinger, 2010, CMAME]
% [Häusler, Lindhorst, Horst, 2011, IJNME]
% [Geniaut, Galenne, 2012, IJSS]
% [Passieux, Rethore, Gravouil, Baietto, 2013, CM]
% [Ambati, Gerasimov, De Lorenzis, 2015, CM]
% [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM]

% clc
clearvars
close all
% myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = false;

setup = 2; % notch geometry setup = 1, 2
filename = ['phasefieldDetLinElasAsymmetricNotchedPlateSetup' num2str(setup) 'Adaptive'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','phasefield',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
formats = {'fig','epsc'};
renderer = 'OpenGL';

gmshoptions = '-v 0';
mmgoptions = '-nomove -v 0';
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
    % cl = clD;
    % cl = 0.01*unit; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2018, AAM]
    cl = 0.025*unit/2; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME]
    % cl = 0.01*unit/2; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME]
    % cl = 0.01*unit/5; % [Mesgarnejad, Bourdin, Khonsari, 2015, CMAME]
    clC = cl; % characteristic length for edge crack/notch
    clH = cl; % characteristic length for circular holes
    % S_phase = gmshasymmetricnotchedplate(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_domain_asymmetric_notched_plate'),2,'gmshoptions',gmshoptions);
    S_phase = gmshasymmetricnotchedplatewithedgesmearedcrack(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_domain_asymmetric_notched_plate'),2,'gmshoptions',gmshoptions);
    
    CL = LIGNE([-b-clC/2,-h],[-b,-h+a]);
    CR = LIGNE([-b+clC/2,-h],[-b,-h+a]);
    
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
    S_phase = addcl(S_phase,CL,'T',1);
    S_phase = addcl(S_phase,CR,'T',1);
    
    d = calc_init_dirichlet(S_phase);
    cl = sizemap(d);
    S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_domain_asymmetric_notched_plate'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    S = S_phase;
    
    S_phase = setmaterial(S_phase,mat_phase);
    S_phase = final(S_phase,'duplicate');
    % S_phase = addcl(S_phase,C,'T',1);
    S_phase = addcl(S_phase,CL,'T',1);
    S_phase = addcl(S_phase,CR,'T',1);
    
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
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    PU = POINT([0.0,h]);
    PLeft = POINT([-9*unit,-h]);
    PRight = POINT([9*unit,-h]);
    
    S = final(S,'duplicate');
    
    ud = 0;
    S = addcl(S,PU,'UY',ud);
    S = addcl(S,PLeft,{'UX','UY'});
    S = addcl(S,PRight,'UY');
    
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
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','CL','CR','C','PU','PLeft','PRight','gc','l','E','g');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','sizemap','CL','CR','C','PU','PLeft','PRight','gc','l','E','g');
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
    
    fprintf('\n+----------+-----------+----------+----------+------------+------------+------------+\n');
    fprintf('|   Iter   |  u (mm)   | Nb nodes | Nb elems |  norm(H)   |  norm(d)   |  norm(u)   |\n');
    fprintf('+----------+-----------+----------+----------+------------+------------+------------+\n');
    fprintf('| %8d | %6.3e | %8d | %8d | %9.4e | %9.4e | %9.4e |\n',0,0,getnbnode(S),getnbelem(S),0,0,0);
    
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
        S_phase = adaptmesh(S_phase,cl,fullfile(pathname,'gmsh_domain_asymmetric_notched_plate'),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
        S = S_phase;
        
        for m=1:length(mats_phase)
            S_phase = setmaterial(S_phase,mats_phase{m},m);
        end
        % S_phase = actualisematerials(S_phase,mats_phase);
        S_phase = final(S_phase,'duplicate');
        S_phase = addcl(S_phase,CL,'T',1);
        S_phase = addcl(S_phase,CR,'T',1);
        
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
        S = addcl(S,PU,'UY',ud);
        S = addcl(S,PLeft,{'UX','UY'});
        S = addcl(S,PRight,'UY');
        
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
        
        fprintf('| %8d | %6.3e | %8d | %8d | %9.4e | %9.4e | %9.4e |\n',i,t(i)*1e3,getnbnode(S),getnbelem(S),norm(Ht{i}),norm(dt{i}),norm(ut{i}));
        
    end
    
    fprintf('+----------+-----------+----------+----------+------------+------------+------------+\n');
    
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
    
    tmax = 900;
    [t,rep] = gettevol(T);
    t = t(1:tmax);
    rep = rep(1:tmax);
    T = TIMEMODEL(t);
    
    Ht_new = cell(1,tamx);
    dt_new = cell(1,tamx);
    ut_new = cell(1,tamx);
    St_new = cell(1,tamx);
    St_phase_new = cell(1,tamx);
    for i=1:tamx
        Ht_new{i} = Ht{i};
        dt_new{i} = dt{i};
        ut_new{i} = ut{i};
        St_new{i} = St{i};
        St_phase_new{i} = St_phase{i};
    end
    Ht = Ht_new;
    dt = dt_new;
    ut = ut_new;
    St = St_new;
    St_phase = St_phase_new;
    
    % DO NOT WORK WITH MESH ADAPTATION
    % u = getmatrixatstep(ut,rep(end));
    u = ut{rep(end)};
    S = St{end};
    
    %% Display domains, boundary conditions and meshes
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
    
    evolModel(T,St,'FrameRate',80,'filename','mesh','pathname',pathname,options{:});
    
%     evolSolutionCell(T,St_phase,Ht,'filename','internal_energy','pathname',pathname,options{:});
    
    evolSolutionCell(T,St_phase,dt,'FrameRate',80,'filename','damage','pathname',pathname,options{:});
    for i=1:2
        evolSolutionCell(T,St,ut,'displ',i,'FrameRate',80,'ampl',ampl,'FrameRate',60,'filename',['displacement_' num2str(i)],'pathname',pathname,options{:});
    end
    
%     for i=1:3
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
        
        for i=1:2
            plotSolution(Sj,uj,'displ',i,'ampl',ampl);
            mysaveas(pathname,['displacement_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        end
        
%         for i=1:3
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
