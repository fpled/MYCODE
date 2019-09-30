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
filename = ['phasefieldDetLinElasAsymmetricNotchedPlateSetup' num2str(setup)];
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
    S_phase = gmshasymmetricnotchedplatewithedgecrack(a,b,clD,clC,clH,unit,fullfile(pathname,'gmsh_domain_asymmetric_notched_plate'));
    S = S_phase;
    
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
    BU = CIRCLE(0.0,h,2.5*unit);
    BL = CIRCLE(-9*unit,-h,2.5*unit);
    BR = CIRCLE(9*unit,-h,2.5*unit);
    
    S_phase = final(S_phase,'duplicate');
    S_phase = addcl(S_phase,C,'T',1);
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
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
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
    % nt = 3000;
    % t = linspace(dt,nt*dt,nt);
    
    % [Ambati, Gerasimov, De Lorenzis, 2015, CM]
    % du = 1e-3 mm during the first 200 time steps (up to u = 0.2 mm)
    % du = 1e-4 mm during the last  1000 time steps (up to u = 0.3 mm)
    dt0 = 1e-3*unit;
    nt0 = 200;
    t0 = linspace(dt0,nt0*dt0,nt0);
    dt1 = 1e-4*unit;
    nt1 = 1000;
    t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
    t = [t0,t1];
    
    T = TIMEMODEL(t);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','C','B','BU','BL','BR','PU','PL','PR','gc','l','E','g');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','C','B','BU','BL','BR','PU','PL','PR','gc','l','E','g');
end

%% Solution
if solveProblem
    
    tTotal = tic;
    
    t = gett(T);
    
    Ht = cell(1,length(T));
    dt = cell(1,length(T));
    ut = cell(1,length(T));
    ft = zeros(1,length(T));
    
    sz_phase = getnbddl(S_phase);
    sz = getnbddl(S);
    H = zeros(sz_phase,1);
    u = zeros(sz,1);
    
    fprintf('\n+----------+-----------+-----------+------------+------------+------------+\n');
    fprintf('|   Iter   |  u [mm]   | f [kN/mm] |  norm(H)   |  norm(d)   |  norm(u)   |\n');
    fprintf('+----------+-----------+-----------+------------+------------+------------+\n');
    
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
        
        % Displacement field
        for m=1:length(mats)
            mats{m} = setparam(mats{m},'E',FENODEFIELD(E.*(g(d)+k)));
        end
        S = actualisematerials(S,mats);
        S = removebc(S);
        ud = -t(i);
        S = addcl(S,PU,'UY',ud);
        S = addcl(S,PL,{'UX','UY'});
        S = addcl(S,PR,'UY');
        
        [A,b] = calc_rigi(S,'nofree');
        b = -b;
        
        u = freematrix(S,A)\b;
        u = unfreevector(S,u);
        
        numddl = findddl(S,'UY',PU);
        f = -A(numddl,:)*u;
        % f = sum(f);
        
        % Update fields
        Ht{i} = double(H);
        dt{i} = d;
        ut{i} = u;
        ft(i) = f;
        
        fprintf('| %8d | %6.3e | %6.3e | %9.4e | %9.4e | %9.4e |\n',i,t(i)*1e3,ft(i)*1e-6,norm(Ht{i}),norm(dt{i}),norm(ut{i}));
        
    end
    
    fprintf('+----------+-----------+-----------+------------+------------+------------+\n');
    
    Ht = TIMEMATRIX(Ht,T,[sz_phase,1]);
    dt = TIMEMATRIX(dt,T,[sz_phase,1]);
    ut = TIMEMATRIX(ut,T,[sz,1]);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'Ht','dt','ut','ft','time');
else
    load(fullfile(pathname,'solution.mat'),'Ht','dt','ut','ft','time');
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
    u = getmatrixatstep(ut,rep(end));
    
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
    
    %% Display force-displacement curve
    figure('Name','Force-displacement')
    clf
    plot(t*1e3,ft*1e-6,'-b','Linewidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN/mm]','Interpreter',interpreter)
    mysaveas(pathname,'force_displacement',formats);
    
    %% Display evolution of solutions
    ampl = 0;
    % ampl = getsize(S)/max(max(abs(getvalue(ut))))/20;
    
    options = {'plotiter',true,'plottime',false};
    framerate = 80;
    
%     evolSolution(S_phase,Ht,'FrameRate',framerate,'filename','internal_energy','pathname',pathname,options{:});
    
    evolSolution(S_phase,dt,'FrameRate',framerate,'filename','damage','pathname',pathname,options{:});
    for i=1:2
        evolSolution(S,ut,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i)],'pathname',pathname,options{:});
    end
    
%     for i=1:3
%         evolSolution(S,ut,'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['epsilon_' num2str(i)],'pathname',pathname,options{:});
%         evolSolution(S,ut,'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['sigma_' num2str(i)],'pathname',pathname,options{:});
%     end
%     
%     evolSolution(S,ut,'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename','epsilon_von_mises','pathname',pathname,options{:});
%     evolSolution(S,ut,'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename','sigma_von_mises','pathname',pathname,options{:});
    
    %% Display solutions at differents instants
    rep = find(abs(t-0.223*unit)<eps | abs(t-0.225*unit)<eps | abs(t-0.227*unit)<eps | abs(t-0.230*unit)<eps)
    for j=1:length(rep)
        close all
        Hj = getmatrixatstep(Ht,rep(j));
        dj = getmatrixatstep(dt,rep(j));
        uj = getmatrixatstep(ut,rep(j));
        
%         plotSolution(S_phase,Hj);
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
