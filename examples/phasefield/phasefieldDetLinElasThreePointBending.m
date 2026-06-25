%% Phase-field fracture model - deterministic linear elasticity problem %%
%  Notched plate under three-point bending                              %%
%%----------------------------------------------------------------------%%

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

test = true; % coarse mesh
% test = false; % fine mesh

Dim = 2; % space dimension Dim = 2, 3
support = 'FlatPunch'; % 'Roller' or 'FlatPunch'
PFmodel = 'Miehe'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'HistoryFieldElem'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 1; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
critConv = 'Energy'; % 'Solution', 'Residual', 'Energy'
initialCrack = 'GeometricCrack'; % 'GeometricCrack', 'GeometricNotch', 'InitialPhaseField'
FEmesh = 'Optim'; % 'Unif' or 'Optim'
structMesh = false; % true or false
optionMesh = []; % [] or 'recombine'
% optionMesh = 'recombine'; % [] or 'recombine'
w = 5e-3; % flat punch / support width
% w = 10e-3; % flat punch / support width

% ws = [5,10]*1e-3; % flat punch / support widths
% for iw=1:length(ws)
% w = ws(iw);
% close all

suffix = '';

foldername = ['threePointBending' support];
if strcmpi(support,'flatpunch')
    foldername = [foldername '_w' num2str(w*1e3) 'mm'];
end
foldername = [foldername '_' num2str(Dim) 'D'];
filename = ['linElas' PFmodel PFsplit PFregularization PFsolver initialCrack...
    'MaxIter' num2str(maxIter)];
if maxIter>1
    filename = [filename 'Tol' num2str(tolConv) num2str(critConv)];
end
filename = [filename 'Mesh' FEmesh];
if structMesh
    filename = [filename 'Structured'];
else
    filename = [filename 'Unstructured'];
end
if ~isempty(optionMesh)
    filename = [filename 'Recombine'];
end
filename = [filename suffix];

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','phasefieldDet',foldername,filename);
if test
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','phasefieldDet_test',foldername,filename);
end
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    L = 160e-3; % length
    H = 40e-3; % height
    ls = L/4; % location of the supports from the outer edges
    % w = 5e-3; % flat punch / support width
    % w = 10e-3; % flat punch / support width
    a = 10e-3; % crack length
    b = L/2; % crack offset from the left outer edge
    e = 40e-3; % thickness
    % if Dim==2
    %     D = DOMAIN(2,[0.0,0.0],[L,H]);
    %     C = LINE([b,0.0],[b,a]);
    % elseif Dim==3
    %     D = DOMAIN(3,[0.0,0.0,0.0],[L,H,e]);
    %     C = QUADRANGLE([b,0.0,0.0],[b,a,0.0],[b,a,e],[b,0.0,e]);
    % end
    
    cl = 0.5e-3;
    if test
        cl = 1e-3;
    end
    switch lower(FEmesh)
        case 'unif' % uniform mesh
            clD = cl;
            B = [];
        case 'optim' % optimized mesh
            % clD = 2e-3;
            % if test
            %     clD = 4e-3;
            % end
            clD = 4*cl;
            VIn = cl; VOut = clD;
            XMin = b-10e-3; XMax = b+10e-3;
            YMin = 0; YMax = H;
            Thickness = min(XMin-ls,L-ls-XMax);
            % Thickness = 0;
            if Dim==2
                B = struct('VIn',VIn,'VOut',VOut,'XMin',XMin,'XMax',XMax,'YMin',YMin,'YMax',YMax,'Thickness',Thickness);
            elseif Dim==3
                ZMin = 0; ZMax = e;
                B = struct('VIn',VIn,'VOut',VOut,'XMin',XMin,'XMax',XMax,'YMin',YMin,'YMax',YMax,'ZMin',ZMin,'ZMax',ZMax,'Thickness',Thickness);
            end
        otherwise
            error('Wrong FE mesh')
    end
    clC = cl;
    clS = cl;
    if structMesh
        % Structured mesh
        switch lower(initialCrack)
            case 'geometriccrack'
                S_phase = gmshThreePointBendingWithSingleEdgeCrackStructured(L,H,ls,w,a,b,e,cl,fullfile(pathname,'gmsh_three_point_bending_single_edge_crack'),Dim,optionMesh);
            case 'geometricnotch'
                c = 2*clC; % crack width
                S_phase = gmshThreePointBendingWithSingleEdgeRectangularNotchStructured(L,H,ls,w,a,b,c,e,cl,fullfile(pathname,'gmsh_three_point_bending_single_edge_crack'),Dim,optionMesh);
            case 'initialphasefield'
                S_phase = gmshThreePointBendingWithSingleEdgeCrackStructured(L,H,ls,w,a,b,e,cl,fullfile(pathname,'gmsh_three_point_bending_single_edge_crack'),Dim,'noduplicate','refinecrack',optionMesh);
            otherwise
                error('Wrong model for initial crack');
        end
        S_phase = concatgroupelem(S_phase);
    else
        % Unstructured mesh
        switch lower(initialCrack)
            case 'geometriccrack'
                S_phase = gmshThreePointBendingWithSingleEdgeCrack(L,H,ls,w,a,b,e,clD,clC,clS,fullfile(pathname,'gmsh_three_point_bending_single_edge_crack'),Dim,optionMesh,'Box',B);
            case 'geometricnotch'
                c = 2*clC; % crack width
                S_phase = gmshThreePointBendingWithSingleEdgeNotch(L,H,ls,w,a,b,c,e,clD,clC,clS,fullfile(pathname,'gmsh_three_point_bending_single_edge_crack'),Dim,optionMesh,'Box',B,'rectangular');
            case 'initialphasefield'
                S_phase = gmshThreePointBendingWithSingleEdgeCrack(L,H,ls,w,a,b,e,clD,clC,clS,fullfile(pathname,'gmsh_three_point_bending_single_edge_crack'),Dim,'noduplicate','refinecrack',optionMesh,'Box',B);
            otherwise
                error('Wrong model for initial crack');
        end
    end
    S = S_phase;
    
    %% Phase-field problem
    %% Material
    % Critical energy release rate (or fracture toughness)
    gc = 10;
    % Regularization parameter (width of the smeared crack)
    l = 2e-3;
    % Small artificial residual stiffness
    % k = 1e-12;
    k = 0;
    
    % Material
    [K,R,Qn] = setphasefieldparam(gc,l,PFregularization);
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'DIM3',e,'PFregularization',PFregularization);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    switch lower(initialCrack)
        case 'geometriccrack'
            C = POINT([b,a]); % crack tip
        case 'geometricnotch'
            C = CIRCLE(b,a-c/2,c/2); % circular notch
            % C = LINE([b-c/2,a],[b+c/2,a]); % rectangular notch
            % C = POINT([b,a]); % V notch
        case 'initialphasefield'
            C = LINE([b,0.0],[b,a]); % crack line
        otherwise
            error('Wrong model for initial crack');
    end
    LU = LINE([0.0,H],[L,H]);
    
    findddlboundary = @(S_phase) findddl(S_phase,'T',LU);
    
    if strcmpi(initialCrack,'geometriccrack')
        S_phase = final(S_phase,'duplicate');
    else
        S_phase = final(S_phase);
    end
    
    S_phase = addbcdamageThreePointBending(S_phase,C,initialCrack);
    
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
    option = 'DEFO'; % plane strain
    % option = 'CONT'; % plane stress
    % Young modulus and Poisson ratio
    E = 20e9;
    NU = 0.2;
    
    % Energetic degradation function
    g = @(d) (1-d).^2;
    % Density
    RHO = 1;
    
    % Material
    d = calc_init_dirichlet(S_phase);
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    P0 = POINT([0.0,0.0]);
    BU = LINE([L/2-w/2,H],[L/2+w/2,H]);
    BL = LINE([ls-w/2,0.0],[ls+w/2,0.0]);
    BR = LINE([L-ls-w/2,0.0],[L-ls+w/2,0.0]);
    
    addbc = @(S,ud) addbcThreePointBending(S,ud,BU,BL,BR,P0);
    findddlforce = @(S) findddl(S,'UY',BU);
    
    if strcmpi(initialCrack,'geometriccrack')
        S = final(S,'duplicate');
    else
        S = final(S);
    end
    
    ud = 0;
    S = addbc(S,ud);
    
    u = calc_init_dirichlet(S);
    mats = MATERIALS(S);
    for m=1:length(mats)
        mats{m} = setparam(mats{m},'u',u);
    end
    S = actualisematerials(S,mats);
    
    %% Stiffness matrices and sollicitation vectors
    % [A,b] = calc_rigi(S);
    % b = -b;
    
    %% Time scheme
    % du = 5e-5 mm during the first 200 time steps (up to u = 10e-3 mm)
    % du = 1e-5 mm during the last 1500 time steps (up to u = 25e-3 mm)
    % dt0 = 5e-8;
    % nt0 = 200;
    % dt1 = 1e-8;
    % nt1 = 1500;
    % if test
    %     dt0 = 5e-7;
    %     nt0 = 20;
    %     dt1 = 1e-7;
    %     nt1 = 150;
    % end
    % t0 = linspace(dt0,nt0*dt0,nt0);
    % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
    % t = [t0,t1];
    % T = TIMEMODEL(t);
    
    % du = 5e-5 mm during the first stage (until the phase-field reaches the threshold value)
    % du = 1e-5 mm during the last stage (as soon as the phase-field exceeds the threshold value, up to u = 0.025 mm)
    dt0 = 5e-8;
    dt1 = 1e-8;
    % dt0 = 1e-7;
    % dt1 = 2e-8;
    if test
        dt0 = 5e-7;
        dt1 = 1e-7;
    end
    tf = 0.025e-3;
    dth = 0.2;
    T = struct('dt0',dt0,'dt1',dt1,'tf',tf,'dth',dth);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','addbc','findddlforce','findddlboundary');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','addbc','findddlforce','findddlboundary');
end

%% Solution
if solveProblem
    tTotal = tic;
    
    displayIter  = true;
    displaySol   = false;
    displayForce = false;
    
    if isa(T,'TIMEMODEL')
        fun = @solvePFDetLinElas;
    else
        fun = @solvePFDetLinElasThreshold;
    end
    switch lower(PFsolver)
        case {'historyfieldelem','historyfieldnode'}
            [dt,ut,ft,Ht,Edt,Eut,output] = fun(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,...
                'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
                'displayiter',displayIter,'displaysol',displaySol,'displayforce',displayForce);
        otherwise
            [dt,ut,ft,~,Edt,Eut,output] = fun(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,...
                'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
                'displayiter',displayIter,'displaysol',displaySol,'displayforce',displayForce);
    end
    % if isa(T,'TIMEMODEL')
    %     fun = @solvePFDetLinElasThreePointBending;
    % else
    %     fun = @solvePFDetLinElasThreePointBendingThreshold;
    % end
    % switch lower(PFsolver)
    %     case {'historyfieldelem','historyfieldnode'}
    %         [dt,ut,ft,Ht,Edt,Eut,output] = fun(S_phase,S,T,PFsolver,BU,BL,BR,P0,LU,...
    %             'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
    %             'displayiter',displayIter,'displaysol',displaySol,'displayforce',displayForce);
    %     otherwise
    %         [dt,ut,ft,~,Edt,Eut,output] = fun(S_phase,S,T,PFsolver,BU,BL,BR,P0,LU,...
    %             'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
    %             'displayiter',displayIter,'displaysol',displaySol,'displayforce',displayForce);
    % end
    
    if isa(T,'TIMEMODEL')
        saveT = false;
    else
        T = gettimemodel(dt);
        saveT = true;
    end
    t = gettevol(T);
    dt_val = getvalue(dt);
    dmaxt = max(dt_val);
    idc = find(dmaxt>=min(0.75,max(dmaxt)),1);
    fc = ft(idc);
    udc = t(idc);
    [fmax,idmax] = max(ft,[],2);
    udmax = t(idmax);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'dt','ut','ft','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        save(fullfile(pathname,'solution.mat'),'Ht','-append');
    end
    if saveT
        save(fullfile(pathname,'solution.mat'),'T','-append');
    end
else
    load(fullfile(pathname,'solution.mat'),'dt','ut','ft','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        load(fullfile(pathname,'solution.mat'),'Ht');
    end
    if ~isa(T,'TIMEMODEL')
        load(fullfile(pathname,'solution.mat'),'T');
    end
end

%% Outputs
filenameResults = fullfile(pathname,'results.txt');
if solveProblem
    fid = fopen(filenameResults,'w');
    fprintf(fid,'Three-point bending with single edge crack\n');
    fprintf(fid,'\n');
    fprintf(fid,'dim      = %d\n',Dim);
    fprintf(fid,'support  = %s\n',support);
    if strcmpi(support,'flatpunch')
        fprintf(fid,'support width = %g mm\n',w*1e3);
    end
    fprintf(fid,'PF model = %s\n',PFmodel);
    fprintf(fid,'PF split = %s\n',PFsplit);
    fprintf(fid,'PF regularization = %s\n',PFregularization);
    fprintf(fid,'PF solver = %s\n',PFsolver);
    fprintf(fid,'nb elements = %g\n',getnbelem(S));
    fprintf(fid,'nb nodes    = %g\n',getnbnode(S));
    fprintf(fid,'nb dofs     = %g\n',getnbddl(S));
    fprintf(fid,'nb time dofs = %g\n',getnbtimedof(T));
    fprintf(fid,'elapsed time = %f s\n',time);
    
    fprintf(fid,'\n');
    fprintf(fid,'fmax  = %g kN\n',fmax*1e-3);
    fprintf(fid,'fc    = %g kN\n',fc*1e-3);
    fprintf(fid,'udmax = %g mm\n',udmax*1e3);
    fprintf(fid,'udc   = %g mm\n',udc*1e3);
    fclose(fid);
end
type(filenameResults) % fprintf('%s', fileread(filenameResults))

%% Display
if Dim==2
    facealpha = 0.1;
    facecolor = 'k';
    facecolordef = 'b';
elseif Dim==3
    facealpha = 1;
    facecolor = 'w';
    facecolordef = 'w';
end

if displayModel
    [t,rep] = gettevol(T);
    
    %% Display domains, boundary conditions and meshes
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 0.5;
    v = calc_init_dirichlet(S);
    [hN,legN] = vectorplot(S,'U',v,ampl,'r','LineWidth',1);
    % legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions_displacement',formats,renderer);
    
    [hD_phase,legD_phase] = plotBoundaryConditions(S_phase,'legend',false);
    % legend(hD_phase,legD_phase,'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions_damage',formats,renderer);
    
    % plotModel(S,'legend',false);
    % mysaveas(pathname,'mesh',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor',facecolor,'FaceAlpha',facealpha,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    u = getmatrixatstep(ut,rep(end));
    ampl = getsize(S)/max(abs(u))/20;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor',facecolordef,'FaceAlpha',facealpha,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor',facecolor,'FaceAlpha',facealpha);
    plot(S+ampl*unfreevector(S,u),'Color','b','FaceColor',facecolordef,'FaceAlpha',facealpha);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
end

%% Display solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display force-displacement curve
    figure('Name','Force vs displacement')
    clf
    plot([0,t]*1e3,[0,ft]*1e-3,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlim tight
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    %% Display maximum damage-displacement curve
    figure('Name','Maximum damage vs displacement')
    clf
    plot([0,t]*1e3,[0,dmaxt],'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlim tight
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Maximum damage','Interpreter',interpreter)
    mysaveas(pathname,'max_damage_displacement',formats);
    mymatlab2tikz(pathname,'max_damage_displacement.tex');
    
    %% Display energy-displacement curves
    figure('Name','Energies vs displacement')
    clf
    plot([0,t]*1e3,[0,Eut],'-b','LineWidth',linewidth)
    hold on
    plot([0,t]*1e3,[0,Edt],'-r','LineWidth',linewidth)
    plot([0,t]*1e3,[0,Eut+Edt],'-k','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlim tight
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Energy [J]','Interpreter',interpreter)
    legend('elastic','fracture','total',...
        'Location','NorthWest','Interpreter',interpreter)
    mysaveas(pathname,'energies_displacement',formats);
    mymatlab2tikz(pathname,'energies_displacement.tex');
    
    %% Display outputs of iterative resolution
    figure('Name','Number of iterations vs displacement')
    clf
    plot(t*1e3,output.iteration,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlim tight
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Number of iterations','Interpreter',interpreter)
    mysaveas(pathname,'nb_iterations_displacement',formats);
    mymatlab2tikz(pathname,'nb_iterations_displacement.tex');
    
    figure('Name','Computing time vs displacement')
    clf
    plot(t*1e3,output.time,'-r','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlim tight
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Computing time [s]','Interpreter',interpreter)
    mysaveas(pathname,'cpu_time_displacement',formats);
    mymatlab2tikz(pathname,'cpu_time_displacement.tex');
    
    figure('Name','Error vs displacement')
    clf
    plot(t*1e3,output.error,'-k','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlim tight
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Error','Interpreter',interpreter)
    mysaveas(pathname,'error_displacement',formats);
    mymatlab2tikz(pathname,'error_displacement.tex');
    
    %% Display solutions at different instants
    ampl = 0;
    if abs(w-5e-3)<eps
        tSnapshots = [12 13 13.5 14 14.5 15 16 18 20 22.5]*1e-6;
    elseif abs(w-10e-3)<eps
        tSnapshots = [11 12 12.5 13 13.5 14 16 18 20 22.5]*1e-6;
    end
    rep = arrayfun(@(x) find(t>x-eps,1),tSnapshots);
    rep = [rep,length(T)];
    % tSnapshots = [tSnapshots,gett1(T)];
    % rep = arrayfun(@(x) find(t>x-eps,1),tSnapshots);
    
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
        % plotSolution(S,uj,'energyint','local','ampl',ampl);
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
    duration = 10; % [s]
    framecount = getnbtimedof(T);
    framerate = framecount/duration;
    
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
    % evolSolution(S,ut,'energyint','local','ampl',ampl,'FrameRate',framerate,'filename','internal_energy_density','pathname',pathname,options{:});
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

% end

% myparallel('stop');
