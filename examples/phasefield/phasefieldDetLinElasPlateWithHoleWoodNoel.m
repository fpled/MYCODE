%% Phase-field fracture model - deterministic linear elasticity problem %%
%  Plate with a central circular hole under compression                 %%
%%----------------------------------------------------------------------%%
% [Noel, Pled, Chevalier, Wilquin, EFM, 2025] (anisotropic phase-field model of He et al.)

% clc
clearvars
close all
% myparallel('start');

%% Experimental data from [Noel, Pled, Chevalier, Wilquin, EFM, 2025]
filenameElas = '_elastic';
filenameGc = '_gc';
pathnameExp = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','phasefield','dataPlateWithHoleWoodNoel');
filenameElas = fullfile(pathnameExp,filenameElas);
filenameGc = fullfile(pathnameExp,filenameGc);
optsElas = detectImportOptions(filenameElas);
optsGc = detectImportOptions(filenameGc);
optsElas.SelectedVariableNames = {'EL','ET','GL','vL'};
optsGc.SelectedVariableNames = {'f_exp_kN_','f_kN_','Gc_mJ_mm2_','l0_mm_'};
T_Elas = readtable(filenameElas,optsElas);
T_Gc = readtable(filenameGc,optsGc);
EL_data = T_Elas.EL; % [MPa]
ET_data = T_Elas.ET; % [MPa]
GL_data = T_Elas.GL; % [MPa]
NUL_data = T_Elas.vL;
f_exp_data = T_Gc.f_exp_kN_; % [kN]
f_data = T_Gc.f_kN_; % [kN]
gc_data = T_Gc.Gc_mJ_mm2_; % [mJ/mm^2]=[kN/m]
l_data = T_Gc.l0_mm_; % [mm]

samplesA = [0,5,7:8,12:16]+1; % samples from crossbeam A
samplesB = [1:4,6,9:11,17]+1; % samples from crossbeam B
samplesSet = union(samplesA,samplesB); % set of samples
samplesRemoved = [2,8,13,17]+1; % removed samples
samplesKept = setdiff(samplesSet,samplesRemoved); % kept samples
samplesKeptA = setdiff(samplesA,samplesRemoved); % selected samples from crossbeam A
samplesKeptB = setdiff(samplesB,samplesRemoved); % selected samples from crossbeam B
numSample = 4+1; % sample number

%% Statistics of elastic properties: mean value and biased standard deviation
fprintf('\n');
fprintf('EL:  mean = %g GPa, std = %g MPa, cv = %g %%\n',mean(EL_data(samplesKept))*1e-3,std(EL_data(samplesKept),1),std(EL_data(samplesKept),1)/mean(EL_data(samplesKept))*1e2);
fprintf('ET:  mean = %g MPa, std = %g MPa, cv = %g %%\n',mean(ET_data(samplesKept)),std(ET_data(samplesKept),1),std(ET_data(samplesKept),1)/mean(ET_data(samplesKept))*1e2);
fprintf('GL:  mean = %g MPa, std = %g MPa, cv = %g %%\n',mean(GL_data(samplesKept)),std(GL_data(samplesKept),1),std(GL_data(samplesKept),1)/mean(GL_data(samplesKept))*1e2);
fprintf('NUL: mean = %g, std = %g, cv = %g %%\n',mean(NUL_data(samplesKept)),std(NUL_data(samplesKept),1),std(NUL_data(samplesKept),1)/mean(NUL_data(samplesKept))*1e2);

fprintf('\n');
fprintf('Group A\n');
fprintf('EL  :  mean = %g GPa, std = %g MPa, cv = %g %%\n',mean(EL_data(samplesKeptA))*1e-3,std(EL_data(samplesKeptA),1),std(EL_data(samplesKeptA),1)/mean(EL_data(samplesKeptA))*1e2);
fprintf('ET  :  mean = %g MPa, std = %g MPa, cv = %g %%\n',mean(ET_data(samplesKeptA)),std(ET_data(samplesKeptA),1),std(ET_data(samplesKeptA),1)/mean(ET_data(samplesKeptA))*1e2);
fprintf('GL  :  mean = %g MPa, std = %g MPa, cv = %g %%\n',mean(GL_data(samplesKeptA)),std(GL_data(samplesKeptA),1),std(GL_data(samplesKeptA),1)/mean(GL_data(samplesKeptA))*1e2);
fprintf('NUL : mean = %g, std = %g, cv = %g %%\n',mean(NUL_data(samplesKeptA)),std(NUL_data(samplesKeptA),1),std(NUL_data(samplesKeptA),1)/mean(NUL_data(samplesKeptA))*1e2);

fprintf('\n');
fprintf('Group B\n');
fprintf('EL  :  mean = %g GPa, std = %g MPa, cv = %g %%\n',mean(EL_data(samplesKeptB))*1e-3,std(EL_data(samplesKeptB),1),std(EL_data(samplesKeptB),1)/mean(EL_data(samplesKeptB))*1e2);
fprintf('ET  :  mean = %g MPa, std = %g MPa, cv = %g %%\n',mean(ET_data(samplesKeptB)),std(ET_data(samplesKeptB),1),std(ET_data(samplesKeptB),1)/mean(ET_data(samplesKeptB))*1e2);
fprintf('GL  :  mean = %g MPa, std = %g MPa, cv = %g %%\n',mean(GL_data(samplesKeptB)),std(GL_data(samplesKeptB),1),std(GL_data(samplesKeptB),1)/mean(GL_data(samplesKeptB))*1e2);
fprintf('NUL : mean = %g, std = %g, cv = %g %%\n',mean(NUL_data(samplesKeptB)),std(NUL_data(samplesKeptB),1),std(NUL_data(samplesKeptB),1)/mean(NUL_data(samplesKeptB))*1e2);

fprintf('\n');
fprintf('Sample #%d\n',numSample);
fprintf('EL  = %g GPa\n',EL_data(numSample)*1e-3);
fprintf('ET  = %g MPa\n',ET_data(numSample));
fprintf('GL  = %g MPa\n',GL_data(numSample));
fprintf('NUL = %g\n',NUL_data(numSample));
fprintf('\n');

%% Input data
setProblem = true;
solveProblem = true;
displayModel = false;
displaySolution = false;
makeMovie = false;
saveParaview = false;

% test = true; % coarse mesh
test = false; % fine mesh

Dim = 2; % space dimension Dim = 2, 3
symmetry = 'Anisot'; % 'Isot' or 'Anisot'. Material symmetry
PFmodel = 'HeFreddi'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT1'; % 'AT1' or 'AT2'
PFsolver = 'HistoryFieldElem'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = Inf; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
critConv = 'Energy'; % 'Solution', 'Residual', 'Energy'
FEmesh = 'Optim'; % 'Unif' or 'Optim'

% PFmodels = {'Bourdin','Amor','Miehe','HeAmor','HeFreddi','Zhang'};
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
% for imaxIter=1:length(maxIters)
% maxIter = maxIters(imaxIter);
% close all

suffix = '';

foldername = ['plateWithHoleWoodNoel_' num2str(Dim) 'D'];
filename = ['linElas' symmetry];
filename = [filename PFmodel PFsplit PFregularization PFsolver...
    'MaxIter' num2str(maxIter)];
if maxIter>1
    filename = [filename 'Tol' num2str(tolConv) num2str(critConv)];
end
filename = [filename 'Mesh' FEmesh suffix];

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
    % [Noel, Pled, Chevalier, Wilquin, EFM, 2025]
    L = 45e-3; % length
    h = 2*L; % height
    e = 20e-3; % thickness
    r = 5e-3; % radius of the hole
    l = L/100; % regularization length (width of the smeared crack)
    if Dim==2
        D = DOMAIN(2,[0.0,0.0],[L,h]);
        C = CIRCLE(L/2,h/2,r);
    elseif Dim==3
        D = DOMAIN(3,[0.0,0.0,0.0],[L,h,e]);
        C = CIRCLE(L/2,h/2,0.0,r);
        C = CYLINDER(C,e); % CYLINDER(L/2,h/2,0.0,r,e);
    end
    
    switch lower(FEmesh)
        case 'unif' % uniform mesh
            cl = l/2;
            if test
                cl = l;
            end
            clD = cl; % characteristic length for domain
            clC = cl; % characteristic length for circular hole
            B = [];
        case 'optim' % optimized mesh
            % [Noel, Pled, Chevalier, Wilquin, EFM, 2025]
            clD = 3*l/2; % characteristic length for domain
            clC = l/2; % characteristic length for circular hole
            if test
                clC = l;
            end
            VIn = clC; VOut = clD;
            switch lower(PFmodel)
                case {'bourdin','amor','heamor'}
                    XMin = 0; XMax = L;
                    YMin = h/2-2*r; YMax = h/2+2*r;
                    % Thickness = h/2-2*r;
                otherwise
                    XMin = L/2-2*r; XMax = L/2+2*r;
                    YMin = 0; YMax = h;
                    % Thickness = L/2-2*r;
            end
            Thickness = 0;
            if Dim==2
                B = struct('VIn',VIn,'VOut',VOut,'XMin',XMin,'XMax',XMax,'YMin',YMin,'YMax',YMax,'Thickness',Thickness);
            elseif Dim==3
                ZMin = 0; ZMax = e;
                B = struct('VIn',VIn,'VOut',VOut,'XMin',XMin,'XMax',XMax,'YMin',YMin,'YMax',YMax,'ZMin',ZMin,'ZMax',ZMax,'Thickness',Thickness);
            end
        otherwise
            error('Wrong FE mesh')
    end
    if Dim==2
        S_phase = gmshDomainWithHole(D,C,clD,clC,fullfile(pathname,'gmsh_plate_with_hole_wood_Noel'),Dim,'Box',B);
        % S_phase = gmsh2femobject(2,fullfile(getfemobjectoptions('path'),'MYCODE','examples','phasefield','gmsh_plate_with_hole_wood_Noel_unif.msh'),2);
        % S_phase = gmsh2femobject(2,fullfile(getfemobjectoptions('path'),'MYCODE','examples','phasefield','gmsh_plate_with_hole_wood_Noel_unif_test.msh'),2);
        % S_phase = gmsh2femobject(2,fullfile(getfemobjectoptions('path'),'MYCODE','examples','phasefield','gmsh_plate_with_hole_wood_Noel_optim.msh'),2);
        % S_phase = gmsh2femobject(2,fullfile(getfemobjectoptions('path'),'MYCODE','examples','phasefield','gmsh_plate_with_hole_wood_Noel_optim_test.msh'),2);
    elseif Dim==3
        S_phase = gmshDomainWithHole(D,C,clD,clC,fullfile(pathname,'gmsh_plate_with_hole_wood_Noel'),Dim,'Box',B,'extrude');
    end
    S = S_phase;
    
    %% Phase-field problem
    %% Material
    % Critical energy release rate (or fracture toughness)
    gc = gc_data(numSample)*1e3; % [J/m^2]=[N/m]
    % Small artificial residual stiffness
    % k = 1e-12;
    k = 0;
    
    % Material
    [K,R,Qn] = setphasefieldparam(gc,l,PFregularization);
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'DIM3',e,'PFregularization',PFregularization);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BRight = LINE([L,0.0],[L,h]);
        BLeft = LINE([0.0,0.0],[0.0,h]);
    elseif Dim==3
        BRight = PLANE([L,0.0,0.0],[L,h,0.0],[L,0.0,e]);
        BLeft = PLANE([0.0,0.0,0.0],[0.0,h,0.0],[0.0,0.0,e]);
    end
    
    findddlboundary = @(S_phase) union(findddl(S_phase,'T',BRight),findddl(S_phase,'T',BLeft));
    
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
    % option = 'DEFO'; % plane strain [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS], [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM]
    option = 'CONT'; % plane stress [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Noel, Pled, Chevalier, Wilquin, EFM, 2025]
    switch lower(symmetry)
        case 'isot' % isotropic material
            % Lame coefficients
            % Young modulus and Poisson ratio
            E = EL_data(numSample)*1e6; % [Pa]
            NU = 0.44;
        case 'anisot' % anisotropic material
            % Longitudinal Young modulus
            EL = EL_data(numSample)*1e6; % [Pa]
            % Transverse Young modulus
            ET = ET_data(numSample)*1e6; % [Pa]
            % Longitudinal shear modulus
            GL = GL_data(numSample)*1e6; % [Pa]
            % Longitudinal Poisson ratio
            NUL = NUL_data(numSample);
            % Transverse Poisson ratio
            NUT = 0.44;
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
        case 'isot' % isotropic material
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
        case 'anisot' % anisotropic material
            mat = ELAS_ISOT_TRANS('AXISL',[0;1;0],'AXIST',[1;0;0],'EL',EL,'ET',ET,'NUL',NUL,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',e);
            mat = setnumber(mat,1);
            S = setoption(S,option);
            S = setmaterial(S,mat);
            Cmat = double(double(calc_opmat(S)));
            mat = ELAS_ANISOT('C',Cmat,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
        otherwise
            error('Wrong material symmetry class');
    end
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BU = LINE([0.0,h],[L,h]);
        BL = LINE([0.0,0.0],[L,0.0]);
    elseif Dim==3
        BU = PLANE([0.0,h,0.0],[L,h,0.0],[0.0,h,e]);
        BL = PLANE([0.0,0.0,0.0],[L,0.0,0.0],[0.0,0.0,e]);
    end
    PD = getvertices(D);
    P1 = POINT(PD{1});
    P2 = POINT(PD{2});

    addbc = @(S,ud) addbcPlateWithHole(S,ud,BU,BL,P1,P2);
    findddlforce = @(S) findddl(S,'UY',BU);
    
    S = final(S);
    
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
    if Dim==2
        % [Noel, Pled, Chevalier, Wilquin, EFM, 2025]
        % du = 24e-5 mm during the first stage (until the phase-field reaches the threshold value)
        % du = 2e-5 mm during the last stage (as soon as the phase-field exceeds the threshold value)
        dt0 = 24e-8;
        dt1 = 2e-8;
        tf = 0.31e-3; % tf = 0.31 mm
        dthreshold = 0.2;
    elseif Dim==3
        % du = 24e-5 mm during the first stage (until the phase-field reaches the threshold value)
        % du = 2e-5 mm during the last stage (as soon as the phase-field exceeds the threshold value)
        dt0 = 24e-8;
        dt1 = 2e-8;
        tf = 0.31e-3; % tf = 0.31 mm
        dthreshold = 0.2;
    end
    T = struct('dt0',dt0,'dt1',dt1,'tf',tf,'dthreshold',dthreshold);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','addbc','findddlforce','findddlboundary');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','addbc','findddlforce','findddlboundary');
end

%% Solution
if solveProblem
    tTotal = tic;
    
    switch lower(PFsolver)
        case {'historyfieldelem','historyfieldnode'}
            [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElasThreshold(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true);
        otherwise
            [dt,ut,ft,~,Edt,Eut,output] = solvePFDetLinElasThreshold(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true);
    end
    % switch lower(PFsolver)
    %     case {'historyfieldelem','historyfieldnode'}
    %         [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElasPlateWithHoleThreshold(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P1,P2,'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true);
    %     otherwise
    %         [dt,ut,ft,~,Edt,Eut,output] = solvePFDetLinElasPlateWithHoleThreshold(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P1,P2,'maxiter',maxIter,'tol',tolConv,'crit',critConv,'displayiter',true);
    % end
    
    T = gettimemodel(dt);
    t = gettevol(T);
    dt_val = getvalue(dt);
    dmaxt = max(dt_val);
    idc = find(dmaxt>=min(0.75,max(dmaxt)),1);
    fc = ft(idc);
    udc = t(idc);
    [fmax,idmax] = max(ft,[],2);
    udmax = t(idmax);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'dt','ut','ft','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','T','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        save(fullfile(pathname,'solution.mat'),'Ht','-append');
    end
else
    load(fullfile(pathname,'solution.mat'),'dt','ut','ft','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','T','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        load(fullfile(pathname,'solution.mat'),'Ht');
    end
end

%% Outputs
if solveProblem
    fid = fopen(fullfile(pathname,'results.txt'),'w');
    fprintf(fid,'Plate with hole\n');
    fprintf(fid,'\n');
    fprintf(fid,'dim      = %d\n',Dim);
    fprintf(fid,'mat sym  = %s\n',symmetry);
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
    
    if Dim==2
        fprintf(fid,'fmax  = %g kN/m\n',fmax*1e-3);
        fprintf(fid,'fc    = %g kN/m\n',fc*1e-3);
    elseif Dim==3
        fprintf(fid,'fmax  = %g kN\n',fmax*1e-3);
        fprintf(fid,'fc    = %g kN\n',fc*1e-3);
    end
    fprintf(fid,'udmax = %g mm\n',udmax*1e3);
    fprintf(fid,'udc   = %g mm\n',udc*1e3);
    fclose(fid);
end

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
    plot(t*1e3,ft*((Dim==2)*1e-6+(Dim==3)*1e-3),'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    %% Display maximum damage-displacement curve
    figure('Name','Maximum damage vs displacement')
    clf
    plot(t*1e3,dmaxt,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Maximum damage','Interpreter',interpreter)
    mysaveas(pathname,'max_damage_displacement',formats);
    mymatlab2tikz(pathname,'max_damage_displacement.tex');
    
    %% Display energy-displacement curves
    figure('Name','Energies vs displacement')
    clf
    plot(t*1e3,Eut,'-b','LineWidth',linewidth)
    hold on
    plot(t*1e3,Edt,'-r','LineWidth',linewidth)
    plot(t*1e3,Eut+Edt,'-k','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Energy [J]','Interpreter',interpreter)
    legend('$\Psi_u$','$\Psi_c$','$\Psi_{\mathrm{tot}}$',...
        'Location','NorthWest','Interpreter','latex')
    mysaveas(pathname,'energies_displacement',formats);
    mymatlab2tikz(pathname,'energies_displacement.tex');
    
    %% Display outputs of iterative resolution
    figure('Name','Number of iterations vs displacement')
    clf
    plot(t*1e3,output.iteration,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
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
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Error','Interpreter',interpreter)
    mysaveas(pathname,'error_displacement',formats);
    mymatlab2tikz(pathname,'error_displacement.tex');
    
    %% Display solutions at different instants
    ampl = 0;
    tSnapshots = [0.1 0.2 0.3 0.4]*1e-3;
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
