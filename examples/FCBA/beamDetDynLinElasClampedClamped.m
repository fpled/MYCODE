%% Clamped-clamped beam deterministic dynamic linear elasticity %%
%%--------------------------------------------------------------%%

% clc
clearvars
close all

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = false;

filename = 'beamDetDynLinElasClampedClamped';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','FCBA',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    % Beam dimensions
    L = 1.0; % [m]
    b = 0.06;
    h = 0.025;
    
    % Points
    P1 = POINT([0.0,0.0]);
    P2 = POINT([L,0.0]);
    P_load = POINT([L/2,0.0]);
    Line = LINE(P1,P2);
    
    elemtype = 'BEAM';
    nbelem = 10;
    S = build_model(Line,'nbelem',nbelem,'elemtype',elemtype);
%     cl = 0.1;
%     S = build_model(Line,'cl',cl,'elemtype',elemtype,'filename',fullfile(pathname,'gmsh_domain'));
    
    %% Materials
    % Gravitational acceleration
    g = 9.81; % [m/s2]
    
    % Cross-section area
    Sec = b*h;
    % Planar second moment of area (or Planar area moment of inertia)
    IY = h*b^3/12;
    IZ = b*h^3/12;
    % Polar second moment of area (or Polar area moment of inertia)
    IX = IY+IZ;
    
    % Young modulus
    E = 12e9; % [Pa]
    % Poisson ratio
    NU = 0.3;
    % Density
    RHO = 700; % [kg/m3]
    
    % Material
    mat = ELAS_BEAM('E',E,'NU',NU,'S',Sec,'IZ',IZ,'IY',IY,'IX',IX,'RHO',RHO);
    mat = setnumber(mat,1);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    S = final(S);
    S = addcl(S,P1,'UY');
    S = addcl(S,P2,'UY');
    S = addcl(S,Line,'UX');
    % S = addcl(S,P1,{'UY','RZ'});
    % S = addcl(S,P2,{'UY','RZ'});
    
    %% Initial conditions
    u0 = zeros(getnbddlfree(S),1);
    v0 = zeros(getnbddlfree(S),1);
    
    %% Time scheme
    t0 = 0;
    t1 = 1;
    nt = 50;
    T = TIMEMODEL(t0,t1,nt);
    
    N = NEWMARKSOLVER(T,'alpha',0,'gamma',1/2,'beta',1/4,'display',false);
    
    fmax = 10;
    tf = 0.3;
    % omega = pi/(t1-t0);
    % loadFunction = @(N) fmax * sin(N,omega); % from t0 to t1
    loadFunction = @(N) fmax * dirac(N,t0,tf,'sin'); % from t0 to tf
    
    %% Mass, stiffness and damping matrices and sollicitation vectors
    M = calc_mass(S);
    K = calc_rigi(S);
    % b = nodalload(S,P_load,'FY',-1);
    delta = L/100;
    fun = @(x) -exp(-((2*x(:,1)-L)/(2*delta)).^2);
    b = bodyload(S,[],'FY',fun);
    b = b*loadFunction(N);
    
    save(fullfile(pathname,'problem.mat'),'S','elemtype','N','M','K','b','u0','v0');
else
    load(fullfile(pathname,'problem.mat'),'S','elemtype','N','M','K','b','u0','v0');
end

%% Solution
if solveProblem
    t = tic;
    [ut,result,vt,at] = ddsolve(N,b,M,K,[],u0,v0);
    time = toc(t);
    
    et = calc_epsilon(S,ut,'node');
    st = calc_sigma(S,ut,'node');
    
    save(fullfile(pathname,'solution.mat'),'ut','result','vt','et','st','time');
else
    load(fullfile(pathname,'solution.mat'),'ut','result','vt','et','st','time');
end

ut = unfreevector(S,ut);
ut_val = getvalue(ut);
Ut = ut_val(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
Uxt = ut_val(findddl(S,'UX'),:);
Uyt = ut_val(findddl(S,'UY'),:);
Rzt = ut_val(findddl(S,'RZ'),:);

Epsxt = et(1);
Gamzt = et(2);
Nt = st(1);
Mzt = st(2);

%% Outputs
fprintf('\n');
fprintf(['spatial mesh : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('time solver : %s\n',class(N));
fprintf('nb time steps = %g\n',getnt(N));
fprintf('nb time dofs  = %g\n',getnbtimedof(N));
fprintf('elapsed time = %f s\n',time);

%% Display
if displaySolution
    %% Display domains and meshes  
    plotDomain(S,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    plotModel(S,'Color','k','FaceColor','k','node',true,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    %% Display evolution of solution
    % ampl = 0;
    ampl = getsize(S)/max(max(abs(Ut)))/5;
    
    i = 2;
    % for i=1:2
        evolSolution(S,ut,'displ',i,'ampl',ampl,'filename',['solution_' num2str(i)],'pathname',pathname);
        % evolSolution(S,vt,'displ',i,'ampl',ampl,'filename',['velocity_' num2str(i)],'pathname',pathname);
        % evolSolution(S,at,'displ',i,'ampl',ampl,'filename',['acceleration_' num2str(i)],'pathname',pathname);
    % end
    
    % DO NOT WORK WITH BEAM ELEMENTS
    % evolSolution(S,ut,'epsilon',1,'ampl',ampl,'filename','Epsx','pathname',pathname);
    % evolSolution(S,ut,'epsilon',2,'ampl',ampl,'filename','Gamz','pathname',pathname);
    % evolSolution(S,ut,'sigma',1,'ampl',ampl,'filename','N','pathname',pathname);
    % evolSolution(S,ut,'sigma',2,'ampl',ampl,'filename','Mz','pathname',pathname);
    
    % DO NOT WORK WITH BEAM ELEMENTS
    % evolSolution(S,ut,'epsilon','mises','ampl',ampl,'filename','EpsVM','pathname',pathname);
    % evolSolution(S,ut,'sigma','mises','ampl',ampl,'filename','SigVM','pathname',pathname);
    
    N = setevolparam(N,'colorbar',true,'FontSize',fontsize);
    
%     figure('Name','Solution Epsx')
%     clf
%     set(gcf,'Color','w')
%     % frame = evol(N,et,S,'compo','EPSX','rescale',true);
%     frame = evol(N,Epsxt,S,'rescale',true);
%     saveMovie(frame,'filename','Epsx','pathname',pathname);
%     
%     figure('Name','Solution Gamz')
%     clf
%     set(gcf,'Color','w')
%     % frame = evol(N,et,S,'compo','GAMZ','rescale',true);
%     frame = evol(N,Gamzt,S,'rescale',true);
%     saveMovie(frame,'filename','Gamz','pathname',pathname);
    
    figure('Name','Solution N')
    clf
    set(gcf,'Color','w')
    % frame = evol(N,st,S,'compo','EFFX','rescale',true);
    frame = evol(N,Nt,S,'rescale',true);
    saveMovie(frame,'filename','N','pathname',pathname);
    
    figure('Name','Solution Mz')
    clf
    set(gcf,'Color','w')
    % frame = evol(N,st,S,'compo','MOMZ','rescale',true);
    frame = evol(N,Mzt,S,'rescale',true);
    saveMovie(frame,'filename','Mz','pathname',pathname);
    
    %% Display solution at different instants
    % ampl = 0;
    ampl = getsize(S)/max(max(abs(Ut)))/5;
    [t,rep] = gettevol(N);
    for k=1:floor(length(rep)/5):length(rep)
        close all
        uk = getmatrixatstep(ut,rep(k));
        % vk = getmatrixatstep(vt,rep(k));
        % ak = getmatrixatstep(at,rep(k));
        ek = getmatrixatstep(et,rep(k));
        sk = getmatrixatstep(st,rep(k));
        
        i = 2;
        % for i=1:2
            plotSolution(S,uk,'displ',i,'ampl',ampl);
            mysaveas(pathname,['solution_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            % plotSolution(S,vk,'displ',i,'ampl',ampl);
            % mysaveas(pathname,['velocity_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            % plotSolution(S,ak,'displ',i,'ampl',ampl);
            % mysaveas(pathname,['acceleration_' num2str(i) '_t' num2str(k-1)],formats,renderer);
        % end
        
%         figure('Name','Solution Epsx')
%         clf
%         plot(ek,S+ampl*uk,'compo','EPSX')
%         colorbar
%         set(gca,'FontSize',fontsize)
%         mysaveas(pathname,['Epsx_t' num2str(k-1)],formats,renderer);
%         
%         figure('Name','Solution Gamz')
%         clf
%         plot(ek,S+ampl*uk,'compo','GAMZ')
%         colorbar
%         set(gca,'FontSize',fontsize)
%         mysaveas(pathname,['Gamz_t' num2str(k-1)],formats,renderer);
        
        figure('Name','Solution N')
        clf
        plot(sk,S+ampl*uk,'compo','EFFX')
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['N_t' num2str(k-1)],formats,renderer);
        
        figure('Name','Solution Mz')
        clf
        plot(sk,S+ampl*uk,'compo','MOMZ')
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['Mz_t' num2str(k-1)],formats,renderer);
    end
end

for t=0:getnt(T)
    uk = Ut(:,t+1);
    rzk = Rzt(:,t+1);
    fields = {uk,rzk};
    fieldnames = {'displacement','rotation'};
    write_vtk_mesh(S,fields,[],fieldnames,[],pathname,filename,1,t);
end
make_pvd_file(pathname,filename,1,getnt(T)+1);
