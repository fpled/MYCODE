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
    L = 1.0;
    b = 0.06;
    h = 0.025;
    
    % Points
    P1 = POINT([0.0,0.0]);
    P2 = POINT([L,0.0]);
    P_load = POINT([L/2,0.0]);
    Line = LIGNE(P1,P2);
    
    elemtype = 'BEAM';
    nbelem = 10;
    pb.S = build_model(Line,'nbelem',nbelem,'elemtype',elemtype);
%     cl = 0.1;
%     pb.S = build_model(Line,'cl',cl,'elemtype',elemtype,'filename',fullfile(pathname,'gmsh_domain'));
    
    %% Materials
    % Gravitational acceleration
    g = 9.81;
    
    % Cross-section area
    Sec = b*h;
    % Planar second moment of area (or Planar area moment of inertia)
    IY = h*b^3/12;
    IZ = b*h^3/12;
    % Polar second moment of area (or Polar area moment of inertia)
    IX = IY+IZ;
    
    % Young modulus
    E = 12e9;
    % Poisson ratio
    NU = 0.3;
    % Density
    RHO = 700;
    
    % Material
    mat = ELAS_BEAM('E',E,'NU',NU,'S',Sec,'IZ',IZ,'IY',IY,'IX',IX,'RHO',RHO);
    mat = setnumber(mat,1);
    pb.S = setmaterial(pb.S,mat);
    
    %% Dirichlet boundary conditions
    pb.S = final(pb.S);
    pb.S = addcl(pb.S,P1,'UY');
    pb.S = addcl(pb.S,P2,'UY');
    pb.S = addcl(pb.S,Line,'UX');
    % pb.S = addcl(pb.S,P1,{'UY','RZ'});
    % pb.S = addcl(pb.S,P2,{'UY','RZ'});
    
    %% Initial conditions
    pb.u0 = zeros(getnbddlfree(pb.S),1);
    pb.v0 = zeros(getnbddlfree(pb.S),1);
    
    %% Time scheme
    t0 = 0;
    t1 = 1;
    nt = 50;
    T = TIMEMODEL(t0,t1,nt);
    
    pb.N = NEWMARKSOLVER(T,'alpha',0,'gamma',1/2,'beta',1/4,'display',false);
    
    fmax = 10;
    tf = 0.3;
    % omega = pi/(t1-t0);
    % pb.loadFunction = @(N) fmax * sin(N,omega); % from t0 to t1
    pb.loadFunction = @(N) fmax * dirac(N,t0,tf,'sin'); % from t0 to tf
    
    %% Mass, stiffness and damping matrices and sollicitation vectors
    pb.M = calc_mass(pb.S);
    pb.K = calc_rigi(pb.S);
    % b = nodalload(pb.S,P_load,'FY',-1);
    delta = L/100;
    fun = @(x) -exp(-((2*x(:,1)-L)/(2*delta)).^2);
    b = bodyload(pb.S,[],'FY',fun);
    pb.b = b*pb.loadFunction(pb.N);
    
    save(fullfile(pathname,'problem.mat'),'pb','elemtype','Line');
else
    load(fullfile(pathname,'problem.mat'),'pb','elemtype','Line');
end

%% Solution
if solveProblem
    t = tic;
    [ut,result,vt,at] = ddsolve(pb.N,pb.b,pb.M,pb.K,[],pb.u0,pb.v0);
    time = toc(t);
    
    et = calc_epsilon(pb.S,ut,'node');
    st = calc_sigma(pb.S,ut,'node');
    
    save(fullfile(pathname,'solution.mat'),'ut','result','vt','et','st','time');
else
    load(fullfile(pathname,'solution.mat'),'ut','result','vt','et','st','time');
end

ut = unfreevector(pb.S,ut);
ut_val = getvalue(ut);
Ut = ut_val(findddl(pb.S,DDL(DDLVECT('U',pb.S.syscoord,'TRANS'))),:);
Uxt = ut_val(findddl(pb.S,'UX'),:);
Uyt = ut_val(findddl(pb.S,'UY'),:);
Rzt = ut_val(findddl(pb.S,'RZ'),:);

Epsxt = et(1);
Gamzt = et(2);
Nt = st(1);
Mzt = st(2);

%% Outputs
fprintf('\n');
fprintf(['spatial mesh : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(pb.S));
fprintf('nb nodes    = %g\n',getnbnode(pb.S));
fprintf('nb dofs     = %g\n',getnbddl(pb.S));
fprintf('time solver : %s\n',class(pb.N));
fprintf('nb time steps = %g\n',getnt(pb.N));
fprintf('nb time dofs  = %g\n',getnbtimedof(pb.N));
fprintf('elapsed time = %f s\n',time);

%% Display
if displaySolution
    %% Display domains and meshes  
    plotDomain(pb.S,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    plotModel(pb.S,'Color','k','FaceColor','k','node',true,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    %% Display evolution of solution
    % ampl = 0;
    ampl = getsize(pb.S)/max(max(abs(Ut)))/5;
    
    i = 2;
    % for i=1:2
        evolSolution(pb.S,ut,'displ',i,'ampl',ampl,'filename',['solution_' num2str(i)],'pathname',pathname);
        % evolSolution(pb.S,vt,'displ',i,'ampl',ampl,'filename',['velocity_' num2str(i)],'pathname',pathname);
        % evolSolution(pb.S,at,'displ',i,'ampl',ampl,'filename',['acceleration_' num2str(i)],'pathname',pathname);
    % end
    
    % DO NOT WORK WITH BEAM ELEMENTS
    % evolSolution(pb.S,ut,'epsilon',1,'ampl',ampl,'filename','Epsx','pathname',pathname);
    % evolSolution(pb.S,ut,'epsilon',2,'ampl',ampl,'filename','Gamz','pathname',pathname);
    % evolSolution(pb.S,ut,'sigma',1,'ampl',ampl,'filename','N','pathname',pathname);
    % evolSolution(pb.S,ut,'sigma',2,'ampl',ampl,'filename','Mz','pathname',pathname);
    
    % DO NOT WORK WITH BEAM ELEMENTS
    % evolSolution(pb.S,ut,'epsilon','mises','ampl',ampl,'filename','EpsVM','pathname',pathname);
    % evolSolution(pb.S,ut,'sigma','mises','ampl',ampl,'filename','SigVM','pathname',pathname);
    
    pb.N = setevolparam(pb.N,'colorbar',true,'FontSize',fontsize);
    
    figure('Name','Solution Epsx')
    clf
    set(gcf,'color','w')
    % frame = evol(pb.N,et,pb.S,'compo','EPSX','rescale',true);
    frame = evol(pb.N,Epsxt,pb.S,'rescale',true);
    saveMovie(frame,'filename','Epsx','pathname',pathname);
    
    figure('Name','Solution Gamz')
    clf
    set(gcf,'color','w')
    % frame = evol(pb.N,et,pb.S,'compo','GAMZ','rescale',true);
    frame = evol(pb.N,Gamzt,pb.S,'rescale',true);
    saveMovie(frame,'filename','Gamz','pathname',pathname);
    
    figure('Name','Solution N')
    clf
    set(gcf,'color','w')
    % frame = evol(pb.N,st,pb.S,'compo','EFFX','rescale',true);
    frame = evol(pb.N,Nt,pb.S,'rescale',true);
    saveMovie(frame,'filename','N','pathname',pathname);
    
    figure('Name','Solution Mz')
    clf
    set(gcf,'color','w')
    % frame = evol(pb.N,st,pb.S,'compo','MOMZ','rescale',true);
    frame = evol(pb.N,Mzt,pb.S,'rescale',true);
    saveMovie(frame,'filename','Mz','pathname',pathname);
    
    %% Display solution at differents instants
    % ampl = 0;
    ampl = getsize(pb.S)/max(max(abs(Ut)))/5;
    [t,rep] = gettevol(pb.N);
    for k=1:floor(length(rep)/5):length(rep)
        close all
        uk = getmatrixatstep(ut,rep(k));
        % vk = getmatrixatstep(vt,rep(k));
        % ak = getmatrixatstep(at,rep(k));
        ek = getmatrixatstep(et,rep(k));
        sk = getmatrixatstep(st,rep(k));
        
        i = 2;
        % for i=1:2
            plotSolution(pb.S,uk,'displ',i,'ampl',ampl);
            mysaveas(pathname,['solution_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            % plotSolution(pb.S,vk,'displ',i,'ampl',ampl);
            % mysaveas(pathname,['velocity_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            % plotSolution(pb.S,ak,'displ',i,'ampl',ampl);
            % mysaveas(pathname,['acceleration_' num2str(i) '_t' num2str(k-1)],formats,renderer);
        % end
        
        figure('Name','Solution Epsx')
        clf
        plot(ek,pb.S+ampl*uk,'compo','EPSX')
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['Epsx_t' num2str(k-1)],formats,renderer);
        
        figure('Name','Solution Gamz')
        clf
        plot(ek,pb.S+ampl*uk,'compo','GAMZ')
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['Gamz_t' num2str(k-1)],formats,renderer);
        
        figure('Name','Solution N')
        clf
        plot(sk,pb.S+ampl*uk,'compo','EFFX')
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['N_t' num2str(k-1)],formats,renderer);
        
        figure('Name','Solution Mz')
        clf
        plot(sk,pb.S+ampl*uk,'compo','MOMZ')
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
    write_vtk_mesh(pb.S,fields,[],fieldnames,[],pathname,filename,1,t);
end
make_pvd_file(pathname,filename,1,getnt(T)+1);
