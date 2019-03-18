%% Monoscale deterministic linear elasticity dynamic beam problem %%
%%----------------------------------------------------------------%%

% clc
clearvars
close all

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

filename = 'dynLinElasBeam';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','monoscaleDet',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    % Beam dimensions
    L = 1.0;
    b = 0.1;
    h = 0.1;
    
    % Points
    P1 = POINT([0.0,0.0]);
    P2 = POINT([L,0.0]);
    P_load = POINT([L/2,0.0]);
    Line = LIGNE(P1,P2);
    
    elemtype = 'BEAM';
    nbelem = 5;
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
    IX = IY+IZ;
    
    % Young modulus
    E = 1;
    % Poisson ratio
    NU = 0.3;
    % Density
    RHO = 1;
    
    % Material
    mat = ELAS_BEAM('E',E,'NU',NU,'S',Sec,'IZ',IZ,'IY',IY,'IX',IX,'RHO',RHO);
    mat = setnumber(mat,1);
    pb.S = setmaterial(pb.S,mat);
    
    %% Dirichlet boundary conditions
    pb.S = final(pb.S);
    pb.S = addcl(pb.S,P1,'UY',0);
    pb.S = addcl(pb.S,P2,'UY',0);
    pb.S = addcl(pb.S,Line,'UX',0);
    % pb.S = addcl(pb.S,P1,{'UY','RZ'},0);
    % pb.S = addcl(pb.S,P2,{'UY','RZ'},0);
    
    %% Initial conditions
    pb.u0 = zeros(getnbddlfree(pb.S),1);
    pb.v0 = zeros(getnbddlfree(pb.S),1);
    
    %% Time scheme
    t0 = 0;
    t1 = 1;
    nt = 10;
    T = TIMEMODEL(t0,t1,nt);
    
    pb.N = NEWMARKSOLVER(T,'alpha',0,'gamma',1/2,'beta',1/4,'display',false);
    
    fmax = 1;
    omega = pi/t1;
    pb.loadFunction = @(N) fmax * sin(N,omega);
    
    %% Mass, stiffness and damping matrices and sollicitation vectors
    pb.M = calc_mass(pb.S);
    pb.A = calc_rigi(pb.S);
    % b = nodalload(pb.S,P_load,'FY',-1);
    delta = L/10;
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
    [ut,result,vt,at] = ddsolve(pb.N,pb.b,pb.M,pb.A,[],pb.u0,pb.v0);
    time = toc(t);
    
    et = calc_epsilon(pb.S,ut);
    st = calc_sigma(pb.S,ut);
    
    save(fullfile(pathname,'solution.mat'),'ut','result','vt','et','st','time');
else
    load(fullfile(pathname,'solution.mat'),'ut','result','vt','et','st','time');
end

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
    
    ut = unfreevector(pb.S,ut);
    ut_val = getvalue(ut);
    Ut = ut_val(findddl(pb.S,DDL(DDLVECT('U',pb.S.syscoord,'TRANS'))),:);
    ampl = getsize(pb.S)/max(max(abs(Ut)))/5;
    
    %% Display evolution of solution
    i = 2;
    % for i=1:2
        evolSolution(pb.S,ut,'displ',i,'ampl',ampl,'filename',['solution_' num2str(i)],'pathname',pathname);
        % evolSolution(pb.S,vt,'displ',i,'ampl',ampl,'filename',['velocity_' num2str(i)],'pathname',pathname);
        % evolSolution(pb.S,at,'displ',i,'ampl',ampl,'filename',['acceleration_' num2str(i)],'pathname',pathname);
    % end
    
    % for i=1:3
    %     evolSolution(pb.S,ut,'epsilon',i,'ampl',ampl,'filename',['epsilon_' num2str(i)],'pathname',pathname);
    %     evolSolution(pb.S,ut,'sigma',i,'ampl',ampl,'filename',['sigma_' num2str(i)],'pathname',pathname);
    % end
    
    % evolSolution(pb.S,ut,'epsilon','mises','ampl',ampl,'filename','epsilon_von_mises','pathname',pathname);
    % evolSolution(pb.S,ut,'sigma','mises','ampl',ampl,'filename','sigma_von_mises','pathname',pathname);
    
    %% Display solution at differents instants
    i = 2;
    % for i=1:2
        [t,rep] = gettevol(pb.N);
        for k=1:floor(length(rep)/5):length(rep)
            close all
            uk = getmatrixatstep(ut,rep(k));
            % vk = getmatrixatstep(vt,rep(k));
            % ak = getmatrixatstep(at,rep(k));
            
            plotSolution(pb.S,uk,'displ',i,'ampl',ampl);
            mysaveas(pathname,['solution_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            % plotSolution(pb.S,vk,'displ',i,'ampl',ampl);
            % mysaveas(pathname,['velocity_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            % plotSolution(pb.S,ak,'displ',i,'ampl',ampl);
            % mysaveas(pathname,['acceleration_' num2str(i) '_t' num2str(k-1)],formats,renderer);
        end
    % end
end
