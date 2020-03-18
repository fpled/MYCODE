%% Clamped-free beam deterministic dynamic linear elasticity %%
%%-----------------------------------------------------------%%

% clc
clearvars
% close all

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = false;

filename = 'beamDetDynLinElasClampedFree';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','FCBA',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

filenameCamera = 'test_3_C001H001S0001.csv';
pathnameCamera = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','materialWoodSlab','resultsCamera');

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Experimental data
filenameExp = fullfile(pathnameCamera,filenameCamera);
opts = detectImportOptions(filenameExp);
opts.SelectedVariableNames = {'time','Point1_Y_'};
T = readtable(filenameExp,opts);
t = T.time;
uyexp = T.Point1_Y_;
nanInd = find(isnan(uyexp));
t(nanInd) = [];
uyexp(nanInd) = [];
uyexp = uyexp-uyexp(end);
startInd = find(uyexp>=uyexp(1),1,'last');
t(1:startInd) = [];
uyexp(1:startInd) = [];
t = t-t(1);

%% Identification
% initial guess
E0 = 10; % Young modulus [GPa]
NU0 = 0.3; % Poisson ratio
delta0 = uyexp(1)*1e2; % initial static displacement [cm]
alpha0 = 5; % mass proportional Rayleigh (viscous) damping coefficient [1e-4]
beta0 = 5; % stiffness proportional Rayleigh (viscous) damping coefficient [1e-4]

param0 = [E0 NU0 delta0 alpha0 beta0];
lb = [0 0 0 0 0];
ub = [Inf 0.5 10 Inf Inf];

optimFun = 'lsqnonlin'; % optimization function
% optimFun = 'fminsearch';
% optimFun = 'fminunc';
% optimFun = 'fmincon';
display = 'off';
tolX = 1e-4; % tolerance on the parameter value
tolFun = 1e-4; % tolerance on the function value

switch optimFun
    case {'lsqnonlin','fminunc','fmincon'}
        options  = optimoptions(optimFun,'Display',display,'TolX',tolX,'TolFun',tolFun);
    case 'fminsearch'
        options = optimset('Display',display,'TolX',tolX,'TolFun',tolFun);
    otherwise
        error(['Wrong optimization function' optimFun])
end

%% Problem
if setProblem
    %% Domains and meshes
    % Beam dimensions
    Ltot = 890e-3; % [m]
    L = 700e-3;
    b = 65e-3;
    h = 15e-3;
    
    % Points
    P1 = POINT([0.0,0.0]);
    P2 = POINT([L,0.0]);
    Line = LIGNE(P1,P2);
    
    elemtype = 'BEAM';
    nbelem = 100;
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
    E = E0*1e9; % [Pa]
    % Poisson ratio
    NU = NU0;
    % Density
    Vol = Sec*Ltot;
    Mass = 430e-3; % [kg]
    RHO = Mass/Vol; % [kg/m3]
    
    % Material
    mat = ELAS_BEAM('E',E,'NU',NU,'S',Sec,'IZ',IZ,'IY',IY,'IX',IX,'RHO',RHO);
    mat = setnumber(mat,1);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    S = final(S);
    S = addcl(S,P1);
    
    %% Initial conditions
    % delta = delta0*1e-2; % [m]
    x = getcoord(getnode(S));
    ux0 = zeros(getnbnode(S),1);
    funuy0 = @(delta) delta/(2*L^3)*(x(:,1).^2).*(3*L-x(:,1));
    funrz0 = @(delta) 3*delta/(2*L^3)*x(:,1).*(2*L-x(:,1));
    funu0 = @(delta) [ux0 funuy0(delta) funrz0(delta)]';
    u0 = funu0(delta);
    u0 = freevector(S,u0(:));
    v0 = zeros(getnbddlfree(S),1);
    
    %% Time scheme
    T = TIMEMODEL(t);
    
    N = NEWMARKSOLVER(T,'alpha',0,'gamma',1/2,'beta',1/4,'display',false);
    
    loadFunction = @(N) zero(N);
    
    %% Mass, stiffness and damping matrices and sollicitation vectors
    M = calc_mass(S);
    K = calc_rigi(S);
    % stiffness proportional Rayleigh (viscous) damping coefficient
    alpha = alpha0*1e-4;
    % mass proportional Rayleigh (viscous) damping coefficient
    beta = beta0*1e-4;
    C = alpha*K + beta*M;
    b = zeros(getnbddlfree(S),1);
    b = b*loadFunction(N);
    
    save(fullfile(pathname,'problem.mat'),'S','elemtype','N','M','b','funu0','v0','uyexp');
else
    load(fullfile(pathname,'problem.mat'),'S','elemtype','N','M','b','funu0','v0','uyexp');
end

%% Solution
if solveProblem
    t = tic;
    
    switch optimFun
        case 'lsqnonlin'
            fun = @(param) funlsqnonlin(param,uyexp,S,N,M,b,funu0,v0);
            [param,err,~,exitflag,output] = lsqnonlin(fun,param0,lb,ub,options);
        case 'fminsearch'
            fun = @(param) funoptim(param,uyexp,S,N,M,b,funu0,v0);
            [param,err,exitflag,output] = fminsearch(fun,param0,options);
        case 'fminunc'
            fun = @(param) funoptim(param,uyexp,S,N,M,b,funu0,v0);
            [param,err,exitflag,output] = fminunc(fun,param0,options);
        case 'fmincon'
            fun = @(param) funoptim(param,uyexp,S,N,M,b,funu0,v0);
            [param,err,exitflag,output] = fmincon(fun,param0,[],[],[],[],lb,ub,[],options);
    end
    
    E = param(1); % [GPa]
    NU = param(2);
    delta = param(3); % [cm]
    alpha = param(4);
    beta = param(5);
    err = sqrt(err)./norm(uyexp);
    
    disp('Optimal parameters');
    disp('------------------');
    fprintf('E  = %g GPa\n',E);
    fprintf('NU = %g\n',NU);
    fprintf('delta = %g cm\n',delta);
    fprintf('alpha = %g\n',alpha);
    fprintf('beta  = %g\n',beta);
    fprintf('err = %g\n',err);
    % fprintf('exitflag = %g\n',exitflag);
    % disp(output);
    
    time = toc(t);
    
    %% Numerical solution
    S = setmaterial();
    u0 = funu0(delta*1e-2);
    K = calc_rigi(S);
    C = alpha*1e-4*K + beta*1e-4*M;
    [ut,result,vt,at] = ddsolve(N,b,M,K,C,u0,v0);
    
    [ut,S] = solveTractionIsotTrans(x,S);
    u = unfreevector(S,u_in);

    et = calc_epsilon(S,ut,'node');
    st = calc_sigma(S,ut,'node');
    
    save(fullfile(pathname,'solution.mat'),'ut','result','vt','et','st','time');
else
    load(fullfile(pathname,'solution.mat'),'ut','result','vt','et','st','time');
end

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
    
    figure('Name','Solution Epsx')
    clf
    set(gcf,'color','w')
    % frame = evol(N,et,S,'compo','EPSX','rescale',true);
    frame = evol(N,Epsxt,S,'rescale',true);
    saveMovie(frame,'filename','Epsx','pathname',pathname);
    
    figure('Name','Solution Gamz')
    clf
    set(gcf,'color','w')
    % frame = evol(N,et,S,'compo','GAMZ','rescale',true);
    frame = evol(N,Gamzt,S,'rescale',true);
    saveMovie(frame,'filename','Gamz','pathname',pathname);
    
    figure('Name','Solution N')
    clf
    set(gcf,'color','w')
    % frame = evol(N,st,S,'compo','EFFX','rescale',true);
    frame = evol(N,Nt,S,'rescale',true);
    saveMovie(frame,'filename','N','pathname',pathname);
    
    figure('Name','Solution Mz')
    clf
    set(gcf,'color','w')
    % frame = evol(N,st,S,'compo','MOMZ','rescale',true);
    frame = evol(N,Mzt,S,'rescale',true);
    saveMovie(frame,'filename','Mz','pathname',pathname);
    
    %% Display solution at differents instants
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
        
        figure('Name','Solution Epsx')
        clf
        plot(ek,S+ampl*uk,'compo','EPSX')
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['Epsx_t' num2str(k-1)],formats,renderer);
        
        figure('Name','Solution Gamz')
        clf
        plot(ek,S+ampl*uk,'compo','GAMZ')
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(pathname,['Gamz_t' num2str(k-1)],formats,renderer);
        
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

%% Display quantity of interest
% uyt: vertical displacement at end point as a function of time
uyt = Uyt(end,:);
t = gett(T);

figure('Name','Quantity of interest : vertical displacement at end point')
clf
plot(t,uyexp*1e2,'-r','LineWidth',1);
hold on
plot(t,uyt*1e2,'-b','LineWidth',1);
hold off
grid on
box on
set(gca,'FontSize',16)
xlabel('Time [s]')
ylabel('Vertical displacement [cm]')
legend('Experimental','Numerical')
mysaveas(pathname,'quantity_of_interest',formats,renderer);
mymatlab2tikz(pathname,'quantity_of_interest.tex');

for t=0:getnt(T)
    uk = Ut(:,t+1);
    rzk = Rzt(:,t+1);
    fields = {uk,rzk};
    fieldnames = {'displacement','rotation'};
    write_vtk_mesh(S,fields,[],fieldnames,[],pathname,filename,1,t);
end
make_pvd_file(pathname,filename,1,getnt(T)+1);
