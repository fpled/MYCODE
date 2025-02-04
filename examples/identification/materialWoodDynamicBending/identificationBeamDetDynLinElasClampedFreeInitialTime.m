%% Clamped-free beam deterministic dynamic linear elasticity %%
%%-----------------------------------------------------------%%

% clc
clearvars
close all

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = false;

junction = false; % junction modeling

% filenameCamera = 'test_3_C001H001S0001';
filenameCamera = 'PoutreConsole4_C001H001S0001';
pathnameCamera = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','materialWoodDynamicBending','resultsCamera');

filename = 'materialWoodDynamicBendingInitialTime';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification',filename,filenameCamera);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc','png'};
renderer = 'OpenGL';

%% Experimental data
filenameExp = fullfile(pathnameCamera,filenameCamera);
opts = detectImportOptions(filenameExp);
opts.SelectedVariableNames = {'time','Point1_Y_'};
T_exp = readtable(filenameExp,opts);
t = T_exp.time;
uy_exp = T_exp.Point1_Y_;
nanInd = find(isnan(uy_exp));
t(nanInd) = [];
uy_exp(nanInd) = [];
if strcmp(filenameCamera,'PoutreConsole4_C001H001S0001')
    uy_exp = uy_exp*1e-3; % conversion from [mm] to [m]
end

duy_exp = diff(uy_exp);
indmax = min(find(duy_exp<0,1,'last'),find(duy_exp>0,1,'last'))+1;
tmax = t(indmax);
indmin = find(duy_exp>0 & t(1:end-1)<tmax,1,'last')+1;
tmin = t(indmin);
offsetmax = uy_exp(indmax);
offsetmin = uy_exp(indmin);
offset = (min(offsetmax,offsetmin)+max(offsetmax,offsetmin))/2; % vertical offset position [m]
uy_exp = uy_exp-offset;

% if uy_exp(1)>0
    [uymax,Imax] = max(-uy_exp);
% else
%     [uymax,Imax] = max(uy_exp);
% end
Imax = find(abs(uy_exp(1:Imax))>abs(uymax),1,'last');
% if uy_exp(1)>0
    Imin = find(duy_exp(1:Imax)>=0,1,'last');
% else
%     Imin = find(duy_exp(1:Imax)<=0,1,'last');
% end
tinitmin = t(Imin);
tinitmax = t(Imax);

fundelta = @(tinit) uy_exp(find(t>=tinit,1));
funuy_exp = @(tinit) uy_exp(t>=tinit);

%% Identification
% initial guess
tinit0 = t(Imax-1); % initial time [s]
delta0 = fundelta(tinit0)*1e2; % initial vertical displacement [cm]
E0 = 13; % Young modulus [GPa]
alpha0 = eps; % mass proportional Rayleigh (viscous) damping coefficient
beta0 = eps; % stiffness proportional Rayleigh (viscous) damping coefficient
if junction
    c0 = 10; % junction rotational stiffness [kN.m/rad]
    J0 = 15; % moment of inertia [kg.m2/rad]=[N.m.s2/rad]
end

disp('Initial parameters');
disp('------------------');
fprintf('tinit = %g s\n',tinit0);
fprintf('delta = %g cm\n',delta0);
fprintf('E     = %g GPa\n',E0);
fprintf('alpha = %g\n',alpha0);
fprintf('beta  = %g\n',beta0);
if junction
    fprintf('c     = %g kN.m/rad\n',c0);
    fprintf('J     = %g kg.m2/rad\n',J0);
end

param0 = [tinit0 E0 alpha0 beta0];
lb = [t(Imin) 10 0 0];
ub = [t(Imax) 15 Inf Inf];
if junction
    param0 = [param0 c0 J0];
    lb = [lb 0 0];
    ub = [ub Inf Inf];
end

% optimFun = 'lsqnonlin'; % optimization function
% optimFun = 'fminsearch';
% optimFun = 'fminunc';
optimFun = 'fmincon';

% display = 'off';
% display = 'iter';
display = 'iter-detailed';
% display = 'final';
% display = 'final-detailed';

% tolX = 1e-5; % tolerance on the parameter value
% tolFun = 1e-5; % tolerance on the function value

switch optimFun
    case {'lsqnonlin','fminunc','fmincon'}
        options  = optimoptions(optimFun,'Display',display);
        % options  = optimoptions(optimFun,'Display',display,'TolX',tolX,'TolFun',tolFun);
    case 'fminsearch'
        options = optimset('Display',display);
        % options = optimset('Display',display,'TolX',tolX,'TolFun',tolFun);
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
    NU = 0.3;
    % Density
    Vol = Sec*Ltot;
    Mass = 430e-3; % [kg]
    RHO = Mass/Vol; % [kg/m3]
    
    % Material
    mat = ELAS_BEAM('E',E,'NU',NU,'S',Sec,'IZ',IZ,'IY',IY,'IX',IX,'RHO',RHO);
    mat = setnumber(mat,1);
    S = setmaterial(S,mat);
    
    if junction
        c = c0*1e3; % [N.m/rad]
        J = J0; % [kg.m2/rad]=[N.m.s2/rad]
    end
    
    %% Dirichlet boundary conditions
    S = final(S);
    if junction
        S = addcl(S,P1,'U');
    else
        S = addcl(S,P1);
    end
    
    %% Initial conditions
    tinit = tinit0; % [s]
    x = getcoord(getnode(S));
    ux0 = zeros(getnbnode(S),1);
    if junction
        lambda = @(E,c) 3*E*IZ/(c*L);
        funuy0 = @(tinit,E,c) fundelta(tinit)/(1+lambda(E,c))*((x(:,1).^2).*(3*L-x(:,1))/(2*L^3) + lambda(E,c)*x(:,1)/L);
        funrz0 = @(tinit,E,c) fundelta(tinit)/(1+lambda(E,c))*(x(:,1).*(2*L-x(:,1))*3/(2*L^3) + lambda(E,c)/L);
        funu0 = @(tinit,E,c) [ux0 funuy0(tinit,E,c) funrz0(tinit,E,c)]';
        % u0 = funu0(tinit,E,c);
    else
        funuy0 = @(tinit) fundelta(tinit)*(x(:,1).^2).*(3*L-x(:,1))/(2*L^3);
        funrz0 = @(tinit) fundelta(tinit)*x(:,1).*(2*L-x(:,1))*3/(2*L^3);
        funu0 = @(tinit) [ux0 funuy0(tinit) funrz0(tinit)]';
        % u0 = funu0(tinit);
    end
    % u0 = freevector(S,u0(:));
    v0 = zeros(getnbddlfree(S),1);
    
    %% Time scheme
    funt = @(tinit) t(t>=tinit)-t(find(t>=tinit,1));
    funT = @(tinit) TIMEMODEL(funt(tinit));
    funN = @(tinit) NEWMARKSOLVER(funT(tinit),'alpha',0,'gamma',1/2,'beta',1/4,'display',false);
    loadFunction = @(N) zero(N);
    
    %% Mass, stiffness and damping matrices and sollicitation vectors
    M = calc_mass(S);
%     K = calc_rigi(S);
%     M0 = M;
%     K0 = K;
%     if junction
%         % [~,numnode,~] = intersect(S,P1,'strict',false);
%         numnode = find(S.node==P1);
%         numddl = findddl(S,'RZ',numnode,'free');
%         K0(numddl,numddl) = K0(numddl,numddl) + c;
%         M0(numddl,numddl) = M0(numddl,numddl) + J;
%     end
%     C0 = alpha0*K0 + beta0*M0;
    b0 = zeros(getnbddlfree(S),1);
    funb = @(tinit) b0*loadFunction(funN(tinit));
    
    save(fullfile(pathname,'problem.mat'),'S','elemtype','funN','M','funb','funu0','v0','P1');
else
    load(fullfile(pathname,'problem.mat'),'S','elemtype','funN','M','funb','funu0','v0','P1');
end

%% Solution
if solveProblem
    t = tic;
    
    switch optimFun
        case 'lsqnonlin'
            fun = @(param) funlsqnonlinBeamDetDynLinElasClampedFreeInitialTime(param,funuy_exp,S,funN,M,funb,funu0,v0,P1);
            [param,err,~,exitflag,output] = lsqnonlin(fun,param0,lb,ub,options);
        case 'fminsearch'
            fun = @(param) funoptimBeamDetDynLinElasClampedFreeInitialTime(param,funuy_exp,S,funN,M,funb,funu0,v0,P1);
            [param,err,exitflag,output] = fminsearch(fun,param0,options);
        case 'fminunc'
            fun = @(param) funoptimBeamDetDynLinElasClampedFreeInitialTime(param,funuy_exp,S,funN,M,funb,funu0,v0,P1);
            [param,err,exitflag,output] = fminunc(fun,param0,options);
        case 'fmincon'
            fun = @(param) funoptimBeamDetDynLinElasClampedFreeInitialTime(param,funuy_exp,S,funN,M,funb,funu0,v0,P1);
            [param,err,exitflag,output] = fmincon(fun,param0,[],[],[],[],lb,ub,[],options);
    end
    
    tinit = param(1); % [s]
    E = param(2); % [GPa]
    alpha = param(3);
    beta = param(4);
    if junction
        c = param(5); % [kN.m/rad]
        J = param(6); % [kg.m2/rad]=[N.m.s2/rad]
    end
    
%     tinit = tinit0; % initial time [s]
%     E = param(1); % [GPa]
%     alpha = param(2);
%     beta = param(3);
%     if junction
%         c = param(4); % [kN.m/rad]
%         J = param(5); % [kg.m2/rad]=[N.m.s2/rad]
%     end
    uy_exp = funuy_exp(tinit);
    delta = fundelta(tinit)*1e2;
    err = sqrt(err)./norm(uy_exp);
    
    disp('Optimal parameters');
    disp('------------------');
    fprintf('tinit = %g s\n',tinit);
    fprintf('delta = %g cm\n',delta);
    fprintf('E     = %g GPa\n',E);
    fprintf('alpha = %g\n',alpha);
    fprintf('beta  = %g\n',beta);
    if junction
        fprintf('c     = %g kN.m/rad\n',c);
        fprintf('J     = %g kg.m2/rad\n',J);
    end
    fprintf('err = %g\n',err);
    % fprintf('exitflag = %g\n',exitflag);
    % disp(output);
    
    timeIdentification = toc(t);
    
    %% Numerical solution
    t = tic;
    [ut,result,vt,at] = solveBeamDetDynLinElasClampedFreeInitialTime(param,S,funN,M,funb,funu0,v0,P1);
    timeSolution = toc(t);
    
    N = funN(tinit);
    
    et = calc_epsilon(S,ut,'node');
    st = calc_sigma(S,ut,'node');
    
    save(fullfile(pathname,'solution.mat'),'ut','result','vt','at','et','st','N','uy_exp','timeIdentification','timeSolution');
else
    load(fullfile(pathname,'solution.mat'),'ut','result','vt','at','et','st','N','uy_exp','timeIdentification','timeSolution');
end

%% Outputs
fprintf('\n');
fprintf(['data file : ' filenameCamera '\n']);
fprintf(['spatial mesh : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('time solver : %s\n',class(N));
fprintf('nb time steps = %g\n',getnt(ut));
fprintf('nb time dofs  = %g\n',getnbtimedof(ut));
fprintf('elapsed time = %f s for identification\n',timeIdentification);
fprintf('elapsed time = %f s for solution\n',timeSolution);

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
    
%     N = setevolparam(N,'colorbar',true,'FontSize',fontsize);
%     
%     figure('Name','Solution Epsx')
%     clf
%     set(gcf,'color','w')
%     % frame = evol(N,et,S,'compo','EPSX','rescale',true);
%     frame = evol(N,Epsxt,S,'rescale',true);
%     saveMovie(frame,'filename','Epsx','pathname',pathname);
%     
%     figure('Name','Solution Gamz')
%     clf
%     set(gcf,'color','w')
%     % frame = evol(N,et,S,'compo','GAMZ','rescale',true);
%     frame = evol(N,Gamzt,S,'rescale',true);
%     saveMovie(frame,'filename','Gamz','pathname',pathname);
%     
%     figure('Name','Solution N')
%     clf
%     set(gcf,'color','w')
%     % frame = evol(N,st,S,'compo','EFFX','rescale',true);
%     frame = evol(N,Nt,S,'rescale',true);
%     saveMovie(frame,'filename','N','pathname',pathname);
%     
%     figure('Name','Solution Mz')
%     clf
%     set(gcf,'color','w')
%     % frame = evol(N,st,S,'compo','MOMZ','rescale',true);
%     frame = evol(N,Mzt,S,'rescale',true);
%     saveMovie(frame,'filename','Mz','pathname',pathname);
    
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
%         
%         figure('Name','Solution N')
%         clf
%         plot(sk,S+ampl*uk,'compo','EFFX')
%         colorbar
%         set(gca,'FontSize',fontsize)
%         mysaveas(pathname,['N_t' num2str(k-1)],formats,renderer);
%         
%         figure('Name','Solution Mz')
%         clf
%         plot(sk,S+ampl*uk,'compo','MOMZ')
%         colorbar
%         set(gca,'FontSize',fontsize)
%         mysaveas(pathname,['Mz_t' num2str(k-1)],formats,renderer);
    end
end

%% Display quantity of interest
% uyt: vertical displacement at end point as a function of time
uyt = Uyt(end,:);
t = gett(ut);

figure('Name','Quantity of interest : vertical displacement at end point')
clf
plot(t,uy_exp*1e2,'-r','LineWidth',1);
hold on
plot(t,uyt*1e2,'-b','LineWidth',1);
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Time [s]')
ylabel('Vertical displacement [cm]')
legend('Experimental','Numerical')
if junction
    mysaveas(pathname,'quantity_of_interest_junction',formats,renderer);
    mymatlab2tikz(pathname,'quantity_of_interest_junction.tex');
else
    mysaveas(pathname,'quantity_of_interest',formats,renderer);
    mymatlab2tikz(pathname,'quantity_of_interest.tex');
end

for t=0:getnt(ut)
    uk = Ut(:,t+1);
    rzk = Rzt(:,t+1);
    fields = {uk,rzk};
    fieldnames = {'displacement','rotation'};
    write_vtk_mesh(S,fields,[],fieldnames,[],pathname,'solution',1,t);
end
make_pvd_file(pathname,'solution',1,getnt(ut)+1);
