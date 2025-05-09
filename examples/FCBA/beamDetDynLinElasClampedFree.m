%% Clamped-free beam deterministic dynamic linear elasticity %%
%%-----------------------------------------------------------%%

% clc
clearvars
% close all

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = false;

junction = false; % junction modeling

% filenameCamera = 'test_3_C001H001S0001';
filenameCamera = 'PoutreConsole4_C001H001S0001';
pathnameCamera = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','materialWoodDynamicBending','resultsCamera');

filename = 'beamDetDynLinElasClampedFree';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','FCBA',filename,filenameCamera);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc','png'};
renderer = 'OpenGL';

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
    E = 13e9; % [Pa]
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
        c = 10e3; % junction bending stiffness [N.m/rad]
        J = 15; % moment of inertia [kg.m2/rad]=[N.m.s2/rad]
    end
    
    %% Dirichlet boundary conditions
    S = final(S);
    if junction
        S = addcl(S,P1,'U');
    else
        S = addcl(S,P1);
    end
    
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
    
    switch filenameCamera
        case 'test_3_C001H001S0001'
            delta_exp = 2.3e-2; % initial vertical displacement [m]
        case 'PoutreConsole4_C001H001S0001'
            delta_exp = 3.2e-2; % initial vertical displacement [m]
    end
    ind = find(uy_exp>delta_exp,1,'last');
    t(1:ind) = [];
    uy_exp(1:ind) = [];
    t = t-t(1);
    
    %% Initial conditions
    delta = uy_exp(1); % initial vertical displacement [m]
    x = getcoord(getnode(S));
    ux0 = zeros(getnbnode(S),1);
    if junction
        lambda = @(E,c) 3*E*IZ/(c*L);
        funuy0 = @(E,c) delta/(1+lambda(E,c))*((x(:,1).^2).*(3*L-x(:,1))/(2*L^3) + lambda(E,c)*x(:,1)/L);
        funrz0 = @(E,c) delta/(1+lambda(E,c))*(x(:,1).*(2*L-x(:,1))*3/(2*L^3) + lambda(E,c)/L);
        funu0 = @(E,c) [ux0 funuy0(E,c) funrz0(E,c)]';
        u0 = funu0(E,c);
    else
        uy0 = delta*(x(:,1).^2).*(3*L-x(:,1))/(2*L^3);
        rz0 = delta*x(:,1).*(2*L-x(:,1))*3/(2*L^3);
        u0 = [ux0 uy0 rz0]';
    end
    u0 = freevector(S,u0(:));
    v0 = zeros(getnbddlfree(S),1);
    
    %% Time scheme
    T = TIMEMODEL(t);
    
    N = NEWMARKSOLVER(T,'alpha',0,'gamma',1/2,'beta',1/4,'display',false);
    
    loadFunction = @(N) zero(N);
    
    %% Mass, stiffness and damping matrices and sollicitation vectors
    M = calc_mass(S);
    K = calc_rigi(S);
    if junction
        % [~,numnode,~] = intersect(S,P1,'strict',false);
        numnode = find(S.node==P1);
        numddl = findddl(S,'RZ',numnode,'free');
        K(numddl,numddl) = K(numddl,numddl) + c;
        M(numddl,numddl) = M(numddl,numddl) + J;
    end
    % stiffness proportional Rayleigh (viscous) damping coefficient
    % alpha = 0;
    alpha = 1e-5;
    % mass proportional Rayleigh (viscous) damping coefficient
    % beta = 0;
    beta = 3;
    C = alpha*K + beta*M;
    % C = zeros(size(K));
    % pl = RHO*g*Sec; % line load (body load for beams) [N/m]
    % b = bodyload(S,[],'FY',pl);
    b = zeros(getnbddlfree(S),1);
    b = b*loadFunction(N);
    
    save(fullfile(pathname,'problem.mat'),'S','elemtype','N','M','K','C','b','u0','v0','P1','uy_exp');
else
    load(fullfile(pathname,'problem.mat'),'S','elemtype','N','M','K','C','b','u0','v0','P1','uy_exp');
end

%% Solution
if solveProblem
    t = tic;
    [ut,result,vt,at] = ddsolve(N,b,M,K,C,u0,v0);
    time = toc(t);
    
    et = calc_epsilon(S,ut,'node');
    st = calc_sigma(S,ut,'node');
    
    save(fullfile(pathname,'solution.mat'),'ut','result','vt','at','et','st','time');
else
    load(fullfile(pathname,'solution.mat'),'ut','result','vt','at','et','st','time');
end

%% Outputs
fprintf('\n');
fprintf(['data file : ' filenameCamera '\n']);
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
    
%     N = setevolparam(N,'colorbar',true,'FontSize',fontsize);
%     
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
%     
%     figure('Name','Solution N')
%     clf
%     set(gcf,'Color','w')
%     % frame = evol(N,st,S,'compo','EFFX','rescale',true);
%     frame = evol(N,Nt,S,'rescale',true);
%     saveMovie(frame,'filename','N','pathname',pathname);
%     
%     figure('Name','Solution Mz')
%     clf
%     set(gcf,'Color','w')
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
mysaveas(pathname,'quantity_of_interest',formats,renderer);
mymatlab2tikz(pathname,'quantity_of_interest.tex');

for t=0:getnt(ut)
    uk = Ut(:,t+1);
    rzk = Rzt(:,t+1);
    fields = {uk,rzk};
    fieldnames = {'displacement','rotation'};
    write_vtk_mesh(S,fields,[],fieldnames,[],pathname,filename,1,t);
end
make_pvd_file(pathname,filename,1,getnt(ut)+1);
