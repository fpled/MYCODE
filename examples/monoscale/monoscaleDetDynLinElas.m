%% Monoscale deterministic linear elasticity dynamic problem %%
%%-----------------------------------------------------------%%

% clc
clearvars
close all

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

filename = 'dynLinElas';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','monoscaleDet',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    D = DOMAIN(2,[0.0,0.0],[2.0,0.4]);
    
    elemtype = 'TRI3';
    % elemtype = 'QUA4';
    % option = 'DEFO'; % plane strain
    option = 'CONT'; % plane stress
    nbelem = [60,12];
    pb.S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'option',option);
%     cl = 0.03;
%     pb.S = build_model(D,'cl',cl,'elemtype',elemtype,'option',option,'filename',fullfile(pathname,'gmsh_domain'));
    
    %% Materials
    % Poisson ratio
    NU = 0;
    % Thickness
    DIM3 = 1;
    % Density
    RHO = 1;
    % Young modulus
    E = 1;
    % Material
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
    mat = setnumber(mat,1);
    pb.S = setmaterial(pb.S,mat);
    
    %% Dirichlet boundary conditions
    L1 = LIGNE(getvertex(D,1),getvertex(D,4));
    L2 = LIGNE(getvertex(D,2),getvertex(D,3));
    
    pb.S = final(pb.S);
    pb.S = addcl(pb.S,L1);
    
    %% Initial conditions
    pb.u0 = zeros(getnbddlfree(pb.S),1);
    pb.v0 = zeros(getnbddlfree(pb.S),1);
    
    %% Time scheme
    t0 = 0;
    t1 = 2;
    nt = 50;
    T = TIMEMODEL(t0,t1,nt);
    
    pb.N = NEWMARKSOLVER(T,'alpha',0.05,'display',false);
    % pb.N = DGTIMESOLVER(T,1,'outputsplit',true,'display',false,'lu',true);
    
    tc = get(T,'t1')/6;
    pb.loadFunction = @(N) rampe(N,t0,tc);
    % pb.loadFunction = @(N) dirac(N,t0,tc);
    % pb.loadFunction = @(N) one(N);
    
    %% Mass, stiffness and damping matrices and sollicitation vectors
    pb.M = calc_mass(pb.S);
    pb.A = calc_rigi(pb.S);
    b = surfload(pb.S,L2,'FX',-1);
    pb.b = b*pb.loadFunction(pb.N);
    
    save(fullfile(pathname,'problem.mat'),'pb','D','L1','L2');
else
    load(fullfile(pathname,'problem.mat'),'pb','D','L1','L2');
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
    figure('Name','Domain')
    clf
    h1 = plot(D,'FaceColor',getfacecolor(1));
    hold on
    h2 = plot(L1,'EdgeColor',getfacecolor(5));
    h3 = plot(L2,'EdgeColor',getfacecolor(6));
    hold off
    set(gca,'FontSize',16)
    l = legend([h1(1),h2(1),h3(1)],'$\Omega$','$\Gamma_D$','$\Gamma_N$');
    set(l,'Interpreter','latex')
    axis image
    axis off
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    figure('Name','Mesh')
    clf
    h1 = plot(pb.S,'FaceColor',getfacecolor(1));
    hold on
    h2 = plot(L1,'EdgeColor',getfacecolor(5));
    h3 = plot(L2,'EdgeColor',getfacecolor(6));
    hold off
    set(gca,'FontSize',16)
    % l = legend([h1(1),h2(1),h3(1)],'$\Omega$','$\Omega_2$','$\Omega_3$','$\Gamma_D^1$','$\Gamma_D^2$');
    % set(l,'Interpreter','latex')
    mysaveas(pathname,'mesh',formats,renderer);
    
    %% Display evolution of solution
    i = 1;
    % for i=1:2
        evolSolution(pb.S,ut,'displ',i,'filename',['evol_solution_' num2str(i)],'pathname',pathname);
        evolSolution(pb.S,ut,'displ',i,'view3',true,'filename',['evol_solution_' num2str(i) '_view3'],'pathname',pathname);
        
        evolSolution(pb.S,vt,'displ',i,'filename',['evol_velocity_' num2str(i)],'pathname',pathname);
        evolSolution(pb.S,vt,'displ',i,'view3',true,'filename',['evol_velocity_' num2str(i) '_view3'],'pathname',pathname);
        
        evolSolution(pb.S,at,'displ',i,'filename',['evol_acceleration_' num2str(i)],'pathname',pathname);
        evolSolution(pb.S,at,'displ',i,'view3',true,'filename',['evol_acceleration_' num2str(i) '_view3'],'pathname',pathname);
    % end
    
    % for i=1:3
    %     evolSolution(pb.S,ut,'epsilon',i,'filename',['evol_epsilon_' num2str(i)],'pathname',pathname);
    %     evolSolution(pb.S,ut,'sigma',i,'filename',['evol_sigma_' num2str(i)],'pathname',pathname);
    % end
    
    % evolSolution(pb.S,ut,'epsilon','mises','filename','evol_epsilon_von_mises','pathname',pathname);
    % evolSolution(pb.S,ut,'sigma','mises','filename','evol_sigma_von_mises','pathname',pathname);
    
    %% Display solution at differents instants
    i = 1;
    % for i=1:2
        [t,rep] = gettevol(pb.N);
        for k=1:floor(length(rep)/5):length(rep)
            close all
            uk = getmatrixatstep(ut,rep(k));
            vk = getmatrixatstep(vt,rep(k));
            ak = getmatrixatstep(at,rep(k));
            
            plotSolution(pb.S,uk,'displ',i);
            mysaveas(pathname,['solution_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            plotSolution(pb.S,uk,'displ',i,'view3',true);
            mysaveas(pathname,['solution_' num2str(i) '_t' num2str(k-1) '_view3'],formats);
            
            plotSolution(pb.S,vk,'displ',i);
            mysaveas(pathname,['velocity_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            plotSolution(pb.S,vk,'displ',i,'view3',true);
            mysaveas(pathname,['velocity_' num2str(i) '_t' num2str(k-1) '_view3'],formats);
            
            plotSolution(pb.S,ak,'displ',i);
            mysaveas(pathname,['acceleration_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            plotSolution(pb.S,ak,'displ',i,'view3',true);
            mysaveas(pathname,['acceleration_' num2str(i) '_t' num2str(k-1) '_view3'],formats);
        end
    % end
end
