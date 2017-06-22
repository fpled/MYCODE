%% Monoscale stochastic transient linear advection-diffusion-reaction problem %%
%%----------------------------------------------------------------------------%%
% [Pares, Diez, Huerta, 2008], [Nouy, 2010]

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');
myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;
testSolution = true;

filename = 'transientLinAdvDiffReac';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','monoscaleSto',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    cl1 = 0.02;
    cl2 = 0.04;
    cl0 = 0.02;
    cltip = 0.01;
    pb.S = gmshcanister(cl1,cl2,cl0,cltip,fullfile(pathname,'gmsh_canister'));
    
    %% Random variables
    d = 3; % parametric dimension
    v = UniformRandomVariable(0,1);
    rv = RandomVector(v,d);
    
    %% Materials
    % Linear diffusion coefficient
    % K1(xi) = 0.01 * (1 + xi)
    % K2 = 0.01
    p = 1;
    basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
    bases = FunctionalBases.duplicate(basis,d);
    rvb = getRandomVector(bases);
    H = FullTensorProductFunctionalBasis(bases);
    I = gaussIntegrationRule(v,2);
    I = I.tensorize(d);
    
    fun = @(xi) 0.01 * (1 + xi(:,1));
    funtr = @(xi) fun(transfer(rvb,rv,xi));
    fun = MultiVariateFunction(funtr,d);
    fun.evaluationAtMultiplePoints = true;
    
    K1 = H.projection(fun,I);
    K2 = 0.01;
    
    % Thermal capacity
    % c1(xi) = 1 + xi
    % c2 = 1
    fun = @(x) 1 + x(:,2);
    funtr = @(x) fun(transfer(rvb,rv,x));
    fun = MultiVariateFunction(funtr,d);
    fun.evaluationAtMultiplePoints = true;
    
    c1 = H.projection(fun,I);
    c2 = 1;
    
    % Advection velocity
    Sadv = pb.S;
    mat = FOUR_ISOT('k',1);
    mat = setnumber(mat,1);
    Sadv = setmaterial(Sadv,mat);
    P = @(i) POINT(getnode(getridge(Sadv,i)));
    L1 = LIGNE(P(5),P(6));
    L2 = LIGNE(P(15),P(16));
    Sadv = final(Sadv);
    Sadv = addcl(Sadv,P(1),'T',0);
    A = calc_rigi(Sadv);
    b1 = surfload(Sadv,L1,'QN',-1);
    b2 = surfload(Sadv,L2,'QN',1);
    b = b1+b2;
    phi = A\b;
    v = 2*FENODEFIELD(calc_sigma(Sadv,phi,'node'));
    V = getvalue(v);
    V = {{FENODEFIELD(V(:,1)),FENODEFIELD(V(:,2))}};
    % Linear reaction parameter
    % R1(xi) = 0.1 * (1 + xi)
    % R2 = 10
    fun = @(x) 0.1 * (1 + x(:,3));
    funtr = @(x) fun(transfer(rvb,rv,x));
    fun = MultiVariateFunction(funtr,d);
    fun.evaluationAtMultiplePoints = true;
    
    R1 = H.projection(fun,I);
    R2 = 10;
    
    % Materials
    mat = MATERIALS();
    mat{1} = FOUR_ISOT('k',K1,'c',c1,'b',V,'r',R1);
    mat{2} = FOUR_ISOT('k',K2,'c',c2,'b',V,'r',R2);
    mat{1} = setnumber(mat{1},1);
    mat{2} = setnumber(mat{2},2);
    pb.S = setmaterial(pb.S,mat{1},1);
    pb.S = setmaterial(pb.S,mat{2},2:3);
    
    %% Dirichlet boundary conditions
    numnode = getnumber(getnode(create_boundary(pb.S)));
    [~,numnode1] = intersect(pb.S,L1);
    [~,numnode2] = intersect(pb.S,L2);
    numnoderest = setdiff(setdiff(numnode,numnode1),numnode2);
    
    pb.S = final(pb.S);
    pb.S = addcl(pb.S,numnode1,'T',1);
    pb.S = addcl(pb.S,numnode2,'T',0);
    permeability = true; % if false, the walls of the canister are considered to be impermeable
    if permeability
        pb.S = addcl(pb.S,numnoderest,'T',0);
    end
    
    %% Initial conditions
    % already taken into account by Dirichlet boundary conditions
    % pb.u0 = calc_init_dirichlet(pb.S);
    
    %% Time scheme
    t0 = 0;
    t1 = 2;
    nt = 100;
    pb.timeModel = TIMEMODEL(t0,t1,nt);
    
    % pb.timeSolver = EULERTIMESOLVER(pb.timeModel,'eulertype','explicit','display',false);
    pb.timeSolver = EULERTIMESOLVER(pb.timeModel,'eulertype','implicit','display',false);
    % pb.timeSolver = DGTIMESOLVER(pb.timeModel,1,'outputsplit',true,'display',false,'lu',true);
    
    pb.loadFunction = @(N) one(N);
    
    %% Mass and stifness matrices and sollicitation vectors
    if ~israndom(pb.S)
        pb.M = calc_mass(pb.S);
        [pb.A,pb.b0] = calc_rigi(pb.S);
        pb.b0 = -pb.b0;
        pb.b = pb.b0*pb.loadFunction(pb.timeSolver);
    end
    
    save(fullfile(pathname,'problem.mat'),'pb','Sadv','v','phi');
else
    load(fullfile(pathname,'problem.mat'),'pb','Sadv','v','phi');
end

%% Solution
if solveProblem
    p = 50;
    basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
    bases = FunctionalBases.duplicate(basis,d);
    rv = getRandomVector(bases);
    
    s = AdaptiveSparseTensorAlgorithm();
    % s.nbSamples = 1;
    % s.addSamplesFactor = 0.1;
    s.tol = 1e-6;
    s.tolStagnation = 1e-1;
    % s.tolOverfit = 1.1;
    % s.bulkParameter = 0.5;
    % s.adaptiveSampling = true;
    % s.adaptationRule = 'reducedmargin';
    s.maxIndex = p;
    % s.display = true;
    % s.displayIterations = true;
    
    ls = LeastSquaresSolver();
    ls.regularization = false;
    % ls.regularizationType = 'l1';
    ls.errorEstimation = true;
    % ls.errorEstimationType = 'leaveout';
    % ls.errorEstimationOptions.correction = true;
    
    % Stationary solution
    pbt = pb;
    pb.timeModel = [];
    pb.timeSolver = [];
    fun = @(xi) solveSystem(calcOperator(funEval(pb,xi)));
    fun = MultiVariateFunction(fun,d,getnbddlfree(pb.S));
    fun.evaluationAtMultiplePoints = false;
    
    t = tic;
    [u,err,~,y] = s.leastSquares(fun,bases,ls,rv);
    time = toc(t);
    
    % Transient solution
    funt = @(xi) solveSystem(calcOperator(funEval(pbt,xi)));
    funt = MultiVariateFunction(funt,d,2);
    funt.evaluationAtMultiplePoints = false;
    tt = tic;
    [ut,errt,~,yt] = s.leastSquaresCell(funt,bases,ls,rv);
    timet = toc(tt);
    
    save(fullfile(pathname,'solution.mat'),'u','err','y','fun','time');
    save(fullfile(pathname,'solution_transient.mat'),'ut','errt','yt','funt','timet');
else
    load(fullfile(pathname,'solution.mat'),'u','err','y','fun','time');
    load(fullfile(pathname,'solution_transient.mat'),'ut','errt','yt','funt','timet');
end

%% Outputs
fprintf('\n');
fprintf('nb elements = %g\n',getnbelem(pbt.S));
fprintf('nb nodes    = %g\n',getnbnode(pbt.S));
fprintf('nb dofs     = %g\n',getnbddl(pbt.S));
fprintf('time solver : %s\n',class(pbt.timeSolver));
fprintf('nb time steps = %g\n',getnt(pbt.timeSolver));
fprintf('nb time dofs  = %g\n',getnbtimedof(pbt.timeSolver));

fprintf('\n');
fprintf('Stationary solution');
fprintf('parametric dimension = %d\n',ndims(u.basis))
fprintf('basis dimension = %d\n',numel(u.basis))
fprintf('order = [ %s ]\n',num2str(max(u.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(u.basis.indices.array)
fprintf('nb samples = %d\n',size(y,1))
fprintf('CV error = %d\n',norm(err))
fprintf('elapsed time = %f s\n',time)

fprintf('\n');
fprintf('Transient solution');
fprintf('parametric dimension = %d\n',ndims(ut.basis))
fprintf('basis dimension = %d\n',numel(ut.basis))
fprintf('order = [ %s ]\n',num2str(max(ut.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(ut.basis.indices.array)
fprintf('nb samples = %d\n',size(yt,1))
fprintf('CV error = %d\n',norm(errt))
fprintf('elapsed time = %f s\n',timet)
fprintf('\n');

%% Display
if displaySolution
    %% Display domains and meshes
    figure('Name','Domain')
    clf
    plot(create_boundary(pb.S));
    hold on
    h1 = plot(pb.S,'selgroup',1,'FaceColor',getfacecolor(1),'EdgeColor','none');
    h2 = plot(pb.S,'selgroup',2,'FaceColor',getfacecolor(2),'EdgeColor','none');
    h3 = plot(pb.S,'selgroup',3,'FaceColor',getfacecolor(3),'EdgeColor','none');
    h4 = plotfacets(pb.S,5,'FaceColor',getfacecolor(4),'EdgeColor','none');
    h5 = plotfacets(pb.S,17,'FaceColor',getfacecolor(5),'EdgeColor','none');
    hold off
    set(gca,'FontSize',16)
    l = legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'$\Omega_1$','$\Omega_2$','$\Omega_3$','$\Gamma_D^1$','$\Gamma_D^2$');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    figure('Name','Advection velocity')
    clf
    plot(phi,Sadv);
    colorbar
    set(gca,'FontSize',16)
    hold on
    ampl = 6;
    quiver(v,Sadv,ampl,'k');
    hold off
    ylim([0,1.7])
    mysaveas(pathname,'advection_velocity',formats,renderer);

    % plotModel(pb.S,'legend',false);
    figure('Name','Mesh')
    clf
    h1 = plot(pb.S,'selgroup',1,'FaceColor',getfacecolor(1));
    hold on
    h2 = plot(pb.S,'selgroup',2,'FaceColor',getfacecolor(2));
    h3 = plot(pb.S,'selgroup',3,'FaceColor',getfacecolor(3));
    h4 = plotfacets(pb.S,5,'FaceColor',getfacecolor(4));
    h5 = plotfacets(pb.S,17,'FaceColor',getfacecolor(5));
    hold off
    set(gca,'FontSize',16)
    % l = legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'$\Omega_1$','$\Omega_2$','$\Omega_3$','$\Gamma_D^1$','$\Gamma_D^2$');
    % set(l,'Interpreter','latex')
    mysaveas(pathname,'mesh',formats,renderer);
    
    %% Display stationary solution
    plotSolution(pb.S,u);
    mysaveas(pathname,'solution',formats,renderer);
    plotSolution(pb.S,u,'surface',true);
    mysaveas(pathname,'solution_surface',formats);
    
    %% Display evolution of transient solution
    evolSolution(pb.S,ut,'filename','evol_solution','pathname',pathname);
    evolSolution(pb.S,ut,'surface',true,'filename','evol_solution_surface','pathname',pathname);
    
    evolSolution(pb.S,vt,'rescale',false,'filename','evol_velocity','pathname',pathname);
    evolSolution(pb.S,vt,'rescale',false,'surface',true,'filename','evol_velocity_surface','pathname',pathname);
    
%     for i=1:2
%         evolSolution(pb.S,ut,'epsilon',i,'filename',['evol_eps_' num2str(i)],'pathname',pathname);
%         evolSolution(pb.S,ut,'sigma',i,'filename',['evol_sig_' num2str(i)],'pathname',pathname);
%     end
    
    %% Display quantity of interest
    % boutput: concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) as a function of time
    % Ioutput: total concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) along the complete time evolution,
    %          corresponding to all the pollutant that the actual filter
    %          (group #1 in mesh) is not able to retain
    ut = unfreevector(pb.S,ut);
    foutput = bodyload(keepgroupelem(pb.S,2),[],'QN',1,'nofree');
    boutput = foutput'*ut;
    
    figure('Name','Quantity of interest')
    clf
    plot(boutput,'-b','LineWidth',1);
    grid on
    box on
    set(gca,'FontSize',16)
    xlabel('Time (s)')
    ylabel('Quantity of interest')
    mysaveas(pathname,'quantity_of_interest',formats,renderer);
    mymatlab2tikz(pathname,'quantity_of_interest.tex');
    
    Ioutput = integrate(boutput);
    fprintf('quantity of interest = %e\n',Ioutput);
end
