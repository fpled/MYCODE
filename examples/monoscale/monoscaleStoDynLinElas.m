%% Monoscale stochastic linear elasticity dynamic problem %%
%%--------------------------------------------------------%%

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

filename = 'dynLinElas';
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
    D = DOMAIN(2,[0.0,0.0],[2.0,0.4]);
    
    elemtype = 'TRI3';
    % elemtype = 'QUA4';
    % option = 'DEFO'; % plane strain
    option = 'CONT'; % plane stress
    nbelem = [60,12];
    pb.S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'option',option);
%     cl = 0.03;
%     pb.S = build_model(D,'cl',cl,'elemtype',elemtype,'option',option,'filename',fullfile(pathname,'gmsh_domain'));
    
    %% Random variables
    d = 2; % parametric dimension
    v = UniformRandomVariable(0,1);
    rv = RandomVector(v,d);
    
    %% Materials
    % Poisson ratio
    NU = 0;
    % Thickness
    DIM3 = 1;
    
    % Density
    % RHO(xi) = 1 + xi
    p = 1;
    basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
    bases = FunctionalBases.duplicate(basis,d);
    rvb = getRandomVector(bases);
    H = FullTensorProductFunctionalBasis(bases);
    I = gaussIntegrationRule(v,2);
    I = I.tensorize(d);
    
    fun = @(xi) 1 + xi(:,1);
    funtr = @(xi) fun(transfer(rvb,rv,xi));
    fun = MultiVariateFunction(funtr,d);
    fun.evaluationAtMultiplePoints = true;
    
    RHO = H.projection(fun,I);
    
    % Young modulus
    % E(xi) = 1 + xi
    fun = @(xi) 1 + xi(:,2);
    funtr = @(xi) fun(transfer(rvb,rv,xi));
    fun = MultiVariateFunction(funtr,d);
    fun.evaluationAtMultiplePoints = true;
    
    E = H.projection(fun,I);
    
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
    
    pb.timeSolver = NEWMARKSOLVER(T,'alpha',0.05,'display',false);
    % pb.timeSolver = DGTIMESOLVER(T,1,'outputsplit',true,'display',false,'lu',true);
    
    tc = get(T,'t1')/6;
    pb.loadFunction = @(N) rampe(N,t0,tc);
    % pb.loadFunction = @(N) dirac(N,t0,tc);
    % pb.loadFunction = @(N) one(N);
    
    %% Mass, stiffness and damping matrices and sollicitation vectors
    if ~israndom(pb.S)
        pb.M = calc_mass(pb.S);
        pb.A = calc_rigi(pb.S);
    end
    
    b = surfload(pb.S,L2,'FX',-1);
    pb.b = b*pb.loadFunction(pb.timeSolver);
    
    save(fullfile(pathname,'problem.mat'),'pb','D','L1','L2');
else
    load(fullfile(pathname,'problem.mat'),'pb','D','L1','L2');
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
    
    fun = @(xi) solveSystem(calcOperator(funEval(pb,xi)));
    fun = MultiVariateFunction(fun,d,3);
    fun.evaluationAtMultiplePoints = false;
    t = tic;
    [Ut,errt,~,yt] = s.leastSquaresCell(fun,bases,ls,rv);
    time = toc(t);
    
    save(fullfile(pathname,'solution.mat'),'Ut','errt','yt','fun','time');
else
    load(fullfile(pathname,'solution.mat'),'Ut','errt','yt','fun','time');
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

ut = Ut{1}; errut = errt{1}; yut = yt{1};
vt = Ut{2}; errvt = errt{2}; yvt = yt{2};
at = Ut{3}; errat = errt{3}; yat = yt{3};

fprintf('\n');
fprintf('parametric dimension = %d for u\n',ndims(ut.basis))
fprintf('                     = %d for v\n',ndims(vt.basis))
fprintf('                     = %d for a\n',ndims(at.basis))
fprintf('basis dimension = %d for u\n',numel(ut.basis))
fprintf('                = %d for v\n',numel(vt.basis))
fprintf('                = %d for a\n',numel(at.basis))
fprintf('order = [ %s ] for u\n',num2str(max(ut.basis.indices.array)))
fprintf('      = [ %s ] for v\n',num2str(max(vt.basis.indices.array)))
fprintf('      = [ %s ] for a\n',num2str(max(at.basis.indices.array)))
% fprintf('multi-index set for u = \n')
% disp(ut.basis.indices.array)
% fprintf('multi-index set for v = \n')
% disp(vt.basis.indices.array)
% fprintf('multi-index set for a = \n')
% disp(at.basis.indices.array)
fprintf('nb samples = %d\n',size(yut,1))
fprintf('CV error = %d for u\n',norm(errut))
fprintf('         = %d for v\n',norm(errvt))
fprintf('         = %d for a\n',norm(errat))
fprintf('elapsed time = %f s\n',time)

%% Test
if testSolution
    Ntest = 100;
    [errttest,xttest,Uttest,yttest] = computeTestErrorCell(Ut,fun,Ntest);
    save(fullfile(pathname,'test.mat'),'utest','errtest','xtest','ytest');
    save(fullfile(pathname,'testt.mat'),'Uttest','errttest','xttest','yttest');
else
    load(fullfile(pathname,'test.mat'),'utest','errtest','xtest','ytest');
    load(fullfile(pathname,'testt.mat'),'Uttest','errttest','xttest','yttest');
end
fprintf('\n');
fprintf('Stationary solution\n');
fprintf('test error = %d\n',errtest)

uttest = Uttest{1}; erruttest = errttest{1}; yuttest = yttest{1};
vttest = Uttest{2}; errvttest = errttest{2}; yvttest = yttest{2};
attest = Uttest{3}; errattest = errttest{3}; yattest = yttest{3};

fprintf('\n');
fprintf('Transient solution\n');
fprintf('test error = %d for u\n',erruttest)
fprintf('test error = %d for v\n',errvttest)
fprintf('test error = %d for a\n',errattest)

%% Display
if displaySolution
    %% Display domains and meshes
%     plotDomain(D,'legend',false);
%     mysaveas(pathname,'domain',formats,renderer);
%     mymatlab2tikz(pathname,'domain.tex');
    
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
    
%     plotModel(pb.S,'legend',false);
%     mysaveas(pathname,'mesh',formats,renderer);
    
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
        evolMean(pb.S,ut,'displ',i,'filename',['evol_solution_' num2str(i)],'pathname',pathname);
        evolMean(pb.S,ut,'displ',i,'view3',true,'filename',['evol_solution_' num2str(i) '_view3'],'pathname',pathname);
        
        evolMean(pb.S,vt,'displ',i,'filename',['evol_velocity_' num2str(i)],'pathname',pathname);
        evolMean(pb.S,vt,'displ',i,'view3',true,'filename',['evol_velocity_' num2str(i) '_view3'],'pathname',pathname);
        
        evolMean(pb.S,at,'displ',i,'filename',['evol_acceleration_' num2str(i)],'pathname',pathname);
        evolMean(pb.S,at,'displ',i,'view3',true,'filename',['evol_acceleration_' num2str(i) '_view3'],'pathname',pathname);
    % end
    
    % for i=1:3
    %     evolMean(pb.S,ut,'epsilon',i,'filename',['evol_epsilon_' num2str(i)],'pathname',pathname);
    %     evolMean(pb.S,ut,'sigma',i,'filename',['evol_sigma_' num2str(i)],'pathname',pathname);
    % end
    
    % evolMean(pb.S,ut,'epsilon','mises','filename','evol_epsilon_von_mises','pathname',pathname);
    % evolMean(pb.S,ut,'sigma','mises','filename','evol_sigma_von_mises','pathname',pathname);
end
