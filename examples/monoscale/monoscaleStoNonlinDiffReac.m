%% Monoscale stochastic non-linear diffusion-reaction problem %%
%%------------------------------------------------------------%%

% clc
clearvars
close all
% rng('default');
myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;
testSolution = true;

filename = 'nonlinDiffReac';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','monoscaleSto',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    D = DOMAIN(2,[0.0,0.0],[1.0,1.0]);
    
    nbelem = [20,20];
    pb.S = build_model(D,'nbelem',nbelem);
    % cl = 0.05;
    % problem.S = build_model(D,'cl',cl,'filename',fullfile(pathname,'gmsh_domain'));
    
    %% Random variables
    d = 2; % parametric dimension
    v = UniformRandomVariable(0,1);
    rv = RandomVector(v,d);
    
    %% Materials
    % IntegrationRule
    p = 1;
    basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
    bases = FunctionalBases.duplicate(basis,d);
    vb = basis.basis.randomVariable;
    rvb = getRandomVector(bases);
    H = FullTensorProductFunctionalBasis(bases);
    I = gaussIntegrationRule(vb,2);
    I = I.tensorize(d);
    
    % Linear diffusion coefficient
    % K(xi) = 1 + xi
    fun = @(xi) 1 + xi(:,1);
    funtr = @(xi) fun(transfer(rvb,rv,xi));
    fun = MultiVariateFunction(funtr,d);
    fun.evaluationAtMultiplePoints = true;
    
    K = H.projection(fun,I);
    
    % Nonlinear reaction parameter
    % R(xi) = xi
    fun = @(x) x(:,2);
    funtr = @(x) fun(transfer(rvb,rv,x));
    fun = MultiVariateFunction(funtr,d);
    fun.evaluationAtMultiplePoints = true;
    
    R = H.projection(fun,I);
    
    % Material
    mat = FOUR_ISOT('k',K,'r3',R);
    mat = setnumber(mat,1);
    pb.S = setmaterial(pb.S,mat);
    
    %% Dirichlet boundary conditions
    pb.S = final(pb.S);
    pb.S = addcl(pb.S,[]);
    
    %% Stiffness matrices and sollicitation vectors
    if ~israndom(pb.S)
        pb.A = @(u) calc_fint(pb.S,u);
        pb.Atang = @(u) calc_rigitang(pb.S,u);
    end
    
    % Source term
    f = 100;
    if ~israndom(f)
        pb.b = bodyload(pb.S,[],'QN',f);
    end
    
    % Solver
    pb.solver = NEWTONSOLVER('type','tangent','increment',true,...
        'maxiter',100,'tol',1e-12,'display',false,'stopini',true);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'pb','d','D');
else
    load(fullfile(pathname,'problem.mat'),'pb','d','D');
end

%% Adaptive sparse least-squares approximation
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
    s.display = true;
    s.displayIterations = true;
    
    ls = LeastSquaresSolver();
    ls.regularization = false;
    % ls.regularizationType = 'l1';
    ls.errorEstimation = true;
    % ls.errorEstimationType = 'leaveout';
    % ls.errorEstimationOptions.correction = true;
    
    % if isanlsolver(pb.solver)
    %     u0 = solveSystem(calcOperator(funEval(pb,mean(rv))));
    %     fun = @(xi) solveSystem(calcOperator(funEval(pb,xi)),'inittype',u0);
    % else
    fun = @(xi) solveSystem(calcOperator(funEval(pb,xi)));
    % end
    fun = MultiVariateFunction(fun,d,getnbddlfree(pb.S));
    fun.evaluationAtMultiplePoints = false;
    
    t = tic;
    [u,err,~,y] = s.leastSquares(fun,bases,ls,rv);
    time = toc(t);
    save(fullfile(pathname,'solution.mat'),'u','err','y','fun','time');
else
    load(fullfile(pathname,'solution.mat'),'u','err','y','fun','time');
end

%% Outputs
fprintf('\n')
fprintf('spatial dimension = %d\n',u.sz)
fprintf('parametric dimension = %d\n',ndims(u.basis))
fprintf('basis dimension = %d\n',numel(u.basis))
fprintf('order = [ %s ]\n',num2str(max(u.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(u.basis.indices.array)
fprintf('nb samples = %d\n',size(y,1))
fprintf('CV error = %d\n',norm(err))
fprintf('elapsed time = %f s\n',time)

%% Test
if testSolution
    Ntest = 100;
    [errtest,xtest,utest,ytest] = computeTestError(u,fun,Ntest);
    save(fullfile(pathname,'test.mat'),'utest','errtest','xtest','ytest');
else
    load(fullfile(pathname,'test.mat'),'utest','errtest','xtest','ytest');
end
fprintf('test error = %d\n',errtest)

%% Display
if displaySolution
    %% Display domains and meshes
    plotDomain(D);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    % plotPartition(pb.S,'legend',false);
    % mysaveas(pathname,'mesh_partition',formats,renderer);
    
    plotModel(pb.S,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);

    %% Display multi-index set
    plotMultiIndexSet(u,'legend',false);
    mysaveas(pathname,'multi_index_set','fig');
    mymatlab2tikz(pathname,'multi_index_set.tex');
    
    %% Display statistical outputs
    % plotStats(pb.S,u);
    
    plotMean(pb.S,u);
    mysaveas(pathname,'mean_solution',formats,renderer);
    
    plotVariance(pb.S,u);
    mysaveas(pathname,'var_solution',formats,renderer);
    
    plotStd(pb.S,u);
    mysaveas(pathname,'std_solution',formats,renderer);
    
    d = ndims(u.basis);
    for i=1:d
        % plotSobolIndices(pb.S,u,i);
        % mysaveas(pathname,['sobol_indices_solution_var_' num2str(i)],formats,renderer);
        
        plotSensitivityIndices(pb.S,u,i);
        mysaveas(pathname,['sensitivity_indices_solution_var_' num2str(i)],formats,renderer);
    end
    
    %% Display random evaluations
    % nbSamples = 3;
    % for i=1:nbSamples
    %     Stest = randomeval(pb.S,xtest(i,:)');
    %     plotSolution(Stest,ytest(i,:)');
    %     mysaveas(pathname,['solution_ref_sample_' num2str(i)],formats,renderer);
    %     plotSolution(Stest,utest(i,:)');
    %     mysaveas(pathname,['solution_sample_' num2str(i)],formats,renderer);
    % end
end

myparallel('stop');
