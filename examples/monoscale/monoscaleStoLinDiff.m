%% Monoscale stochastic linear diffusion problem %%
%%-----------------------------------------------%%

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

Dim = 2; % space dimension Dim = 2, 3
filename = ['linDiff_' num2str(Dim) 'D'];
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
    if Dim==2
        D = DOMAIN(2,[0.0,0.0],[1.0,1.0]);
    elseif Dim==3
        D = DOMAIN(3,[0.0,0.0,0.0],[1.0,1.0,1.0]);
    end
    
    nbelem = repmat(20,1,Dim);
    pb.S = build_model(D,'nbelem',nbelem);
    % cl = 0.05;
    % pb.S = build_model(D,'cl',cl,'filename',fullfile(pathname,'gmsh_domain'));
    
    %% Random variables
    d = 1; % parametric dimension
    v = UniformRandomVariable(0,1);
    rv = RandomVector(v,d);
    
    %% Materials
    % IntegrationRule
    p = 1;
    % basis = PolynomialFunctionalBasis(orthonormalPolynomials(v),0:p);
    % bases = FunctionalBases.duplicate(basis,d);
    bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
    bases = FunctionalBases(bases);
    H = FullTensorProductFunctionalBasis(bases);
    I = gaussIntegrationRule(v,2);
    I = I.tensorize(d);
    
    % Linear diffusion coefficient
    % K(xi) = 1 + xi
    fun = @(xi) 1 + xi(:,1);
    fun = UserDefinedFunction(fun,d);
    fun.evaluationAtMultiplePoints = true;
    
    K = H.projection(fun,I);
    
    % Material
    mat = FOUR_ISOT('k',K);
    mat = setnumber(mat,1);
    pb.S = setmaterial(pb.S,mat);
    
    %% Dirichlet boundary conditions
    pb.S = final(pb.S);
    pb.S = addcl(pb.S,[]);
    
    %% Stiffness matrices and sollicitation vectors
    if ~israndom(pb.S)
        pb.A = calc_rigi(pb.S);
    end
    
    % Source term
    f = 100;
    if ~israndom(f)
        pb.b = bodyload(pb.S,[],'QN',f);
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'pb','d','D','rv');
else
    load(fullfile(pathname,'problem.mat'),'pb','d','D','rv');
end

%% Adaptive sparse least-squares approximation
if solveProblem
    p = 50;
    bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
    bases = FunctionalBases(bases);
    
    s = AdaptiveSparseTensorAlgorithm();
    s.tol = 1e-12;
    s.tolStagnation = 1e-1;
    s.display = true;
    s.displayIterations = true;
    
    ls = LeastSquaresSolver();
    ls.errorEstimation = true;
    
    % if isanlsolver(pb.solver)
    %     u0 = solveSystem(calcOperator(funEval(pb,mean(rv))));
    %     fun = @(xi) solveSystem(calcOperator(funEval(pb,xi)),'inittype',u0);
    % else
    fun = @(xi) solveSystem(calcOperator(funEval(pb,xi)));
    % end
    fun = UserDefinedFunction(fun,d,getnbddlfree(pb.S));
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
fprintf('basis dimension = %d\n',cardinal(u.basis))
fprintf('order = [ %s ]\n',num2str(max(u.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(u.basis.indices.array)
fprintf('nb samples = %d\n',size(y,1))
fprintf('CV error = %d\n',norm(err))
fprintf('elapsed time = %f s\n',time)

%% Test
if testSolution
    N = 100;
    errL2 = testError(u,fun,N,rv);
    save(fullfile(pathname,'test.mat'),'errL2');
else
    load(fullfile(pathname,'test.mat'),'errL2');
end
fprintf('mean squared error = %d\n',errL2)

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
    % xtest = random(rv,nbSamples);
    % ytest = fun(xtest);
    % utest = u(xtest);
    % for i=1:nbSamples
    %     Stest = randomeval(pb.S,xtest(i,:));
    %     plotSolution(Stest,ytest(i,:)');
    %     mysaveas(pathname,['solution_ref_sample_' num2str(i)],formats,renderer);
    %     plotSolution(Stest,utest(i,:)');
    %     mysaveas(pathname,['solution_sample_' num2str(i)],formats,renderer);
    % end
end

myparallel('stop');
