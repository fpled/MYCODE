%% Monoscale stochastic linear diffusion problem with circular inclusions - Isotropic case %%
%%-----------------------------------------------------------------------------------------%%
% [Beck, Nobile, Tamellini, Tempone, 2011,2014], [Chkifa, Cohen, Migliorati, Tempone, 2015]

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

filename = 'linDiffCircInclusionsIso';
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
    
    r = 0.13;
    B = cell(1,9);
    B{1} = CIRCLE(0.2,0.2,r);
    B{2} = CIRCLE(0.2,0.5,r);
    B{3} = CIRCLE(0.2,0.8,r);
    B{4} = CIRCLE(0.5,0.8,r);
    B{5} = CIRCLE(0.8,0.8,r);
    B{6} = CIRCLE(0.8,0.5,r);
    B{7} = CIRCLE(0.8,0.2,r);
    B{8} = CIRCLE(0.5,0.2,r);
    B{9} = DOMAIN(2,[0.4,0.4],[0.6,0.6]);
    
    cl = 0.02;
    pb.S = gmshdomainwithinclusion(D,B,cl,cl,fullfile(pathname,'gmsh_circular_inclusions'));
    
    %% Random variables
    d = 8; % parametric dimension d = 2, 4, 8
    v = UniformRandomVariable(-0.99,-0.2);
    rv = RandomVector(v,d);
    
    %% Materials
    % Deterministic subdomains
    % Linear diffusion coefficient
    K_det = 1;
    mat_det = FOUR_ISOT('k',K_det);
    mat_det = setnumber(mat_det,0);
    k = [0 9];
    pb.S = setmaterial(pb.S,mat_det,k+1);
    
    % Stochastic subdomains
    % IntegrationRule
    p = 1;
    basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
    bases = FunctionalBases.duplicate(basis,d);
    vb = basis.basis.randomVariable;
    rvb = getRandomVector(bases);
    H = FullTensorProductFunctionalBasis(bases);
    I = gaussIntegrationRule(vb,2);
    I = I.tensorize(d);
    
    mat_sto = MATERIALS();
    for i=1:d
        % Linear diffusion coefficient
        % K_sto(xi) = 1 + xi
        fun = @(xi) 1 + xi(:,i);
        funtr = @(xi) fun(transfer(rvb,rv,xi));
        fun = MultiVariateFunction(funtr,d);
        fun.evaluationAtMultiplePoints = true;
        
        K_sto = H.projection(fun,I);
        
        mat_sto{i} = FOUR_ISOT('k',K_sto);
        mat_sto{i} = setnumber(mat_sto{i},i);
    end
    switch d
        case 2
            k = [2 4 6 8];
            pb.S = setmaterial(pb.S,mat_sto{1},k+1);
            k = [1 3 5 7];
            pb.S = setmaterial(pb.S,mat_sto{2},k+1);
        case 4
            k = [1 5];
            pb.S = setmaterial(pb.S,mat_sto{1},k+1);
            k = [2 6];
            pb.S = setmaterial(pb.S,mat_sto{2},k+1);
            k = [3 7];
            pb.S = setmaterial(pb.S,mat_sto{3},k+1);
            k = [4 8];
            pb.S = setmaterial(pb.S,mat_sto{4},k+1);
        case 8
            for i=1:d
                pb.S = setmaterial(pb.S,mat_sto{i},i+1);
            end
        otherwise
            error('Wrong number of variables')
    end
    
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
        pb.b = bodyload(keepgroupelem(pb.S,10),[],'QN',f);
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'pb','d','D','B');
else
    load(fullfile(pathname,'problem.mat'),'pb','d','D','B');
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
    s.tol = 1e-2;
    s.tolStagnation = 5e-2;
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
    plotDomain(D,B);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    plotPartition(pb.S,'legend',false);
    mysaveas(pathname,'mesh_partition',formats,renderer);
    
    plotModel(pb.S,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    %% Display multi-index set
    for i=1:2:d
        plotMultiIndexSet(u,'dim',[i i+1],'legend',false);
        mysaveas(pathname,['multi_index_set_dim_' num2str(i) '_' num2str(i+1)],'fig');
        mymatlab2tikz(pathname,['multi_index_set_dim_' num2str(i) '_' num2str(i+1) '.tex']);
    end
    
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
        plotSobolIndices(pb.S,u,i);
        mysaveas(pathname,['sobol_indices_solution_var_' num2str(i)],formats,renderer);
        
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
