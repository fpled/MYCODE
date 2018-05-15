%% Multiscale stochastic linear diffusion problem with n circular inclusions - Anisotropic case %%
%%----------------------------------------------------------------------------------------------%%
% [Beck, Nobile, Tamellini, Tempone, 2011,2014]

% clc
clearvars
close all
% rng('default');
myparallel('start');

%% Input data
setProblem = true;
directSolver = true;
iterativeSolver = true;
displaySolution = true;

n = 4; % number of inclusions n = 4
filename = ['linDiff' num2str(n) 'CircInclusionsAniso'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','multiscaleSto',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    % Global
    glob = Global();
    globOut = GlobalOutside();
    
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
    glob.S = gmshdomainwithinclusion(D,B,cl,cl,fullfile(pathname,['gmsh_circular_' num2str(n) '_inclusions']));
    
    % Patches
    patches = Patches(n);
    
    D_patch = cell(1,n);
    for k=1:n
        D_patch{k} = B{k};
        patches.patches{k}.S = keepgroupelem(glob.S,2*k);
        patches.patches{k}.S = removenodewithoutelem(patches.patches{k}.S);
        patches.patches{k}.S = keepeleminnode(patches.patches{k}.S);
    end
    
    % Partition of global mesh
    glob = partition(glob,patches);
    
    %% Random variables
    d = n; % parametric dimension
    v = UniformRandomVariable(-0.99,0);
    rv = RandomVector(v,d);
    
    %% Materials
    % Linear diffusion coefficient
    K_out = 1;
    K_patch = cell(1,n);
    K_in = cell(1,n);
    
    % IntegrationRule
    p = 1;
    basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
    bases = FunctionalBases.duplicate(basis,d);
    vb = basis.basis.randomVariable;
    rvb = getRandomVector(bases);
    H = FullTensorProductFunctionalBasis(bases);
    I = gaussIntegrationRule(vb,2);
    I = I.tensorize(d);
    
    g = [1 0.9 0.75 0.6];
    for k=1:n
        % K_patch(xi) = 1 + g * xi
        % K_in(x)     = 1
        fun = @(xi) 1 + g(k) * xi(:,k);
        funtr = @(xi) fun(transfer(rvb,rv,xi));
        fun = MultiVariateFunction(funtr,d);
        fun.evaluationAtMultiplePoints = true;
        
        K_patch{k} = H.projection(fun,I);
        K_in{k} = 1;
    end
    
    % Complementary subdomain
    mat_out = FOUR_ISOT('k',K_out);
    mat_out = setnumber(mat_out,0);
    glob.S = setmaterial(glob.S,mat_out,getnumgroupelemwithparam(glob.S,'partition',0));
    
    % Patches
    mat_patch = MATERIALS();
    for k=1:n
        mat_patch{k} = FOUR_ISOT('k',K_patch{k});
        mat_patch{k} = setnumber(mat_patch{k},k);
        patches.patches{k}.S = setmaterial(patches.patches{k}.S,mat_patch{k});
    end
    
    % Fictitious patches
    mat_in = MATERIALS();
    for k=1:n
        mat_in{k} = FOUR_ISOT('k',K_in{k});
        mat_in{k} = setnumber(mat_in{k},k);
        glob.S = setmaterial(glob.S,mat_in{k},getnumgroupelemwithparam(glob.S,'partition',k));
    end
    
    %% Dirichlet boundary conditions
    % Global
    glob.S = final(glob.S);
    glob.S = addcl(glob.S,[]);
    glob.S_out = getfinalmodelpart(glob.S,0);
    % S_in = cell(1,n);
    % for k=1:n
    %     S_in{k} = getfinalmodelpart(glob.S,k);
    % end
    
    % Complementary subdomain
    globOut.S = glob.S_out;
    
    % Patches
    for k=1:n
        patches.patches{k}.S = final(patches.patches{k}.S);
    end
    
    % Interfaces
    interfaces = Interfaces(patches);
    
    %% Stiffness matrices and sollicitation vectors
    % Source term
    f = 100;
        
    % Global
    glob.A = calc_rigi(glob.S);
    for k=1:n
        glob.A_in{k} = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
    end
    glob.b_out = bodyload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),[],'QN',f);
    
    % Complementary subdomain
    globOut.A = calc_rigi(globOut.S);
    globOut.b = bodyload(globOut.S,[],'QN',f);
    
    % Patches
    for k=1:n
        if ~israndom(patches.patches{k}.S)
            patches.patches{k}.A = calc_rigi(patches.patches{k}.S);
        end
        if ~israndom(f)
            patches.patches{k}.b = bodyload(patches.patches{k}.S,[],'QN',f);
        end
    end
    
    %% Mass matrices
    for k=1:n
        interfaces.interfaces{k}.M = calc_massgeom(interfaces.interfaces{k}.S);
    end
    
    %% Projection operators
    glob.P_out = calcProjection(glob);
    for k=1:n
        interfaces.interfaces{k}.P_glob = calcProjection(glob,interfaces.interfaces{k});
        interfaces.interfaces{k}.P_globOut = interfaces.interfaces{k}.P_glob*glob.P_out';
        interfaces.interfaces{k}.P_patch = calcProjection(patches.patches{k},interfaces.interfaces{k});
    end
    
    %% Parameters for global and local problems
    % Global problem
    glob.increment = true;
    
    % Local problems
    for k=1:n
        patches.patches{k}.changeOfVariable = false;
        patches.patches{k}.increment = true;
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'glob','globOut','patches','interfaces','D','D_patch');
else
    load(fullfile(pathname,'problem.mat'),'glob','globOut','patches','interfaces','D','D_patch');
end

%% Direct solver
if directSolver
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
    
    DS = DirectSolver();
    DS.changeOfVariable = false;
    DS.display = true;
    
    [U_ref,w_ref,lambda_ref,output_ref] = DS.solveRandom(globOut,patches,interfaces,s,bases,ls,rv);
    save(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','output_ref','s','bases','ls','rv');
else
    load(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','output_ref','s','bases','ls','rv');
end

%% Outputs
fprintf('\n')
fprintf('spatial dimension = %d for U_ref\n',U_ref.sz)
for k=1:n
    fprintf('                  = %d for w_ref{%u}\n',w_ref{k}.sz,k)
    fprintf('                  = %d for lambda_ref{%u}\n',lambda_ref{k}.sz,k)
end
fprintf('parametric dimension = %d\n',ndims(U_ref.basis))
fprintf('basis dimension = %d for U_ref\n',numel(U_ref.basis))
for k=1:n
    fprintf('                = %d for w_ref{%u}\n',numel(w_ref{k}.basis),k)
    fprintf('                = %d for lambda_ref{%u}\n',numel(lambda_ref{k}.basis),k)
end
fprintf('order = [ %s ] for U_ref\n',num2str(max(U_ref.basis.indices.array)))
for k=1:n
    fprintf('      = [ %s ] for w_ref{%u}\n',num2str(max(w_ref{k}.basis.indices.array)),k)
    fprintf('      = [ %s ] for lambda_ref{%u}\n',num2str(max(lambda_ref{k}.basis.indices.array)),k)
end
% fprintf('multi-index set for U_ref = \n')
% disp(num2str(U_ref.basis.indices.array))
% for k=1:n
%     fprintf('multi-index set for w_ref{%u} = \n',k)
%     disp(num2str(w_ref{k}.basis.indices.array))
%     fprintf('multi-index set for lambda_ref{%u} = \n',k)
%     disp(num2str(lambda_ref{k}.basis.indices.array))
% end
fprintf('nb samples = %d\n',output_ref.nbSamples)
fprintf('CV error = %d for U_ref\n',norm(output_ref.CVErrorGlobalSolution))
for k=1:n
    fprintf('         = %d for w_ref{%u}\n',norm(output_ref.CVErrorLocalSolution{k}),k)
    fprintf('         = %d for lambda_ref{%u}\n',norm(output_ref.CVErrorLagrangeMultiplier{k}),k)
end
fprintf('elapsed time = %f s\n',output_ref.time)

%% Global-local Iterative solver
if iterativeSolver
    s.tol = 1e-2;
    s.tolStagnation = 1e-1;
    s.display = true;
    s.displayIterations = false;
    
    IS = IterativeSolver();
    IS.maxIterations = 20;
    IS.tolerance = eps;
    IS.relaxation = 'Aitken';
    IS.updateRelaxationParameter = true;
    IS.errorCriterion = 'reference';
    IS.referenceSolution = {U_ref,w_ref,lambda_ref};
    IS.display = true;
    IS.displayIterations = true;

    [U,w,lambda,output] = IS.solveRandom(glob,patches,interfaces,s,bases,ls,rv);
    save(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
else
    load(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
end

%% Outputs
fprintf('\n')
fprintf('spatial dimension = %d for U\n',U.sz)
for k=1:n
    fprintf('                  = %d for w{%u}\n',w{k}.sz,k)
    fprintf('                  = %d for lambda{%u}\n',lambda{k}.sz,k)
end
fprintf('parametric dimension = %d\n',ndims(U.basis))
fprintf('basis dimension = %d for U\n',numel(U.basis))
for k=1:n
    fprintf('                = %d for w{%u}\n',numel(w{k}.basis),k)
    fprintf('                = %d for lambda{%u}\n',numel(lambda{k}.basis),k)
end
fprintf('order = [ %s ] for U\n',num2str(max(U.basis.indices.array)))
for k=1:n
    fprintf('      = [ %s ] for w{%u}\n',num2str(max(w{k}.basis.indices.array)),k)
    fprintf('      = [ %s ] for lambda{%u}\n',num2str(max(lambda{k}.basis.indices.array)),k)
end
% fprintf('multi-index set for U = \n')
% disp(num2str(U.basis.indices.array))
% for k=1:n
%     fprintf('multi-index set for w{%u} = \n',k)
%     disp(num2str(w{k}.basis.indices.array))
%     fprintf('multi-index set for lambda{%u} = \n',k)
%     disp(num2str(lambda{k}.basis.indices.array))
% end
fprintf('elapsed time = %f s\n',output.totalTime)

%% Display
if displaySolution
    %% Display domains and meshes
    plotDomain(D,D_patch);
    mysaveas(pathname,'domain_global_patches',formats,renderer);
    mymatlab2tikz(pathname,'domain_global_patches.tex');
    
    % plotPartition(glob,'legend',false);
    % mysaveas(pathname,'mesh_partition',formats,renderer);
    
    plotModel(glob,patches,'legend',false);
    mysaveas(pathname,'mesh_global_patches',formats,renderer);
    
    % plotModel(glob);
    % plotModel(patches);
    % plotModel(interfaces);
    
    %% Display evolutions of error indicator, stagnation indicator, CPU time, relaxation parameter w.r.t. number of iterations
    plotError(output);
    mysaveas(pathname,'error','fig');
    mymatlab2tikz(pathname,'error.tex');
    
    plotStagnation(output);
    mysaveas(pathname,'stagnation','fig');
    mymatlab2tikz(pathname,'stagnation.tex');
    
    plotErrorGlobalSolution(output);
    mysaveas(pathname,'error_global_solution','fig');
    mymatlab2tikz(pathname,'error_global_solution.tex');
    
    plotStagnationGlobalSolution(output);
    mysaveas(pathname,'stagnation_global_solution','fig');
    mymatlab2tikz(pathname,'stagnation_global_solution.tex');
    
    plotCPUTime(output,'legend',false);
    mysaveas(pathname,'cpu_time','fig');
    mymatlab2tikz(pathname,'cpu_time.tex');
    
    plotRelaxationParameter(output,'legend',false);
    mysaveas(pathname,'relaxation_parameter','fig');
    mymatlab2tikz(pathname,'relaxation_parameter.tex');
    
    plotNbSamples(output);
    mysaveas(pathname,'nb_samples','fig');
    mymatlab2tikz(pathname,'nb_samples.tex');
    
    plotDimStochasticBasis(output);
    mysaveas(pathname,'dim_stochastic_basis','fig');
    mymatlab2tikz(pathname,'dim_stochastic_basis.tex');
    
    plotCVError(output);
    mysaveas(pathname,'cv_error','fig');
    mymatlab2tikz(pathname,'cv_error.tex');
    
    %% Display multi-index sets
    for i=1:2:d
        plotMultiIndexSet(U,'dim',[i i+1],'legend',false);
        mysaveas(pathname,['multi_index_set_global_solution_dim_' num2str(i) '_' num2str(i+1)],'fig');
        mymatlab2tikz(pathname,['multi_index_set_global_solution_dim_' num2str(i) '_' num2str(i+1) '.tex']);
    end
    
    for k=1:n
        for i=1:2:d
            plotMultiIndexSet(w{k},'dim',[i i+1],'legend',false);
            mysaveas(pathname,['multi_index_set_local_solution_' num2str(k) '_dim_' num2str(i) '_' num2str(i+1)],'fig');
            mymatlab2tikz(pathname,['multi_index_set_local_solution_' num2str(k) '_dim_' num2str(i) '_' num2str(i+1) '.tex']);
            
            plotMultiIndexSet(lambda{k},'dim',[i i+1],'legend',false);
            mysaveas(pathname,['multi_index_set_Lagrange_multiplier_' num2str(k) '_dim_' num2str(i) '_' num2str(i+1)],'fig');
            mymatlab2tikz(pathname,['multi_index_set_Lagrange_multiplier_' num2str(k) '_dim_' num2str(i) '_' num2str(i+1) '.tex']);
        end
    end
    
    %% Display statistical outputs
    % plotStatsAllSolutions(glob,patches,interfaces,U,w,lambda);
    
    plotMeanGlobalSolution(glob,U);
    mysaveas(pathname,'mean_global_solution',formats,renderer);
    
    % plotMeanLocalSolution(patches,w);
    % mysaveas(pathname,'mean_local_solution',formats,renderer);
    
    % plotMeanLagrangeMultiplier(interfaces,lambda);
    % mysaveas(pathname,'mean_Lagrange_multiplier',formats,renderer);
    
    plotMeanMultiscaleSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'mean_multiscale_solution',formats,renderer);
    
    plotMeanGlobalLocalSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'mean_global_local_solution',formats,renderer);
    plotMeanGlobalLocalSolution(glob,patches,interfaces,U,w,'view3',true);
    mysaveas(pathname,'mean_global_local_solution_view3',formats,renderer);
    
    plotVarianceGlobalSolution(glob,U);
    mysaveas(pathname,'var_global_solution',formats,renderer);
    
    % plotVarianceLocalSolution(patches,w);
    % mysaveas(pathname,'var_local_solution',formats,renderer);
    
    % plotVarianceLagrangeMultiplier(interfaces,lambda);
    % mysaveas(pathname,'var_Lagrange_multiplier',formats,renderer);
    
    plotVarianceMultiscaleSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'var_multiscale_solution',formats,renderer);
    
    plotVarianceGlobalLocalSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'var_global_local_solution',formats,renderer);
    plotVarianceGlobalLocalSolution(glob,patches,interfaces,U,w,'view3',true);
    mysaveas(pathname,'var_global_local_solution_view3',formats,renderer);
    
    plotStdGlobalSolution(glob,U);
    mysaveas(pathname,'std_global_solution',formats,renderer);
    
    % plotStdLocalSolution(patches,w);
    % mysaveas(pathname,'std_local_solution',formats,renderer);
    
    % plotStdLagrangeMultiplier(interfaces,lambda);
    % mysaveas(pathname,'std_Lagrange_multiplier',formats,renderer);
    
    plotStdMultiscaleSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'std_multiscale_solution',formats,renderer);
    
    plotStdGlobalLocalSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'std_global_local_solution',formats,renderer);
    plotStdGlobalLocalSolution(glob,patches,interfaces,U,w,'view3',true);
    mysaveas(pathname,'std_global_local_solution_view3',formats,renderer);
    
    d = ndims(U.basis);
    for i=1:d
        plotSobolIndicesMultiscaleSolution(glob,patches,interfaces,U,w,i);
        mysaveas(pathname,['sobol_indices_multiscale_solution_var_' num2str(i)],formats,renderer);
        
        plotSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,U,w,i);
        mysaveas(pathname,['sensitivity_indices_multiscale_solution_var_' num2str(i)],formats,renderer);
    end
    
    %% Quantities of interest
    % I1: mean value of U over square subdomain B{9}
    % I2: mean value of u over domain D
    % I3: mean value of the gradient of u over domain D
    
    %% Display random evaluations
    % nbsamples = 3;
    % for i=1:nbsamples
    %     xi = random(rv,1);
    %     U_xi = U(xi);
    %     w_xi = cellfun(@(x) x(xi),w,'UniformOutput',false);
    %     lambda_xi = cellfun(@(x) x(xi),lambda,'UniformOutput',false);
    %     % plotAllSolutions(glob,patches.eval(xi),interfaces,U_xi',cellfun(@(x) x',w_xi,'UniformOutput',false),cellfun(@(x) x',lambda_xi,'UniformOutput',false));
    %     plotGlobalSolution(glob,U_xi');
    %     % plotLocalSolution(patches,cellfun(@(x) x',w_xi,'UniformOutput',false));
    %     % plotLagrangeMultiplier(interfaces,cellfun(@(x) x',lambda_xi,'UniformOutput',false));
    %     plotMultiscaleSolution(glob,patches.eval(xi),interfaces,U_xi',cellfun(@(x) x',w_xi,'UniformOutput',false));
    % end
end

myparallel('stop');
