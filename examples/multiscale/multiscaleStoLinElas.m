%% Multiscale stochastic linear elasticity problem %%
%%-------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');
myparallel('start');

%% Input data
setProblem = true;
directSolver = true;
iterativeSolver = true;
displaySolution = true;

n = 4; % number of patches n = 1, 2, 4
filename = ['linElas' num2str(n) 'Patches'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','multiscaleSto',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    % Global
    glob = Global();
    glob_out = GlobalOutside();
    
    D = DOMAIN(2,[0.0,0.0],[1.0,1.0]);
    
    option = 'DEFO'; % plane strain
    nbelem = [20,20];
    glob.S = build_model(D,'nbelem',nbelem,'option',option);
    % cl = 0.05;
    % glob.S = build_model(D,'cl',cl,'option',option,'filename',fullfile(pathname,'gmsh_domain'));
    
    % Patches
    patches = Patches(n);
    
    D_patch = cell(1,n);
    switch n
        case 1
            D_patch{1} = DOMAIN(2,[0.4,0.4],[0.6,0.6]);
        case 2
            D_patch{1} = DOMAIN(2,[0.1,0.1],[0.3,0.3]);
            D_patch{2} = DOMAIN(2,[0.7,0.7],[0.9,0.9]);
        case 4
            D_patch{1} = DOMAIN(2,[0.1,0.1],[0.3,0.3]);
            D_patch{2} = DOMAIN(2,[0.1,0.7],[0.3,0.9]);
            D_patch{3} = DOMAIN(2,[0.7,0.7],[0.9,0.9]);
            D_patch{4} = DOMAIN(2,[0.7,0.1],[0.9,0.3]);
        otherwise
            error('Wrong number of patches')
    end
    
    nbelem_patch = [40,40];
    for k=1:n
        patches.patches{k}.S = build_model(D_patch{k},'nbelem',nbelem_patch);
        patches.patches{k}.S = setoption(patches.patches{k}.S,option);
    end
    % cl_patch = 0.005;
    % for k=1:n
    %     patches.patches{k}.S = build_model(D_patch{k},'cl',cl_patch,'filename',fullfile(pathname,['gmsh_patch_' num2str(k)]));
    %     patches.patches{k}.S = setoption(patches.patches{k}.S,option);
    % end
    
    % Partition of global mesh
    glob = partition(glob,patches);
    
    %% Random variables
    d = n; % parametric dimension
    v = UniformRandomVariable(0,1);
    rv = RandomVector(v,d);
    
    %% Materials
    % Poisson ratio
    NU = 0.3;
    % Thickness
    DIM3 = 1;
    % Density
    RHO = 1;
    % Young modulus
    E_out = 1;
    E_patch = cell(1,n);
    E_in = cell(1,n);
    
    p = 1;
    basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
    bases = FunctionalBases.duplicate(basis,d);
    vb = basis.basis.randomVariable;
    rvb = getRandomVector(bases);
    H = FullTensorProductFunctionalBasis(bases);
    I = gaussIntegrationRule(vb,2);
    I = I.tensorize(d);
    
    for k=1:n
        patch = patches.patches{k};
        % E_patch(x,xi) = 1 + f(x) * xi
        % E_in(x)       = 1
        % with f(x) = 1 if ||x-c||_Inf < L
        %           = 0 if ||x-c||_Inf >= L
        L = norm(getsize(D_patch{k}),Inf)/4;
        c = getcenter(D_patch{k});
        f = @(x) distance(x,c,Inf)<L;
        fun = @(xi) ones(size(xi,1),patch.S.nbnode) + xi(:,k) * double(squeeze(f(patch.S.node)))';
        funtr = @(xi) fun(transfer(rvb,rv,xi));
        fun = MultiVariateFunction(funtr,d,patch.S.nbnode);
        fun.evaluationAtMultiplePoints = true;
        
        E_patch{k} = FENODEFIELD(H.projection(fun,I));
        E_in{k} = 1;
    end
    
    % Complementary subdomain
    mat_out = ELAS_ISOT('E',E_out,'NU',NU,'RHO',RHO,'DIM3',DIM3);
    mat_out = setnumber(mat_out,0);
    glob.S = setmaterial(glob.S,mat_out,getnumgroupelemwithparam(glob.S,'partition',0));
    
    % Patches
    mat_patch = MATERIALS();
    for k=1:n
        mat_patch{k} = ELAS_ISOT('E',E_patch{k},'NU',NU,'RHO',RHO,'DIM3',DIM3);
        mat_patch{k} = setnumber(mat_patch{k},k);
        patches.patches{k}.S = setmaterial(patches.patches{k}.S,mat_patch{k});
    end
    
    % Fictitious patches
    mat_in = MATERIALS();
    for k=1:n
        mat_in{k} = ELAS_ISOT('E',E_in{k},'NU',NU,'RHO',RHO,'DIM3',DIM3);
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
    glob_out.S = glob.S_out;
    
    % Patches
    for k=1:n
        patches.patches{k}.S = final(patches.patches{k}.S);
    end
    
    % Interfaces
    interfaces = Interfaces(patches);
    
    %% Stiffness matrices and sollicitation vectors
    % Traction force density
    f = [0;-100];
    
    % Global
    glob.A = calc_rigi(glob.S);
    for k=1:n
        glob.A_in{k} = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
    end
    glob.b_out = bodyload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),[],{'FX','FY'},f);
    
    % Complementary subdomain
    glob_out.A = calc_rigi(glob_out.S);
    glob_out.b = bodyload(glob_out.S,[],{'FX','FY'},f);
    
    % Patches
    for k=1:n
        if ~israndom(patches.patches{k}.S)
            patches.patches{k}.A = calc_rigi(patches.patches{k}.S);
        end
        if ~israndom(f)
            patches.patches{k}.b = bodyload(patches.patches{k}.S,[],{'FX','FY'},f);
        end
    end
    
    %% Mass matrices
    for k=1:n
        interfaces.interfaces{k}.M = calc_massgeom(interfaces.interfaces{k}.S);
    end
    
    %% Projection operators
    glob.P_out = calcProjection(glob);
    for k=1:n
        [interfaces.interfaces{k}.P_glob] = calcProjection(interfaces.interfaces{k},glob);
        [interfaces.interfaces{k}.P_glob_out,numnode] = calcProjection(interfaces.interfaces{k},glob_out);
        interfaces.interfaces{k}.P_patch = calcProjection(patches.patches{k},interfaces.interfaces{k});
        % plotProjectionOperator(glob,patches.patches{k},numnode);
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
    save(fullfile(pathname,'problem.mat'),'glob','glob_out','patches','interfaces','D','D_patch');
else
    load(fullfile(pathname,'problem.mat'),'glob','glob_out','patches','interfaces','D','D_patch');
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
    
    [U_ref,w_ref,lambda_ref,output_ref] = DS.solveRandom(glob_out,patches,interfaces,s,bases,ls,rv);
    save(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','output_ref');
else
    load(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','output_ref');
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
    s.tol = 1e-3;
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
    plotMultiIndexSet(U,'legend',false)
    mysaveas(pathname,'multi_index_set_global_solution','fig');
    mymatlab2tikz(pathname,'multi_index_set_global_solution.tex');
    
    for k=1:n
        plotMultiIndexSet(w{k},'legend',false)
        mysaveas(pathname,['multi_index_set_local_solution_' num2str(k)],'fig');
        mymatlab2tikz(pathname,['multi_index_set_local_solution_' num2str(k) '.tex']);
        
        plotMultiIndexSet(lambda{k},'legend',false)
        mysaveas(pathname,['multi_index_set_Lagrange_multiplier_' num2str(k)],'fig');
        mymatlab2tikz(pathname,['multi_index_set_Lagrange_multiplier_' num2str(k) '.tex']);
    end
    
    %% Display statistical outputs
    for i=1:2
        % plotStatsAllSolutions(glob,patches,interfaces,U,w,lambda,'displ',i);
        
        plotMeanGlobalSolution(glob,U,'displ',i);
        mysaveas(pathname,['mean_global_solution_' num2str(i)],formats,renderer);
        
        % plotMeanLocalSolution(patches,w,'displ',i);
        % mysaveas(pathname,['mean_local_solution_' num2str(i)],formats,renderer);
        
        % plotMeanLagrangeMultiplier(interfaces,lambda,'displ',i);
        % mysaveas(pathname,['mean_Lagrange_multiplier_' num2str(i)],formats,renderer);
        
        plotMeanMultiscaleSolution(glob,patches,interfaces,U,w,'displ',i);
        mysaveas(pathname,['mean_multiscale_solution_' num2str(i)],formats,renderer);
        
        plotMeanGlobalLocalSolution(glob,patches,interfaces,U,w,'displ',i);
        mysaveas(pathname,['mean_global_local_solution_' num2str(i)],formats,renderer);
        
        plotMeanGlobalLocalSolution(glob,patches,interfaces,U,w,'displ',i,'view3',true);
        mysaveas(pathname,['mean_global_local_solution_' num2str(i) '_surf'],formats,renderer);
        
        plotVarianceGlobalSolution(glob,U,'displ',i);
        mysaveas(pathname,['var_global_solution_' num2str(i)],formats,renderer);
        
        % plotVarianceLocalSolution(patches,w,'displ',i);
        % mysaveas(pathname,['var_local_solution_' num2str(i)],formats,renderer);
        
        % plotVarianceLagrangeMultiplier(interfaces,lambda,'displ',i);
        % mysaveas(pathname,['var_Lagrange_multiplier_' num2str(i)],formats,renderer);
        
        plotVarianceMultiscaleSolution(glob,patches,interfaces,U,w,'displ',i);
        mysaveas(pathname,['var_multiscale_solution_' num2str(i)],formats,renderer);
        
        plotVarianceGlobalLocalSolution(glob,patches,interfaces,U,w,'displ',i);
        mysaveas(pathname,['var_global_local_solution_' num2str(i)],formats,renderer);
        
        plotVarianceGlobalLocalSolution(glob,patches,interfaces,U,w,'displ',i,'view3',true);
        mysaveas(pathname,['var_global_local_solution_' num2str(i) '_surf'],formats,renderer);
        
        plotStdGlobalSolution(glob,U,'displ',i);
        mysaveas(pathname,['std_global_solution_' num2str(i)],formats,renderer);
        
        % plotStdLocalSolution(patches,w,'displ',i);
        % mysaveas(pathname,['std_local_solution_' num2str(i)],formats,renderer);
        
        % plotStdLagrangeMultiplier(interfaces,lambda,'displ',i);
        % mysaveas(pathname,['std_Lagrange_multiplier_' num2str(i)],formats,renderer);
        
        plotStdMultiscaleSolution(glob,patches,interfaces,U,w,'displ',i);
        mysaveas(pathname,['std_multiscale_solution_' num2str(i)],formats,renderer);
        
        plotStdGlobalLocalSolution(glob,patches,interfaces,U,w,'displ',i);
        mysaveas(pathname,['std_global_local_solution_' num2str(i)],formats,renderer);
        
        plotStdGlobalLocalSolution(glob,patches,interfaces,U,w,'displ',i,'view3',true);
        mysaveas(pathname,['std_global_local_solution_' num2str(i) '_surf'],formats,renderer);
        
        for j=1:d
            plotSobolIndicesMultiscaleSolution(glob,patches,interfaces,U,w,j,'displ',i);
            mysaveas(pathname,['sobol_indices_multiscale_solution_' num2str(i) '_var_' num2str(j)],formats,renderer);
            
            plotSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,U,w,j,'displ',i);
            mysaveas(pathname,['sensitivity_indices_multiscale_solution_' num2str(i) '_var_' num2str(j)],formats,renderer);
        end
    end
    
    %% Display random evaluations
    % nbsamples = 3;
    % for i=1:nbsamples
    %     xi = random(rv,1);
    %     U_xi = U(xi);
    %     w_xi = cellfun(@(x) x(xi),w,'UniformOutput',false);
    %     lambda_xi = cellfun(@(x) x(xi),lambda,'UniformOutput',false);
    %     for j=1:2
    %         % plotAllSolutions(glob,patches.eval(xi),interfaces,U_xi',cellfun(@(x) x',w_xi,'UniformOutput',false),cellfun(@(x) x',lambda_xi,'UniformOutput',false),'displ',j);
    %         plotGlobalSolution(glob,U_xi','displ',j);
    %         % plotLocalSolution(patches,cellfun(@(x) x',w_xi,'UniformOutput',false),'displ',j);
    %         % plotLagrangeMultiplier(interfaces,cellfun(@(x) x',lambda_xi,'UniformOutput',false),'displ',j);
    %         plotMultiscaleSolution(glob,patches.eval(xi),interfaces,U_xi',cellfun(@(x) x',w_xi,'UniformOutput',false),'displ',j);
    %     end
    % end
end

myparallel('stop');
