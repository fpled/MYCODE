%% Multiscale deterministic linear diffusion formulation problem %%
%%---------------------------------------------------------------%%

% clc
clearvars
close all
% myparallel('start');

%% Input data
setProblem = true;
directSolver = true;
iterativeSolver = true;
displaySolution = true;

n = 4; % number of patches n = 1, 2, 4
filename = ['linDiff' num2str(n) 'Patches'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','multiscaleDet',filename);
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
    
    nbelem = [20,20];
    glob.S = build_model(D,'nbelem',nbelem);
    % cl = 0.05;
    % glob.S = build_model(D,'cl',cl,'filename',fullfile(pathname,'gmsh_domain'));
    glob.S = final(glob.S);
    glob.S = addcl(glob.S,[]);
    
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
        patches.patches{k}.S = final(patches.patches{k}.S);
    end
    % cl_patch = 0.005;
    % for k=1:n
    %     patches.patches{k}.S = build_model(D_patch{k},'cl',cl_patch,'filename',fullfile(pathname,['gmsh_patch_' num2str(k)]));
    %     patches.patches{k}.S = final(patches.patches{k}.S);
    % end
    
    % Partition of global mesh
    glob = partition(glob,D_patch);
    glob.S_out = getfinalmodelpart(glob.S,0);
    % S_in = cell(1,n);
    % for k=1:n
    %     S_in{k} = getfinalmodelpart(glob.S,k);
    % end
    
    % Complementary subdomain
    globOut.S = glob.S_out;
    
    % Interfaces
    interfaces = Interfaces(patches);
    
    %% Bilinear and linear forms
    % Linear diffusion coefficient
    K_out = 1;
    K_patch = cell(1,n);
    K_in = cell(1,n);
    for k=1:n
        patch = patches.patches{k};
        % K_patch(x) = 1 + beta_patch * f(x)
        % K_in(x)    = 1 + beta_in * f(x)
        % with f(x) = alpha*exp( -ampl*||x-c||_2^2/L^2 ) if ||x-c||_Inf < L
        %           = 0                                 if ||x-c||_Inf >= L
        % alpha = 10;
        % ampl = 2;
        % L = norm(getsize(D_patch{k}),Inf)/4;
        % c = getcenter(D_patch{k});
        % f = @(x) (distance(x,c,Inf)<L) * alpha * exp(-ampl*distance(x,c,2).^2/L^2);
        % beta_patch = 1;
        % K_patch{k} = 1 + beta_patch * squeeze(f(patch.S.node));
        % beta_in = 0;
        % K_in{k} = 1 + beta_in * squeeze(f(glob.S.node));
        
        % K_patch(x) = 1 + f(x)
        % K_in(x)    = 1
        % with f(x) = 1 if ||x-c||_Inf < L
        %           = 0 if ||x-c||_Inf >= L
        L = norm(getsize(D_patch{k}),Inf)/4;
        c = getcenter(D_patch{k});
        f = @(x) distance(x,c,Inf)<L;
        K_patch{k} = ones(patch.S.nbnode,1) + double(squeeze(f(patch.S.node)));
        K_in{k} = 1;
    end
    
    % Complementary subdomain
    % a(u,v) = int( K.grad(u).grad(v) )
    % a_out = BILINFORM(1,1,K_out);
    a_out = DIFFUSIONFORM(K_out);
    a_out = setfree(a_out,1);
    
    % Patches
    % a(u,v) = int( K.grad(u).grad(v) )
    a_patch = cell(1,n);
    for k=1:n
        % a_patch{k} = BILINFORM(1,1,K_patch{k}); % uniform value
        % a_patch{k} = BILINFORM(1,1,K_patch{k},0); % nodal values
        a_patch{k} = DIFFUSIONFORM(K_patch{k});
        a_patch{k} = setfree(a_patch{k},1);
    end
    
    % Fictitious patches
    % a(u,v) = int( K.grad(u).grad(v) )
    a_in = cell(1,n);
    for k=1:n
        % a_in{k} = BILINFORM(1,1,K_in{k}); % uniform value
        % a_in{k} = BILINFORM(1,1,K_in{k},0); % nodal values
        a_in{k} = DIFFUSIONFORM(K_in{k});
        a_in{k} = setfree(a_in{k},1);
    end
    
    % Source term
    f = 100;
    
    % Global
    % l(v) = int( f.v )
    l = LINFORM(0,f);
    l = setfree(l,1);
    
    % Patches
    % l(v) = int( f.v )
    l_patch = l;
    
    %% Stiffness matrices and sollicitation vectors
    % Global
    a_out = setselgroup(a_out,getnumgroupelemwithparam(glob.S,'partition',0));
    A_out = calc_matrix(a_out,glob.S);
    glob.A = A_out;
    for k=1:n
        a_in{k} = setselgroup(a_in{k},getnumgroupelemwithparam(glob.S,'partition',k));
        glob.A_in{k} = calc_matrix(a_in{k},glob.S);
        glob.A = glob.A + glob.A_in{k};
    end
    l_out = setselgroup(l,getnumgroupelemwithparam(glob.S,'partition',0));
    glob.b_out = calc_vector(l_out,glob.S);
    
    % Complementary subdomain
    globOut.A = calc_matrix(a_out,globOut.S);
    globOut.b = calc_vector(l,globOut.S);
    
    % Patches
    for k=1:n
        patches.patches{k}.A = calc_matrix(a_patch{k},patches.patches{k}.S);
        patches.patches{k}.b = calc_vector(l_patch,patches.patches{k}.S);
    end
    
    %% Mass matrices
    % a = BILINFORM(0,0,1);
    % a = setfree(a,1);
    for k=1:n
        % interfaces.interfaces{k}.M = calc_matrix(a,interfaces.interfaces{k}.S);
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
    DS = DirectSolver();
    DS.changeOfVariable = false;
    DS.display = true;
    
    [U_ref,w_ref,lambda_ref,output_ref] = DS.solve(globOut,patches,interfaces);
    save(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','output_ref');
else
    load(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','output_ref');
end

%% Outputs
fprintf('\n')
fprintf('spatial dimension = %d for U_ref\n',length(U_ref))
for k=1:n
    fprintf('                  = %d for w_ref{%u}\n',length(w_ref{k}),k)
    fprintf('                  = %d for lambda_ref{%u}\n',length(lambda_ref{k}),k)
end
fprintf('elapsed time = %f s\n',output_ref.time)

%% Global-local Iterative solver
if iterativeSolver
    IS = IterativeSolver();
    IS.maxIterations = 50;
    IS.tolerance = eps;
    IS.relaxation = 'Aitken';
    IS.updateRelaxationParameter = true;
    IS.errorCriterion = 'reference';
    IS.referenceSolution = {U_ref,w_ref,lambda_ref};
    IS.display = true;
    IS.displayIterations = true;
    
    [U,w,lambda,output] = IS.solve(glob,patches,interfaces);
    save(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
else
    load(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
end

%% Outputs
fprintf('\n')
fprintf('spatial dimension = %d for U\n',length(U))
for k=1:n
    fprintf('                  = %d for w{%u}\n',length(w{k}),k)
    fprintf('                  = %d for lambda{%u}\n',length(lambda{k}),k)
end
fprintf('elapsed time = %f s\n',output.totalTime)

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
    
    %% Display evolutions of error indicator, stagnation indicator, CPU time w.r.t. number of iterations
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
    
    %% Display solutions
    % plotAllSolutions(glob,patches,interfaces,U,w,lambda);
    % mysaveas(pathname,'all_solutions',formats,renderer);
    
    plotGlobalSolution(glob,U);
    mysaveas(pathname,'global_solution',formats,renderer);
    
    % plotLocalSolution(patches,w);
    % mysaveas(pathname,'local_solution',formats,renderer);
    
    % plotLagrangeMultiplier(interfaces,lambda);
    % mysaveas(pathname,'Lagrange_multiplier',formats,renderer);
    
    plotMultiscaleSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'multiscale_solution',formats,renderer);
    
    plotGlobalLocalSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'global_local_solution',formats,renderer);
end

% myparallel('stop');
