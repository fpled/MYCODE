%% Multiscale deterministic transient linear advection-diffusion-reaction problem %%
%%--------------------------------------------------------------------------------%%
% [Pares, Diez, Huerta, 2008], [Nouy, 2010]

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% myparallel('start');

%% Input data
setProblem = true;
directSolver = true;
iterativeSolver = true;
displaySolution = true;

n = 3; % number of patches
filename = ['transientLinAdvDiffReac' num2str(n) 'Patches'];
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','multiscaleDet',filename);
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
    
    % Patches
    patches = Patches(n);
    
    D_patch = cell(1,n);
    D_patch{1} = DOMAIN(2,[0.9,0.45],[1.0,0.55]);
    D_patch{2} = DOMAIN(2,[0.51,0.45],[0.61,0.55]);
    D_patch{3} = DOMAIN(2,[0.1,0.45],[0.2,0.55]);
    
    cl1 = 0.02;
    cl2 = 0.04;
    cl0 = 0.02;
    cltip = 0.01;
    glob.S = gmshcanistermulti(D_patch,cl1,cl2,cl0,cltip,cl1,fullfile(pathname,'gmsh_canister_multi'));
    
    nbelem_patch = [20,20];
    for k=1:n
        patches.patches{k}.S = build_model(D_patch{k},'nbelem',nbelem_patch);
    end
    % cl_patch = 0.005;
    % for k=1:n
    %     patches.patches{k}.S = build_model(D_patch{k},'cl',cl_patch,'filename',fullfile(pathname,['gmsh_patch_' num2str(k)]));
    % end
    
    % Partition of global mesh
    glob = partition(glob,D_patch);
    
    %% Materials
    % Linear diffusion coefficient
    K_out = 0.01;
    K_patch = cell(1,n);
    K_in = cell(1,n);
    % Thermal capacity
    c_out = 1;
    c_patch = cell(1,n);
    c_in = cell(1,n);
    % Advection velocity
    Sc = glob.S;
    mat = FOUR_ISOT('k',1);
    mat = setnumber(mat,1);
    Sc = setmaterial(Sc,mat);
    P = @(i) POINT(getnode(getridge(Sc,i)));
    L1 = LIGNE(P(5),P(6));
    L2 = LIGNE(P(15),P(16));
    Sc = final(Sc);
    Sc = addcl(Sc,P(1),'T',0);
    A = calc_rigi(Sc);
    b1 = surfload(Sc,L1,'QN',-1);
    b2 = surfload(Sc,L2,'QN',1);
    b = b1+b2;
    phi = A\b;
    v = FENODEFIELD(calc_sigma(Sc,phi,'node'));
    V = 2*getvalue(v);
    V = {{FENODEFIELD(V(:,1)),FENODEFIELD(V(:,2))}};
    % Linear reaction parameter
    R1_out = 0.1;
    R2_out = 10;
    R_patch = cell(1,n);
    R_in = cell(1,n);
    for k=1:n
        patch = patches.patches{k};
        % K_patch(x)  = K_out * (1 + f(x))
        % K_in(x)     = K_out
        % c_patch(x)  = c_out * (1 + f(x))
        % c_in(x)     = c_out
        % R_patch(x)  = R1_out * (1 + f(x))
        % R_in(x)     = R1_out
        % with f(x) = 1 if ||x-c||_Inf < L
        %           = 0 if ||x-c||_Inf >= L
        % L = norm(getsize(D_patch{k}),Inf)/4;
        % c = getcenter(D_patch{k});
        % f = @(x) distance(x,c,Inf)<L;
        % K_patch{k} = K_out * FENODEFIELD(ones(patch.S.nbnode,1) + double(squeeze(f(patch.S.node))));
        % c_patch{k} = c_out * FENODEFIELD(ones(patch.S.nbnode,1) + double(squeeze(f(patch.S.node))));
        % R_patch{k} = R1_out * FENODEFIELD(ones(patch.S.nbnode,1) + double(squeeze(f(patch.S.node))));
        K_patch{k} = K_out;
        c_patch{k} = c_out;
        R_patch{k} = R1_out;
        K_in{k} = K_out;
        c_in{k} = c_out;
        R_in{k} = R1_out;
    end
    
    % Complementary subdomain
    mat_out = MATERIALS();
    mat_out{1} = FOUR_ISOT('k',K_out,'c',c_out,'b',V,'r',R1_out);
    mat_out{2} = FOUR_ISOT('k',K_out,'c',c_out,'b',V,'r',R2_out);
    mat_out{1} = setnumber(mat_out{1},1);
    mat_out{2} = setnumber(mat_out{2},2);
    glob.S = setmaterial(glob.S,mat_out{1},1);
    glob.S = setmaterial(glob.S,mat_out{2},2:3);
    
    % Patches
    mat_patch = MATERIALS();
    for k=1:n
        mat_patch{k} = FOUR_ISOT('k',K_patch{k},'c',c_patch{k},'b',V,'r',R_patch{k});
        mat_patch{k} = setnumber(mat_patch{k},2+k);
        patches.patches{k}.S = setmaterial(patches.patches{k}.S,mat_patch{k});
    end
    
    % Fictitious patches
    mat_in = MATERIALS();
    for k=1:n
        mat_in{k} = FOUR_ISOT('k',K_in{k},'c',c_in{k},'b',V,'r',R_in{k});
        mat_in{k} = setnumber(mat_in{k},2+k);
        glob.S = setmaterial(glob.S,mat_in{k},getnumgroupelemwithparam(glob.S,'partition',k));
    end
    
    %% Dirichlet boundary conditions
    % Global
    numnode = getnumber(getnode(create_boundary(glob.S)));
    [~,numnode1] = intersect(glob.S,L1);
    [~,numnode2] = intersect(glob.S,L2);
    numnoderest = setdiff(setdiff(numnode,numnode1),numnode2);
    
    glob.S = final(glob.S);
    glob.S = addcl(glob.S,numnode1,'T',1);
    glob.S = addcl(glob.S,numnode2,'T',0);
    permeability = true; % if false, the walls of the canister are considered to be impermeable
    if permeability
        glob.S = addcl(glob.S,numnoderest,'T',0);
    end
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
    
    %% Initial conditions
    % Global
    glob.u0 = calc_init_dirichlet(glob.S);
    
    % Complementary subdomain
    glob_out.u0 = calc_init_dirichlet(glob_out.S);
    
    % Patches
    for k=1:n
        patches.patches{k}.u0 = calc_init_dirichlet(patches.patches{k}.S);
    end
    
    %% Time scheme
    t0 = 0;
    t1 = 2;
    nt = 100;
    T = TIMEMODEL(t0,t1,nt);
    
    % N = EULERTIMESOLVER(T,'eulertype','explicit','display',false);
    % N = EULERTIMESOLVER(T,'eulertype','implicit','display',false);
    N = DGTIMESOLVER(T,1,'outputsplit',true,'display',false,'lu',true);
    
    % Global
    glob_out.timeSolver = N;
    glob_out.timeOrder = 1;
    
    % Complementary subdomain
    glob.timeSolver = N;
    glob.timeOrder = 1;
    
    % Patches
    for k=1:n
        patches.patches{k}.timeSolver = N;
        patches.patches{k}.timeOrder = 1;
    end
    
    %% Mass and stifness matrices and sollicitation vectors
    % Global
    glob.M = calc_mass(glob.S);
    glob.A = calc_rigi(glob.S);
    [~,glob.b0_out] = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',0));
    glob.b0_out = -glob.b0_out;
    for k=1:n
        glob.M_in{k} = calc_mass(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
        glob.A_in{k} = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
    end
    glob.b_out = glob.b0_out*one(glob.timeSolver);
    
    % Complementary subdomain
    glob_out.M = calc_mass(glob_out.S);
    [glob_out.A,glob_out.b0] = calc_rigi(glob_out.S);
    glob_out.b0 = -glob_out.b0;
    glob_out.b = glob_out.b0*one(glob_out.timeSolver);
    
    % Patches
    for k=1:n
        patches.patches{k}.M = calc_mass(patches.patches{k}.S);
        [patches.patches{k}.A,patches.patches{k}.b0] = calc_rigi(patches.patches{k}.S);
        patches.patches{k}.b0 = -patches.patches{k}.b0;
        patches.patches{k}.b = patches.patches{k}.b0*one(patches.patches{k}.timeSolver);
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
    save(fullfile(pathname,'problem.mat'),'glob','glob_out','patches','interfaces','N','D_patch','Sc','v','phi');
else
    load(fullfile(pathname,'problem.mat'),'glob','glob_out','patches','interfaces','N','D_patch','Sc','v','phi');
end

%% Direct solver
if directSolver
    DS = DirectSolver();
    DS.changeOfVariable = false;
    DS.timeSolver = N;
    DS.timeOrder = 1;
    DS.display = true;
    
    [U_ref,w_ref,lambda_ref,output_ref] = DS.solve(glob_out,patches,interfaces);
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

%% Display
if displaySolution
    %% Display domains and meshes
    figure('Name','Domain')
    clf
    plot(create_boundary(glob.S));
    h1 = plot(glob.S,'selgroup',[1,3+(1:n)],'FaceColor',getfacecolor(1),'EdgeColor','none');
    h2 = plot(glob.S,'selgroup',2,'FaceColor',getfacecolor(2),'EdgeColor','none');
    h3 = plot(glob.S,'selgroup',3,'FaceColor',getfacecolor(3),'EdgeColor','none');
    h4 = plotfacets(glob.S,5,'FaceColor',getfacecolor(4),'EdgeColor','none');
    h5 = plotfacets(glob.S,17,'FaceColor',getfacecolor(5),'EdgeColor','none');
    set(gca,'FontSize',16)
    l = legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'$\Omega_1$','$\Omega_2$','$\Omega_3$','$\Gamma_1$','$\Gamma_2$');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    figure('Name','Advection velocity')
    clf
    plot(phi,Sc);
    colorbar
    set(gca,'FontSize',16)
    hold on
    ampl = 6;
    quiver(v,Sc,ampl,'k');
    ylim([0,1.7])
    mysaveas(pathname,'advection_velocity',formats,renderer);
    
    plotDomain(glob.S,D_patch);
%     figure('Name','Domain')
%     clf
%     h1 = plot(glob.S,'selgroup',1:3,'FaceColor',getfacecolor(1),'EdgeColor','none');
%     h2 = plot(glob.S,'selgroup',2,'FaceColor',getfacecolor(2),'EdgeColor','none');
%     h3 = plot(glob.S,'selgroup',3,'FaceColor',getfacecolor(3),'EdgeColor','none');
%     h4 = plotfacets(pb.S,5,'FaceColor',getfacecolor(4),'EdgeColor','none');
%     h5 = plotfacets(pb.S,17,'FaceColor',getfacecolor(5),'EdgeColor','none');
%     set(gca,'FontSize',16)
%     l = legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'$\Omega_1$','$\Omega_2$','$\Omega_0$','$\Gamma_1$','$\Gamma_2$');
%     set(l,'Interpreter','latex')
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
    % evolAllSolutions(glob,patches,interfaces,U,w,lambda);
    % mysaveas(pathname,'all_solutions',formats,renderer);
    
    evolGlobalSolution(glob,U);
    mysaveas(pathname,'global_solution',formats,renderer);
    
    % evolLocalSolution(patches,w);
    % mysaveas(pathname,'local_solution',formats,renderer);
    %
    % evolLagrangeMultiplier(interfaces,lambda);
    % mysaveas(pathname,'Lagrange_multiplier',formats,renderer);
    
    evolMultiscaleSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'multiscale_solution',formats,renderer);
    
    evolGlobalLocalSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'global_local_solution',formats,renderer);
    
    evolGlobalLocalSolution(glob,patches,interfaces,U,w,'view3',true);
    mysaveas(pathname,'global_local_solution_surf',formats,renderer);
    
    %% Display quantity of interest
    % boutput: concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) as a function of time
    % Ioutput: total concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) along the complete time evolution,
    %          corresponding to all the pollutant that the actual filter
    %          (group #1 in mesh) is not able to retain
    foutput = bodyload(keepgroupelem(glob.S,2),[],'QN',1,'nofree');
    foutput_ref = bodyload(keepgroupelem(glob_out.S,2),[],'QN',1,'nofree');
    boutput = foutput'*unfreevector(glob.S,U);
    boutput_ref = foutput_ref'*unfreevector(glob_out.S,U_ref);
    
    figure('Name','Quantity of interest')
    clf
    plot(boutput,'-b','LineWidth',1);
    hold on
    plot(boutput_ref,'-r','LineWidth',1);
    grid on
    box on
    set(gca,'FontSize',16)
    xlabel('Time (s)')
    ylabel('Quantity of interest')
    legend('multiscale','monoscale')
    mysaveas(pathname,'quantity_of_interest',formats,renderer);
    
    Ioutput = integrate(boutput);
    Ioutput_ref = integrate(boutput_ref);
    errOutput = norm(Ioutput-Ioutput_ref)/Ioutput_ref;
    fprintf('quantity of interest           = %e\n',Ioutput);
    fprintf('reference quantity of interest = %e\n',Ioutput_ref);
    fprintf('error in quantity of interest  = %e\n',erroutput);
end

% myparallel('stop');
