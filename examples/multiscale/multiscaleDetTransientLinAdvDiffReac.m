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
    globOut = GlobalOutside();
    
    % Patches
    patches = Patches(n);
    
    D_patch = cell(1,n);
    D_patch{1} = DOMAIN(2,[0.85,0.40],[1.05,0.60]);
    D_patch{2} = DOMAIN(2,[0.45,0.40],[0.65,0.60]);
    D_patch{3} = DOMAIN(2,[0.05,0.40],[0.25,0.60]);
    
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
    Sadv = glob.S;
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
    V_out = V;
    V_in = cell(1,n);
    for k=1:n
        V_in{k} = V;
    end
    Sadv_patch = cell(1,n);
    for k=1:n
        Sadv_patch{k} = patches.patches{k}.S;
        Sadv_patch{k} = setmaterial(Sadv_patch{k},mat);
        Sadv_patch{k} = final(Sadv_patch{k});
    end
    phi_patch = cell(1,n);
    v_patch = cell(1,n);
    V_patch = cell(1,n);
    for k=1:n
        phi_patch{k} = calcProjection(Sadv_patch{k},Sadv)'*phi;
        v_patch{k} = 2*FENODEFIELD(calc_sigma(Sadv_patch{k},phi_patch{k},'node'));
        V_patch{k} = getvalue(v_patch{k});
        V_patch{k} = {{FENODEFIELD(V_patch{k}(:,1)),FENODEFIELD(V_patch{k}(:,2))}};
    end
    
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
    mat_out{1} = FOUR_ISOT('k',K_out,'c',c_out,'b',V_out,'r',R1_out);
    mat_out{2} = FOUR_ISOT('k',K_out,'c',c_out,'b',V_out,'r',R2_out);
    mat_out{1} = setnumber(mat_out{1},1);
    mat_out{2} = setnumber(mat_out{2},2);
    glob.S = setmaterial(glob.S,mat_out{1},1);
    glob.S = setmaterial(glob.S,mat_out{2},2:3);
    
    % Patches
    mat_patch = MATERIALS();
    for k=1:n
        mat_patch{k} = FOUR_ISOT('k',K_patch{k},'c',c_patch{k},'b',V_patch{k},'r',R_patch{k});
        mat_patch{k} = setnumber(mat_patch{k},2+k);
        patches.patches{k}.S = setmaterial(patches.patches{k}.S,mat_patch{k});
    end
    
    % Fictitious patches
    mat_in = MATERIALS();
    for k=1:n
        mat_in{k} = FOUR_ISOT('k',K_in{k},'c',c_in{k},'b',V_in{k},'r',R_in{k});
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
    globOut.S = glob.S_out;
    
    % Patches
    for k=1:n
        patches.patches{k}.S = final(patches.patches{k}.S);
    end
    
    % Interfaces
    interfaces = Interfaces(patches);
    
    %% Initial conditions
    % already taken into account by Dirichlet boundary conditions
    % Global
    % glob.u0 = calc_init_dirichlet(glob.S);
    
    % Complementary subdomain
    % globOut.u0 = calc_init_dirichlet(globOut.S);
    
    % Patches
    % for k=1:n
    %     patches.patches{k}.u0 = calc_init_dirichlet(patches.patches{k}.S);
    % end
    
    %% Time scheme
    t0 = 0;
    t1 = 2;
    nt = 100;
    T = TIMEMODEL(t0,t1,nt);
    
    % N = EULERTIMESOLVER(T,'eulertype','explicit','display',false);
    N = EULERTIMESOLVER(T,'eulertype','implicit','display',false);
    % N = DGTIMESOLVER(T,1,'outputsplit',true,'display',false,'lu',true);
    
    % Global
    globOut.timeSolver = N;
    globOut.timeOrder = 1;
    
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
    globOut.M = calc_mass(globOut.S);
    [globOut.A,globOut.b0] = calc_rigi(globOut.S);
    globOut.b0 = -globOut.b0;
    globOut.b = globOut.b0*one(globOut.timeSolver);
    
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
        interfaces.interfaces{k}.P_glob = calcProjection(glob,interfaces.interfaces{k});
        interfaces.interfaces{k}.P_globOut = calcProjection(globOut,interfaces.interfaces{k});
        interfaces.interfaces{k}.P_patch = calcProjection(patches.patches{k},interfaces.interfaces{k});
    end
    
    %% Parameters for global and local problems
    % Global problem
    glob.increment = false;
    
    % Local problems
    for k=1:n
        patches.patches{k}.changeOfVariable = false;
        patches.patches{k}.increment = false;
    end
    
    %% Static solution
    glob_static = glob;
    glob_static.timeSolver = [];
    glob_static.timeOrder = [];
    glob_static.b_out = glob.b0_out;
    
    globOut_static = globOut;
    globOut_static.timeSolver = [];
    globOut_static.timeOrder = [];
    globOut_static.b = globOut.b0;
    
    patches_static = patches;
    interfaces_static = interfaces;
    for k=1:n
        patches_static.patches{k}.timeSolver = [];
        patches_static.patches{k}.timeOrder = [];
        patches_static.patches{k}.b = patches.patches{k}.b0;
    end
    
    %% Save variables
    save(fullfile(pathname,'problem_static.mat'),'glob_static','globOut_static','patches_static','interfaces_static');
    save(fullfile(pathname,'problem_dynamic.mat'),'glob','globOut','patches','interfaces','N','D_patch','Sadv','Sadv_patch','v','v_patch','phi','phi_patch');
else
    load(fullfile(pathname,'problem_static.mat'),'glob_static','globOut_static','patches_static','interfaces_static');
    load(fullfile(pathname,'problem_dynamic.mat'),'glob','globOut','patches','interfaces','N','D_patch','Sadv','Sadv_patch','v','v_patch','phi','phi_patch');
end

%% Direct solver
if directSolver
    DSt = DirectSolver();
    DSt.changeOfVariable = false;
    DSt.timeSolver = N;
    DSt.timeOrder = 1;
    DSt.display = true;
    
    DS = DSt;
    DS.timeSolver = [];
    DS.timeOrder = [];
    
    [U_ref,w_ref,lambda_ref,output_ref] = DS.solve(globOut_static,patches_static,interfaces_static);
    % [Ut_ref,wt_ref,lambdat_ref,outputt_ref] = DSt.solve(globOut,patches,interfaces);
    save(fullfile(pathname,'reference_solution_static.mat'),'U_ref','w_ref','lambda_ref','output_ref');
    % save(fullfile(pathname,'reference_solution_dynamic.mat'),'Ut_ref','wt_ref','lambdat_ref','outputt_ref');
else
    load(fullfile(pathname,'reference_solution_static.mat'),'U_ref','w_ref','lambda_ref','output_ref');
    % load(fullfile(pathname,'reference_solution_dynamic.mat'),'Ut_ref','wt_ref','lambdat_ref','outputt_ref');
end

%% Outputs
fprintf('\n')
fprintf('spatial dimension = %d for U_ref\n',length(U_ref))
for k=1:n
    fprintf('                  = %d for w_ref{%u}\n',length(w_ref{k}),k)
    fprintf('                  = %d for lambda_ref{%u}\n',length(lambda_ref{k}),k)
end
fprintf('elapsed time = %f s\n',output_ref.time)

% fprintf('\n')
% fprintf('spatial dimension = %d for U_ref\n',size(Ut_ref,1))
% for k=1:n
%     fprintf('                  = %d for w_ref{%u}\n',size(wt_ref{k},1),k)
%     fprintf('                  = %d for lambda_ref{%u}\n',size(lambdat_ref{k},1),k)
% end
% fprintf('nb time steps = %g\n',getnt(N))
% fprintf('nb time dofs  = %g\n',getnbtimedof(N))
% fprintf('elapsed time = %f s\n',outputt_ref.time)

%% Global-local Iterative solver
if iterativeSolver
    ISt = IterativeSolver();
    ISt.maxIterations = 10;
    ISt.tolerance = eps;
    ISt.relaxation = 'Aitken';
    ISt.updateRelaxationParameter = true;
    ISt.errorCriterion = 'reference';
    % ISt.referenceSolution = {Ut_ref,wt_ref,lambdat_ref};
    ISt.display = true;
    ISt.displayIterations = true;
    
    IS = ISt;
    IS.referenceSolution = {U_ref,w_ref,lambda_ref};
    
    [U,w,lambda,output] = IS.solve(glob_static,patches_static,interfaces);
    % [Ut,wt,lambdat,outputt] = ISt.solve(glob,patches,interfaces);
    save(fullfile(pathname,'solution_static.mat'),'U','w','lambda','output');
    % save(fullfile(pathname,'solution_dynamic.mat'),'Ut','wt','lambdat','outputt');
else
    load(fullfile(pathname,'solution_static.mat'),'U','w','lambda','output');
    % load(fullfile(pathname,'solution_dynamic.mat'),'Ut','wt','lambdat','outputt');
end

%% Outputs
fprintf('\n')
fprintf('Static solution\n')
fprintf('spatial dimension = %d for U\n',length(U))
for k=1:n
    fprintf('                  = %d for w{%u}\n',length(w{k}),k)
    fprintf('                  = %d for lambda{%u}\n',length(lambda{k}),k)
end
fprintf('elapsed time = %f s\n',output.totalTime)

% fprintf('\n')
% fprintf('Dynamic solution\n')
% fprintf('spatial dimension = %d for U\n',size(Ut,1))
% for k=1:n
%     fprintf('                  = %d for w{%u}\n',size(wt{k},1),k)
%     fprintf('                  = %d for lambda{%u}\n',size(lambdat{k},1),k)
% end
% fprintf('elapsed time = %f s\n',outputt.totalTime)

%% Display
if displaySolution
    %% Display domains and meshes
    figure('Name','Domain')
    clf
    plot(create_boundary(glob.S));
    hold on
    h1 = plot(glob.S,'selgroup',[1,3+(1:n)],'FaceColor',getfacecolor(1),'EdgeColor','none');
    h2 = plot(glob.S,'selgroup',2,'FaceColor',[0 0 1],'EdgeColor','none');
    h3 = plot(glob.S,'selgroup',3,'FaceColor',[0 1 0],'EdgeColor','none');
    h4 = plotfacets(glob.S,5,'FaceColor',[0.63,0.13,0.94],'EdgeColor','none');
    h5 = plotfacets(glob.S,17,'FaceColor',[1 0.5 0],'EdgeColor','none');
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
    mysaveas(pathname,'advection_velocity_global',formats,renderer);
    
    figure('Name','Advection velocity')
    clf
    Sadv_out = getfinalmodelpart(Sadv,0);
    phi_out = calc_P_free(Sadv,Sadv_out)*phi;
    plot(phi_out,Sadv_out);
    hold on
    for k=1:n
        plot(phi_patch{k},Sadv_patch{k});
    end
    colorbar
    set(gca,'FontSize',16)
    ampl = 6;
    v_out = calc_P(Sadv,Sadv_out)*v;
    quiver(v_out,Sadv_out,ampl,'k');
    for k=1:n
        quiver(v_patch{k},Sadv_patch{k},ampl,'k');
    end
    hold off
    ylim([0,1.7])
    mysaveas(pathname,'advection_velocity_global_patches',formats,renderer);
    
    % plotDomain(glob.S,D_patch);
    % mysaveas(pathname,'domain_global_patches',formats,renderer);
    % mymatlab2tikz(pathname,'domain_global_patches.tex');
    
    figure('Name','Domain')
    clf
    plot(create_boundary(glob.S));
    hold on
    h1 = plot(glob.S,'selgroup',1:3,'FaceColor',getfacecolor(1),'EdgeColor','none');
    h2 = plot(glob.S,'selgroup',2,'FaceColor',[0 0 1],'EdgeColor','none');
    h3 = plot(glob.S,'selgroup',3,'FaceColor',[0 1 0],'EdgeColor','none');
    h4 = plotfacets(glob.S,5,'FaceColor',[0.63,0.13,0.94],'EdgeColor','none');
    h5 = plotfacets(glob.S,17,'FaceColor',[1 0.5 0],'EdgeColor','none');
    h_patch = cell(1,n);
    leg_patch = cell(1,n);
    for k=1:n
        plot(create_boundary(patches.patches{k}.S));
        h_patch{k} = plot(patches.patches{k}.S,'FaceColor',getfacecolor(1+k),'EdgeColor','none');
        leg_patch{k} = ['$\Lambda_{' num2str(k) '}$'];
    end
    hold off
    set(gca,'FontSize',16)
    l = legend([h1(1),h2(1),h3(1),h4(1),h5(1),h_patch{:}],...
        '$\Omega_1 \setminus \Lambda$','$\Omega_2$','$\Omega_0$','$\Gamma_D^1$','$\Gamma_D^2$',leg_patch{:});
    set(l,'Interpreter','latex')
    mysaveas(pathname,'domain_global_patches',formats,renderer);
    mymatlab2tikz(pathname,'domain_global_patches.tex');
    
    % plotPartition(glob,'legend',false);
    % mysaveas(pathname,'mesh_partition',formats,renderer);
    
    % plotModel(glob,patches,'legend',false);
    % mysaveas(pathname,'mesh_global_patches',formats,renderer);
    
    figure('Name','Mesh')
    clf
    h1 = plot(glob.S,'selgroup',1:3,'FaceColor',getfacecolor(1));
    hold on
    h2 = plot(glob.S,'selgroup',2,'FaceColor',[0 0 1]);
    h3 = plot(glob.S,'selgroup',3,'FaceColor',[0 1 0]);
    h4 = plotfacets(glob.S,5,'FaceColor',[0.63,0.13,0.94]);
    h5 = plotfacets(glob.S,17,'FaceColor',[1 0.5 0]);
    h_patch = cell(1,n);
    leg_patch = cell(1,n);
    for k=1:n
        h_patch{k} = plot(patches.patches{k}.S,'FaceColor',getfacecolor(1+k));
        leg_patch{k} = ['$\Lambda_{' num2str(k) '}$'];
    end
    hold off
    set(gca,'FontSize',16)
    l = legend([h1(1),h2(1),h3(1),h4(1),h5(1),h_patch{:}],...
        '$\Omega_1 \setminus \Lambda$','$\Omega_2$','$\Omega_0$','$\Gamma_D^1$','$\Gamma_D^2$',leg_patch{:});
    set(l,'Interpreter','latex')
    mysaveas(pathname,'mesh_global_patches',formats,renderer);
    
    %% Display evolutions of error indicator, stagnation indicator, CPU time w.r.t. number of iterations for static solution
    plotError(output);
    mysaveas(pathname,'error_static','fig');
    mymatlab2tikz(pathname,'error_static.tex');
    
    plotStagnation(output);
    mysaveas(pathname,'stagnation_static','fig');
    mymatlab2tikz(pathname,'stagnation_static.tex');
    
    plotErrorGlobalSolution(output);
    mysaveas(pathname,'error_global_solution_static','fig');
    mymatlab2tikz(pathname,'error_global_solution_static.tex');
    
    plotStagnationGlobalSolution(output);
    mysaveas(pathname,'stagnation_global_solution_static','fig');
    mymatlab2tikz(pathname,'stagnation_global_solution_static.tex');
    
    plotCPUTime(output,'legend',false);
    mysaveas(pathname,'cpu_time_static','fig');
    mymatlab2tikz(pathname,'cpu_time_static.tex');
    
    plotRelaxationParameter(output,'legend',false);
    mysaveas(pathname,'relaxation_parameter_static','fig');
    mymatlab2tikz(pathname,'relaxation_parameter_static.tex');
    
    %% Display static solutions
    % plotAllSolutions(glob,patches,interfaces,U,w,lambda);
    % mysaveas(pathname,'all_solutions_static',formats,renderer);
    
    plotGlobalSolution(glob,U);
    mysaveas(pathname,'global_solution_static',formats,renderer);
    
    % plotLocalSolution(patches,w);
    % mysaveas(pathname,'local_solution_static',formats,renderer);
    %
    % plotLagrangeMultiplier(interfaces,lambda);
    % mysaveas(pathname,'Lagrange_multiplier_static',formats,renderer);
    
    plotMultiscaleSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'multiscale_solution_static',formats,renderer);
    
    plotGlobalLocalSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'global_local_solution_static',formats,renderer);
    
    plotGlobalLocalSolution(glob,patches,interfaces,U,w,'view3',true);
    mysaveas(pathname,'global_local_solution_surf_static',formats,renderer);
    
    %% Display evolutions of error indicator, stagnation indicator, CPU time w.r.t. number of iterations for dynamic solution
%     plotError(outputt);
%     mysaveas(pathname,'error_dynamic','fig');
%     mymatlab2tikz(pathname,'error_dynamic.tex');
%     
%     plotStagnation(outputt);
%     mysaveas(pathname,'stagnation_dynamic','fig');
%     mymatlab2tikz(pathname,'stagnation_dynamic.tex');
%     
%     plotErrorGlobalSolution(outputt);
%     mysaveas(pathname,'error_global_solution_dynamic','fig');
%     mymatlab2tikz(pathname,'error_global_solution_dynamic.tex');
%     
%     plotStagnationGlobalSolution(outputt);
%     mysaveas(pathname,'stagnation_global_solution_dynamic','fig');
%     mymatlab2tikz(pathname,'stagnation_global_solution_dynamic.tex');
%     
%     plotCPUTime(outputt,'legend',false);
%     mysaveas(pathname,'cpu_time_dynamic','fig');
%     mymatlab2tikz(pathname,'cpu_time_dynamic.tex');
%     
%     plotRelaxationParameter(outputt,'legend',false);
%     mysaveas(pathname,'relaxation_parameter_dynamic','fig');
%     mymatlab2tikz(pathname,'relaxation_parameter_dynamic.tex');
    
    %% Display dynamic solutions
%     % evolAllSolutions(glob,patches,interfaces,Ut,wt,lambdat);
%     % mysaveas(pathname,'all_solutions_dynamic',formats,renderer);
%     
%     evolGlobalSolution(glob,Ut);
%     mysaveas(pathname,'global_solution_dynamic',formats,renderer);
%     
%     % evolLocalSolution(patches,wt);
%     % mysaveas(pathname,'local_solution_dynamic',formats,renderer);
%     %
%     % evolLagrangeMultiplier(interfaces,lambdat);
%     % mysaveas(pathname,'Lagrange_multiplier_dynamic',formats,renderer);
%     
%     evolMultiscaleSolution(glob,patches,interfaces,Ut,wt);
%     mysaveas(pathname,'multiscale_solution_dynamic',formats,renderer);
%     
%     evolGlobalLocalSolution(glob,patches,interfaces,Ut,wt);
%     mysaveas(pathname,'global_local_solution_dynamic',formats,renderer);
%     
%     evolGlobalLocalSolution(glob,patches,interfaces,Ut,wt,'view3',true);
%     mysaveas(pathname,'global_local_solution_surf_dynamic',formats,renderer);
    
    %% Display quantity of interest
    % boutput: concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) as a function of time
    % Ioutput: total concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) along the complete time evolution,
    %          corresponding to all the pollutant that the actual filter
    %          (group #1 in mesh) is not able to retain
    foutput = bodyload(keepgroupelem(glob.S,2),[],'QN',1,'nofree');
    foutput_ref = bodyload(keepgroupelem(globOut.S,2),[],'QN',1,'nofree');
    boutput = foutput'*unfreevector(glob.S,Ut);
    boutput_ref = foutput_ref'*unfreevector(globOut.S,Ut_ref);
    
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
