%% Monoscale stochastic transient linear advection-diffusion-reaction problem %%
%%----------------------------------------------------------------------------%%
% [Pares, Diez, Huerta, 2008, CMAME]
% [Nouy, 2010, CMAME]

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

filename = 'transientLinAdvDiffReac';
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
    cl1 = 0.02;
    cl2 = 0.04;
    cl0 = 0.02;
    cltip = 0.01;
    pb.S = gmshcanister(cl1,cl2,cl0,cltip,fullfile(pathname,'gmsh_canister'));
    
    %% Random variables
    d = 3; % parametric dimension
    r = UniformRandomVariable(0,1);
    rv = RandomVector(r,d);
    
    %% Materials
    % IntegrationRule
    p = 1;
    bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
    bases = FunctionalBases(bases);
    H = FullTensorProductFunctionalBasis(bases);
    I = gaussIntegrationRule(r,2);
    I = I.tensorize(d);
    
    % Linear diffusion coefficient
    % K1(xi) = 0.01 * (1 + 0.25 * (2 * xi - 1))
    % K2 = 0.01
    fun = @(xi) 0.01 * (1 + 0.25 * (2 * xi(:,1) - 1));
    fun = UserDefinedFunction(fun,d);
    fun.evaluationAtMultiplePoints = true;
    
    K1 = H.projection(fun,I);
    K2 = 0.01;
    
    % Thermal capacity
    % c1(xi) = 1 * (1 + 0.1 * (2 * xi - 1))
    % c2 = 1
    fun = @(xi) 1 * (1 + 0.1 * (2 * xi(:,2) - 1));
    fun = UserDefinedFunction(fun,d);
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
    % R1(xi) = 0.1 * (1 + 0.25 * (2 * xi - 1))
    % R2 = 10
    fun = @(xi) 0.1 * (1 + 0.25 * (2 * xi(:,3) - 1));
    fun = UserDefinedFunction(fun,d);
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
    T = TIMEMODEL(t0,t1,nt);
    
    % pb.timeSolver = EULERTIMESOLVER(T,'eulertype','explicit','display',false);
    pb.timeSolver = EULERTIMESOLVER(T,'eulertype','implicit','display',false);
    % pb.timeSolver = DGTIMESOLVER(T,1,'outputsplit',true,'display',false,'lu',true);
    pb.timeOrder = 1;
    
    pb.loadFunction = @(N) one(N);
    
    %% Mass and stifness matrices and sollicitation vectors
    if ~israndom(pb.S)
        pb.M = calc_mass(pb.S);
        [pb.A,pb.b0] = calc_rigi(pb.S);
        pb.b0 = -pb.b0;
        pb.b = pb.b0*pb.loadFunction(pb.timeSolver);
    end
    
    save(fullfile(pathname,'problem.mat'),'pb','Sadv','v','phi','rv');
else
    load(fullfile(pathname,'problem.mat'),'pb','Sadv','v','phi','rv');
end

%% Solution
if solveProblem
    p = 50;
    bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
    bases = FunctionalBases(bases);
    
    s = AdaptiveSparseTensorAlgorithm();
    s.tolStagnation = 1e-1;
    s.display = true;
    s.displayIterations = true;
    
    ls = LinearModelLearningSquareLoss();
    ls.errorEstimation = true;
    ls.sharedCoefficients = false;
    
    pbt = pb;
    
    % Stationary solution
    s.tol = 1e-12;
    pb.timeSolver = [];
    pb.loadFunction = [];
    fun = @(xi) solveSystem(calcOperator(funEval(pb,xi)));
    fun = UserDefinedFunction(fun,d,getnbddlfree(pb.S));
    fun.evaluationAtMultiplePoints = false;
    
    t = tic;
    [u,err,~,y] = s.leastSquares(fun,bases,ls,rv);
    time = toc(t);
    
    % Transient solution
    s.tol = 1e-1;
    pb = pbt;
    funt = @(xi) solveSystem(calcOperator(funEval(pb,xi)));
    funt = UserDefinedFunction(funt,d,[getnbddlfree(pb.S),getnbtimedof(pb.timeSolver)]);
    funt.evaluationAtMultiplePoints = false;
    
    funtCell = @(xi) solveSystemCell(calcOperator(funEval(pb,xi)));
    funtCell = UserDefinedFunction(funtCell,d,2);
    funtCell.evaluationAtMultiplePoints = false;
    
    t = tic;
    [ut,errt,~,yt] = s.leastSquares(funt,bases,ls,rv);
    timet = toc(t);
    
    save(fullfile(pathname,'solution.mat'),'u','err','y','fun','time');
    save(fullfile(pathname,'solution_transient.mat'),'ut','errt','yt','funt','funtCell','timet');
else
    load(fullfile(pathname,'solution.mat'),'u','err','y','fun','time');
    load(fullfile(pathname,'solution_transient.mat'),'ut','errt','yt','funt','funtCell','timet');
end

%% Outputs
fprintf('\n');
fprintf('nb elements = %g\n',getnbelem(pb.S));
fprintf('nb nodes    = %g\n',getnbnode(pb.S));
fprintf('nb dofs     = %g\n',getnbddl(pb.S));
fprintf('time solver : %s\n',class(pb.timeSolver));
fprintf('nb time steps = %g\n',getnt(pb.timeSolver));
fprintf('nb time dofs  = %g\n',getnbtimedof(pb.timeSolver));

fprintf('\n');
fprintf('Stationary solution\n');
fprintf('-------------------\n');
fprintf('spatial dimension = %d\n',u.sz)
fprintf('parametric dimension = %d\n',ndims(u.basis))
fprintf('basis dimension = %d\n',cardinal(u.basis))
fprintf('order = [ %s ]\n',num2str(max(u.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(u.basis.indices.array)
fprintf('nb samples = %d\n',size(y,1))
fprintf('CV error = %d\n',norm(err))
fprintf('elapsed time = %f s\n',time)

fprintf('\n');
fprintf('Transient solution\n');
fprintf('------------------\n');
fprintf('spatial dimension = %d\n',ut.sz(1))
fprintf('time dimension = %d\n',ut.sz(2))
fprintf('parametric dimension = %d\n',ndims(ut.basis))
fprintf('basis dimension = %d\n',cardinal(ut.basis))
fprintf('order = [ %s ]\n',num2str(max(ut.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(ut.basis.indices.array)
fprintf('nb samples = %d\n',size(yt,1))
fprintf('CV error = %d\n',norm(errt))
fprintf('elapsed time = %f s\n',timet)

T = gettimemodel(pb.timeSolver);
sz = [getnbddl(pb.S),getnbtimedof(T)];

vt_data = zeros([cardinal(ut.basis),ut.sz]);
vt_data_tot = zeros([cardinal(ut.basis),sz]);
for k=1:cardinal(ut.basis)
    vtk = diff(pb.timeSolver,TIMEMATRIX(reshape(ut.data(k,:,:),ut.sz),T));
    vt_data(k,:,:) = getvalue(vtk);
    vt_data_tot(k,:,:) = getvalue(unfreevector(pb.S,vtk)-calc_init_dirichlet(pb.S)*one(T));
end
vt = FunctionalBasisArray(vt_data,ut.basis,ut.sz);
vt_tot = FunctionalBasisArray(vt_data_tot,ut.basis,sz);

%% Test
if testSolution
    N = 100;
    errL2 = testError(u,fun,N,rv);
    errL2t = testError(ut,funt,N,rv);
    save(fullfile(pathname,'test.mat'),'errL2','errL2t');
else
    load(fullfile(pathname,'test.mat'),'errL2','errL2t');
end
fprintf('\n');
fprintf('Stationary solution\n');
fprintf('-------------------\n');
fprintf('mean squared error = %d\n',errL2)

fprintf('\n');
fprintf('Transient solution\n');
fprintf('------------------\n');
fprintf('mean squared error = %d\n',errL2t)

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
    l = legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'$\Omega_1$','$\Omega_2$','$\Omega_3$','$\Gamma_D^1$','$\Gamma_D^2$','Location','NorthEast');
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
    l = legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'$\Omega_1$','$\Omega_2$','$\Omega_3$','$\Gamma_D^1$','$\Gamma_D^2$','Location','NorthEast');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'mesh',formats,renderer);
    
    %% Display multi-index set for stationary solution
    plotMultiIndexSet(u,'legend',false);
    mysaveas(pathname,'multi_index_set_solution_stationary','fig');
    mymatlab2tikz(pathname,'multi_index_set_solution_stationary.tex');
    
    %% Display multi-index set for transient solution
    plotMultiIndexSet(ut,'legend',false);
    mysaveas(pathname,'multi_index_set_solution_transient','fig');
    mymatlab2tikz(pathname,'multi_index_set_solution_transient.tex');
    
    %% Display statistical outputs for stationary solution
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
    
    %% Display evolution of statistical outputs for transient solution
    T = gettimemodel(pb.timeSolver);
    
    evolMean(pb.S,T,ut,'filename','mean_solution','pathname',pathname);
    evolMean(pb.S,T,ut,'surface',true,'filename','mean_solution_surface','pathname',pathname);
    
    evolMean(pb.S,T,vt_tot,'rescale',false,'filename','mean_velocity','pathname',pathname);
    evolMean(pb.S,T,vt_tot,'rescale',false,'surface',true,'filename','mean_velocity_surface','pathname',pathname);
    
    evolVariance(pb.S,T,ut,'filename','var_solution','pathname',pathname);
    evolVariance(pb.S,T,ut,'surface',true,'filename','var_solution_surface','pathname',pathname);
    
    evolVariance(pb.S,T,vt,'rescale',false,'filename','var_velocity','pathname',pathname);
    evolVariance(pb.S,T,vt,'rescale',false,'surface',true,'filename','var_velocity_surface','pathname',pathname);
    
    evolStd(pb.S,T,ut,'filename','std_solution','pathname',pathname);
    evolStd(pb.S,T,ut,'surface',true,'filename','std_solution_surface','pathname',pathname);
    
    evolStd(pb.S,T,vt,'rescale',false,'filename','std_velocity','pathname',pathname);
    evolStd(pb.S,T,vt,'rescale',false,'surface',true,'filename','std_velocity_surface','pathname',pathname);
    
    d = ndims(ut.basis);
    for i=1:d
        % evolSobolIndices(pb.S,T,ut,i,'filename',['sobol_indices_solution_var_' num2str(i)],'pathname',pathname);
        % evolSobolIndices(pb.S,T,ut,i,'surface',true,'filename',['sobol_indices_solution_var_' num2str(i) '_surface'],'pathname',pathname);
        
        % evolSobolIndices(pb.S,T,vt,i,'filename',['sobol_indices_velocity_var_' num2str(i)],'pathname',pathname);
        % evolSobolIndices(pb.S,T,vt,i,'surface',true,'filename',['sobol_indices_velocity_var_' num2str(i) '_surface'],'pathname',pathname);
        
        evolSensitivityIndices(pb.S,T,ut,i,'filename',['sensitivity_indices_solution_var_' num2str(i)],'pathname',pathname);
        evolSensitivityIndices(pb.S,T,ut,i,'surface',true,'filename',['sensitivity_indices_solution_var_' num2str(i) '_surface'],'pathname',pathname);
        
        evolSensitivityIndices(pb.S,T,vt,i,'filename',['sensitivity_indices_velocity_var_' num2str(i)],'pathname',pathname);
        evolSensitivityIndices(pb.S,T,vt,i,'surface',true,'filename',['sensitivity_indices_velocity_var_' num2str(i) '_surface'],'pathname',pathname);
    end
    
    %% Display statistical outputs at different instants for transient solution
%     [t,rep] = gettevol(pb.timeSolver);
%     for k=1:floor(length(rep)/4):length(rep)
%         close all
%         uk_data = ut.data(:,:,rep(k));
%         vk_data = vt.data(:,:,rep(k));
%         vk_data_tot = vt_tot.data(:,:,rep(k));
%         uk = FunctionalBasisArray(uk_data,ut.basis,ut.sz(1));
%         vk = FunctionalBasisArray(vk_data,vt.basis,vt.sz(1));
%         vk_tot = FunctionalBasisArray(vk_data_tot,vt_tot.basis,vt_tot.sz(1));
%         
%         plotMean(pb.S,uk);
%         mysaveas(pathname,['mean_solution_t' num2str(k-1)],formats,renderer);
%         plotMean(pb.S,uk,'surface',true);
%         mysaveas(pathname,['mean_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotMean(pb.S,vk_tot);
%         mysaveas(pathname,['mean_velocity_t' num2str(k-1)],formats,renderer);
%         plotMean(pb.S,vk_tot,'surface',true);
%         mysaveas(pathname,['mean_velocity_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotVariance(pb.S,uk);
%         mysaveas(pathname,['var_solution_t' num2str(k-1)],formats,renderer);
%         plotVariance(pb.S,uk,'surface',true);
%         mysaveas(pathname,['var_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotVariance(pb.S,vk);
%         mysaveas(pathname,['var_velocity_t' num2str(k-1)],formats,renderer);
%         plotVariance(pb.S,vk,'surface',true);
%         mysaveas(pathname,['var_velocity_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotStd(pb.S,uk);
%         mysaveas(pathname,['std_solution_t' num2str(k-1)],formats,renderer);
%         plotStd(pb.S,uk,'surface',true);
%         mysaveas(pathname,['std_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotStd(pb.S,vk);
%         mysaveas(pathname,['std_velocity_t' num2str(k-1)],formats,renderer);
%         plotStd(pb.S,vk,'surface',true);
%         mysaveas(pathname,['std_velocity_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         d = ndims(uk.basis);
%         for i=1:d
%             % plotSobolIndices(pb.S,uk,i);
%             % mysaveas(pathname,['sobol_indices_solution_var_' num2str(i) '_t' num2str(k-1)],formats,renderer);
%             % plotSobolIndices(pb.S,uk,i,'surface',true);
%             % mysaveas(pathname,['sobol_indices_solution_var_' num2str(i) '_t' num2str(k-1) '_surface'],formats,renderer);
%             
%             % plotSobolIndices(pb.S,vk,i);
%             % mysaveas(pathname,['sobol_indices_velocity_var_' num2str(i) '_t' num2str(k-1)],formats,renderer);
%             % plotSobolIndices(pb.S,vk,i,'surface',true);
%             % mysaveas(pathname,['sobol_indices_velocity_var_' num2str(i) '_t' num2str(k-1) '_surface'],formats,renderer);
%             
%             plotSensitivityIndices(pb.S,uk,i);
%             mysaveas(pathname,['sensitivity_indices_solution_var_' num2str(i) '_t' num2str(k-1)],formats,renderer);
%             plotSensitivityIndices(pb.S,uk,i,'surface',true);
%             mysaveas(pathname,['sensitivity_indices_solution_var_' num2str(i) '_t' num2str(k-1) '_surface'],formats,renderer);
%             
%             plotSensitivityIndices(pb.S,vk,i);
%             mysaveas(pathname,['sensitivity_indices_velocity_var_' num2str(i) '_t' num2str(k-1)],formats,renderer);
%             plotSensitivityIndices(pb.S,vk,i,'surface',true);
%             mysaveas(pathname,['sensitivity_indices_velocity_var_' num2str(i) '_t' num2str(k-1) '_surface'],formats,renderer);
%         end
%     end
    
    %% Display quantity of interest
    % boutput: mean of concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) as a function of time
    % Ioutput: mean of total concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) along the complete time evolution,
    %          corresponding to all the pollutant that the actual filter
    %          (group #1 in mesh) is not able to retain
    foutput = bodyload(keepgroupelem(pb.S,2),[],'QN',1,'nofree');
    
    mean_ut = mean(ut);
    mean_ut = reshape(mean_ut,ut.sz);
    mean_ut = TIMEMATRIX(mean_ut,T);
    mean_ut = unfreevector(pb.S,mean_ut);
    
    std_ut = std(ut);
    std_ut = reshape(std_ut,ut.sz);
    std_ut = TIMEMATRIX(std_ut,T);
    std_ut = unfreevector(pb.S,std_ut)-calc_init_dirichlet(pb.S)*one(T);
    
    mean_boutput = foutput'*mean_ut;
    std_boutput = foutput'*std_ut;
    
    probs = [0.025 0.975];
    nbSamples = 100;
    x = random(rv,nbSamples);
    samples_ut = ut(x);
    samples_boutput = zeros(nbSamples,ut.sz(2));
    for i=1:nbSamples
        sample_ut = reshape(samples_ut(i,:,:),ut.sz);
        sample_ut = TIMEMATRIX(sample_ut,T);
        sample_ut = unfreevector(pb.S,sample_ut)-calc_init_dirichlet(pb.S)*one(T);
        samples_boutput(i,:) = foutput'*sample_ut;
    end
    ci_boutput = quantile(samples_boutput,probs);
    
    figure('Name','Quantity of interest')
    clf
    [t,rep] = gettplot(pb.timeSolver);
    ciplot(ci_boutput(1,:),ci_boutput(2,:),t,'b');
    alpha(0.2)
    hold on
    plot(mean_boutput,'-b','LineWidth',1);
    hold off
    grid on
    box on
    set(gca,'FontSize',16)
    xlabel('Time [s]')
    ylabel('Concentration of pollutant')
    l = legend({['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value'});
    set(l,'Interpreter','latex')
    mysaveas(pathname,'quantity_of_interest',formats,renderer);
    mymatlab2tikz(pathname,'quantity_of_interest.tex');
    
    mean_Ioutput = integrate(mean_boutput);
    std_Ioutput = integrate(std_boutput);
    lowerci_Ioutput = integrate(TIMEMATRIX(ci_boutput(1,:),T));
    upperci_Ioutput = integrate(TIMEMATRIX(ci_boutput(2,:),T));
    fprintf('mean of quantity of interest = %e\n',mean_Ioutput);
    fprintf('std  of quantity of interest = %e\n',std_Ioutput);
    fprintf('%d%% confidence interval of quantity of interest = [%e,%e]\n',(probs(2)-probs(1))*100,lowerci_Ioutput,upperci_Ioutput);
end

myparallel('stop');
