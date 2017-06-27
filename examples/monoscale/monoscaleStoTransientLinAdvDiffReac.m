%% Monoscale stochastic transient linear advection-diffusion-reaction problem %%
%%----------------------------------------------------------------------------%%
% [Pares, Diez, Huerta, 2008], [Nouy, 2010]

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');
myparallel('start');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = false;
testSolution = true;

filename = 'transientLinAdvDiffReac';
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
    cl1 = 0.02;
    cl2 = 0.04;
    cl0 = 0.02;
    cltip = 0.01;
    pb.S = gmshcanister(cl1,cl2,cl0,cltip,fullfile(pathname,'gmsh_canister'));
    
    %% Random variables
    d = 3; % parametric dimension
    v = UniformRandomVariable(0,1);
    rv = RandomVector(v,d);
    
    %% Materials
    % Linear diffusion coefficient
    % K1(xi) = 0.01 * (1 + xi)
    % K2 = 0.01
    p = 1;
    basis = PolynomialFunctionalBasis(LegendrePolynomials(),0:p);
    bases = FunctionalBases.duplicate(basis,d);
    rvb = getRandomVector(bases);
    H = FullTensorProductFunctionalBasis(bases);
    I = gaussIntegrationRule(v,2);
    I = I.tensorize(d);
    
    fun = @(xi) 0.01 * (1 + xi(:,1));
    funtr = @(xi) fun(transfer(rvb,rv,xi));
    fun = MultiVariateFunction(funtr,d);
    fun.evaluationAtMultiplePoints = true;
    
    K1 = H.projection(fun,I);
    K2 = 0.01;
    
    % Thermal capacity
    % c1(xi) = 1 + 0.1 * xi
    % c2 = 1
    fun = @(x) 1 + 0.1 * x(:,2);
    funtr = @(x) fun(transfer(rvb,rv,x));
    fun = MultiVariateFunction(funtr,d);
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
    % R1(xi) = 0.1 * (1 + xi)
    % R2 = 10
    fun = @(x) 0.1 * (1 + x(:,3));
    funtr = @(x) fun(transfer(rvb,rv,x));
    fun = MultiVariateFunction(funtr,d);
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
    
    pb.loadFunction = @(N) one(N);
    
    %% Mass and stifness matrices and sollicitation vectors
    if ~israndom(pb.S)
        pb.M = calc_mass(pb.S);
        [pb.A,pb.b0] = calc_rigi(pb.S);
        pb.b0 = -pb.b0;
        pb.b = pb.b0*pb.loadFunction(pb.timeSolver);
    end
    
    save(fullfile(pathname,'problem.mat'),'pb','Sadv','v','phi');
else
    load(fullfile(pathname,'problem.mat'),'pb','Sadv','v','phi');
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
    s.tol = 1e-3;
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
    
    pbt = pb;
    
    % Stationary solution
    pb.timeSolver = [];
    pb.loadFunction = [];
    fun = @(xi) solveSystem(calcOperator(funEval(pb,xi)));
    fun = MultiVariateFunction(fun,d,getnbddlfree(pb.S));
    fun.evaluationAtMultiplePoints = false;
    
    t = tic;
    [u,err,~,y] = s.leastSquares(fun,bases,ls,rv);
    time = toc(t);
    
    % Transient solution
    pb = pbt;
    funt = @(xi) solveSystem(calcOperator(funEval(pb,xi)));
    funt = MultiVariateFunction(funt,d,getnbddlfree(pb.S)*getnbtimedof(pb.timeSolver));
    funt.evaluationAtMultiplePoints = false;
    funtCell = @(xi) solveSystemCell(calcOperator(funEval(pb,xi)));
    funtCell = MultiVariateFunction(funtCell,d,3);
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
fprintf('spatial dimension = %d\n',u.sz)
fprintf('parametric dimension = %d\n',ndims(u.basis))
fprintf('basis dimension = %d\n',numel(u.basis))
fprintf('order = [ %s ]\n',num2str(max(u.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(u.basis.indices.array)
fprintf('nb samples = %d\n',size(y,1))
fprintf('CV error = %d\n',norm(err))
fprintf('elapsed time = %f s\n',time)

T = gettimemodel(pb.timeSolver);
sz_free = [getnbddlfree(pb.S),getnbtimedof(T)];
sz = [getnbddl(pb.S),getnbtimedof(T)];

utoutput = ut;
ut_data = ut.data;
ut_data = reshape(ut_data,[numel(ut.basis),sz_free]);
ut = FunctionalBasisArray(ut_data,ut.basis,sz_free);

% vt_data = vt.data;
% vt_data = reshape(vt_data,[numel(vt.basis),sz_free]);
% vt = FunctionalBasisArray(vt_data,vt.basis,sz_free);
% v0 = calc_init_dirichlet(pb.S)*one(T);
% vt_data_unfree = zeros([numel(vt.basis),sz]);
% for k=1:numel(vt.basis)
%     vtk = TIMEMATRIX(reshape(vt_data(k,:,:),sz_free),T);
%     vtk = getvalue(unfreevector(pb.S,vtk)-v0);
%     vt_data_unfree(k,:,:) = reshape(vtk,[1,sz]);
% end
% vt_unfree = FunctionalBasisArray(vt_data_unfree,vt.basis,sz);

fprintf('\n');
fprintf('Transient solution\n');
fprintf('spatial dimension = %d\n',ut.sz(1))
fprintf('time dimension = %d\n',ut.sz(2))
fprintf('parametric dimension = %d\n',ndims(ut.basis))
fprintf('basis dimension = %d\n',numel(ut.basis))
fprintf('order = [ %s ]\n',num2str(max(ut.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(ut.basis.indices.array)
fprintf('nb samples = %d\n',size(yt,1))
fprintf('CV error = %d\n',norm(errt))
fprintf('elapsed time = %f s\n',timet)

%% Test
if testSolution
    Ntest = 100;
    [errtest,xtest,utest,ytest] = computeTestError(u,fun,Ntest);
    [errttest,xttest,uttest,yttest] = computeTestError(utoutput,funt,Ntest);
    [Errttest,Xttest,Uttest,Yttest] = computeTestErrorCell({utoutput,vt},funtCell,Ntest);
    save(fullfile(pathname,'test.mat'),'utest','errtest','xtest','ytest');
    save(fullfile(pathname,'testt.mat'),'uttest','errttest','xttest','yttest');
    save(fullfile(pathname,'testt.mat'),'Uttest','Errttest','Xttest','Yttest');
else
    load(fullfile(pathname,'test.mat'),'utest','errtest','xtest','ytest');
    load(fullfile(pathname,'testt.mat'),'uttest','errttest','xttest','yttest');
    load(fullfile(pathname,'testt.mat'),'Uttest','Errttest','Xttest','Yttest');
end
fprintf('\n');
fprintf('Stationary solution\n');
fprintf('test error = %d\n',errtest)

fprintf('\n');
fprintf('Transient solution\n');
fprintf('test error = %d\n',errttest)

erruttest = errttest{1};
errvttest = errttest{2};
fprintf('test error = %d for u\n',erruttest)
fprintf('test error = %d for v\n',errvttest)

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
    mysaveas(pathname,'advection_velocity',formats,renderer);

    % plotModel(pb.S,'legend',false);
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
    % l = legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'$\Omega_1$','$\Omega_2$','$\Omega_3$','$\Gamma_D^1$','$\Gamma_D^2$');
    % set(l,'Interpreter','latex')
    mysaveas(pathname,'mesh',formats,renderer);
    
    %% Display multi-index set for stationary solution
    plotMultiIndexSet(u,'legend',false);
    mysaveas(pathname,'multi_index_set_solution_stationary','fig');
    mymatlab2tikz(pathname,'multi_index_set_solution_stationary.tex');
    
    %% Display multi-index set for transient solution
    plotMultiIndexSet(ut,'legend',false);
    mysaveas(pathname,'multi_index_set_solution_transient','fig');
    mymatlab2tikz(pathname,'multi_index_set_solution_transient.tex');
    
%     plotMultiIndexSet(vt,'legend',false);
%     mysaveas(pathname,'multi_index_set_velocity','fig');
%     mymatlab2tikz(pathname,'multi_index_set_velocity.tex');
    
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
        plotSobolIndices(pb.S,u,i);
        mysaveas(pathname,['sobol_indices_solution_var_' num2str(i)],formats,renderer);
        
        plotSensitivityIndices(pb.S,u,i);
        mysaveas(pathname,['sensitivity_indices_solution_var_' num2str(i)],formats,renderer);
    end
    
    %% Display evolution of statistical outputs for transient solution
    evolMean(pb.S,T,ut,'filename','evol_mean_solution','pathname',pathname);
    evolMean(pb.S,T,ut,'surface',true,'filename','evol_mean_solution_surface','pathname',pathname);
    
%     evolMean(pb.S,T,vt_unfree,'rescale',false,'filename','evol_mean_velocity','pathname',pathname);
%     evolMean(pb.S,T,vt_unfree,'rescale',false,'surface',true,'filename','evol_mean_velocity_surface','pathname',pathname);
    
%     for i=1:2
%         evolMean(pb.S,T,ut,'epsilon',i,'filename',['evol_mean_eps_' num2str(i)],'pathname',pathname);
%         evolMean(pb.S,T,ut,'sigma',i,'filename',['evol_mean_sig_' num2str(i)],'pathname',pathname);
%     end
    
    evolVariance(pb.S,T,ut,'filename','evol_variance_solution','pathname',pathname);
    evolVariance(pb.S,T,ut,'surface',true,'filename','evol_variance_solution_surface','pathname',pathname);
    
%     evolVariance(pb.S,T,vt,'rescale',false,'filename','evol_variance_velocity','pathname',pathname);
%     evolVariance(pb.S,T,vt,'rescale',false,'surface',true,'filename','evol_variance_velocity_surface','pathname',pathname);
    
%     for i=1:2
%         evolVariance(pb.S,T,ut,'epsilon',i,'filename',['evol_variance_eps_' num2str(i)],'pathname',pathname);
%         evolVariance(pb.S,T,ut,'sigma',i,'filename',['evol_variance_sig_' num2str(i)],'pathname',pathname);
%     end
    
    evolStd(pb.S,T,ut,'filename','evol_std_solution','pathname',pathname);
    evolStd(pb.S,T,ut,'surface',true,'filename','evol_std_solution_surface','pathname',pathname);
    
%     evolStd(pb.S,T,vt,'rescale',false,'filename','evol_std_velocity','pathname',pathname);
%     evolStd(pb.S,T,vt,'rescale',false,'surface',true,'filename','evol_std_velocity_surface','pathname',pathname);
    
%     for i=1:2
%         evolStd(pb.S,T,ut,'epsilon',i,'filename',['evol_std_eps_' num2str(i)],'pathname',pathname);
%         evolStd(pb.S,T,ut,'sigma',i,'filename',['evol_std_sig_' num2str(i)],'pathname',pathname);
%     end
    
    d = ndims(ut.basis);
    for i=1:d
        evolSobolIndices(pb.S,T,ut,i,'filename',['evol_sobol_indices_solution_var_' num2str(i)],'pathname',pathname);
%         evolSobolIndices(pb.S,T,vt,i,'filename',['evol_sobol_indices_velocity_var_' num2str(i)],'pathname',pathname);
        
%        for j=1:2
%            evolSobolIndices(pb.S,T,ut,i,'epsilon',j,'filename',['evol_sobol_indices_eps_' num2str(j) '_var_' num2str(i)],'pathname',pathname);
%            evolSobolIndices(pb.S,T,ut,i,'sigma',j,'filename',['evol_sobol_indices_sig_' num2str(j) '_var_' num2str(i)],'pathname',pathname);
%        end
        
        evolSensitivityIndices(pb.S,T,ut,i,'filename',['evol_sensitivity_indices_solution_var_' num2str(i)],'pathname',pathname);
%         evolSensitivityIndices(pb.S,T,vt,i,'filename',['evol_sensitivity_indices_velocity_var_' num2str(i)],'pathname',pathname);
        
%        for j=1:2
%            evolSensitivityIndices(pb.S,T,ut,i,'epsilon',j,'filename',['evol_sensitivity_indices_eps_' num2str(j) '_var_' num2str(i)],'pathname',pathname);
%            evolSensitivityIndices(pb.S,T,ut,i,'sigma',j,'filename',['evol_sensitivity_indices_sig_' num2str(j) '_var_' num2str(i)],'pathname',pathname);
%        end
    end
    
    %% Display quantity of interest
    % boutput: concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) as a function of time
    % Ioutput: total concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) along the complete time evolution,
    %          corresponding to all the pollutant that the actual filter
    %          (group #1 in mesh) is not able to retain
%     ut = unfreevector(pb.S,ut);
%     foutput = bodyload(keepgroupelem(pb.S,2),[],'QN',1,'nofree');
%     boutput = foutput'*ut;
%     
%     figure('Name','Quantity of interest')
%     clf
%     plot(boutput,'-b','LineWidth',1);
%     grid on
%     box on
%     set(gca,'FontSize',16)
%     xlabel('Time (s)')
%     ylabel('Quantity of interest')
%     mysaveas(pathname,'quantity_of_interest',formats,renderer);
%     mymatlab2tikz(pathname,'quantity_of_interest.tex');
%     
%     Ioutput = integrate(boutput);
%     fprintf('quantity of interest = %e\n',Ioutput);
end
