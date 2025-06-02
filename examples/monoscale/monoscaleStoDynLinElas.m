%% Monoscale stochastic linear elasticity dynamic problem %%
%%--------------------------------------------------------%%

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

filename = 'dynLinElas';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','monoscaleSto',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
fontsize = 16;
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    D = DOMAIN(2,[0.0,0.0],[2.0,0.4]);
    
    elemtype = 'TRI3';
    % elemtype = 'QUA4';
    % option = 'DEFO'; % plane strain
    option = 'CONT'; % plane stress
    nbelem = [60,12];
    pb.S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'option',option);
%     cl = 0.03;
%     pb.S = build_model(D,'cl',cl,'elemtype',elemtype,'option',option,'filename',fullfile(pathname,'gmsh_domain'));
    
    %% Random variables
    d = 2; % parametric dimension
    v = UniformRandomVariable(0,1);
    rv = RandomVector(v,d);
    
    %% Materials
    % Poisson ratio
    NU = 0;
    % Thickness
    DIM3 = 1;
    
    % IntegrationRule
    p = 1;
    bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
    bases = FunctionalBases(bases);
    H = FullTensorProductFunctionalBasis(bases);
    I = gaussIntegrationRule(v,2);
    I = I.tensorize(d);
    
    % Density
    % RHO(xi) = 1 + 0.05 * (2 * xi - 1)
    fun = @(xi) 1 + 0.05 * (2 * xi(:,1) - 1);
    fun = UserDefinedFunction(fun,d);
    fun.evaluationAtMultiplePoints = true;
    
    RHO = H.projection(fun,I);
    
    % Young modulus
    % E(xi) = 1 + 0.5 * (2 * xi - 1) = xi + 0.5
    fun = @(xi) xi(:,2) + 0.5;
    fun = UserDefinedFunction(fun,d);
    fun.evaluationAtMultiplePoints = true;
    
    E = H.projection(fun,I);
    
    % Material
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
    mat = setnumber(mat,1);
    pb.S = setmaterial(pb.S,mat);
    
    %% Dirichlet boundary conditions
    L1 = LIGNE(getvertex(D,1),getvertex(D,4));
    L2 = LIGNE(getvertex(D,2),getvertex(D,3));
    
    pb.S = final(pb.S);
    pb.S = addcl(pb.S,L1);
    
    %% Initial conditions
    pb.u0 = zeros(getnbddlfree(pb.S),1);
    pb.v0 = zeros(getnbddlfree(pb.S),1);
    
    %% Time scheme
    t0 = 0;
    t1 = 2;
    nt = 50;
    T = TIMEMODEL(t0,t1,nt);
    
    pb.timeSolver = NEWMARKSOLVER(T,'alpha',0.05,'display',false);
    % pb.timeSolver = DGTIMESOLVER(T,1,'outputsplit',true,'display',false,'lu',true);
    pb.timeOrder = 2;
    
    tc = get(T,'t1')/6;
    pb.loadFunction = @(N) rampe(N,t0,tc);
    % pb.loadFunction = @(N) dirac(N,t0,tc);
    % pb.loadFunction = @(N) one(N);
    
    %% Mass, stiffness and damping matrices and sollicitation vectors
    if ~israndom(pb.S)
        pb.M = calc_mass(pb.S);
        pb.A = calc_rigi(pb.S);
    end
    
    b = surfload(pb.S,L2,'FX',-1);
    pb.b = b*pb.loadFunction(pb.timeSolver);
    
    save(fullfile(pathname,'problem.mat'),'pb','D','L1','L2','rv');
else
    load(fullfile(pathname,'problem.mat'),'pb','D','L1','L2','rv');
end

%% Solution
if solveProblem
    p = 50;
    bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
    bases = FunctionalBases(bases);
    
    s = AdaptiveSparseTensorAlgorithm();
    s.tol = 1e-2;
    s.tolStagnation = 1e-1;
    s.display = true;
    s.displayIterations = true;
    
    ls = LinearModelLearningSquareLoss();
    ls.errorEstimation = true;
    ls.sharedCoefficients = false;
    
    funt = @(xi) solveSystem(calcOperator(funEval(pb,xi)));
    funt = UserDefinedFunction(funt,d,getnbddlfree(pb.S)*getnbtimedof(pb.timeSolver));
    funt.evaluationAtMultiplePoints = false;
    
    t = tic;
    [ut,errt,~,yt] = s.leastSquares(funt,bases,ls,rv);
    time = toc(t);
    
    save(fullfile(pathname,'solution.mat'),'ut','errt','yt','funt','time');
else
    load(fullfile(pathname,'solution.mat'),'ut','errt','yt','funt','time');
end

%% Outputs
fprintf('\n');
fprintf('nb elements = %g\n',getnbelem(pb.S));
fprintf('nb nodes    = %g\n',getnbnode(pb.S));
fprintf('nb dofs     = %g\n',getnbddl(pb.S));
fprintf('time solver : %s\n',class(pb.timeSolver));
fprintf('nb time steps = %g\n',getnt(pb.timeSolver));
fprintf('nb time dofs  = %g\n',getnbtimedof(pb.timeSolver));

T = gettimemodel(pb.timeSolver);
sz = [getnbddlfree(pb.S),getnbtimedof(T)];

utoutput = ut;
ut.data = reshape(ut.data,[cardinal(ut.basis),sz]);
ut.sz = sz;

fprintf('\n');
fprintf('spatial dimension = %d\n',ut.sz(1))
fprintf('time dimension = %d\n',ut.sz(2))
fprintf('parametric dimension = %d\n',ndims(ut.basis))
fprintf('basis dimension = %d\n',cardinal(ut.basis))
fprintf('order = [ %s ]\n',num2str(max(ut.basis.indices.array)))
% fprintf('multi-index set = \n')
% disp(ut.basis.indices.array)
fprintf('nb samples = %d\n',size(yt,1))
fprintf('CV error = %d\n',norm(errt))
fprintf('elapsed time = %f s\n',time)

%% Test
if testSolution
    N = 100;
    errL2t = testError(utoutput,funt,N,rv);
    save(fullfile(pathname,'test.mat'),'errL2t');
else
    load(fullfile(pathname,'test.mat'),'errL2t');
end
fprintf('mean squared error = %d\n',errL2t)

%% Display
if displaySolution
    %% Display domains and meshes
    figure('Name','Domain')
    clf
    h1 = plot(D,'FaceColor',getfacecolor(1));
    hold on
    h2 = plot(L1,'EdgeColor',getfacecolor(5));
    h3 = plot(L2,'EdgeColor',getfacecolor(6));
    hold off
    set(gca,'FontSize',fontsize)
    legend([h1(1),h2(1),h3(1)],'$\Omega$','$\Gamma_D$','$\Gamma_N$','Interpreter','latex')
    axis image
    axis off
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    figure('Name','Mesh')
    clf
    h1 = plot(pb.S,'FaceColor',getfacecolor(1));
    hold on
    h2 = plot(L1,'EdgeColor',getfacecolor(5));
    h3 = plot(L2,'EdgeColor',getfacecolor(6));
    hold off
    set(gca,'FontSize',fontsize)
    legend([h1(1),h2(1),h3(1)],'$\Omega$','$\Omega_2$','$\Omega_3$','$\Gamma_D^1$','$\Gamma_D^2$','Interpreter','latex')
    mysaveas(pathname,'mesh',formats,renderer);
    
    %% Display multi-index set
    plotMultiIndexSet(ut,'legend',false);
    mysaveas(pathname,'multi_index_set','fig');
    mymatlab2tikz(pathname,'multi_index_set.tex');
    
    %% Display evolution of statistical outputs
    T = gettimemodel(pb.timeSolver);
    
    i = 1;
%     for i=1:2
        evolMean(pb.S,T,ut,'displ',i,'filename',['mean_solution_' num2str(i)],'pathname',pathname);
%         evolMean(pb.S,T,vt,'displ',i,'filename',['mean_velocity_' num2str(i)],'pathname',pathname);
%         evolMean(pb.S,T,at,'displ',i,'filename',['mean_acceleration_' num2str(i)],'pathname',pathname);
        
        evolVariance(pb.S,T,ut,'displ',i,'filename',['var_solution_' num2str(i)],'pathname',pathname);
%         evolVariance(pb.S,T,vt,'displ',i,'filename',['var_velocity_' num2str(i)],'pathname',pathname);
%         evolVariance(pb.S,T,at,'displ',i,'filename',['var_acceleration_' num2str(i)],'pathname',pathname);
        
        evolStd(pb.S,T,ut,'displ',i,'filename',['std_solution_' num2str(i)],'pathname',pathname);
%         evolStd(pb.S,T,vt,'displ',i,'filename',['std_velocity_' num2str(i)],'pathname',pathname);
%         evolStd(pb.S,T,at,'displ',i,'filename',['std_acceleration_' num2str(i)],'pathname',pathname);
        
        d = ndims(ut.basis);
        for j=1:d
%             evolSobolIndices(pb.S,T,ut,j,'displ',i,'filename',['sobol_indices_solution_' num2str(i) '_var_' num2str(j)],'pathname',pathname);
%             evolSobolIndices(pb.S,T,vt,j,'displ',i,'filename',['sobol_indices_velocity_' num2str(i) '_var_' num2str(j)],'pathname',pathname);
%             evolSobolIndices(pb.S,T,at,j,'displ',i,'filename',['sobol_indices_acceleration_' num2str(i) '_var_' num2str(j)],'pathname',pathname);
            
            evolSensitivityIndices(pb.S,T,ut,j,'displ',i,'filename',['sensitivity_indices_solution_' num2str(i) '_var_' num2str(j)],'pathname',pathname);
%             evolSensitivityIndices(pb.S,T,vt,j,'displ',i,'filename',['sensitivity_indices_velocity_' num2str(i) '_var_' num2str(j)],'pathname',pathname);
%             evolSensitivityIndices(pb.S,T,at,j,'displ',i,'filename',['sensitivity_indices_acceleration_' num2str(i) '_var_' num2str(j)],'pathname',pathname);
        end
%     end
    
    %% Display statistical outputs at different instants
    [t,rep] = gettevol(pb.timeSolver);
    for k=1:floor(length(rep)/5):length(rep)
        close all
        uk_data = ut.data(:,:,rep(k));
%         vk_data = vt.data(:,:,rep(k));
%         ak_data = at_tot.data(:,:,rep(k));
        uk = FunctionalBasisArray(uk_data,ut.basis,ut.sz(1));
%         vk = FunctionalBasisArray(vk_data,vt.basis,vt.sz(1));
%         ak = FunctionalBasisArray(ak_data,at.basis,at.sz(1));
        
        i = 1;
        % for i=1:2
            plotMean(pb.S,uk,'displ',i);
            mysaveas(pathname,['mean_solution_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
%             plotMean(pb.S,vk,'displ',i);
%             mysaveas(pathname,['mean_velocity_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
%             plotMean(pb.S,ak,'displ',i);
%             mysaveas(pathname,['mean_acceleration_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            plotVariance(pb.S,uk,'displ',i);
            mysaveas(pathname,['var_solution_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
%             plotVariance(pb.S,vk,'displ',i);
%             mysaveas(pathname,['var_velocity_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
%             plotVariance(pb.S,ak,'displ',i);
%             mysaveas(pathname,['var_acceleration_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            plotStd(pb.S,uk,'displ',i);
            mysaveas(pathname,['std_solution_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
%             plotStd(pb.S,vk,'displ',i);
%             mysaveas(pathname,['std_velocity_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
%             plotStd(pb.S,ak,'displ',i);
%             mysaveas(pathname,['std_acceleration_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            d = ndims(uk.basis);
            for j=1:d
%                 plotSobolIndices(pb.S,uk,j,'displ',i);
%                 mysaveas(pathname,['sobol_indices_solution_' num2str(i) '_var_' num2str(j) '_t' num2str(k-1)],formats,renderer);
%                 plotSobolIndices(pb.S,vk,j,'displ',i);
%                 mysaveas(pathname,['sobol_indices_velocity_' num2str(i) '_var_' num2str(j) '_t' num2str(k-1)],formats,renderer);
%                 plotSobolIndices(pb.S,ak,j,'displ',i);
%                 mysaveas(pathname,['sobol_indices_acceleration_' num2str(i) '_var_' num2str(j) '_t' num2str(k-1)],formats,renderer);
                
                plotSensitivityIndices(pb.S,uk,j,'displ',i);
                mysaveas(pathname,['sensitivity_indices_solution_' num2str(i) '_var_' num2str(j) '_t' num2str(k-1)],formats,renderer);
%                 plotSensitivityIndices(pb.S,vk,j,'displ',i);
%                 mysaveas(pathname,['sensitivity_indices_velocity_' num2str(i) '_var_' num2str(j) '_t' num2str(k-1)],formats,renderer);
%                 plotSensitivityIndices(pb.S,ak,j,'displ',i);
%                 mysaveas(pathname,['sensitivity_indices_acceleration_' num2str(i) '_var_' num2str(j) '_t' num2str(k-1)],formats,renderer);
            end
        % end
    end
end

myparallel('stop');
