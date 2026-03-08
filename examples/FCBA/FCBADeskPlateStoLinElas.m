%% FCBA desk plate stochastic linear elasticity %%
%%----------------------------------------------%%

% clc
clearvars
close all
rng('default');
s = rng;
myparallel('start');

%% Input data
estimateParam = false;
solveProblem = true;
displaySolution = false;
displayCv = false;

% tests = {'StaticHori1'};     % strength test under static horizontal load 1
% tests = {'StaticHori2'};     % strength test under static horizontal load 2
% tests = {'StaticHori3'};     % strength test under static horizontal load 3 (lifting)
% tests = {'StaticHori4'};     % strength test under static horizontal load 4 (lifting)
% tests = {'StaticVert'};      % strength test under static vertical load
% tests = {'DurabilityHori1'}; % durability test under horizontal load 1
% tests = {'DurabilityHori2'}; % durability test under horizontal load 2
% tests = {'DurabilityHori3'}; % durability test under horizontal load 3 (lifting)
% tests = {'DurabilityHori4'}; % durability test under horizontal load 4 (lifting)
% tests = {'StabilityVert'};   % stability test under vertical load
% tests = {'Impact'};          % vertical impact test
% tests = {'Drop'};            % drop test

tests = {'StaticHori1','StaticHori2',...
    'StaticVert',...
    'DurabilityHori1','DurabilityHori2',...
    'StabilityVert'};

% elemtypes = {'DKT'};       % Kirchhoff-Love (classical) plate theory
% elemtypes = {'DST'};       % Reissner-Mindlin (first-order shear) plate theory
elemtypes = {'DKT','DST'}; % Both plate theories

pointLoad = false; % point load

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'epsc','png'};
renderer = 'OpenGL';

for it=1:length(tests)
    test = tests{it};
    switch lower(test)
        case {'statichori1','statichori2'}
            loads = [100, 200]; % point load, F1=F2=100, 200 [N]
        case {'statichori3','statichori4'}
            loads = 100; % point load, F3=F4=100 [N]
        case 'staticvert'
            loads = [300, 400, 500]; % point load, F=300, 400, 500 [N]
        case {'durabilityhori1','durabilityhori2','durabilityhori3','durabilityhori4'}
            loads = 100; % point load, F=100 [N]
        case 'stabilityvert'
            loads = 400; % point load, V=400 [N]
        case 'impact'
            loads = 180e-3; % height [m]
        case 'drop'
            loads = 100e-3; % height [m]
    end
    for il=1:length(loads)
        loading = loads(il);
        if pointLoad
            filename = ['FCBADeskPlateStoLinElas' test '_' num2str(loading) 'N_PointLoad'];
        else
            filename = ['FCBADeskPlateStoLinElas' test '_' num2str(loading) 'N'];
        end

for ie=1:length(elemtypes)
    elemtype = elemtypes{ie};
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','FCBA',filename,elemtype);
    if ~exist(pathname,'dir')
        mkdir(pathname);
    end

%% Problem
if solveProblem
    %% Domains and meshes
    % Plates dimensions
    a12 = 750e-3; % [m]
    b12 = 396e-3;
    a3 = 1006e-3;
    b3 = 501e-3;
    a5 = 940e-3;
    b5 = 113e-3;
    % Plates 1 and 2
    c = 74e-3;
    d = 147e-3;
    e = 63e-3;
    f = 30e-3;
    % Plate 5
    a = 20e-3;
    b = 30e-3;
    % Thickness
    % Same thickness for all the plates
    h = 15e-3;
    %
    x1 = (a5+h)/2;
    y1_14 = -b12+c;
    y1_23 = c;
    z1_12 = 0;
    z1_34 = a12+h/2;
    Q1 = QUADRANGLE([x1,y1_14,z1_12],[x1,y1_23,z1_12],...
                    [x1,y1_23,z1_34],[x1,y1_14,z1_34]);
    x2 = -(a5+h)/2;
    y2_14 = -b12+c;
    y2_23 = c;
    z2_12 = 0;
    z2_34 = a12+h/2;
    Q2 = QUADRANGLE([x2,y2_14,z2_12],[x2,y2_23,z2_12],...
                    [x2,y2_23,z2_34],[x2,y2_14,z2_34]);
    x3_14 = -a3/2;
    x3_23 = a3/2;
    y3_12 = c-(b3+b12)/2;
    y3_34 = c+(b3-b12)/2;
    z3 = a12+h/2;
    Q3 = QUADRANGLE([x3_14,y3_12,z3],[x3_23,y3_12,z3],...
                    [x3_23,y3_34,z3],[x3_14,y3_34,z3]);
    x5a_14 = -(a5+h)/2;
    x5a_23 = (a5+h)/2;
    y5a = 0;
    z5a_12 = a12-d-e-b;
    z5a_34 = a12-d+a;
    Q5a = QUADRANGLE([x5a_14,y5a,z5a_12],[x5a_23,y5a,z5a_12],...
                     [x5a_23,y5a,z5a_34],[x5a_14,y5a,z5a_34]);
    x5b_14 = -(a5+h)/2;
    x5b_23 = (a5+h)/2;
    y5b = 0;
    z5b_12 = f-b;
    z5b_34 = f-b+b5;
    Q5b = QUADRANGLE([x5b_14,y5b,z5b_12],[x5b_23,y5b,z5b_12],...
                     [x5b_23,y5b,z5b_34],[x5b_14,y5b,z5b_34]);
    
    % Points
    L3 = getedges(Q3);
    x_hori = {double(getcenter(L3{2})),double(getcenter(L3{4})),...
              double(getcenter(L3{3})),double(getcenter(L3{1}))};
    x_vert = double(getcenter(Q3));
    x_dura = {[x3_23,y3_12+50e-3,z3],[x3_14,y3_12+50e-3,z3],...
              [x3_23-50e-3,y3_12,z3],[x3_23-50e-3,y3_34,z3]};
    x_stab = double(getcenter(L3{1}))+[0.0,50e-3,0.0];
    x_meas = [x_hori,x_vert,x_dura,x_stab];
    P_hori = cellfun(@(x) POINT(x),x_hori,'UniformOutput',false);
    P_vert = POINT(x_vert);
    P_dura = cellfun(@(x) POINT(x),x_dura,'UniformOutput',false);
    P_stab = POINT(x_stab);
    P_meas = cellfun(@(x) POINT(x),x_meas,'UniformOutput',false);
    
    % Plates meshes
    % elemtype = 'DST';
    cl = h/3;
    cl_12 = cl;
    cl_3 = cl;
    cl_5 = cl;
    r_load = 40e-3;
    r_masse = 100e-3;
    C_masse = CIRCLE(0.0,y3_12+b3/2,z3,r_masse);
    x_masse = double(getcenter(C_masse));
    if pointLoad
        PbQ3 = {x_hori{4},x_dura{3},x_dura{1},x_hori{1},...
                x_dura{4},x_hori{3},x_hori{2},x_dura{2}};
        if ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4')
            S = gmshFCBAdeskpointload(Q1,Q2,Q3,Q5a,Q5b,C_masse,PbQ3,x_stab,x_masse,...
                cl_12,cl_12,cl_3,cl_5,cl_5,cl_3,cl_3,cl_3,cl_3,...
                fullfile(pathname,['gmsh_desk_' elemtype]),3);
        else
            S = gmshFCBAdeskpointload(Q1,Q2,Q3,Q5a,Q5b,C_masse,PbQ3,x_stab,x_masse,...
                cl_12,cl_12,cl_3,cl_5,cl_5,cl_3,cl_3,cl_3,cl_3,...
                fullfile(pathname,['gmsh_desk_' elemtype]),3,'recombine');
        end
    else
        L_hori{1} = LINE(x_hori{1}+[0,-r_load,0],x_hori{1}+[0,r_load,0]);
        L_hori{2} = LINE(x_hori{2}+[0,r_load,0],x_hori{2}+[0,-r_load,0]);
        L_hori{3} = LINE(x_hori{3}+[r_load,0,0],x_hori{3}+[-r_load,0,0]);
        L_hori{4} = LINE(x_hori{4}+[-r_load,0,0],x_hori{4}+[r_load,0,0]);
        L_dura{1} = LINE(x_dura{1}+[0,-r_load,0],x_dura{1}+[0,r_load,0]);
        L_dura{2} = LINE(x_dura{2}+[0,r_load,0],x_dura{2}+[0,-r_load,0]);
        L_dura{3} = LINE(x_dura{3}+[-r_load,0,0],x_dura{3}+[r_load,0,0]);
        L_dura{4} = LINE(x_dura{4}+[r_load,0,0],x_dura{4}+[-r_load,0,0]);
        LbQ3 = {L_hori{4},L_dura{3},L_dura{1},L_hori{1},...
                L_dura{4},L_hori{3},L_hori{2},L_dura{2}};
        C_vert = CIRCLE(x_vert(1),x_vert(2),x_vert(3),r_load);
        C_stab = CIRCLE(x_stab(1),x_stab(2),x_stab(3),r_load);
        if ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4')
            S = gmshFCBAdesk(Q1,Q2,Q3,Q5a,Q5b,C_masse,LbQ3,C_stab,C_vert,...
                cl_12,cl_12,cl_3,cl_5,cl_5,cl_3,cl_3,cl_3,cl_3,...
                fullfile(pathname,['gmsh_desk_' elemtype]),3);
        else
            S = gmshFCBAdesk(Q1,Q2,Q3,Q5a,Q5b,C_masse,LbQ3,C_stab,C_vert,...
                cl_12,cl_12,cl_3,cl_5,cl_5,cl_3,cl_3,cl_3,cl_3,...
                fullfile(pathname,['gmsh_desk_' elemtype]),3,'recombine');
        end
    end
    S = convertelem(S,elemtype);
    
    %% Random variables
    % Data
    filenameAna = 'data_ET_GL.mat';
    filenameNum = 'data_EL_NUL.mat';
    pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','identification','materialParticleBoard');
    load(fullfile(pathnameIdentification,filenameAna));
    load(fullfile(pathnameIdentification,filenameNum));
    
    % Material symmetry
    materialSym = 'isotTrans';
    
    % Parameterization
    useRedParam = true; % reduced parameterization
    
    % Number of samples
    N = 2e3;
    
    switch lower(materialSym)
        case 'isot'
            % Data
            N_data = length(mean_ET_data);
            E_data = mean_ET_data*1e-3; % [GPa]
            NU_data = 0.1+0.2*rand(N_data,1); % artificial data for NU uniformly distributed from 0.1 to 0.3
            G_data = E_data./(2*(1+NU_data)); % [GPa]
            lambda_data = E_data.*NU_data./((1+NU_data).*(1-2*NU_data)); % [GPa]
            C1_data = lambda_data + 2/3*G_data; % [GPa]
            C2_data = G_data; % [GPa]
            C_data = [C1_data(:) C2_data(:)]; % [GPa]
            
            % Empirical estimates
            mC_data = mean(C_data,1);
            % vC_data = var(C_data,0,1); % vC_data = size(C_data,1)/(size(C_data,1)-1)*moment(C_data,2,1);
            stdC_data = std(C_data,0,1); % stdC_data = sqrt(vC_data);
            cvC_data = stdC_data./mC_data;
            % sC_data = sqrt(norm(vC_data));
            % mCnorm_data = norm(mC_data);
            % dC_data = sC_data/mCnorm_data;
            % phiC_data = log(96*C_data(:,1).*C_data(:,2).^5);
            % nuC_data = mean(phiC_data,1);
            % fC_data = [mC_data nuC_data];
            
            % Initial parameter values
            la0 = mean([1-1/cvC_data(1)^2 (1-1/cvC_data(2)^2)/5]); % la < 1/5
            a01 = 1-la0;   % a1 > 0
            a02 = 1-5*la0; % a2 > 0
            la01 = a01/mC_data(1);   % la1 > 0 (match mathematical expectation with empirical mean value for C1)
            la02 = a02/mC_data(2); % la2 > 0 (match mathematical expectation with empirical mean value for C2)
            % la01 = -la0/mC_data(1);   % la1 > 0 (match mode with empirical mean value for C1)
            % la02 = -5*la0/mC_data(2); % la2 > 0 (match mode with empirical mean value for C2)
            b01 = 1/la01;  % b1 > 0
            b02 = 1/la02;  % b2 > 0
            
            % Initial parameter vector
            if useRedParam
                lambda0 = la0; % reduced parameterization
            else
                lambda0 = [la01 la02 la0]; % full parameterization
            end
            
            % Parameter estimation
            method = 'mle'; % parameter estimation method = 'mle' or 'lse'
            if estimateParam
                switch lower(method)
                    case 'mle'
                        % Maximum likelihood estimation
                        % loglfval0 = -nloglfElasIsot(lambda0,C_data);
                        % lambda = mleStoLinElasIsot(C_data,lambda0,'display','iter');
                        % loglfval = -nloglfElasIsot(lambda,C_data);
                        [lambda,loglfval,loglfval0] = mleStoLinElasIsot(C_data,lambda0,'display','iter');
                    case 'lse'
                        % Least-squares estimation
                        % lambda = lseStoLinElasIsot(C_data,lambda0,'display','iter-detailed');
                        [lambda,resnorm,residual,exitflag,output] = lseStoLinElasIsot(C_data,lambda0,'display','iter-detailed');
                end
            else
                filenameParam = 'param_lambda.mat';
                foldernameParam = 'modelStoLinElasIsot_ElasTensor_';
                if useRedParam
                    foldernameParam = [foldernameParam 'ReducedParam_'];
                else
                    foldernameParam = [foldernameParam 'FullParam_'];
                end
                foldernameParam = [foldernameParam method];
                pathnameParam = fullfile(getfemobjectoptions('path'),'MYCODE',...
                    'results','identification',foldernameParam);
                load(fullfile(pathnameParam,filenameParam));
            end
            
            % Optimal parameter values
            if useRedParam
                la  = lambda(1);     % la < 1/5
                a1  = 1-la;          % a1 > 0
                a2  = 1-5*la;        % a2 > 0
                la1 = a1/mC_data(1); % la1 > 0
                la2 = a2/mC_data(2); % la2 > 0
            else
                la1 = lambda(1); % la1 > 0
                la2 = lambda(2); % la2 > 0
                la  = lambda(3); % la < 1/5
                a1  = 1-la;      % a1 > 0
                a2  = 1-5*la;    % a2 > 0
            end
            b1 = 1/la1; % b1 > 0
            b2 = 1/la2; % b2 > 0
            
            % Sample set
            C1_sample = gamrnd(a1,b1,N,1)*1e9; % [Pa]
            C2_sample = gamrnd(a2,b2,N,1)*1e9; % [Pa]
            lambda_sample = C1_sample-2/3*C2_sample; % [Pa]
            E_sample = (9*C1_sample.*C2_sample)./(3*C1_sample+C2_sample); % [Pa]
            NU_sample = (3*C1_sample-2*C2_sample)./(6*C1_sample+2*C2_sample);
            
        case 'isottrans'
            % Data
            N_data = length(mean_ET_data);
            ET_data = mean_ET_data*1e-3; % [GPa]
            GL_data = mean_GL_data*1e-3; % [GPa]
            EL_data = mean_EL_data*1e-3; % [GPa]
            % NUL_data = mean_NUL_data;
            NUL_data = 0.03+0.03*rand(N_data,1); % artificial data for NUL uniformly distributed from 0.03 to 0.06
            NUT_data = 0.1+0.2*rand(N_data,1);   % artificial data for NUT uniformly distributed from 0.1 to 0.3
            GT_data = ET_data./(2*(1+NUT_data)); % [GPa]
            kT_data = (EL_data.*ET_data)./(2*(1-NUT_data).*EL_data-4*ET_data.*(NUL_data).^2); % [GPa]
            C1_data = EL_data + 4*(NUL_data.^2).*kT_data; % [GPa]
            C2_data = 2*kT_data; % [GPa]
            C3_data = 2*sqrt(2)*kT_data.*NUL_data; % [GPa]
            C4_data = 2*GT_data; % [GPa]
            C5_data = 2*GL_data; % [GPa]
            C_data = [C1_data(:) C2_data(:) C3_data(:) C4_data(:) C5_data(:)];
            
            % Empirical estimates
            mC_data = mean(C_data,1);
            % vC_data = var(C_data,0,1); % vC_data = size(C_data,1)/(size(C_data,1)-1)*moment(C_data,2,1);
            stdC_data = std(C_data,0,1); % stdC_data = sqrt(vC_data);
            cvC_data = stdC_data./mC_data;
            % sC_data = sqrt(norm(vC_data));
            % mCnorm_data = norm(mC_data);
            % dC_data = sC_data/mCnorm_data;
            % phiC_data = log((C_data(:,1).*C_data(:,2)-C_data(:,3).^2).*(C_data(:,4).^2).*(C_data(:,5).^2));
            % nuC_data = mean(phiC_data,1);
            % fC_data = [mC_data nuC_data];
            
            % Initial parameter values
            cvC45 = cvC_data(4:5); % empirical coefficient of variations (cv) for C4 and C5
            la0 = (1-1/mean(cvC45).^2)/2; % la < 1/2, la < 0 for MH sampling (match cv with averaged empirical cv for C4 and C5)
            la01 = -(mC_data(2)*la0)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la1 > 0 (match mode with empirical mean value for C1)
            la02 = -(mC_data(1)*la0)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la2 > 0 (match mode with empirical mean value for C2)
            la03 = (2*mC_data(3)*la0)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la3 in R such that 2*sqrt(la1*la2)-la3 > 0 (match mode with empirical mean value for C3)
            a0 = 1-2*la0; % a > 0
            la04 = a0/mC_data(4); % la4 > 0 (match mathematical expectation with empirical mean value for C4)
            la05 = a0/mC_data(5); % la5 > 0 (match mathematical expectation with empirical mean value for C5)
            % la04 = -2*la0/mC_data(4); % la4 > 0 (match mode with empirical mean value for C4)
            % la05 = -2*la0/mC_data(5); % la5 > 0 (match mode with empirical mean value for C5)
            b04 = 1/la04; % b4 > 0
            b05 = 1/la05; % b5 > 0
            
            % Initial parameter vector
            if useRedParam
                lambda0 = [la01 la02 la03 la0]; % reduced parameterization
            else
                lambda0 = [la01 la02 la03 la04 la05 la0]; % full parameterization
            end
            
            % Markov Chain Monte Carlo (MCMC) method
            MCMCalgo = 'IMH'; % algorithm for MCMC method = 'IMH', 'RWMH' or 'SS'
            useParallel = true; % parallel computing to evaluate gradients of objective function
            
            % Least-squares estimation with MCMC method
            if estimateParam
                Nconv = 1e6; % number of samples to ensure convergence of statistical estimates
                % lambda = lseStoLinElasIsotTrans(C_data,lambda0,Nconv,MCMCalgo,s,'display','iter-detailed','useParallel',useParallel);
                [lambda,fval,exitflag,output] = lseStoLinElasIsotTrans(C_data,lambda0,Nconv,MCMCalgo,s,'display','iter-detailed','useParallel',useParallel);
            else
                filenameParam = 'param_lambda.mat';
                foldernameParam = 'modelStoLinElasIsotTrans_ElasTensor_';
                if useRedParam
                    foldernameParam = [foldernameParam 'ReducedParam_'];
                else
                    foldernameParam = [foldernameParam 'FullParam_'];
                end
                foldernameParam = [foldernameParam MCMCalgo];
                pathnameParam = fullfile(getfemobjectoptions('path'),'MYCODE',...
                    'results','identification',foldernameParam);
                load(fullfile(pathnameParam,filenameParam));
            end
            
            % Optimal parameter values
            la1 = lambda(1); % la1 > 0
            la2 = lambda(2); % la2 > 0
            la3 = lambda(3); % la3 in R such that 2*sqrt(la1*la2)-la3 > 0
            if useRedParam
                la  = lambda(4);    % la < 1/2, la < 0 for MH sampling using a trivariate normal proposal pdf with positive-definite covariance matrix Sigma
                a   = 1-2*la;       % a > 0
                la4 = a/mC_data(4); % la4 > 0
                la5 = a/mC_data(5); % la5 > 0
            else
                la4 = lambda(4); % la4 > 0
                la5 = lambda(5); % la5 > 0
                la  = lambda(6); % la < 1/2, la < 0 for MH sampling using a trivariate normal proposal pdf with positive-definite covariance matrix Sigma
                a   = 1-2*la;    % a > 0
            end
            b4 = 1/la4; % b4 > 0
            b5 = 1/la5; % b5 > 0
            
            % Sample generation
            rng(s); % initialize the random number generator using the default settings contained in s
            switch lower(MCMCalgo)
                case {'imh','rwmh'}
                    [C123_sample,accept] = mhsampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),N,MCMCalgo);
                case 'ss'
                    [C123_sample,neval] = slicesampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),N);
                otherwise
                    error(['MCMC algorithm ' MCMC ' not implemented'])
            end
            C1_sample = C123_sample(:,1)*1e9; % [Pa]
            C2_sample = C123_sample(:,2)*1e9; % [Pa]
            C3_sample = C123_sample(:,3)*1e9; % [Pa]
            C4_sample = gamrnd(a,b4,N,1)*1e9; % [Pa]
            C5_sample = gamrnd(a,b5,N,1)*1e9; % [Pa]
            
            % Sample set
            kT_sample = C2_sample/2; % [Pa]
            NUL_sample = (C3_sample./kT_sample)/(2*sqrt(2));
            EL_sample = C1_sample - 4*(NUL_sample.^2).*kT_sample; % [Pa]
            GT_sample = C4_sample/2; % [Pa]
            GL_sample = C5_sample/2; % [Pa]
            ET_sample = 4./(1./kT_sample+1./GT_sample+4*(NUL_sample.^2)./EL_sample); % [Pa]
            NUT_sample = (ET_sample./GT_sample)/2-1;
            
        otherwise
            error('Wrong material symmetry !')
    end
    
    %% Materials
    % Gravitational acceleration
    g = 9.81; % [m/s2]
    
    % Density
    Mass_total = 13.9; % [kg]
    Vol_total = h*(a12*b12*2+a3*b3+a5*b5*2);
    RHO = Mass_total/(Vol_total); % [kg/m3]
    
    % Material symmetry
    switch lower(materialSym)
        case 'isot'
            % Young modulus
            E = mean(E_sample);
            % Poisson ratio
            NU = mean(NU_sample);
            % Material
            mat = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',h,'k',5/6);
        case 'isottrans'
            % Transverse Young modulus
            ET = mean(ET_sample);
            % Longitudinal shear modulus
            GL = mean(GL_sample);
            % Longitudinal Young modulus
            % EL = mean(EL_sample);
            % Longitudinal Poisson ratio
            % NUL = mean(NUL_sample);
            % Transverse Poisson ratio
            NUT = mean(NUT_sample);
            % Material
            mat = ELAS_SHELL_ISOT_TRANS('ET',ET,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',h,'k',5/6);
        otherwise
            error('Wrong material symmetry !')
    end
    mat = setnumber(mat,1);
    S = setmaterial(S,mat);
    
    %% Neumann boundary conditions
    p_plate = RHO*g*h; % surface load (body load for plates) [N/m2]
    L_hori_dura = 2*r_load;
    Sec_vert_stab = pi*r_load^2;
    switch lower(test)
        case {'statichori1','statichori2','statichori3','statichori4'}
            masse = 50.5; % [kg]
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse; % surface load (body load for plates) [N/m2]
            if pointLoad
                p = loading; % point load [N]
            else
                p = loading/L_hori_dura; % line load (surface load for plates) [N/m]
            end
            slope = 0;
        case 'staticvert'
            if pointLoad
                p = loading; % point load [N]
            else
                p = loading/Sec_vert_stab; % surface load (body load for plates) [N/m2]
            end
        case {'durabilityhori1','durabilityhori2','durabilityhori3','durabilityhori4'}
            masse = 50.5; % [kg]
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse; % surface load (body load for plates) [N/m2]
            if pointLoad
                p = loading; % point load [N]
            else
                p = loading/L_hori_dura; % line load (surface load for plates) [N/m]
            end
        case 'stabilityvert'
            if pointLoad
                p = loading; % point load [N]
            else
                p = loading/Sec_vert_stab; % surface load (body load for plates) [N/m2]
            end
        case 'impact'
            p = loading; % height [m]
        case 'drop'
            p = loading; % height [m]
    end
    
    %% Dirichlet boundary conditions
    S = final(S);
    L1 = getedge(Q1,1);
    L2 = getedge(Q2,1);
    L5b = getedge(Q5b,1);
    [~,numnode1] = intersect(S,L1);
    [~,numnode2] = intersect(S,L2);
    [~,numnode5b] = intersect(S,L5b);
    switch lower(test)
        case {'statichori1','statichori2'}
            S = addcl(S,numnode2);
            S = addcl(S,union(numnode1,numnode5b),{'UY','UZ'});
        case {'statichori3','statichori4'}
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
        case 'staticvert'
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
        case {'durabilityhori1','durabilityhori2','durabilityhori3','durabilityhori4'}
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
        case 'stabilityvert'
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
        case {'impact','drop'}
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
    end
    
    %% Sollicitation vector
    switch lower(test)
        case {'statichori1','statichori2','statichori3','statichori4'}
            if strcmpi(test,'statichori1')
                if pointLoad
                    f = nodalload(S,P_hori{1},{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
                    if isempty(ispointin(P_hori{1},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{1},{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
                end
            elseif strcmpi(test,'statichori2')
                if pointLoad
                    f = nodalload(S,P_hori{2},{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
                    if isempty(ispointin(P_hori{2},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{2},{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
                end
            elseif strcmpi(test,'statichori3')
                if pointLoad
                    f = nodalload(S,P_hori{3},{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
                    if isempty(ispointin(P_hori{3},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{3},{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
                end
            elseif strcmpi(test,'statichori4')
                if pointLoad
                    f = nodalload(S,P_hori{4},{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
                    if isempty(ispointin(P_hori{4},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{4},{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
                end
            end
            if pointLoad
                f = f + bodyload(keepgroupelem(S,4),[],'FZ',-p_masse);
            else
                f = f + bodyload(keepgroupelem(S,[4,5]),[],'FZ',-p_masse);
            end
        case 'staticvert'
            if pointLoad
                f = nodalload(S,P_vert,'FZ',-p);
                if isempty(ispointin(P_vert,POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            else
                f = bodyload(keepgroupelem(S,4),[],'FZ',-p);
            end
        case {'durabilityhori1','durabilityhori2','durabilityhori3','durabilityhori4'}
            if strcmpi(test,'durabilityhori1')
                if pointLoad
                    f = nodalload(S,P_dura{1},'FX',-p);
                    if isempty(ispointin(P_dura{1},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_dura{1},'FX',-p);
                end
            elseif strcmpi(test,'durabilityhori2')
                if pointLoad
                    f = nodalload(S,P_dura{2},'FX',p);
                    if isempty(ispointin(P_dura{2},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_dura{2},'FX',p);
                end
            elseif strcmpi(test,'durabilityhori3')
                if pointLoad
                    f = nodalload(S,P_dura{3},'FY',p);
                    if isempty(ispointin(P_dura{3},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_dura{3},'FY',p);
                end
            elseif strcmpi(test,'durabilityhori4')
                if pointLoad
                    f = nodalload(S,P_dura{4},'FY',-p);
                    if isempty(ispointin(P_dura{4},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_dura{4},'FY',-p);
                end
            end
            if pointLoad
                f = f + bodyload(keepgroupelem(S,4),[],'FZ',-p_masse);
            else
                f = f + bodyload(keepgroupelem(S,[4,5]),[],'FZ',-p_masse);
            end
        case 'stabilityvert'
            if pointLoad
                f = nodalload(S,P_stab,'FZ',-p);
                if isempty(ispointin(P_stab,POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            else
                f = bodyload(keepgroupelem(S,6),[],'FZ',-p);
            end
        case {'impact','drop'}
            error('Not implemented')
    end
    f = f + bodyload(S,[],'FZ',-p_plate);
    
    %% Stiffness matrix and solution
    t = tic;
    u = sparse(getnbddlfree(S),N);
    if ~verLessThan('matlab','9.2') % introduced in R2017a
        q = parallel.pool.DataQueue;
        afterEach(q, @(~) nUpdateProgressBar(N));
    else
        q = [];
    end
    textprogressbar('Solving problem: ');
    parfor i=1:N
        if ~verLessThan('matlab','9.2') % introduced in R2017a
            send(q,1); % send(q,1) is enough, since the callback only needs a signal that one iteration has finished, and not the iteration index i
            % send(q,i);
        end
        mati = mat;
        switch lower(materialSym)
            case 'isot'
                % Young modulus
                Ei = E_sample(i);
                % Poisson ratio
                NUi = NU_sample(i);
                % Material
                mati = setparam(mati,'E',Ei);
                mati = setparam(mati,'NU',NUi);
            case 'isottrans'
                % Transverse Young modulus
                ETi = ET_sample(i);
                % Longitudinal shear modulus
                GLi = GL_sample(i);
                % Longitudinal Young modulus
                % ELi = EL_sample(i);
                % Longitudinal Poisson ratio
                % NULi = NUL_sample(i);
                % Transverse Poisson ratio
                % NUTi = 0.2;
                NUTi = NUT_sample(i);
                % Material
                mati = setparam(mati,'ET',ETi);
                mati = setparam(mati,'GL',GLi);
                mati = setparam(mati,'NUT',NUTi);
        end
        Si = setmaterial(S,mati);
        % Stiffness matrix
        Ai = calc_rigi(Si);
        % Solution
        u(:,i) = Ai\f;
    end
    textprogressbar(' done');
    time = toc(t);
    
    mean_u = mean(u,2);
    mean_u = unfreevector(S,mean_u);
    
    std_u = std(u,0,2);
    std_u = unfreevector(S,std_u);
    
    probs = [0.025 0.975];
    % probs = [0.01 0.99];
    ci_u = quantile(u,probs,2);
    ci_u = unfreevector(S,ci_u);
    
    % mean_U = mean_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    % mean_Ux = mean_u(findddl(S,'UX'),:);
    % mean_Uy = mean_u(findddl(S,'UY'),:);
    % mean_Uz = mean_u(findddl(S,'UZ'),:);
    
    % mean_R = mean_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    % mean_Rx = mean_u(findddl(S,'RX'),:);
    % mean_Ry = mean_u(findddl(S,'RY'),:);
    % mean_Rz = mean_u(findddl(S,'RZ'),:);
    
    % std_U = std_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    % std_Ux = std_u(findddl(S,'UX'),:);
    % std_Uy = std_u(findddl(S,'UY'),:);
    % std_Uz = std_u(findddl(S,'UZ'),:);
    
    % std_R = std_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    % std_Rx = std_u(findddl(S,'RX'),:);
    % std_Ry = std_u(findddl(S,'RY'),:);
    % std_Rz = std_u(findddl(S,'RZ'),:);
    
    % ci_U = ci_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    % ci_Ux = ci_u(findddl(S,'UX'),:);
    % ci_Uy = ci_u(findddl(S,'UY'),:);
    % ci_Uz = ci_u(findddl(S,'UZ'),:);
    
    % ci_R = ci_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    % ci_Rx = ci_u(findddl(S,'RX'),:);
    % ci_Ry = ci_u(findddl(S,'RY'),:);
    % ci_Rz = ci_u(findddl(S,'RZ'),:);
    
    %% Test solution
    switch lower(test)
        case 'statichori1'
            P = P_hori{2};
        case 'statichori2'
            P = P_hori{1};
        case 'statichori3'
            P = P_hori{4};
        case 'statichori4'
            P = P_hori{3};
        case 'staticvert'
            P = P_vert;
        case 'durabilityhori1'
            P = P_dura{2};
        case 'durabilityhori2'
            P = P_dura{1};
        case 'durabilityhori3'
            P = P_dura{4};
        case 'durabilityhori4'
            P = P_dura{3};
        case 'stabilityvert'
            P = P_stab;
    end
    mean_ux = eval_sol(S,mean_u,P,'UX');
    mean_uy = eval_sol(S,mean_u,P,'UY');
    mean_uz = eval_sol(S,mean_u,P,'UZ');
    
    mean_rx = eval_sol(S,mean_u,P,'RX');
    mean_ry = eval_sol(S,mean_u,P,'RY');
    mean_rz = eval_sol(S,mean_u,P,'RZ');
    
    std_ux = eval_sol(S,std_u,P,'UX');
    std_uy = eval_sol(S,std_u,P,'UY');
    std_uz = eval_sol(S,std_u,P,'UZ');
    
    std_rx = eval_sol(S,std_u,P,'RX');
    std_ry = eval_sol(S,std_u,P,'RY');
    std_rz = eval_sol(S,std_u,P,'RZ');
    
    ci_ux = eval_sol(S,ci_u,P,'UX');
    ci_uy = eval_sol(S,ci_u,P,'UY');
    ci_uz = eval_sol(S,ci_u,P,'UZ');
    
    ci_rx = eval_sol(S,ci_u,P,'RX');
    ci_ry = eval_sol(S,ci_u,P,'RY');
    ci_rz = eval_sol(S,ci_u,P,'RZ');
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','elemtype','loading','f','N',...
        'a12','b12','a3','b3','a5','b5','h');
    save(fullfile(pathname,'solution.mat'),'u','mean_u','std_u','ci_u','probs','time');
    save(fullfile(pathname,'test_solution.mat'),'P',...
        'mean_ux','mean_uy','mean_uz',...
        'mean_rx','mean_ry','mean_rz',...
        'std_ux','std_uy','std_uz',...
        'std_rx','std_ry','std_rz',...
        'ci_ux','ci_uy','ci_uz',...
        'ci_rx','ci_ry','ci_rz');
else
    load(fullfile(pathname,'problem.mat'),'S','elemtype','loading','f','N',...
        'a12','b12','a3','b3','a5','b5','h');
    load(fullfile(pathname,'solution.mat'),'u','mean_u','std_u','ci_u','probs','time');
    load(fullfile(pathname,'test_solution.mat'),'P',...
        'mean_ux','mean_uy','mean_uz',...
        'mean_rx','mean_ry','mean_rz',...
        'std_ux','std_uy','std_uz',...
        'std_rx','std_ry','std_rz',...
        'ci_ux','ci_uy','ci_uz',...
        'ci_rx','ci_ry','ci_rz');
end

%% Outputs
filenameResults = fullfile(pathname,'results.txt');
fid = fopen(filenameResults,'w');
fprintf(fid,'2D Plate Desk\n');
fprintf(fid,'\n');
fprintf(fid,'test = %s\n',test);
fprintf(fid,'load = %g N\n',loading);
fprintf(fid,'mesh = %s elements\n',elemtype);
fprintf(fid,'nb elements = %g\n',getnbelem(S));
fprintf(fid,'nb nodes    = %g\n',getnbnode(S));
fprintf(fid,'nb dofs     = %g\n',getnbddl(S));
fprintf(fid,'span-to-thickness ratio of plates 1 and 2 = %g\n',min(a12,b12)/h);
fprintf(fid,'span-to-thickness ratio of plate 3 = %g\n',min(a3,b3)/h);
fprintf(fid,'span-to-thickness ratio of plates 5a and 5b = %g\n',min(a5,b5)/h);
fprintf(fid,'nb samples = %g\n',N);
fprintf(fid,'elapsed time = %f s\n',time);

switch lower(test)
    case 'statichori1'
        if loading==100
            ux_exp_start = -6.88*1e-3;
            ux_exp_end = -[10.5 10.51 10.44 10.8 10.72 10.62 10.67 10.65 10.66 10.87 10.86]*1e-3;
        elseif loading==200
            ux_exp_start = -6.16*1e-3;
            ux_exp_end = -[16.78 16.74 16.72 17.13 17 16.8 16.87 16.78 17.04 16.82 16.71 17.17]*1e-3;
        end
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(mean_ux-ux_exp)/norm(ux_exp);
        env_ux_exp = [min(ux_exp_end - ux_exp_start),max(ux_exp_end - ux_exp_start)];
        
        fprintf(fid,'\n');
        fprintf(fid,'Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf(fid,'mean(ux) = %g m, std(ux) = %g m, ci(ux) = [%g %g] m\n',mean_ux,std_ux,ci_ux(1),ci_ux(2));
        fprintf(fid,'ux_exp   = %g m, error(mean(ux)) = %.3e, env(ux_exp) = [%g %g] m\n',ux_exp,err_ux,env_ux_exp);
        fprintf(fid,'mean(uy) = %g m, std(uy) = %g m, ci(uy) = [%g %g] m\n',mean_uy,std_uy,ci_uy(1),ci_uy(2));
        fprintf(fid,'mean(uz) = %g m, std(uz) = %g m, ci(uz) = [%g %g] m\n',mean_uz,std_uz,ci_uz(1),ci_uz(2));
    case 'statichori2'
        if loading==100
            ux_exp_start = 2.12*1e-3;
            ux_exp_end = [6.22 6.17 6.26 6.31 6.33 6.24 6.26 6.4 6.26 6.49 6.48 6.42 6.36 6.56 6.37 6.39]*1e-3;
        elseif loading==200
            ux_exp_start = 1.91*1e-3;
            ux_exp_end = [12.45 12.68 12.66 12.65 12.71 12.64 12.82 12.73 12.89 12.86 12.79 12.86]*1e-3;
        end
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(mean_ux-ux_exp)/norm(ux_exp);
        env_ux_exp = [min(ux_exp_end - ux_exp_start),max(ux_exp_end - ux_exp_start)];
        
        fprintf(fid,'\n');
        fprintf(fid,'Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf(fid,'mean(ux) = %g m, std(ux) = %g m, ci(ux) = [%g %g] m\n',mean_ux,std_ux,ci_ux(1),ci_ux(2));
        fprintf(fid,'ux_exp   = %g m, error(mean(ux)) = %.3e, env(ux_exp) = [%g %g] m\n',ux_exp,err_ux,env_ux_exp);
        fprintf(fid,'mean(uy) = %g m, std(uy) = %g m, ci(uy) = [%g %g] m\n',mean_uy,std_uy,ci_uy(1),ci_uy(2));
        fprintf(fid,'mean(uz) = %g m, std(uz) = %g m, ci(uz) = [%g %g] m\n',mean_uz,std_uz,ci_uz(1),ci_uz(2));
    case 'statichori3'
        uy_exp_start = -3.77*1e-3;
        uy_exp_end = -[4.71 4.73 4.69 4.56 4.47 4.73]*1e-3;
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(mean_uy-uy_exp)/norm(uy_exp);
        env_uy_exp = [min(uy_exp_end - uy_exp_start),max(uy_exp_end - uy_exp_start)];
        
        fprintf(fid,'\n');
        fprintf(fid,'Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf(fid,'mean(ux) = %g m, std(ux) = %g m, ci(ux) = [%g %g] m\n',mean_ux,std_ux,ci_ux(1),ci_ux(2));
        fprintf(fid,'mean(uy) = %g m, std(uy) = %g m, ci(uy) = [%g %g] m\n',mean_uy,std_uy,ci_uy(1),ci_uy(2));
        fprintf(fid,'uy_exp   = %g m, error(mean(uy)) = %.3e, env(uy_exp) = [%g %g] m\n',uy_exp,err_uy,env_uy_exp);
        fprintf(fid,'mean(uz) = %g m, std(uz) = %g m, ci(uz) = [%g %g] m\n',mean_uz,std_uz,ci_uz(1),ci_uz(2));
    case 'statichori4'
        uy_exp_start = 9.71*1e-3;
        uy_exp_end = [12.21 12.2 12.2 12.23 12.2 12.19 12.21]*1e-3;
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(mean_uy-uy_exp)/norm(uy_exp);
        env_uy_exp = [min(uy_exp_end - uy_exp_start),max(uy_exp_end - uy_exp_start)];
        
        fprintf(fid,'\n');
        fprintf(fid,'Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf(fid,'mean(ux) = %g m, std(ux) = %g m, ci(ux) = [%g %g] m\n',mean_ux,std_ux,ci_ux(1),ci_ux(2));
        fprintf(fid,'mean(uy) = %g m, std(uy) = %g m, ci(uy) = [%g %g] m\n',mean_uy,std_uy,ci_uy(1),ci_uy(2));
        fprintf(fid,'uy_exp   = %g m, error(mean(uy)) = %.3e, env(uy_exp) = [%g %g] m\n',uy_exp,err_uy,env_uy_exp);
        fprintf(fid,'mean(uz) = %g m, std(uz) = %g m, ci(uz) = [%g %g] m\n',mean_uz,std_uz,ci_uz(1),ci_uz(2));
    case 'staticvert'
        if loading==300
            uz_exp_start = -0.69*1e-3;
            uz_exp_end = -[10.10 9.88 9.64 9.88 9.94 9.79 9.92 9.93 9.82 9.95]*1e-3;
        elseif loading==400
            uz_exp_start = -0.75*1e-3;
            uz_exp_end = -[13.45 13.52 13.56 13.64 13.65 13.74 13.75 13.44 13.74 13.53]*1e-3;
        elseif loading==500
            uz_exp_start = -0.78*1e-3;
            uz_exp_end = -[16.66 16.57 16.59 16.78 16.55 16.69 16.75 16.59 16.73 16.76]*1e-3;
        end
        uz_exp = mean(uz_exp_end - uz_exp_start);
        err_uz = norm(mean_uz-uz_exp)/norm(uz_exp);
        env_uz_exp = [min(uz_exp_end - uz_exp_start),max(uz_exp_end - uz_exp_start)];
        
        fprintf(fid,'\n');
        fprintf(fid,'Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf(fid,'mean(ux) = %g m, std(ux) = %g m, ci(ux) = [%g %g] m\n',mean_ux,std_ux,ci_ux(1),ci_ux(2));
        fprintf(fid,'mean(uy) = %g m, std(uy) = %g m, ci(uy) = [%g %g] m\n',mean_uy,std_uy,ci_uy(1),ci_uy(2));
        fprintf(fid,'mean(uz) = %g m, std(uz) = %g m, ci(uz) = [%g %g] m\n',mean_uz,std_uz,ci_uz(1),ci_uz(2));
        fprintf(fid,'uz_exp   = %g m, error(mean(uz)) = %.3e, env(uz_exp) = [%g %g] m\n',uz_exp,err_uz,env_uz_exp);
    case 'durabilityhori1'
        ux_exp_start = -4.42*1e-3;
        ux_exp_end = -[8.4 8.3 8.37 8.41 8.54 8.39 8.56 8.48 8.46 8.49 8.49 8.43 8.55 8.52]*1e-3;
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(mean_ux-ux_exp)/norm(ux_exp);
        env_ux_exp = [min(ux_exp_end - ux_exp_start),max(ux_exp_end - ux_exp_start)];
        
        fprintf(fid,'\n');
        fprintf(fid,'Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf(fid,'mean(ux) = %g m, std(ux) = %g m, ci(ux) = [%g %g] m\n',mean_ux,std_ux,ci_ux(1),ci_ux(2));
        fprintf(fid,'ux_exp   = %g m, error(mean(ux)) = %.3e, env(ux_exp) = [%g %g] m\n',ux_exp,err_ux,env_ux_exp);
        fprintf(fid,'mean(uy) = %g m, std(uy) = %g m, ci(uy) = [%g %g] m\n',mean_uy,std_uy,ci_uy(1),ci_uy(2));
        fprintf(fid,'mean(uz) = %g m, std(uz) = %g m, ci(uz) = [%g %g] m\n',mean_uz,std_uz,ci_uz(1),ci_uz(2));
    case 'durabilityhori2'
        ux_exp_start = 3.48*1e-3;
        ux_exp_end = [7.89 7.85 8.1 8.4 8.36 8.55 8.27 8.27 8.47 8.49 8.64 8.35 8.5 8.63 8.73]*1e-3;
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(mean_ux-ux_exp)/norm(ux_exp);
        env_ux_exp = [min(ux_exp_end - ux_exp_start),max(ux_exp_end - ux_exp_start)];
        
        fprintf(fid,'\n');
        fprintf(fid,'Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf(fid,'mean(ux) = %g m, std(ux) = %g m, ci(ux) = [%g %g] m\n',mean_ux,std_ux,ci_ux(1),ci_ux(2));
        fprintf(fid,'ux_exp   = %g m, error(mean(ux)) = %.3e, env(ux_exp) = [%g %g] m\n',ux_exp,err_ux,env_ux_exp);
        fprintf(fid,'mean(uy) = %g m, std(uy) = %g m, ci(uy) = [%g %g] m\n',mean_uy,std_uy,ci_uy(1),ci_uy(2));
        fprintf(fid,'mean(uz) = %g m, std(uz) = %g m, ci(uz) = [%g %g] m\n',mean_uz,std_uz,ci_uz(1),ci_uz(2));
    case 'durabilityhori3'
        uy_exp_start = 3.35*1e-3;
        uy_exp_end = [6.16 5.76 5.97 5.81 5.84 5.61 5.86 5.64 5.62 5.68]*1e-3;
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(mean_uy-uy_exp)/norm(uy_exp);
        env_uy_exp = [min(uy_exp_end - uy_exp_start),max(uy_exp_end - uy_exp_start)];
        
        fprintf(fid,'\n');
        fprintf(fid,'Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf(fid,'mean(ux) = %g m, std(ux) = %g m, ci(ux) = [%g %g] m\n',mean_ux,std_ux,ci_ux(1),ci_ux(2));
        fprintf(fid,'mean(uy) = %g m, std(uy) = %g m, ci(uy) = [%g %g] m\n',mean_uy,std_uy,ci_uy(1),ci_uy(2));
        fprintf(fid,'uy_exp   = %g m, error(mean(uy)) = %.3e, env(uy_exp) = [%g %g] m\n',uy_exp,err_uy,env_uy_exp);
        fprintf(fid,'mean(uz) = %g m, std(uz) = %g m, ci(uz) = [%g %g] m\n',mean_uz,std_uz,ci_uz(1),ci_uz(2));
    case 'durabilityhori4'
        uy_exp_start = -3.75*1e-3;
        uy_exp_end = -[3.89 3.88 3.89 3.88 3.89]*1e-3;
        uy_exp = mean(uy_exp_end - uy_exp_start);
        err_uy = norm(mean_uy-uy_exp)/norm(uy_exp);
        env_uy_exp = [min(uy_exp_end - uy_exp_start),max(uy_exp_end - uy_exp_start)];
        
        fprintf(fid,'\n');
        fprintf(fid,'Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf(fid,'mean(ux) = %g m, std(ux) = %g m, ci(ux) = [%g %g] m\n',mean_ux,std_ux,ci_ux(1),ci_ux(2));
        fprintf(fid,'mean(uy) = %g m, std(uy) = %g m, ci(uy) = [%g %g] m\n',mean_uy,std_uy,ci_uy(1),ci_uy(2));
        fprintf(fid,'uy_exp   = %g m, error(mean(uy)) = %.3e, env(uy_exp) = [%g %g] m\n',uy_exp,err_uy,env_uy_exp);
        fprintf(fid,'mean(uz) = %g m, std(uz) = %g m, ci(uz) = [%g %g] m\n',mean_uz,std_uz,ci_uz(1),ci_uz(2));
    case 'stabilityvert'
        uz_exp_start = -1.93*1e-3;
        uz_exp_end = -[18.46 18.44 18.53 18.58 18.59 18.7 18.77 18.73 18.85 18.76]*1e-3;
        uz_exp = mean(uz_exp_end - uz_exp_start);
        err_uz = norm(mean_uz-uz_exp)/norm(uz_exp);
        env_uz_exp = [min(uz_exp_end - uz_exp_start),max(uz_exp_end - uz_exp_start)];
        
        fprintf(fid,'\n');
        fprintf(fid,'Displacement u at point (%g,%g,%g) m\n',double(P));
        fprintf(fid,'mean(ux) = %g m, std(ux) = %g m, ci(ux) = [%g %g] m\n',mean_ux,std_ux,ci_ux(1),ci_ux(2));
        fprintf(fid,'mean(uy) = %g m, std(uy) = %g m, ci(uy) = [%g %g] m\n',mean_uy,std_uy,ci_uy(1),ci_uy(2));
        fprintf(fid,'mean(uz) = %g m, std(uz) = %g m, ci(uz) = [%g %g] m\n',mean_uz,std_uz,ci_uz(1),ci_uz(2));
        fprintf(fid,'uz_exp   = %g m, error(mean(uz)) = %.3e, env(uz_exp) = [%g %g] m\n',uz_exp,err_uz,env_uz_exp);
end

fprintf(fid,'\n');
fprintf(fid,'Rotation r at point (%g,%g,%g) m\n',double(P));
fprintf(fid,'mean(rx) = %g rad = %g deg, std(rx) = %g rad = %g deg, ci(rx) = [%g %g] rad = [%g %g] deg\n',mean_rx,rad2deg(mean_rx),std_rx,rad2deg(std_rx),ci_rx(1),ci_rx(2),rad2deg(ci_rx(1)),rad2deg(ci_rx(2)));
fprintf(fid,'mean(ry) = %g rad = %g deg, std(ry) = %g rad = %g deg, ci(ry) = [%g %g] rad = [%g %g] deg\n',mean_ry,rad2deg(mean_ry),std_ry,rad2deg(std_ry),ci_ry(1),ci_ry(2),rad2deg(ci_ry(1)),rad2deg(ci_ry(2)));
fprintf(fid,'mean(rz) = %g rad = %g deg, std(rz) = %g rad = %g deg, ci(rz) = [%g %g] rad = [%g %g] deg\n',mean_rz,rad2deg(mean_rz),std_rz,rad2deg(std_rz),ci_rz(1),ci_rz(2),rad2deg(ci_rz(1)),rad2deg(ci_rz(2)));
fclose(fid);
type(filenameResults) % fprintf('%s', fileread(filenameResults))
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain(S,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 20;
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',linewidth);
    hP = plot(P,'g+');
    legend([hD,hN,hP],[legD,legN,'measure'],'Location','NorthEastOutside')
    %legend([hD,hN,hP],[legD,legN,'mesure'],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    mean_U = mean_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    ampl = getsize(S)/max(abs(mean_U))/20;
    plotModelDeflection(S,mean_u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*mean_u,'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    %% Display solution
    ampl = 0;
    % ampl = getsize(S)/max(abs(mean_U))/20;
    options = {'solid',true};
    % options = {};
    
    mean_u_mm = mean_u*1e3; % [mm]
    std_u_mm = std_u*1e3; % [mm]
    
    switch lower(test)
        case {'statichori1','statichori2','durabilityhori1','durabilityhori2'}
            plotSolution(S,mean_u_mm,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Ux',formats,renderer);
            
            plotSolution(S,std_u_mm,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Ux',formats,renderer);
        case {'statichori3','statichori4','durabilityhori3','durabilityhori4'}
            plotSolution(S,mean_u_mm,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uy',formats,renderer);
            
            plotSolution(S,std_u_mm,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uy',formats,renderer);
        case {'staticvert','stabilityvert','impact','drop'}
            plotSolution(S,mean_u_mm,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uz',formats,renderer);
            
            plotSolution(S,std_u_mm,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uz',formats,renderer);
    end
end

%% Display convergence
if displayCv
    % N = size(u,2);
    switch lower(test)
        case {'statichori1','statichori2','durabilityhori1','durabilityhori2'}
            comp = 'UX';
            mean_u_exp = ux_exp;
            lb_u_exp = env_ux_exp(1);
            ub_u_exp = env_ux_exp(2);
        case {'statichori3','statichori4','durabilityhori3','durabilityhori4'}
            comp = 'UY';
            mean_u_exp = uy_exp;
            lb_u_exp = env_uy_exp(1);
            ub_u_exp = env_uy_exp(2);
        case {'stabilityvert','staticvert','impact','drop'}
            comp = 'UZ';
            mean_u_exp = uz_exp;
            lb_u_exp = env_uz_exp(1);
            ub_u_exp = env_uz_exp(2);
    end
    % Response samples at point P
    n = 1:N;
    uP = arrayfun(@(i) eval_sol(S,u(:,i),P,comp),n);
    % Cumulative mean and standard deviation
    c1 = cumsum(uP);
    c2 = cumsum(uP.^2);
    means_u = c1 ./ n;
    stds_u = sqrt(max((c2 - c1.^2 ./ n) ./ (n-1), 0));
    % Cumulative quantiles
    lbs_u = zeros(1,N);
    ubs_u = zeros(1,N);
    for i=1:N
        qi = quantile(uP(1:i),probs);
        lbs_u(i) = qi(1);
        ubs_u(i) = qi(2);
    end
    % % Cumulative mean, standard deviation, and quantiles
    % means_u = zeros(1,N);
    % stds_u  = zeros(1,N);
    % lbs_u = zeros(1,N);
    % ubs_u = zeros(1,N);
    % for i=1:N
    %     ui = uP(1:i);
    %     means_u(i) = mean(ui);
    %     stds_u(i)  = std(ui,0);
    %     qi = quantile(ui,probs);
    %     lbs_u(i) = qi(1);
    %     ubs_u(i) = qi(2);
    % end
    % means_u = arrayfun(@(i) mean(uP(1:i)),n);
    % stds_u  = arrayfun(@(i) std(uP(1:i),0),n);
    % lbs_u   = arrayfun(@(i) quantile(uP(1:i),probs(1)),n);
    % ubs_u   = arrayfun(@(i) quantile(uP(1:i),probs(2)),n);
    
    figure('Name','Convergence solution')
    clf
    ciplot(lbs_u*1e3,ubs_u*1e3,1:N,'b');
    hold on
    ciplot([lb_u_exp,lb_u_exp]*1e3,[ub_u_exp,ub_u_exp]*1e3,[1,N],'r');
    alpha(0.2)
    plot(1:N,means_u*1e3,'-b','LineWidth',linewidth)
    plot([1,N],[mean_u_exp,mean_u_exp]*1e3,'-r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of samples','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    switch lower(test)
        case {'stabilityvert','staticvert'}
            ylabel('Vertical displacement [mm]','Interpreter',interpreter)
            %ylabel('D\''eplacement vertical [mm]','Interpreter',interpreter)
        case {'statichori1','statichori2','durabilityhori1','durabilityhori2',...
                'statichori3','statichori4','durabilityhori3','durabilityhori4'}
            ylabel('Horizontal displacement [mm]','Interpreter',interpreter)
            %ylabel('D\''eplacement horizontal [mm]','Interpreter',interpreter)
    end
    switch lower(test)
        case {'stabilityvert','staticvert','statichori1','statichori2','durabilityhori1','durabilityhori2',...
                'statichori3','statichori4','durabilityhori3','durabilityhori4'}
            legend({['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
                'experimental envelope',...
                'mean value','experimental mean value'},'Interpreter',interpreter)
            % legend({['intervalle de confiance \`a $' num2str((probs(2)-probs(1))*100) '\%$'],...
            %     'enveloppe exp\''erimentale',...
            %    'valeur moyenne','valeur moyenne exp\''erimentale'},'Interpreter',interpreter)
        case{'impact','drop'}
            
    end
    mysaveas(pathname,'convergence_solution',formats,renderer);
    mymatlab2tikz(pathname,'convergence_solution.tex');
    
    figure('Name','Convergence mean')
    clf
    plot(1:N,means_u*1e3,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Mean value [mm]','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Moyenne [mm]','Interpreter',interpreter)
    mysaveas(pathname,'convergence_mean',formats,renderer);
    mymatlab2tikz(pathname,'convergence_mean.tex');
    
    figure('Name','Convergence standard deviation')
    clf
    plot(1:N,stds_u*1e3,'-r','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of samples','Interpreter',interpreter)
    ylabel('Standard deviation [mm]','Interpreter',interpreter)
    %xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
    %ylabel('Ecart-type [mm]','Interpreter',interpreter)
    mysaveas(pathname,'convergence_std',formats,renderer);
    mymatlab2tikz(pathname,'convergence_std.tex');
end

end
end
end

myparallel('stop');

function nUpdateProgressBar(N)
persistent j
if isempty(j)
    j = 0;
end
j = j+1;
textprogressbar(j/N*100,sprintf('(%d/%d)',j,N));
if j >= N
    j = [];
end
end
