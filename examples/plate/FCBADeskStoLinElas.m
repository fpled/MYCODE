%% FCBA desk stochastic linear elasticity %%
%%----------------------------------------%%

% clc
clearvars
close all
% rng('default');
myparallel('start');

%% Input data
solveProblem = true;
displaySolution = true;
displayCv = true;

% test = 'Stability'; % stability test under vertical load
% test = 'StaticHori1'; % test under static horizontal load 1
% test = 'StaticHori2'; % test under static horizontal load 2
% test = 'StaticHori3'; % test under static horizontal load 3 (lifting)
% test = 'StaticHori4'; % test under static horizontal load 4 (lifting)
test = 'StaticVert'; % test under static vertical load
% test = 'Fatigue1'; % fatigue test under horizontal load 1 
% test = 'Fatigue2'; % fatigue test under horizontal load 2 
% test = 'Fatigue3'; % fatigue test under horizontal load 3 (lifting) 
% test = 'Fatigue4'; % fatigue test under horizontal load 4 (lifting)
% test = 'Impact'; % vertical impact test
% test = 'Drop'; % drop test

pointwiseLoading = 1; % pointwise loading

formats = {'fig','epsc2'};
renderer = 'OpenGL';

if pointwiseLoading
    filename = ['FCBADeskStoLinElasIsot' test 'PointwiseLoading'];
else
    filename = ['FCBADeskStoLinElasIsot' test];
end
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','plate',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Problem
if solveProblem
    %% Domains and meshes
    % Plates Dimensions
    a12 = 750e-3; % m
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
    x_fati = {[x3_23,y3_12+50e-3,z3],[x3_14,y3_12+50e-3,z3],...
                   [x3_23-50e-3,y3_12,z3],[x3_23-50e-3,y3_34,z3]};
    x_stab = double(getcenter(L3{1}))+[0.0,50e-3,0.0];
    x_meas = [x_hori,x_vert,x_fati,x_stab];
    P_hori = cellfun(@(x) POINT(x),x_hori,'UniformOutput',false);
    P_vert = POINT(x_vert);
    P_fati = cellfun(@(x) POINT(x),x_fati,'UniformOutput',false);
    P_stab = POINT(x_stab);
    P_meas = cellfun(@(x) POINT(x),x_meas,'UniformOutput',false);
    
    % Plates meshes
    elemtype = 'DKT';
    cl = h;
    cl_12 = cl;
    cl_3 = cl;
    cl_5 = cl;
    r_load = 40e-3;
    r_masse = 100e-3;
    C_masse = CIRCLE(0.0,y3_12+b3/2,z3,r_masse);
    x_masse = double(getcoord(getcenter(C_masse)));
    %
    L1_a = LIGNE([x5a_23,y5a,z5a_12],[x5a_23,y5a,z5a_34]);
    L1_b = LIGNE([x5b_23,y5b,z5b_12],[x5b_23,y5b,z5b_34]);
    S1 = gmshFCBAdesk12(Q1,L1_a,L1_b,cl_12,cl_5,cl_5,...
        fullfile(pathname,['gmsh_desk_1_' elemtype '_cl_' num2str(cl_12)]),3);
    S1 = convertelem(S1,elemtype);
    %
    L2_a = LIGNE([x5a_14,y5a,z5a_12],[x5a_14,y5a,z5a_34]);
    L2_b = LIGNE([x5b_14,y5b,z5b_12],[x5b_14,y5b,z5b_34]);
    S2 = gmshFCBAdesk12(Q2,L2_a,L2_b,cl_12,cl_5,cl_5,...
        fullfile(pathname,['gmsh_desk_2_' elemtype '_cl_' num2str(cl_12)]),3);
    S2 = convertelem(S2,elemtype);
    %
    L3_1 = LIGNE([x1,y1_23,z1_34],[x1,y1_14,z1_34]);
    L3_2 = LIGNE([x2,y2_23,z2_34],[x2,y2_14,z2_34]);
    if pointwiseLoading
        PbQ3 = {x_hori{4},x_fati{3},x_fati{1},x_hori{1},...
                x_fati{4},x_hori{3},x_hori{2},x_fati{2}};
        S3 = gmshFCBAdesk3simplified(Q3,C_masse,L3_1,L3_2,PbQ3,x_stab,x_masse,...
            cl_3,cl_3,cl_12,cl_12,cl_3,cl_3,cl_3,...
            fullfile(pathname,['gmsh_desk_3_' elemtype '_cl_' num2str(cl_3)]),3);
    else
        L_hori{1} = LIGNE(x_hori{1}+[0,-r_load,0],x_hori{1}+[0,r_load,0]);
        L_hori{2} = LIGNE(x_hori{2}+[0,r_load,0],x_hori{2}+[0,-r_load,0]);
        L_hori{3} = LIGNE(x_hori{3}+[r_load,0,0],x_hori{3}+[-r_load,0,0]);
        L_hori{4} = LIGNE(x_hori{4}+[-r_load,0,0],x_hori{4}+[r_load,0,0]);
        L_fati{1} = LIGNE(x_fati{1}+[0,-r_load,0],x_fati{1}+[0,r_load,0]);
        L_fati{2} = LIGNE(x_fati{2}+[0,r_load,0],x_fati{2}+[0,-r_load,0]);
        L_fati{3} = LIGNE(x_fati{3}+[-r_load,0,0],x_fati{3}+[r_load,0,0]);
        L_fati{4} = LIGNE(x_fati{4}+[r_load,0,0],x_fati{4}+[-r_load,0,0]);
        LbQ3 = {L_hori{4},L_fati{3},L_fati{1},L_hori{1},L_fati{4},L_hori{3},L_hori{2},L_fati{2}};
        C_vert = CIRCLE(x_vert(1),x_vert(2),x_vert(3),r_load);
        C_stab = CIRCLE(x_stab(1),x_stab(2),x_stab(3),r_load);
        S3 = gmshFCBAdesk3(Q3,C_masse,L3_1,L3_2,LbQ3,C_stab,C_vert,...
            cl_3,cl_3,cl_12,cl_12,cl_3,cl_3,cl_3,...
            fullfile(pathname,['gmsh_desk_3_' elemtype '_cl_' num2str(cl_3)]),3);
    end
    S3 = convertelem(S3,elemtype);
    %
    S5a = build_model(Q5a,'cl',cl_5,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_desk_5a_' elemtype '_cl_' num2str(cl_5)]));
    S5a = convertelem(S5a,elemtype);
    %
    S5b = build_model(Q5b,'cl',cl_5,'elemtype',elemtype,...
        'filename',fullfile(pathname,['gmsh_desk_5b_' elemtype '_cl_' num2str(cl_5)]));
    S5b = convertelem(S5b,elemtype);
    
    %% Random variables
    % Data
    filenameAna = 'data_ET_GL.mat';
    filenameNum = 'data_EL_NUL.mat';
    pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','identification','materialParticleBoard');
    load(fullfile(pathnameIdentification,filenameAna));
    load(fullfile(pathnameIdentification,filenameNum));
    
    % Sample number
    sample = 'B';
    
    % Material symmetry
    materialSym = 'isot';
    
    % Number of samples
    N = 1e3;
    
    switch lower(materialSym)
        case 'isot'
            % Data
            E_data = zeros(1,27);
            G_data = zeros(1,27);
            for i=1:27
                sampleNum = [sample num2str(i)];
                E_data(i) = eval(['mean_ET_' sampleNum '_data;']); % GPa
                G_data(i) = eval(['mean_GL_' sampleNum '_data;'])*13*1e-3; % GPa
            end
            NU_data = E_data./(2*G_data)-1;
            lambda_data = E_data.*NU_data./((1+NU_data).*(1-2*NU_data));
            C1_data = lambda_data + 2/3*G_data;
            C2_data = G_data;
            
            % Maximum likelihood estimation
            n_data = length(C1_data);
            m_data = length(C2_data);
            data = [C1_data C2_data];
            
            nloglf = @(lambda,data,cens,freq) -n_data*( (1-lambda(3))*log(lambda(1))...
                - gammaln(1-lambda(3)) ) - m_data*( (1-5*lambda(3))*log(lambda(2))...
                - gammaln(1-5*lambda(3)) ) + lambda(3)*( sum(log(data(1:n_data)))...
                + 5*sum(log(data(n_data+1:end))) ) + lambda(1)*sum(data(1:n_data))...
                + lambda(2)*sum(data(n_data+1:end));
            lambda = mle(data,'nloglf',nloglf,'start',[1 1 0],...
                'lowerbound',[0 0 -Inf],'upperbound',[Inf Inf 1/5]);
            a1 = 1-lambda(3);
            b1 = 1/lambda(1);
            a2 = 1-5*lambda(3);
            b2 = 1/lambda(2);
            
            % Sample set
            C1_sample = gamrnd(a1,b1,N,1)*1e9; % Pa
            C2_sample = gamrnd(a2,b2,N,1)*1e9; % Pa
            lambda_sample = C1_sample-2/3*C2_sample; % Pa
            E_sample = (9*C1_sample.*C2_sample)./(3*C1_sample+C2_sample); % Pa
            NU_sample = (3*C1_sample-2*C2_sample)./(6*C1_sample+2*C2_sample);
        case 'isottrans'
            % Data
            ET_data = zeros(1,27);
            GL_data = zeros(1,27);
            EL_data = zeros(1,27);
            NUL_data = zeros(1,27);
            for i=1:27
                sampleNum = [sample num2str(i)];
                ET_data(i) = eval(['mean_ET_' sampleNum '_data;']); % GPa
                GL_data(i) = eval(['mean_GL_' sampleNum '_data;'])*1e-3; % GPa
                EL_data(i) = eval(['mean_EL_' sampleNum '_data;'])*1e-3; % GPa
                NUL_data(i) = eval(['mean_NUL_' sampleNum '_data;']);
            end
            % NUT_data = ;
            
            % Markov Chain Monte-Carlo (MCMC) method based on 
            % Metropolis-Hasting algorithm
            
            % Sample set
            % ET_sample = ;
            % GL_sample = ;
            % EL_sample = ;
            % NUL_sample = ;
            % NUT_sample = ;
            
        otherwise
            error('Wrong material symmetry !')
    end
    
    %% Materials
    % Gravitational acceleration
    g = 9.81;
    
    % Density
    Mass_total = 13.9; % kg
    Vol_total = h*(a12*b12*2+a3*b3+a5*b5*2);
    RHO = Mass_total/(Vol_total);
    
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
            NUT = 0.25;
            % Material
            mat = ELAS_SHELL_ISOT_TRANS('ET',ET,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',h,'k',5/6);
    end
    mat = setnumber(mat,1);
    S1 = setmaterial(S1,mat);
    S2 = setmaterial(S2,mat);
    S3 = setmaterial(S3,mat);
    S5a = setmaterial(S5a,mat);
    S5b = setmaterial(S5b,mat);
    
    S = union(S1,S2,S3,S5a,S5b);
    
    %% Neumann boundary conditions
    p_plate = RHO*g*h; % surface load (body load for plates)
    Sec_stab_vert = pi*r_load^2;
    L_hori_fati = 2*r_load;
    switch lower(test)
        case 'stability'
            p = 400; % pointwise load
            if ~pointwiseLoading
                p = p/Sec_stab_vert; % surface load (body load for plates)
            end
        case {'statichori1','statichori2','statichori3','statichori4'}
            masse = 50.5;
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse; % surface load (body load for plates)
            p = 100; % pointwise load, F1=F2=100N or 200N, F3=F4=100N
            if ~pointwiseLoading
                p = p/L_hori_fati; % line load (surface load for plates)
            end
            slope = 0;
        case 'staticvert'
            p = 300; % pointwise load, 300N, 400N or 500N
            if ~pointwiseLoading
                p = p/Sec_stab_vert; % surface load (body load for plates)
            end
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            masse = 50.5;
            Sec_masse = pi*r_masse^2;
            p_masse = masse*g/Sec_masse; % surface load (body load for plates)
            p = 100; % pointwise load
            if ~pointwiseLoading
                p = p/L_hori_fati; % line load (surface load for plates)
            end
        case 'impact'
            H = 180e-3;
        case 'drop'
            H = 100e-3;
    end
    
    %% Dirichlet boundary conditions
    L1_1 = getedge(Q1,1);
    L2_1 = getedge(Q2,1);
    L5b_1 = getedge(Q5b,1);
    [~,numnode1] = intersect(S,L1_1);
    [~,numnode2] = intersect(S,L2_1);
    [~,numnode5b] = intersect(S,L5b_1);
    
    S = final(S);
    switch lower(test)
        case 'stability'
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
        case {'statichori1','statichori2'}
            S = addcl(S,numnode2);
            S = addcl(S,union(numnode1,numnode5b),{'UY','UZ'});
        case {'statichori3','statichori4'}
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
        case 'staticvert'
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
        case {'impact','drop'}
            S = addcl(S,union(numnode1,numnode2));
            S = addcl(S,numnode5b,'UZ');
    end
    
    %% Sollicitation vector
    switch lower(test)
        case 'stability'
            if pointwiseLoading
                f = nodalload(S,P_stab,'FZ',-p);
                if isempty(ispointin(P_stab,POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            else
                f = bodyload(S,C_stab,'FZ',-p);
            end
        case {'statichori1','statichori2','statichori3','statichori4'}
            if strcmpi(test,'statichori1')
                if pointwiseLoading
                    f = nodalload(S,P_hori{1},{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
                    if isempty(ispointin(P_hori{1},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{1},{'FX','FZ'},-p*[cosd(slope);sind(slope)]);
                end
            elseif strcmpi(test,'statichori2')
                if pointwiseLoading
                    f = nodalload(S,P_hori{2},{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
                    if isempty(ispointin(P_hori{2},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{2},{'FX','FZ'},p*[cosd(slope);-sind(slope)]);
                end
            elseif strcmpi(test,'statichori3')
                if pointwiseLoading
                    f = nodalload(S,P_hori{3},{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
                    if isempty(ispointin(P_hori{3},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{3},{'FY','FZ'},-p*[cosd(slope);sind(slope)]);
                end
            elseif strcmpi(test,'statichori4')
                if pointwiseLoading
                    f = nodalload(S,P_hori{4},{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
                    if isempty(ispointin(P_hori{4},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_hori{4},{'FY','FZ'},p*[cosd(slope);-sind(slope)]);
                end
            end
            if pointwiseLoading
                f = f + bodyload(keepgroupelem(S,4),[],'FZ',-p_masse);
            else
                f = f + bodyload(keepgroupelem(S,[4,5]),[],'FZ',-p_masse);
            end
        case 'staticvert'
            if pointwiseLoading
                f = nodalload(S,P_vert,'FZ',-p);
                if isempty(ispointin(P_vert,POINT(S.node)))
                    error('Pointwise load must be applied to a node of the mesh')
                end
            else
                f = bodyload(S,C_vert,'FZ',-p);
            end
        case {'fatigue1','fatigue2','fatigue3','fatigue4'}
            if strcmpi(test,'fatigue1')
                if pointwiseLoading
                    f = nodalload(S,P_fati{1},'FX',-p);
                    if isempty(ispointin(P_fati{1},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_fati{1},'FX',-p);
                end
            elseif strcmpi(test,'fatigue2')
                if pointwiseLoading
                    f = nodalload(S,P_fati{2},'FX',p);
                    if isempty(ispointin(P_fati{2},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_fati{2},'FX',p);
                end
            elseif strcmpi(test,'fatigue3')
                if pointwiseLoading
                    f = nodalload(S,P_fati{3},'FY',p);
                    if isempty(ispointin(P_fati{3},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                else
                    f = surfload(S,L_fati{3},'FY',p);
                end
            elseif strcmpi(test,'fatigue4')
                if pointwiseLoading
                    f = nodalload(S,P_fati{4},'FY',-p);
                    if isempty(ispointin(P_fati{4},POINT(S.node)))
                        error('Pointwise load must be applied to a node of the mesh')
                    end
                    
                else
                    f = surfload(S,L_fati{4},'FY',-p);
                end
            end
            if pointwiseLoading
                f = f + bodyload(keepgroupelem(S,4),[],'FZ',-p_masse);
            else
                f = f + bodyload(keepgroupelem(S,[4,5]),[],'FZ',-p_masse);
            end
        case {'impact','drop'}
            error('Not implemented')
    end
    f = f + bodyload(S,[],'FZ',-p_plate);
    
    %% Stiffness matrix and solution
    t = tic;
    u = sparse(getnbddlfree(S),N);
    parfor i=1:N
        switch lower(materialSym)
            case 'isot'
                % Young modulus
                Ei = E_sample(i);
                % Poisson ratio
                NUi = NU_sample(i);
                % Material
                mati = setparam(mat,'E',Ei);
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
                NUTi = 0.25;
                % Material
                mati = setparam(mat,'ET',NUi);
                mati = setparam(mati,'GL',Ei);
                mati = setparam(mati,'NUT',NUTi);
        end
        Si = setmaterial(S,mati);
        % Stiffness matrix
        Ai = calc_rigi(Si);
        % Solution
        u(:,i) = Ai\f;
    end
    time = toc(t);
    
    mean_u = mean(u,2);
    mean_u = unfreevector(S,mean_u);
    
    std_u = std(u,0,2);
    std_u = unfreevector(S,std_u);
    
    mean_U = mean_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    mean_Ux = mean_u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,mean_u,S.node,'UX')),:);
    mean_Uy = mean_u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,mean_u,S.node,'UY')),:);
    mean_Uz = mean_u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,mean_u,S.node,'UZ')),:);
    
    mean_R = mean_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    mean_Rx = mean_u(findddl(S,'RX'),:); % Rx = double(squeeze(eval_sol(S,mean_u,S.node,'RX')),:);
    mean_Ry = mean_u(findddl(S,'RY'),:); % Ry = double(squeeze(eval_sol(S,mean_u,S.node,'RY')),:);
    mean_Rz = mean_u(findddl(S,'RZ'),:); % Rz = double(squeeze(eval_sol(S,mean_u,S.node,'RZ')),:);
    
    std_U = std_u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    std_Ux = std_u(findddl(S,'UX'),:); % Ux = double(squeeze(eval_sol(S,std_u,S.node,'UX')),:);
    std_Uy = std_u(findddl(S,'UY'),:); % Uy = double(squeeze(eval_sol(S,std_u,S.node,'UY')),:);
    std_Uz = std_u(findddl(S,'UZ'),:); % Uz = double(squeeze(eval_sol(S,std_u,S.node,'UZ')),:);
    
    std_R = std_u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    std_Rx = std_u(findddl(S,'RX'),:); % Rx = double(squeeze(eval_sol(S,std_u,S.node,'RX')),:);
    std_Ry = std_u(findddl(S,'RY'),:); % Ry = double(squeeze(eval_sol(S,std_u,S.node,'RY')),:);
    std_Rz = std_u(findddl(S,'RZ'),:); % Rz = double(squeeze(eval_sol(S,std_u,S.node,'RZ')),:);
    
    %% Test solution
    mean_ux_P_vert = eval_sol(S,mean_u,P_vert,'UX');
    mean_uy_P_vert = eval_sol(S,mean_u,P_vert,'UY');
    mean_uz_P_vert = eval_sol(S,mean_u,P_vert,'UZ');
    mean_ux_P_stab = eval_sol(S,mean_u,P_stab,'UX');
    mean_uy_P_stab = eval_sol(S,mean_u,P_stab,'UY');
    mean_uz_P_stab = eval_sol(S,mean_u,P_stab,'UZ');
    mean_ux__P_hori(1) = eval_sol(S,mean_u,P_hori{1},'UX');
    mean_uy__P_hori(1) = eval_sol(S,mean_u,P_hori{1},'UY');
    mean_uz__P_hori(1) = eval_sol(S,mean_u,P_hori{1},'UZ');
    mean_ux_P_hori(2) = eval_sol(S,mean_u,P_hori{2},'UX');
    mean_uy_P_hori(2) = eval_sol(S,mean_u,P_hori{2},'UY');
    mean_uz_P_hori(2) = eval_sol(S,mean_u,P_hori{2},'UZ');
    mean_ux_P_hori(3) = eval_sol(S,mean_u,P_hori{3},'UX');
    mean_uy_P_hori(3) = eval_sol(S,mean_u,P_hori{3},'UY');
    mean_uz_P_hori(3) = eval_sol(S,mean_u,P_hori{3},'UZ');
    mean_ux_P_hori(4) = eval_sol(S,mean_u,P_hori{4},'UX');
    mean_uy_P_hori(4) = eval_sol(S,mean_u,P_hori{4},'UY');
    mean_uz_P_hori(4) = eval_sol(S,mean_u,P_hori{4},'UZ');    
    mean_ux_P_fati(1) = eval_sol(S,mean_u,P_fati{1},'UX');
    mean_uy_P_fati(1) = eval_sol(S,mean_u,P_fati{1},'UY');
    mean_uz_P_fati(1) = eval_sol(S,mean_u,P_fati{1},'UZ');
    mean_ux_P_fati(2) = eval_sol(S,mean_u,P_fati{2},'UX');
    mean_uy_P_fati(2) = eval_sol(S,mean_u,P_fati{2},'UY');
    mean_uz_P_fati(2) = eval_sol(S,mean_u,P_fati{2},'UZ');
    mean_ux_P_fati(3) = eval_sol(S,mean_u,P_fati{3},'UX');
    mean_uy_P_fati(3) = eval_sol(S,mean_u,P_fati{3},'UY');
    mean_uz_P_fati(3) = eval_sol(S,mean_u,P_fati{3},'UZ');
    mean_ux_P_fati(4) = eval_sol(S,mean_u,P_fati{4},'UX');
    mean_uy_P_fati(4) = eval_sol(S,mean_u,P_fati{4},'UY');
    mean_uz_P_fati(4) = eval_sol(S,mean_u,P_fati{4},'UZ');
    
    mean_rx_P_vert = eval_sol(S,mean_u,P_vert,'RX');
    mean_ry_P_vert = eval_sol(S,mean_u,P_vert,'RY');
    mean_rz_P_vert = eval_sol(S,mean_u,P_vert,'RZ');
    mean_rx_P_stab = eval_sol(S,mean_u,P_stab,'RX');
    mean_ry_P_stab = eval_sol(S,mean_u,P_stab,'RY');
    mean_rz_P_stab = eval_sol(S,mean_u,P_stab,'RZ');
    mean_rx__P_hori(1) = eval_sol(S,mean_u,P_hori{1},'RX');
    mean_ry__P_hori(1) = eval_sol(S,mean_u,P_hori{1},'RY');
    mean_rz__P_hori(1) = eval_sol(S,mean_u,P_hori{1},'RZ');
    mean_rx_P_hori(2) = eval_sol(S,mean_u,P_hori{2},'RX');
    mean_ry_P_hori(2) = eval_sol(S,mean_u,P_hori{2},'RY');
    mean_rz_P_hori(2) = eval_sol(S,mean_u,P_hori{2},'RZ');
    mean_rx_P_hori(3) = eval_sol(S,mean_u,P_hori{3},'RX');
    mean_ry_P_hori(3) = eval_sol(S,mean_u,P_hori{3},'RY');
    mean_rz_P_hori(3) = eval_sol(S,mean_u,P_hori{3},'RZ');
    mean_rx_P_hori(4) = eval_sol(S,mean_u,P_hori{4},'RX');
    mean_ry_P_hori(4) = eval_sol(S,mean_u,P_hori{4},'RY');
    mean_rz_P_hori(4) = eval_sol(S,mean_u,P_hori{4},'RZ');   
    mean_rx_P_fati(1) = eval_sol(S,mean_u,P_fati{1},'RX');
    mean_ry_P_fati(1) = eval_sol(S,mean_u,P_fati{1},'RY');
    mean_rz_P_fati(1) = eval_sol(S,mean_u,P_fati{1},'RZ');
    mean_rx_P_fati(2) = eval_sol(S,mean_u,P_fati{2},'RX');
    mean_ry_P_fati(2) = eval_sol(S,mean_u,P_fati{2},'RY');
    mean_rz_P_fati(2) = eval_sol(S,mean_u,P_fati{2},'RZ');
    mean_rx_P_fati(3) = eval_sol(S,mean_u,P_fati{3},'RX');
    mean_ry_P_fati(3) = eval_sol(S,mean_u,P_fati{3},'RY');
    mean_rz_P_fati(3) = eval_sol(S,mean_u,P_fati{3},'RZ');
    mean_rx_P_fati(4) = eval_sol(S,mean_u,P_fati{4},'RX');
    mean_ry_P_fati(4) = eval_sol(S,mean_u,P_fati{4},'RY');
    mean_rz_P_fati(4) = eval_sol(S,mean_u,P_fati{4},'RZ');
    
    std_ux_P_vert = eval_sol(S,std_u,P_vert,'UX');
    std_uy_P_vert = eval_sol(S,std_u,P_vert,'UY');
    std_uz_P_vert = eval_sol(S,std_u,P_vert,'UZ');
    std_ux_P_stab = eval_sol(S,std_u,P_stab,'UX');
    std_uy_P_stab = eval_sol(S,std_u,P_stab,'UY');
    std_uz_P_stab = eval_sol(S,std_u,P_stab,'UZ');
    std_ux__P_hori(1) = eval_sol(S,std_u,P_hori{1},'UX');
    std_uy__P_hori(1) = eval_sol(S,std_u,P_hori{1},'UY');
    std_uz__P_hori(1) = eval_sol(S,std_u,P_hori{1},'UZ');
    std_ux_P_hori(2) = eval_sol(S,std_u,P_hori{2},'UX');
    std_uy_P_hori(2) = eval_sol(S,std_u,P_hori{2},'UY');
    std_uz_P_hori(2) = eval_sol(S,std_u,P_hori{2},'UZ');
    std_ux_P_hori(3) = eval_sol(S,std_u,P_hori{3},'UX');
    std_uy_P_hori(3) = eval_sol(S,std_u,P_hori{3},'UY');
    std_uz_P_hori(3) = eval_sol(S,std_u,P_hori{3},'UZ');
    std_ux_P_hori(4) = eval_sol(S,std_u,P_hori{4},'UX');
    std_uy_P_hori(4) = eval_sol(S,std_u,P_hori{4},'UY');
    std_uz_P_hori(4) = eval_sol(S,std_u,P_hori{4},'UZ');    
    std_ux_P_fati(1) = eval_sol(S,std_u,P_fati{1},'UX');
    std_uy_P_fati(1) = eval_sol(S,std_u,P_fati{1},'UY');
    std_uz_P_fati(1) = eval_sol(S,std_u,P_fati{1},'UZ');
    std_ux_P_fati(2) = eval_sol(S,std_u,P_fati{2},'UX');
    std_uy_P_fati(2) = eval_sol(S,std_u,P_fati{2},'UY');
    std_uz_P_fati(2) = eval_sol(S,std_u,P_fati{2},'UZ');
    std_ux_P_fati(3) = eval_sol(S,std_u,P_fati{3},'UX');
    std_uy_P_fati(3) = eval_sol(S,std_u,P_fati{3},'UY');
    std_uz_P_fati(3) = eval_sol(S,std_u,P_fati{3},'UZ');
    std_ux_P_fati(4) = eval_sol(S,std_u,P_fati{4},'UX');
    std_uy_P_fati(4) = eval_sol(S,std_u,P_fati{4},'UY');
    std_uz_P_fati(4) = eval_sol(S,std_u,P_fati{4},'UZ');
    
    std_rx_P_vert = eval_sol(S,std_u,P_vert,'RX');
    std_ry_P_vert = eval_sol(S,std_u,P_vert,'RY');
    std_rz_P_vert = eval_sol(S,std_u,P_vert,'RZ');
    std_rx_P_stab = eval_sol(S,std_u,P_stab,'RX');
    std_ry_P_stab = eval_sol(S,std_u,P_stab,'RY');
    std_rz_P_stab = eval_sol(S,std_u,P_stab,'RZ');
    std_rx__P_hori(1) = eval_sol(S,std_u,P_hori{1},'RX');
    std_ry__P_hori(1) = eval_sol(S,std_u,P_hori{1},'RY');
    std_rz__P_hori(1) = eval_sol(S,std_u,P_hori{1},'RZ');
    std_rx_P_hori(2) = eval_sol(S,std_u,P_hori{2},'RX');
    std_ry_P_hori(2) = eval_sol(S,std_u,P_hori{2},'RY');
    std_rz_P_hori(2) = eval_sol(S,std_u,P_hori{2},'RZ');
    std_rx_P_hori(3) = eval_sol(S,std_u,P_hori{3},'RX');
    std_ry_P_hori(3) = eval_sol(S,std_u,P_hori{3},'RY');
    std_rz_P_hori(3) = eval_sol(S,std_u,P_hori{3},'RZ');
    std_rx_P_hori(4) = eval_sol(S,std_u,P_hori{4},'RX');
    std_ry_P_hori(4) = eval_sol(S,std_u,P_hori{4},'RY');
    std_rz_P_hori(4) = eval_sol(S,std_u,P_hori{4},'RZ');    
    std_rx_P_fati(1) = eval_sol(S,std_u,P_fati{1},'RX');
    std_ry_P_fati(1) = eval_sol(S,std_u,P_fati{1},'RY');
    std_rz_P_fati(1) = eval_sol(S,std_u,P_fati{1},'RZ');
    std_rx_P_fati(2) = eval_sol(S,std_u,P_fati{2},'RX');
    std_ry_P_fati(2) = eval_sol(S,std_u,P_fati{2},'RY');
    std_rz_P_fati(2) = eval_sol(S,std_u,P_fati{2},'RZ');
    std_rx_P_fati(3) = eval_sol(S,std_u,P_fati{3},'RX');
    std_ry_P_fati(3) = eval_sol(S,std_u,P_fati{3},'RY');
    std_rz_P_fati(3) = eval_sol(S,std_u,P_fati{3},'RZ');
    std_rx_P_fati(4) = eval_sol(S,std_u,P_fati{4},'RX');
    std_ry_P_fati(4) = eval_sol(S,std_u,P_fati{4},'RY');
    std_rz_P_fati(4) = eval_sol(S,std_u,P_fati{4},'RZ');
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S','S1','S2','S3','S5a','S5b',...
        'elemtype','a12','b12','a3','b3','a5','b5','h','Sec_stab_vert','L_hori_fati',...
        'f','p','pointwiseLoading');
    save(fullfile(pathname,'solution.mat'),'u','mean_u','std_u','time',...
        'mean_U','mean_Ux','mean_Uy','mean_Uz',...
        'mean_R','mean_Rx','mean_Ry','mean_Rz',...
        'std_U','std_Ux','std_Uy','std_Uz',...
        'std_R','std_Rx','std_Ry','std_Rz');
    save(fullfile(pathname,'test_solution.mat'),'P_vert','P_stab',...
        'P_hori','P_fati',...
        'mean_ux_P_vert','mean_uy_P_vert','mean_uz_P_vert',...
        'mean_ux_P_stab','mean_uy_P_stab','mean_uz_P_stab',...
        'mean_ux_P_hori','mean_uy_P_hori','mean_uz_P_hori',...        
        'mean_ux_P_fati','mean_uy_P_fati','mean_uz_P_fati',...      
        'mean_rx_P_vert','mean_ry_P_vert','mean_rz_P_vert',...
        'mean_rx_P_stab','mean_ry_P_stab','mean_rz_P_stab',...
        'mean_rx_P_hori','mean_ry_P_hori','mean_rz_P_hori',...
        'mean_rx_P_fati','mean_ry_P_fati','mean_rz_P_fati',...
        'std_ux_P_vert','std_uy_P_vert','std_uz_P_vert',...
        'std_ux_P_stab','std_uy_P_stab','std_uz_P_stab',...
        'std_ux_P_hori','std_uy_P_hori','std_uz_P_hori',...
        'std_ux_P_fati','std_uy_P_fati','std_uz_P_fati',...
        'std_rx_P_vert','std_ry_P_vert','std_rz_P_vert',...
        'std_rx_P_stab','std_ry_P_stab','std_rz_P_stab',...
        'std_rx_P_hori','std_ry_P_hori','std_rz_P_hori',...
        'std_rx_P_fati','std_ry_P_fati','std_rz_P_fati');
else
    load(fullfile(pathname,'problem.mat'),'S','S1','S2','S3','S5a','S5b',...
        'elemtype','a12','b12','a3','b3','a5','b5','h','Sec_stab_vert','L_hori_fati',...
        'f','p','pointwiseLoading');
    load(fullfile(pathname,'solution.mat'),'u','mean_u','std_u','time',...
        'mean_U','mean_Ux','mean_Uy','mean_Uz',...
        'mean_R','mean_Rx','mean_Ry','mean_Rz',...
        'std_U','std_Ux','std_Uy','std_Uz',...
        'std_R','std_Rx','std_Ry','std_Rz');
    load(fullfile(pathname,'test_solution.mat'),'P_vert','P_stab',...
        'P_hori','P_fati',...
        'mean_ux_P_vert','mean_uy_P_vert','mean_uz_P_vert',...
        'mean_ux_P_stab','mean_uy_P_stab','mean_uz_P_stab',...
        'mean_ux_P_hori','mean_uy_P_hori','mean_uz_P_hori',...        
        'mean_ux_P_fati','mean_uy_P_fati','mean_uz_P_fati',...      
        'mean_rx_P_vert','mean_ry_P_vert','mean_rz_P_vert',...
        'mean_rx_P_stab','mean_ry_P_stab','mean_rz_P_stab',...
        'mean_rx_P_hori','mean_ry_P_hori','mean_rz_P_hori',...
        'mean_rx_P_fati','mean_ry_P_fati','mean_rz_P_fati',...
        'std_ux_P_vert','std_uy_P_vert','std_uz_P_vert',...
        'std_ux_P_stab','std_uy_P_stab','std_uz_P_stab',...
        'std_ux_P_hori','std_uy_P_hori','std_uz_P_hori',...
        'std_ux_P_fati','std_uy_P_fati','std_uz_P_fati',...
        'std_rx_P_vert','std_ry_P_vert','std_rz_P_vert',...
        'std_rx_P_stab','std_ry_P_stab','std_rz_P_stab',...
        'std_rx_P_hori','std_ry_P_hori','std_rz_P_hori',...
        'std_rx_P_fati','std_ry_P_fati','std_rz_P_fati');
end

%% Outputs
fprintf('\nDesk\n');
fprintf(['test : ' test '\n']);
fprintf(['mesh : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('span-to-thickness ratio of plates 1 and 2 = %g\n',min(a12,b12)/h);
fprintf('span-to-thickness ratio of plate 3 = %g\n',min(a3,b3)/h);
fprintf('span-to-thickness ratio of plates 5a and 5b = %g\n',min(a5,b5)/h);
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

switch lower(test)
    case 'staticvert'
        disp('Displacement u at point'); disp(P_vert);
        fprintf('mean(ux) = %g, std(ux) = %g\n',mean_ux_P_vert,std_ux_P_vert);
        fprintf('mean(uy) = %g, std(uy) = %g\n',mean_uy_P_vert,std_uy_P_vert);
        fprintf('mean(uz) = %g, std(uz) = %g\n',mean_uz_P_vert,std_uz_P_vert);
        if (pointwiseLoading && p==300) || (~pointwiseLoading && p==300/Sec_stab_vert)
            uz_exp_start = -0.69*1e-3;
            uz_exp_end = -[10.10 9.88 9.64 9.88 9.94 9.79 9.92 9.93 9.82 9.95]*1e-3;
        elseif (pointwiseLoading && p==400) || (~pointwiseLoading && p==400/Sec_stab_vert)
            uz_exp_start = -0.75*1e-3;
            uz_exp_end = -[13.45 13.52 13.56 13.64 13.65 13.74 13.75 13.44 13.74 13.53]*1e-3;
        elseif (pointwiseLoading && p==500) || (~pointwiseLoading && p==500/Sec_stab_vert)
            uz_exp_start = -0.78*1e-3;
            uz_exp_end = -[16.66 16.57 16.59 16.78 16.55 16.69 16.75 16.59 16.73 16.76]*1e-3;
        end
        uz_exp = mean(uz_exp_end - uz_exp_start);
        err_uz = norm(mean_uz_P_vert-uz_exp)/norm(uz_exp);
        fprintf('uz_exp   = %g, error    = %.3e\n',uz_exp,err_uz);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_vert);
        fprintf('mean(rx) = %g, std(rx) = %g\n',mean_rx_P_vert,std_rx_P_vert);
        fprintf('mean(ry) = %g, std(ry) = %g\n',mean_ry_P_vert,std_ry_P_vert);
        fprintf('mean(rz) = %g, std(rz) = %g\n',mean_rz_P_vert,std_rz_P_vert);
        fprintf('\n');
    case 'stability'
        disp('Displacement u at point'); disp(P_stab);
        fprintf('mean(ux) = %g, std(ux) = %g\n',mean_ux_P_stab,std_ux_P_stab);
        fprintf('mean(uy) = %g, std(uy) = %g\n',mean_uy_P_stab,std_uy_P_stab);
        fprintf('mean(uz) = %g, std(uz) = %g\n',mean_uz_P_stab,std_uz_P_stab);
        uz_exp_start = -1.93*1e-3;
        uz_exp_end = -[18.46 18.44 18.53 18.58 18.59 18.7 18.77 18.73 18.85 18.76]*1e-3;
        uz_exp = mean(uz_exp_end - uz_exp_start);
        err_uz = norm(mean_uz_P_stab-uz_exp)/norm(uz_exp);
        fprintf('uz_exp   = %g, error    = %.3e\n',uz_exp,err_uz);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_stab);
        fprintf('mean(rx) = %g, std(rx) = %g\n',mean_rx_P_stab,std_rx_P_stab);
        fprintf('mean(ry) = %g, std(ry) = %g\n',mean_ry_P_stab,std_ry_P_stab);
        fprintf('mean(rz) = %g, std(rz) = %g\n',mean_rz_P_stab,std_rz_P_stab);
        fprintf('\n');
    case 'statichori1'
        disp('Displacement u at point'); disp(P_hori{2});
        fprintf('mean(ux) = %g, std(ux) = %g\n',mean_ux_P_hori(2),std_ux_P_hori(2));
        fprintf('mean(uy) = %g, std(uy) = %g\n',mean_uy_P_hori(2),std_uy_P_hori(2));
        fprintf('mean(uz) = %g, std(uz) = %g\n',mean_uz_P_hori(2),std_uz_P_hori(2));
        if (pointwiseLoading && p==100) || (~pointwiseLoading && p==100/L_hori_fati)
            ux_exp_start = -6.88*1e-3;
            ux_exp_end = -[10.5 10.51 10.44 10.8 10.72 10.62 10.67 10.65 10.66 10.87 10.86]*1e-3;
        elseif (pointwiseLoading && p==200) || (~pointwiseLoading && p==200/L_hori_fati)
            ux_exp_start = -6.16*1e-3;
            ux_exp_end = -[16.78 16.74 16.72 17.13 17 16.8 16.87 16.78 17.04 16.82 16.71 17.17]*1e-3;
        end
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(mean_ux_P_hori(2)-ux_exp)/norm(ux_exp);
        fprintf('ux_exp   = %g, error    = %.3e\n',ux_exp,err_ux);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_hori{2});
        fprintf('mean(rx) = %g, std(rx) = %g\n',mean_rx_P_hori(2),std_rx_P_hori(2));
        fprintf('mean(ry) = %g, std(ry) = %g\n',mean_ry_P_hori(2),std_ry_P_hori(2));
        fprintf('mean(rz) = %g, std(rz) = %g\n',mean_rz_P_hori(2),std_rz_P_hori(2));
        fprintf('\n');
    case 'statichori2'
        disp('Displacement u at point'); disp(P_hori{1});
        fprintf('mean(ux) = %g, std(ux) = %g\n',mean_ux__P_hori(1),std_ux__P_hori(1));
        fprintf('mean(uy) = %g, std(uy) = %g\n',mean_uy__P_hori(1),std_uy__P_hori(1));
        fprintf('mean(uz) = %g, std(uz) = %g\n',mean_uz__P_hori(1),std_uz__P_hori(1));
        if (pointwiseLoading && p==100) || (~pointwiseLoading && p==100/L_hori_fati)
            ux_exp_start = 2.12*1e-3;
            ux_exp_end = [6.22 6.17 6.26 6.31 6.33 6.24 6.26 6.4 6.26 6.49 6.48 6.42 6.36 6.56 6.37 6.39]*1e-3;
        elseif (pointwiseLoading && p==200) || (~pointwiseLoading && p==200/L_hori_fati)
            ux_exp_start = 1.91*1e-3;
            ux_exp_end = [12.45 12.68 12.66 12.65 12.71 12.64 12.82 12.73 12.89 12.86 12.79 12.86]*1e-3;
        end
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(mean_ux__P_hori(1)-ux_exp)/norm(ux_exp);
        fprintf('ux_exp   = %g, error    = %.3e\n',ux_exp,err_ux);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_hori{1});
        fprintf('mean(rx) = %g, std(rx) = %g\n',mean_rx__P_hori(1),std_rx__P_hori(1));
        fprintf('mean(ry) = %g, std(ry) = %g\n',mean_ry__P_hori(1),std_ry__P_hori(1));
        fprintf('mean(rz) = %g, std(rz) = %g\n',mean_rz__P_hori(1),std_rz__P_hori(1));
        fprintf('\n'); 
    case 'fatigue1'
        disp('Displacement u at point'); disp(P_fati{2});
        fprintf('mean(ux) = %g, std(ux) = %g\n',mean_ux_P_fati(2),std_ux_P_fati(2));
        fprintf('mean(uy) = %g, std(uy) = %g\n',mean_uy_P_fati(2),std_uy_P_fati(2));
        fprintf('mean(uz) = %g, std(uz) = %g\n',mean_uz_P_fati(2),std_uz_P_fati(2));
        ux_exp_start = -4.42*1e-3;
        ux_exp_end = -[8.4 8.3 8.37 8.41 8.54 8.39 8.56 8.48 8.46 8.49 8.49 8.43 8.55 8.52]*1e-3;   
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(mean_ux_P_fati(2)-ux_exp)/norm(ux_exp);
        fprintf('ux_exp   = %g, error    = %.3e\n',ux_exp,err_ux);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_fati{2});
        fprintf('mean(rx) = %g, std(rx) = %g\n',mean_rx_P_fati(2),std_rx_P_fati(2));
        fprintf('mean(ry) = %g, std(ry) = %g\n',mean_ry_P_fati(2),std_ry_P_fati(2));
        fprintf('mean(rz) = %g, std(rz) = %g\n',mean_rz_P_fati(2),std_rz_P_fati(2));
        fprintf('\n');
    case 'fatigue2'
        disp('Displacement u at point'); disp(P_fati{1});
        fprintf('mean(ux) = %g, std(ux) = %g\n',mean_ux_P_fati(1),std_ux_P_fati(1));
        fprintf('mean(uy) = %g, std(uy) = %g\n',mean_uy_P_fati(1),std_uy_P_fati(1));
        fprintf('mean(uz) = %g, std(uz) = %g\n',mean_uz_P_fati(1),std_uz_P_fati(1));
        ux_exp_start = 3.48*1e-3;
        ux_exp_end = [7.89 7.85 8.1 8.4 8.36 8.55 8.27 8.27 8.47 8.49 8.64 8.35 8.5 8.63 8.73]*1e-3;   
        ux_exp = mean(ux_exp_end - ux_exp_start);
        err_ux = norm(mean_ux_P_fati(1)-ux_exp)/norm(ux_exp);
        fprintf('ux_exp   = %g, error    = %.3e\n',ux_exp,err_ux);
        fprintf('\n');
        disp('Rotation r at point'); disp(P_fati{1});
        fprintf('mean(rx) = %g, std(rx) = %g\n',mean_rx_P_fati(1),std_rx_P_fati(1));
        fprintf('mean(ry) = %g, std(ry) = %g\n',mean_ry_P_fati(1),std_ry_P_fati(1));
        fprintf('mean(rz) = %g, std(rz) = %g\n',mean_rz_P_fati(1),std_rz_P_fati(1));
        fprintf('\n');   
end

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain(S,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 8;
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
    % legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    ampl = getsize(S)/max(abs(mean_u))/10;
    plotModelDeflection(S,mean_u,'ampl',ampl,'Color','b','FaceColor','b',...
        'FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*mean_u,'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    %% Display solution
    % ampl = 0;
    ampl = getsize(S)/max(abs(mean_u))/10;
    options = {'solid',true};
    % options = {};
    
    switch lower(test)
        case 'stability'
            plotSolution(S,mean_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uz',formats,renderer);
            
            plotSolution(S,std_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uz',formats,renderer);
        case {'statichori1','statichori2'}
            plotSolution(S,mean_u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Ux',formats,renderer);
            
            plotSolution(S,std_u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Ux',formats,renderer);
        case {'statichori3','statichori4'}
            plotSolution(S,mean_u,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uy',formats,renderer);
            
            plotSolution(S,std_u,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uy',formats,renderer);
        case 'staticvert'
            plotSolution(S,mean_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uz',formats,renderer);
            
            plotSolution(S,std_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uz',formats,renderer);
        case {'fatigue1','fatigue2'}
            plotSolution(S,mean_u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Ux',formats,renderer);
            
            plotSolution(S,std_u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Ux',formats,renderer);
        case {'fatigue3','fatigue4'}
            plotSolution(S,mean_u,'displ',2,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uy',formats,renderer);
            
            plotSolution(S,std_u,'displ',1,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uy',formats,renderer);
        case 'impact'
            plotSolution(S,mean_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uz',formats,renderer);
            
            plotSolution(S,std_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uz',formats,renderer);
        case 'drop'
            plotSolution(S,mean_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'mean_Uz',formats,renderer);
            
            plotSolution(S,std_u,'displ',3,'ampl',ampl,options{:});
            mysaveas(pathname,'std_Uz',formats,renderer);
    end
    
    % plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'Rx',formats,renderer);
    %
    % plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'Ry',formats,renderer);
end

%% Display convergence Monte-Carlo
if displayCv
    N = size(u,2);
    means_u = arrayfun(@(x) norm(mean(u(:,1:x),2)),1:N);
    stds_u = arrayfun(@(x) norm(std(u(:,1:x),0,2)),1:N);
    
    figure('Name','Convergence empirical mean')
    clf
    plot(1:N,means_u,'-b','LineWidth',1)
    grid on
    box on
    set(gca,'FontSize',16)
    % xlabel('Nombre de r\''ealisations','Interpreter','latex')
    % ylabel('Moyenne empirique','Interpreter','latex')
    xlabel('Number of samples','Interpreter','latex')
    ylabel('Empirical mean','Interpreter','latex')
    mysaveas(pathname,'convergence_empirical_mean','fig');
    mymatlab2tikz(pathname,'convergence_empirical_mean.tex');
    
    figure('Name','Convergence empirical standard deviation')
    clf
    plot(1:N,stds_u,'-r','LineWidth',1)
    grid on
    box on
    set(gca,'FontSize',16)
    % xlabel('Nombre de r\''ealisations','Interpreter','latex')
    % ylabel('Ecart-type empirique','Interpreter','latex')
    xlabel('Number of samples','Interpreter','latex')
    ylabel('Empirical standard deviation','Interpreter','latex')
    mysaveas(pathname,'convergence_empirical_std','fig');
    mymatlab2tikz(pathname,'convergence_empirical_std.tex');
end

myparallel('stop');
