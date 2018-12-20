%% FCBA desk beam deterministic linear elasticity %%
%%------------------------------------------------%%

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

tests = {'StaticVert'}; % test under static vertical load

filename = 'FCBADeskBeamJunctionDetLinElas';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','plate',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

junction = true; % junction modeling

%% Problem
if solveProblem
    %% Domains and meshes
    % Beams dimensions
    a1 = 750e-3; % m
    b1 = 396e-3;
    a2 = 940e-3;
    b2 = 501e-3;
    b3 = 113e-3;
    % Same thickness for all the beams
    h = 15e-3;
    
    L1 = a1+h/2;
    L2 = (a2+h)/2;
    
    % Points
    P1 = POINT([0.0,0.0]);
    P2 = POINT([0.0,L1]);
    P_load = POINT([L2,L1]);
    P3 = POINT([2*L2,L1]);
    P4 = POINT([2*L2,0.0]);
    
    L_beam{1} = LIGNE(P1,P2);
    L_beam{2} = LIGNE(P2,P_load);
    L_beam{3} = LIGNE(P_load,P3);
    L_beam{4} = LIGNE(P4,P3);
    
    cl_beam = h/20;
    S_beam = cellfun(@(L,n) build_model(L,'cl',cl_beam,'elemtype','BEAM','filename',fullfile(pathname,['gmsh_beam_' num2str(n) '_cl_' num2str(cl_beam)])),L_beam,num2cell(1:length(L_beam)),'UniformOutput',false);
    
    %% Materials
    % Gravitational acceleration
    g = 9.81;
    
    % Density
    RHO = 707.1384; % kg/m3
    Vol_total = 2*(L1*b1+L2*(b2+2*b3))*h;
    Mass_total = Vol_total*RHO; % kg
    
    % Cross-section area
    Sec1 = b1*h;
    Sec2 = b2*h;
    % Planar second moment of area (or Planar area moment of inertia)
    IY1 = h*b1^3/12;
    IY2 = h*b2^3/12;
    IZ1 = b1*h^3/12;
    IZ2 = b2*h^3/12;
    IX1 = IY1+IZ1;
    IX2 = IY2+IZ2;
    
    % Data
    filenameAna = 'data_ET_GL.mat';
    filenameNum = 'data_EL_NUL.mat';
    
    pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','identification','materialParticleBoard');
    load(fullfile(pathnameIdentification,filenameAna));
    load(fullfile(pathnameIdentification,filenameNum));
    
    % Material symmetry
    materialSym = 'isot';
    
    switch lower(materialSym)
        case 'isot'
            % Young modulus
            E = mean(mean_ET_data)*1e6; % Pa
            %E = 2e9; % Pa
            % Shear modulus
            %G = mean(mean_GL_data)*1e6*13; % Pa
            % Poisson ratio
            %NU = E./(2*G)-1;
            NU = 0.25;
            % Material
            mat_1 = ELAS_BEAM('E',E,'NU',NU,'S',Sec1,'IZ',IZ1,'IY',IY1,'IX',IX1,'RHO',RHO);
            mat_1 = setnumber(mat_1,1);
            mat_2 = ELAS_BEAM('E',E,'NU',NU,'S',Sec2,'IZ',IZ2,'IY',IY2,'IX',IX2,'RHO',RHO);
            mat_2 = setnumber(mat_2,2);
            S_beam([1,4]) = cellfun(@(S) setmaterial(S,mat_1),S_beam([1,4]),'UniformOutput',false);
            S_beam([2,3]) = cellfun(@(S) setmaterial(S,mat_2),S_beam([2,3]),'UniformOutput',false);
        case 'isottrans'
            % Transverse Young modulus
            ET = mean(mean_ET_data)*1e6; % Pa
            % Longitudinal shear modulus
            GL = mean(mean_GL_data)*1e6; % Pa
            % Longitudinal Young modulus
            % EL = mean(mean_EL_data)*1e6; % Pa
            % Longitudinal Poisson ratio
            % NUL = mean(mean_NUL_data);
            % Transverse Poisson ratio
            NUT = 0.25;
            % Material
            mat_1 = ELAS_BEAM_ISOT_TRANS('ET',ET,'NUT',NUT,'GL',GL,'S',Sec1,'IZ',IZ1,'IY',IY1,'IX',IX1,'RHO',RHO);
            mat_1 = setnumber(mat_1,1);
            mat_2 = ELAS_BEAM_ISOT_TRANS('ET',ET,'NUT',NUT,'GL',GL,'S',Sec2,'IZ',IZ2,'IY',IY2,'IX',IX2,'RHO',RHO);
            mat_2 = setnumber(mat_2,2);
            S_beam([1,4]) = cellfun(@(S) setmaterial(S,mat_1),S_beam([1,4]),'UniformOutput',false);
            S_beam([2,3]) = cellfun(@(S) setmaterial(S,mat_2),S_beam([2,3]),'UniformOutput',false);
        otherwise
            error('Wrong material symmetry !')
    end
    
    if junction
        S23 = union(S_beam{2:3});
        S = union(S_beam{1},S23,S_beam{4},'duplicate');
    else
        S = union(S_beam{:});
    end
    
    %% Neumann boundary conditions
    p1 = RHO*g*Sec1; % line load (body load for beams)
    p2 = RHO*g*Sec2; % line load (body load for beams)
    p = 300; % pointwise load, 300N, 400N or 500N
    
    %% Dirichlet boundary conditions
    if junction
        S = final(S,'duplicate');
        % [~,numnode2,~] = intersect(S,P2,'strict',false);
        % [~,numnode3,~] = intersect(S,P3,'strict',false);
        numnode2 = find(S.node==P2);
        numnode3 = find(S.node==P3);
        S = addclperiodic(S,numnode2(1),numnode2(2),'U');
        S = addclperiodic(S,numnode3(1),numnode3(2),'U');
    else
        S = final(S);
    end
    P_support = [P1 P4];
    S = addcl(S,P_support);
    
    %% Stiffness matrix and sollicitation vector
    A = calc_rigi(S);
    f = nodalload(S,P_load,'FY',-p);
    f = f + bodyload(keepgroupelem(S,[1,4]),[],'FY',-p1);
    f = f + bodyload(keepgroupelem(S,[2,3]),[],'FY',-p2);
    
    if junction
        k = 1e2; % additonal junction rotational stiffness
        numddl2 = findddl(S,'RZ',numnode2,'free');
        numddl3 = findddl(S,'RZ',numnode3,'free');
        A_add = [k -k;-k k];
        A(numddl2,numddl2) = A(numddl2,numddl2) + A_add;
        A(numddl3,numddl3) = A(numddl3,numddl3) + A_add;
    end
    
    %% Solution
    t = tic;
    u = A\f;
    time = toc(t);
    
    u = unfreevector(S,u);
    
    % U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    Ux = u(findddl(S,'UX'),:);
    Uy = u(findddl(S,'UY'),:);
    % R = u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    Rz = u(findddl(S,'RZ'),:);
    
    x = getcoord(S.node);
    
    e = calc_epsilon(S,u,'smooth');
    s = calc_sigma(S,u,'smooth');
    
    Epsx = e(1);
    Gamz = e(2);
    N = s(1);
    Mz = s(2);
    
    %% Reference solution
    if junction
        XA = (3*p+4*p2*L2)/2*L2^2*IZ1/(L1*(L1*IZ2+4*L2*IZ1*(1+E*IZ2/(k*L2))) + 12*L2/L1*IZ1/Sec2*(IZ2*(1+E*IZ1/(k*L1))+L2/L1*IZ1));
        % MA = (p+4/3*p2*L2)/2*L2^2*IZ1*(-L1+6*L2/L1^2*IZ1/Sec2)/(L1*(L1*IZ2+4*L2*IZ1) + 12*L2/L1*IZ1/Sec2*(IZ2*(1+E*IZ1/(k*L1))+L2/L1*IZ1));
    else
        XA = (3*p+4*p2*L2)/2*L2^2*IZ1/(L1*(L1*IZ2+4*L2*IZ1) + 12*L2/L1*IZ1/Sec2*(IZ2+L2/L1*IZ1));
        % MA = (p+4/3*p2*L2)/2*L2^2*IZ1*(-L1+6*L2/L1^2*IZ1/Sec2)/(L1*(L1*IZ2+4*L2*IZ1) + 12*L2/L1*IZ1/Sec2*(IZ2+L2/L1*IZ1));
    end
    YA = p1*L1+p2*L2+p/2;
    MA = XA*(-L1/3+2*L2/L1^2*IZ1/Sec2);
    
    fun_Ux = cell(2,1);
    fun_Uy = cell(2,1);
    fun_Rz = cell(2,1);
    fun_N  = cell(2,1);
    fun_Mz = cell(2,1);
    fun_Epsx = cell(2,1);
    fun_Gamz = cell(2,1);
    
    % Horizontal beam in local coordinates system
    fun_Ux{1} = @(x) 1/(E*Sec1)*(p1*x/2-YA).*x;
    fun_Uy{1} = @(x) -1/(E*IZ1)*(XA*x.^3/6+MA*x.^2/2);
    fun_Rz{1} = @(x) -1/(E*IZ1)*(XA*x.^2/2+MA*x);
    fun_N{1}  = @(x) p1*x-YA;
    fun_Ty{1} = @(x) XA*ones(size(x));
    fun_Mz{1} = @(x) -XA*x-MA;
    fun_Epsx{1} = @(x) fun_N{1}(x)/(E*Sec1);
    fun_Gamz{1} = @(x) fun_Mz{1}(x)/(E*IZ1);
    
    % Vertical beam in local coordinates system
    fun_Ux{2} = @(x) XA/(E*Sec2)*(L2-x);
    fun_Uy{2} = @(x) 1/(E*IZ2)*(-p2*x.^4/24 + (p2*L2+p/2)*x.^3/6 - (XA*L1+MA)*x.^2/2) - L1/(E*IZ1)*(XA*L1/2+MA)*x - L1/(E*Sec1)*(p1*L1/2+p2*L2+p/2);
    if junction
        fun_Uy{2} = @(x) fun_Uy{2}(x) - 1/k*(XA*L1+MA)*x;
    end
    fun_Rz{2} = @(x) 1/(E*IZ2)*(p2/6*(3*L2^2-(x-L2).^2) + p/4*(x+L2) - (XA*L1+MA)).*(x-L2);
    % fun_Rz{2} = @(x) 1/(E*IZ2)*(-p2*x.^3/6 + (p2*L2+p/2)*x.^2/2 - (XA*L1+MA)*x) - L1/(E*IZ1)*(XA*L1/2+MA);
    % if junction
    %     fun_Rz{2} = @(x) fun_Rz{2}(x) - 1/k*(XA*L1+MA);
    % end
    fun_N{2}  = @(x) -XA*ones(size(x));
    fun_Ty{2} = @(x) p2*(x-L2)-p/2;
    % fun_Ty{2} = @(x) p2*x+p1*L1-YA;
    fun_Mz{2} = @(x) -p2*x.^2/2+(p2*L2+p/2)*x-XA*L1-MA;
    fun_Epsx{2} = @(x) fun_N{2}(x)/(E*Sec2);
    fun_Gamz{2} = @(x) fun_Mz{2}(x)/(E*IZ2);
    
    for i=1:2
        fun_Ux{i} = MultiVariateFunction(fun_Ux{i},1);
        fun_Uy{i} = MultiVariateFunction(fun_Uy{i},1);
        fun_Rz{i} = MultiVariateFunction(fun_Rz{i},1);
        fun_N{i}  = MultiVariateFunction(fun_N{i},1);
        fun_Ty{i} = MultiVariateFunction(fun_Ty{i},1);
        fun_Mz{i} = MultiVariateFunction(fun_Mz{i},1);
        fun_Epsx{i} = MultiVariateFunction(fun_Epsx{i},1);
        fun_Gamz{i} = MultiVariateFunction(fun_Gamz{i},1);
        fun_Ux{i}.evaluationAtMultiplePoints = true;
        fun_Uy{i}.evaluationAtMultiplePoints = true;
        fun_Rz{i}.evaluationAtMultiplePoints = true;
        fun_N{i}.evaluationAtMultiplePoints  = true;
        fun_Ty{i}.evaluationAtMultiplePoints = true;
        fun_Mz{i}.evaluationAtMultiplePoints = true;
        fun_Epsx{i}.evaluationAtMultiplePoints = true;
        fun_Gamz{i}.evaluationAtMultiplePoints = true;
    end
    
    dim = zeros(getnbgroupelem(S),1);
    numnode = cell(getnbgroupelem(S),1);
    Ux_ex = zeros(getnbnode(S),1);
    Uy_ex = zeros(getnbnode(S),1);
    Rz_ex = zeros(getnbnode(S),1);
    N_ex  = cellfun(@(x) zeros(size(x)),getvalue(N),'UniformOutput',false);
    Mz_ex = cellfun(@(x) zeros(size(x)),getvalue(Mz),'UniformOutput',false);
    Epsx_ex = cellfun(@(x) zeros(size(x)),getvalue(Epsx),'UniformOutput',false);
    Gamz_ex = cellfun(@(x) zeros(size(x)),getvalue(Gamz),'UniformOutput',false);
    s_ex = cellfun(@(x) zeros(size(x)),getvalue(s),'UniformOutput',false);
    e_ex = cellfun(@(x) zeros(size(x)),getvalue(e),'UniformOutput',false);
    for i=1:getnbgroupelem(S)
        numnode = getnumnodeingroupelem(S,i);
        if i==1
            xi = x(numnode,2);
            Ux_ex(numnode) = -fun_Uy{1}(xi);
            Uy_ex(numnode) = fun_Ux{1}(xi);
            Rz_ex(numnode) = fun_Rz{1}(xi);
            N_ex{i}(:,:,numnode)  = fun_N{1}(xi);
            Mz_ex{i}(:,:,numnode) = fun_Mz{1}(xi);
            Epsx_ex{i}(:,:,numnode) = fun_Epsx{1}(xi);
            Gamz_ex{i}(:,:,numnode) = fun_Gamz{1}(xi);
        elseif i==2
            xi = x(numnode,1);
            Ux_ex(numnode) = fun_Ux{2}(xi);
            Uy_ex(numnode) = fun_Uy{2}(xi);
            Rz_ex(numnode) = fun_Rz{2}(xi);
            N_ex{i}(:,:,numnode)  = fun_N{2}(xi);
            Mz_ex{i}(:,:,numnode) = fun_Mz{2}(xi);
            Epsx_ex{i}(:,:,numnode) = fun_Epsx{2}(xi);
            Gamz_ex{i}(:,:,numnode) = fun_Gamz{2}(xi);
        elseif i==3
            xi = 2*L2-fliplr(x(numnode,1));
            Ux_ex(numnode) = -fun_Ux{2}(xi);
            Uy_ex(numnode) = fun_Uy{2}(xi);
            Rz_ex(numnode) = -fun_Rz{2}(xi);
            N_ex{i}(:,:,numnode)  = fun_N{2}(xi);
            Mz_ex{i}(:,:,numnode) = fun_Mz{2}(xi);
            Epsx_ex{i}(:,:,numnode) = fun_Epsx{2}(xi);
            Gamz_ex{i}(:,:,numnode) = fun_Gamz{2}(xi);
        elseif i==4
            xi = x(numnode,2);
            Ux_ex(numnode) = fun_Uy{1}(xi);
            Uy_ex(numnode) = fun_Ux{1}(xi);
            Rz_ex(numnode) = -fun_Rz{1}(xi);
            N_ex{i}(:,:,numnode)  = fun_N{1}(xi);
            Mz_ex{i}(:,:,numnode) = -fun_Mz{1}(xi);
            Epsx_ex{i}(:,:,numnode) = fun_Epsx{1}(xi);
            Gamz_ex{i}(:,:,numnode) = -fun_Gamz{1}(xi);
        end
        s_ex{i}(1,:,numnode) = N_ex{i}(1,:,numnode);
        s_ex{i}(2,:,numnode) = Mz_ex{i}(1,:,numnode);
        e_ex{i}(1,:,numnode) = Epsx_ex{i}(1,:,numnode);
        e_ex{i}(2,:,numnode) = Gamz_ex{i}(1,:,numnode);
    end
    N_ex  = FEELEMFIELD(N_ex,'type',gettype(N),'storage',getstorage(N),'ddl',get(N,'ddl'));
    Mz_ex = FEELEMFIELD(Mz_ex,'type',gettype(Mz),'storage',getstorage(Mz),'ddl',get(Mz,'ddl'));
    Epsx_ex = FEELEMFIELD(Epsx_ex,'type',gettype(Epsx),'storage',getstorage(Epsx),'ddl',get(Epsx,'ddl'));
    Gamz_ex = FEELEMFIELD(Gamz_ex,'type',gettype(Gamz),'storage',getstorage(Gamz),'ddl',get(Gamz,'ddl'));
    s_ex = FEELEMFIELD(s_ex,'type',gettype(s),'storage',getstorage(s),'ddl',get(s,'ddl'));
    e_ex = FEELEMFIELD(e_ex,'type',gettype(e),'storage',getstorage(e),'ddl',get(e,'ddl'));
    
    u_ex = [Ux_ex Uy_ex Rz_ex]';
    u_ex = u_ex(:);
    
    err_Ux = norm(Ux-Ux_ex)/norm(Ux_ex);
    err_Uy = norm(Uy-Uy_ex)/norm(Uy_ex);
    err_Rz = norm(Rz-Rz_ex)/norm(Rz_ex);
    err_N  = zeros(getnbgroupelem(S),1);
    err_Mz = zeros(getnbgroupelem(S),1);
    err_Epsx = zeros(getnbgroupelem(S),1);
    err_Gamz = zeros(getnbgroupelem(S),1);
    for i=1:getnbgroupelem(S)
        Ni     = reshape(N{i},[getnbnode(S),1]);
        Ni_ex  = reshape(N_ex{i},[getnbnode(S),1]);
        Mzi    = reshape(Mz{i},[getnbnode(S),1]);
        Mzi_ex = reshape(Mz_ex{i},[getnbnode(S),1]);
        Epsxi    = reshape(Epsx{i},[getnbnode(S),1]);
        Epsxi_ex = reshape(Epsx_ex{i},[getnbnode(S),1]);
        Gamzi    = reshape(Gamz{i},[getnbnode(S),1]);
        Gamzi_ex = reshape(Gamz_ex{i},[getnbnode(S),1]);
        err_N(i)  = norm(Ni-Ni_ex)/norm(Ni_ex);
        err_Mz(i) = norm(Mzi-Mzi_ex)/norm(Mzi_ex);
        err_Epsx(i) = norm(Epsxi-Epsxi_ex)/norm(Epsxi_ex);
        err_Gamz(i) = norm(Gamzi-Gamzi_ex)/norm(Gamzi_ex);
    end
    
    %% Test solution
    % [~,numnode2,~] = intersect(S,P2,'strict',false);
    numnode2 = find(S.node==P2);
    xP2 = x(numnode2,:);
    
    ux = eval_sol(S,u,P2,'UX');
    uy = eval_sol(S,u,P2,'UY');
    if junction
        % ux = eval_sol(S.groupelem{1},S.node,u,P2,'UX'); % ux = eval_sol(S.groupelem{2},S.node,u,P2,'UX');
        % uy = eval_sol(S.groupelem{1},S.node,u,P2,'UY'); % uy = eval_sol(S.groupelem{2},S.node,u,P2,'UY');
        rz = [eval_sol(S.groupelem{1},S.node,u,P2,'RZ') ...
            eval_sol(S.groupelem{2},S.node,u,P2,'RZ')];
        rz = double(rz);
    else
        rz = eval_sol(S,u,P2,'RZ');
    end
    
    n = 0;
    mz = 0;
    epsx = 0;
    gamz = 0;
    for i=1:getnbgroupelem(S)
        Ni  = reshape(abs(N{i}),[getnbnode(S),1]);
        Mzi = reshape(abs(Mz{i}),[getnbnode(S),1]);
        Epsxi = reshape(abs(Epsx{i}),[getnbnode(S),1]);
        Gamzi = reshape(abs(Gamz{i}),[getnbnode(S),1]);
        for j=1:length(numnode2)
            ni = double(Ni(numnode2(j)));
            mzi = double(Mzi(numnode2(j)));
            epsxi = double(Epsxi(numnode2(j)));
            gamzi = double(Gamzi(numnode2(j)));
            n = max(n,ni);
            mz = max(mz,mzi);
            epsx = max(epsx,epsxi);
            gamz = max(gamz,gamzi);
        end
    end
    
    ux_ex = Ux_ex(numnode2(1));
    uy_ex = Uy_ex(numnode2(1));
    if junction
        rz_ex(1) = Rz_ex(numnode2(1));
        rz_ex(2) = Rz_ex(numnode2(2));
    else
        rz_ex = Rz_ex(numnode2(1));
    end
    
    n_ex = 0;
    mz_ex = 0;
    epsx_ex = 0;
    gamz_ex = 0;
    for i=1:getnbgroupelem(S)
        Ni_ex  = reshape(abs(N_ex{i}),[getnbnode(S),1]);
        Mzi_ex = reshape(abs(Mz_ex{i}),[getnbnode(S),1]);
        Epsxi_ex = reshape(abs(Epsx_ex{i}),[getnbnode(S),1]);
        Gamzi_ex = reshape(abs(Gamz_ex{i}),[getnbnode(S),1]);
        for j=1:length(numnode2)
            ni_ex = double(Ni_ex(numnode2(j)));
            mzi_ex = double(Mzi_ex(numnode2(j)));
            epsxi_ex = double(Epsxi_ex(numnode2(j)));
            gamzi_ex = double(Gamzi_ex(numnode2(j)));
            n_ex = max(n_ex,ni_ex);
            mz_ex = max(mz_ex,mzi_ex);
            epsx_ex = max(epsx_ex,epsxi_ex);
            gamz_ex = max(gamz_ex,gamzi_ex);
        end
    end
    
    err_ux = norm(ux-ux_ex)/norm(ux_ex);
    err_uy = norm(uy-uy_ex)/norm(uy_ex);
    if junction
        err_rz(1) = norm(rz(1)-rz_ex(1))/norm(rz_ex(1));
        err_rz(2) = norm(rz(2)-rz_ex(2))/norm(rz_ex(2));
    else
        err_rz = norm(rz-rz_ex)/norm(rz_ex);
    end
    err_n = norm(n-n_ex)/norm(n_ex);
    err_mz = norm(mz-mz_ex)/norm(mz_ex);
    err_epsx = norm(epsx-epsx_ex)/norm(epsx_ex);
    err_gamz = norm(gamz-gamz_ex)/norm(gamz_ex);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S',...
        'L1','L2','b1','b2','b3','h',...
        'f','p','junction');
    save(fullfile(pathname,'solution.mat'),'u','s','e','time',...
        'Ux','Uy','Rz','N','Mz','Epsx','Gamz');
    save(fullfile(pathname,'reference_solution.mat'),'u_ex','s_ex','e_ex',...
        'Ux_ex','Uy_ex','Rz_ex','N_ex','Mz_ex','Epsx_ex','Gamz_ex',...
        'err_Ux','err_Uy','err_Rz','err_N','err_Mz','err_Epsx','err_Gamz');
    save(fullfile(pathname,'test_solution.mat'),'P2',...
        'ux','uy','rz','n','mz','epsx','gamz',...
        'ux_ex','uy_ex','rz_ex','n_ex','mz_ex','epsx_ex','gamz_ex',...
        'err_ux','err_uy','err_rz','err_n','err_mz','err_epsx','err_gamz');
else
    load(fullfile(pathname,'problem.mat'),'S',...
        'L1','L2','b1','b2','b3','h',...
        'f','p','junction');
    load(fullfile(pathname,'solution.mat'),'u','s','e','time',...
        'Ux','Uy','Rz','N','Mz','Epsx','Gamz');
    load(fullfile(pathname,'reference_solution.mat'),'u_ex','s_ex','e_ex',...
        'Ux_ex','Uy_ex','Rz_ex','N_ex','Mz_ex','Epsx_ex','Gamz_ex',...
        'err_Ux','err_Uy','err_Rz','err_N','err_Mz','err_Epsx','err_Gamz');
    load(fullfile(pathname,'test_solution.mat'),'P2',...
        'ux','uy','rz','n','mz','epsx','gamz',...
        'ux_ex','uy_ex','rz_ex','n_ex','mz_ex','epsx_ex','gamz_ex',...
        'err_ux','err_uy','err_rz','err_n','err_mz','err_epsx','err_gamz');
end

%% Outputs
fprintf('\nDesk\n');
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('error = %.3e for Ux\n',err_Ux);
fprintf('      = %.3e for Uy\n',err_Uy);
fprintf('      = %.3e for Rz\n',err_Rz);
for i=1:getnbgroupelem(S)
    fprintf('      = %.3e for N in groupelem #%d\n',err_N(i),i);
end
for i=1:getnbgroupelem(S)
    fprintf('      = %.3e for Mz in groupelem #%d\n',err_Mz(i),i);
end
for i=1:getnbgroupelem(S)
    fprintf('      = %.3e for Epsx in groupelem #%d\n',err_Epsx(i),i);
end
for i=1:getnbgroupelem(S)
    fprintf('      = %.3e for Gamz in groupelem #%d\n',err_Gamz(i),i);
end
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

disp('Displacement u and rotation r at point'); disp(P2);
fprintf('ux    = %g m\n',ux);
fprintf('ux_ex = %g m, error = %g\n',ux_ex,err_ux);
fprintf('uy    = %g m\n',uy);
fprintf('uy_ex = %g m, error = %g\n',uy_ex,err_uy);
if junction
    fprintf('rz(1) = %g rad = %g deg\n',rz(1),rad2deg(rz(1)));
    fprintf('rz(2) = %g rad = %g deg\n',rz(2),rad2deg(rz(2)));
    fprintf('rz_ex(1) = %g rad = %g deg, error = %g\n',rz_ex(1),rad2deg(rz_ex(1)),err_rz(1));
    fprintf('rz_ex(2) = %g rad = %g deg, error = %g\n',rz_ex(2),rad2deg(rz_ex(2)),err_rz(2));
    fprintf('|rz(1) - rz(2)|       = %g rad = %g deg\n',abs(rz(1)-rz(2)),rad2deg(abs(rz(1)-rz(2))));
    fprintf('|rz_ex(1) - rz_ex(2)| = %g rad = %g deg\n',abs(rz_ex(1)-rz_ex(2)),rad2deg(abs(rz_ex(1)-rz_ex(2))));
else
    fprintf('rz    = %g rad = %g deg\n',rz,rad2deg(rz));
    fprintf('rz_ex = %g rad = %g deg, error = %g\n',rz_ex,rad2deg(rz_ex),err_rz);
end
fprintf('\n');

disp('Force N and moment Mz at point'); disp(P2);
fprintf('N     = %g N\n',n);
fprintf('N_ex  = %g N, error = %g\n',n_ex,err_n);
fprintf('Mz    = %g N.m\n',mz);
fprintf('Mz_ex = %g N.m, error = %g\n',mz_ex,err_mz);
fprintf('\n');

disp('Axial deformation Epsx and shear deformation Gamz at point'); disp(P2);
fprintf('Epsx    = %g\n',epsx);
fprintf('Epsx_ex = %g, error = %g\n',epsx_ex,err_epsx);
fprintf('Gamz    = %g\n',gamz);
fprintf('Gamz_ex = %g, error = %g\n',gamz_ex,err_gamz);
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain(S,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 5;
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
    hP = plot(P2,'g+');
    %legend([hD,hN,hP],[legD,legN,'measure'],'Location','NorthEastOutside')
    legend([hD,hN,hP],[legD,legN,'mesure'],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','node',true,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    ampl = getsize(S)/max(abs(U))/10;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','node',true,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    plotModelDeflection(S,u_ex,'ampl',ampl,'Color','r','FaceColor','r','node',true,'legend',false);
    mysaveas(pathname,'mesh_deflected_ex',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','node',true);
    plot(S+ampl*u,'Color','b','FaceColor','b','node',true);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','node',true);
    plot(S+ampl*u_ex,'Color','r','FaceColor','r','node',true);
    mysaveas(pathname,'meshes_deflected_ex',formats,renderer);
    
    %% Display solution
    ampl = 0;
    % ampl = getsize(S)/max(abs(U))/10;
    options = {'solid',true};
    % options = {};
    
    plotSolution(S,u,'displ',1,'ampl',ampl,options{:});
    mysaveas(pathname,'Ux',formats,renderer);
    
    plotSolution(S,u_ex,'displ',1,'ampl',ampl,options{:});
    mysaveas(pathname,'Ux_ex',formats,renderer);
    
    plotSolution(S,u,'displ',2,'ampl',ampl,options{:});
    mysaveas(pathname,'Uy',formats,renderer);
    
    plotSolution(S,u_ex,'displ',2,'ampl',ampl,options{:});
    mysaveas(pathname,'Uy_ex',formats,renderer);
    
    plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
    mysaveas(pathname,'Rz',formats,renderer);
    
    plotSolution(S,u_ex,'rotation',1,'ampl',ampl,options{:});
    mysaveas(pathname,'Rz_ex',formats,renderer);
    
    % DO NOT WORK WITH BEAM ELEMENTS
    % plotSolution(S,u,'epsilon',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'eps_x',formats,renderer);
    %
    % plotSolution(S,u,'epsilon',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'gam_z',formats,renderer);
    %
    % plotSolution(S,u,'sigma',1,'ampl',ampl,options{:});
    % mysaveas(pathname,'eff_x',formats,renderer);
    %
    % plotSolution(S,u,'sigma',2,'ampl',ampl,options{:});
    % mysaveas(pathname,'mom_z',formats,renderer);
    
    u = unfreevector(S,u);
    
    figure('Name','Solution eps_x')
    clf
    plot(e,S+ampl*u,'compo','EPSX')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Epsx',formats,renderer);
    
    figure('Name','Solution eps_x_ex')
    clf
    plot(e_ex,S+ampl*u_ex,'compo','EPSX')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Epsx_ex',formats,renderer);
    
    figure('Name','Solution gam_z')
    clf
    plot(e,S+ampl*u,'compo','GAMZ')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Gamz',formats,renderer);
    
    figure('Name','Solution gam_z_ex')
    clf
    plot(e_ex,S+ampl*u_ex,'compo','GAMZ')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Gamz_ex',formats,renderer);
    
    figure('Name','Solution N')
    clf
    plot(s,S+ampl*u,'compo','EFFX')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'N',formats,renderer);
    
    figure('Name','Solution N_ex')
    clf
    plot(s_ex,S+ampl*u_ex,'compo','EFFX')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'N_ex',formats,renderer);
    
    figure('Name','Solution Mz')
    clf
    plot(s,S+ampl*u,'compo','MOMZ')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Mz',formats,renderer);
    
    figure('Name','Solution Mz_ex')
    clf
    plot(s_ex,S+ampl*u_ex,'compo','MOMZ')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'Mz_ex',formats,renderer);
end
