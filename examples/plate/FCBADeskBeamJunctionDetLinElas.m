%% FCBA desk beam deterministic linear elasticity %%
%%------------------------------------------------%%

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

tests = {'StaticVert'}; % test under static vertical load

filename = 'FCBADeskBeamDetLinElas';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','plate',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

junction = true; % junction

%% Problem
if solveProblem
    %% Domains and meshes
    % Beams dimensions
    L1 = 750e-3; % m
    b1 = 396e-3;
    L2 = 1006e-3;
    b2 = 5001e-3;
    % Same thickness for all the beams
    h = 15e-3;
    
    % Points
    P1 = POINT([0.0,0.0]);
    P2 = POINT([0.0,L1]);
    P_load = POINT([L2/2,L1]);
    P3 = POINT([L2,L1]);
    P4 = POINT([L2,0.0]);
    
    L_beam{1} = LIGNE(P1,P2);
    L_beam{2} = LIGNE(P2,P_load);
    L_beam{3} = LIGNE(P_load,P3);
    L_beam{4} = LIGNE(P4,P3);
    
    cl_beam = h/2;
    S_beam = cellfun(@(L,n) build_model(L,'cl',cl_beam,'elemtype','BEAM','filename',fullfile(pathname,['gmsh_beam_' num2str(n) '_cl_' num2str(cl_beam)])),L_beam,num2cell(1:length(L_beam)),'UniformOutput',false);
    
    %% Materials
    % Gravitational acceleration
    g = 9.81;
    % Young modulus
    E_beam = 1;
    % Poisson ratio
    NU_beam = 0.3;
    % Density
    RHO_beam = 1;
    % Cross-section area
    S1 = b1*h;
    S2 = b2*h;
    % Planar second moment of area (or Planar area moment of inertia)
    IY1 = h*b1^3/12;
    IY2 = h*b1^3/12;
    IZ1 = b1*h^3/12;
    IZ2 = b2*h^3/12;
    IX1 = IY1+IZ1;
    IX2 = IY2+IZ2;
    
    mat_beam{1} = ELAS_BEAM('E',E_beam,'NU',NU_beam,'S',S1,'IZ',IZ1,'IY',IY1,'IX',IX1,'RHO',RHO_beam);
    mat_beam{1} = setnumber(mat_beam{1},1);
    mat_beam{2} = ELAS_BEAM('E',E_beam,'NU',NU_beam,'S',S2,'IZ',IZ2,'IY',IY2,'IX',IX2,'RHO',RHO_beam);
    mat_beam{2} = setnumber(mat_beam{2},2);
    S_beam([1,4]) = cellfun(@(S) setmaterial(S,mat_beam{1}),S_beam([1,4]),'UniformOutput',false);
    S_beam([2,3]) = cellfun(@(S) setmaterial(S,mat_beam{2}),S_beam([2,3]),'UniformOutput',false);
    
    if junction
        S23 = union(S_beam{2:3});
        S = union(S_beam{1},S23,S_beam{4},'duplicate');
    else
        S = union(S_beam{:});
    end
    
    %% Neumann boundary conditions
    p = 1;
    
    %% Dirichlet boundary conditions
    if junction
        P_corner = [P2 P3];
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
    
    if junction
        k = 1e4; % additonal junction rotational stiffness
        nbddlfree = getnbddlfree(S);
        numddl2 = findddl(S,'RZ',numnode2,'free');
        numddl3 = findddl(S,'RZ',numnode3,'free');
        numddl = setdiff(1:nbddlfree,[numddl2,numddl3]);
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
    
    N = s(1);
    Mz = s(2);
    
    %% Reference solution
    if junction
        XA = 3/8*p*L2^2*IZ1/(L1*(L1*IZ2+2*L2*IZ1) + 3/S2*L2/L1*IZ1*(2*IZ2+L2/L1*IZ1));
        % MA = p/8*L2^2*IZ1*(-L1+3*L2/L1^2*IZ1/S2)/(L1*(L1*IZ2+2*L2*IZ1) + 3/S2*L2/L1*IZ1*(2*IZ2+L2/L1*IZ1));
    else
        XA = 3/8*p*L2^2*IZ1/(L1*(L1*IZ2+2*L2*IZ1) + 3/S2*L2/L1*IZ1*(2*IZ2+L2/L1*IZ1) - 6*L1*E_beam*IZ1*IZ2/k);
        % MA = p/8*L2^2*IZ1*(-L1+3*L2/L1^2*IZ1/S2)/(L1*(L1*IZ2+2*L2*IZ1) + 3/S2*L2/L1*IZ1*(2*IZ2+L2/L1*IZ1)  - 6*L1*E_beam*IZ1*IZ2/k);
    end
    YA = p/2;
    MA = XA*(-L1/3+L2/L1^2*IZ1/S2);
    
    fun_N = @(x) -YA*x(:,2);
    N2_ex = @(x) -XA;
    Ty1_ex = @(x) XA;
    Ty2_ex = @(x) -YA;
    Mz1_ex = @(x) -XA*x-MA;
    Mz2_ex = @(x) p/2*x-XA*L1-MA;
    
    ux1_ex = @(x) -p/(2*E_beam*S1)*x;
    ux2_ex = @(x) XA/(E_beam*S2)*(-x+L2/2);
    rz1_ex = @(x) -1/(E_beam*IZ1)*(XA*x.^2/2+MA*x);
    rz2_ex = @(x) 1/(E_beam*IZ2)*(p/4*(x+L2/2) - (XA*L1+MA)).*(x-L2/2);
    % rz2_ex = @(x) 1/(E_beam*IZ2)*(p/4*x.^2 - (XA*L1+MA)*x) - 1/(E_beam*IZ1)*(XA*L1^2/2+MA*L1);
    % if junction
    %     rz2_ex = @(x) rz2_ex(x) + 1/k*(XA*L1+MA);
    % end
    uy1_ex = @(x) -1/(E_beam*IZ1)*(XA*x.^3/6+MA*x.^2/2);
    uy2_ex = @(x) 1/(E_beam*IZ2)*(p/2*x.^3/6 - (XA*L1+MA)*x.^2/2) - L1/(E_beam*IZ1)*(XA*L1/2+MA)*x - p/(2*E_beam*S1)*L1;
    if junction
        uy2_ex = @(x) uy2_ex(x) + 1/k*(XA*L1+MA)*x;
    end
    
    fun_Ux = MultiVariateFunction(fun_Ux,3);
    fun_Uy = MultiVariateFunction(fun_Uy,3);
    fun_Rz = MultiVariateFunction(fun_Rz,3);
    fun_N  = MultiVariateFunction(fun_N,3);
    fun_Mz = MultiVariateFunction(fun_Mz,3);
    fun_Ux.evaluationAtMultiplePoints = true;
    fun_Uy.evaluationAtMultiplePoints = true;
    fun_Rz.evaluationAtMultiplePoints = true;
    fun_N.evaluationAtMultiplePoints  = true;
    fun_Mz.evaluationAtMultiplePoints = true;
    
    Ux_ex = fun_Ux(x);
    Uy_ex = fun_Uy(x);
    Rz_ex = fun_Rz(x);
    N_ex  = fun_N(x);
    Mz_ex = fun_Mz(x);
    
%     ind_Ux = find(~isinf(Ux_ex));
%     err_Ux = norm(Ux(ind_Ux)-Ux_ex(ind_Ux))/norm(Ux_ex(ind_Ux));
    err_Ux = norm(Ux-Ux_ex)/norm(Ux_ex);
    err_Uy = norm(Uy-Uy_ex)/norm(Uy_ex);
    err_Rz = norm(Rz-Rz_ex)/norm(Rz_ex);
    err_N = norm(N-N_ex)/norm(N_ex);
    err_Mz = norm(Mz-Mz_ex)/norm(Mz_ex);
    
    %% Test solution
    xP2 = double(getcoord(P2));
    
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
    
    [~,numnode,~] = intersect(S,P2,'strict',false);
    n = 0;
    mz = 0;
    for i=1:getnbgroupelem(S)
        Ni = reshape(abs(N{i}),[getnbnode(S),1]);
        Mzi = reshape(abs(Mz{i}),[getnbnode(S),1]);
        ni = double(Ni(numnode));
        mzi = double(Mzi(numnode));
        n = max(n,ni);
        mz = max(mz,mzi);
    end
    
%     ux_ex = fun_Ux(xP);
%     uy_ex = fun_Uy(xP);
%     rz_ex = fun_Rz(xP);
%     n_ex  = fun_N(xP);
%     mz_ex = fun_Mz(xP);
    
    ux_ex = ux2_ex(0); % ux_ex = -uy1_ex(L1);
    uy_ex = uy2_ex(0); % uy_ex = ux1_ex(L1);
    if junction
        rz_ex(1) = rz1_ex(L1);
        rz_ex(2) = rz2_ex(0);
    else
        rz_ex = rz2_ex(0); % rz_ex = rz1_ex(L1);
    end
    n_ex = max(abs(N1_ex(L1)),abs(N2_ex(0)));
    mz_ex = max(abs(Mz1_ex(L1)),abs(Mz2_ex(0)));
    
    err_ux = norm(ux-ux_ex)/norm(ux_ex);
    err_uy = norm(uy-uy_ex)/norm(uy_ex);
    if junction
        err_rz(1) = norm(rz(1)-rz_ex(1))/norm(rz_ex(1));
        err_rz(2) = norm(rz(2)-rz_ex(2))/norm(rz_ex(2));
    else
        err_rz = norm(rz-rz_ex)/norm(rz_ex);
    end
    err_n = norm(N-n_ex)/norm(n_ex);
    err_mz = norm(mz-mz_ex)/norm(mz_ex);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'S',...
        'L1','L2','b1','b2','h',...
        'f','p','junction');
    save(fullfile(pathname,'solution.mat'),'u','time',...
        'Ux','Uy','Rz','N','Mz');
    save(fullfile(pathname,'reference_solution.mat'),...
        'Ux_ex','Uy_ex','Rz_ex','N_ex','Mz_ex',...
        'err_Ux','err_Uy','err_Rz','err_N','err_Mz');
    save(fullfile(pathname,'test_solution.mat'),'P2',...
        'ux','uy','rz','n','mz',...
        'ux_ex','uy_ex','rz_ex','n_ex','mz_ex',...
        'err_ux','err_uy','err_rz','err_n','err_mz');
else
    load(fullfile(pathname,'problem.mat'),'S',...
        'L1','L2','b1','b2','h',...
        'f','p','junction');
    load(fullfile(pathname,'solution.mat'),'u','time',...
        'Ux','Uy','Rz','N','Mz');
    load(fullfile(pathname,'reference_solution.mat'),...
        'Ux_ex','Uy_ex','Rz_ex','N_ex','Mz_ex',...
        'err_Ux','err_Uy','err_Rz','err_N','err_Mz');
    load(fullfile(pathname,'test_solution.mat'),'P2',...
        'ux','uy','rz','n','mz',...
        'ux_ex','uy_ex','rz_ex','n_ex','mz_ex',...
        'err_ux','err_uy','err_rz','err_n','err_mz');
end

%% Outputs
fprintf('\nDesk\n');
fprintf('nb elements = %g\n',getnbelem(S));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('elapsed time = %f s\n',time);
fprintf('\n');

disp('Displacement u and rotation r at point'); disp(P2);
fprintf('ux    = %g mm\n',ux);
fprintf('ux_ex = %g mm, error = %g\n',ux_ex,err_ux);
fprintf('uy    = %g mm\n',uy);
fprintf('uy_ex = %g mm, error = %g\n',uy_ex,err_uy);
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
fprintf('N     = %g N\n',N);
fprintf('N_ex  = %g N, error = %g\n',N_ex,err_N);
fprintf('Mz    = %g N.mm\n',Mz);
fprintf('Mz_ex = %g N.mm, error = %g\n',Mz_ex,err_Mz);
fprintf('\n');

disp('Maximum force N and moment Mz');
fprintf('N_max     = %g N\n',N_max);
fprintf('N_max_ex  = %g N, error = %g\n',n_max_ex,err_N_max);
fprintf('Mz_max    = %g N.mm\n',Mz_max);
fprintf('Mz_max_ex = %g N.mm, error = %g\n',mz_max_ex,err_Mz_max);
fprintf('\n');

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain(S,'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 0.1;
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',1);
    hP = plot(P2,'g+');
    legend([hD,hN,hP],[legD,legN,'measure'],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','node',true,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    % ampl = getsize(S)/max(abs(u))/5;
    ampl = 2;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','node',true,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','node',true);
    plot(S+ampl*u,'Color','b','FaceColor','b','node',true);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
    
    ampl = 0;
    % ampl = getsize(S)/max(abs(u))/5;
    
    plotSolution(S,u,'displ',1,'ampl',ampl);
    mysaveas(pathname,'u_x',formats,renderer);
    
    plotSolution(S,u,'displ',2,'ampl',ampl);
    mysaveas(pathname,'u_y',formats,renderer);
    
    plotSolution(S,u,'displ',3,'ampl',ampl);
    mysaveas(pathname,'r_z',formats,renderer)
    
    % DO NOT WORK WITH BEAM ELEMENTS
    % plotSolution(S,u,'epsilon',1,'ampl',ampl);
    % mysaveas(pathname,'eps_x',formats,renderer);
    %
    % plotSolution(S,u,'epsilon',2,'ampl',ampl);
    % mysaveas(pathname,'gam_z',formats,renderer);
    %
    % plotSolution(S,u,'sigma',1,'ampl',ampl);
    % mysaveas(pathname,'eff_x',formats,renderer);
    %
    % plotSolution(S,u,'sigma',2,'ampl',ampl);
    % mysaveas(pathname,'mom_z',formats,renderer);
    
    u = unfreevector(S,u);
    
    figure('Name','Solution eps_x')
    clf
    plot(e,S+ampl*u,'compo','EPSX')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'eps_x',formats,renderer);
    
    figure('Name','Solution gam_z')
    clf
    plot(e,S+ampl*u,'compo','GAMZ')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'gam_z',formats,renderer);
    
    figure('Name','Solution eff_x')
    clf
    plot(s,S+ampl*u,'compo','EFFX')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'eff_x',formats,renderer);
    
    figure('Name','Solution mom_z')
    clf
    plot(s,S+ampl*u,'compo','MOMZ')
    colorbar
    set(gca,'FontSize',fontsize)
    mysaveas(pathname,'mom_z',formats,renderer);
end
