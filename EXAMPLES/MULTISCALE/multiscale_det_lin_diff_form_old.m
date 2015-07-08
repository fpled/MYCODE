%% Multiscale deterministic linear diffusion form old %%
%%----------------------------------------------------%%

% clc
clear all
close all

%% Input data
FileName = 'multiscale_det_lin_diff_form_old';
PathName = [getfemobjectoptions('path') 'MYCODE/RESULTS/' FileName '/'];
% PathName = '/Users/Op/Dropbox/ANTHONY-FLORENT-MATHILDE/PATCH_NONLINEAR/figures/';
if ~exist(PathName,'dir')
    dos(['mkdir ' PathName]);
end
set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object
fontsize = 16;

dim = 2; % dimension

% Data of domain D and mesh S
% P = POINT([0.0,0.0;1.0,1.0]); % size of domain D defined from point P
P = POINT([0.0,0.0;5.0,5.0]);
% nbelem = [50,50]; % number of elements for mesh S
nbelem = [20,20];
cl = 0.25; % characteristic length for mesh S
regular_mesh = 1;

% Data of domain D_patch{k} and mesh S_patch{k}
nbpatch = 4;
P_patch = cell(1,nbpatch);
nbelem_patch = cell(1,nbpatch);
cl_patch = cell(1,nbpatch);
for k=1:nbpatch
    if k == 1
%         P_patch{k} = POINT([0.5,0.5;0.9,0.9]); % size of domain D_patch{k} defined from point P_patch{k}
%         nbelem_patch{k} = [20,20]; % number of elements for mesh S_patch{k}
%         cl_patch{k} = 0.1; % characteristic length for mesh s_patch{k}
        P_patch{k} = POINT([3.0,3.0;4.0,4.0]);
        nbelem_patch{k} = [40,40];
        cl_patch{k} = 0.025;
    elseif k ==2
        P_patch{k} = POINT([1.0,0.5;2.0,1.5]);
        nbelem_patch{k} = [40,40];
        cl_patch{k} = 0.025;
    elseif k ==3
        P_patch{k} = POINT([1.0,2.0;2.0,3.0]);
        nbelem_patch{k} = [40,40];
        cl_patch{k} = 0.025;
    elseif k ==4
        P_patch{k} = POINT([3.0,0.5;4.0,1.5]);
        nbelem_patch{k} = [40,40];
        cl_patch{k} = 0.025;
    end
end
regular_mesh_patch = 1;

%% Parameters

full_projection = 0; % compute full projection operator

% Initial problem
solve_reference = true; % solve initial reference problem
save_reference = false; % save reference solution
load_reference = false; % load reference solution
change_of_variable_initial_problem = 0; % reformulation of initial problem by introducing a change of variable w=w_tilde+z
initial_solver = 'direct'; % 'direct' : direct resolution of initial problem
                           % 'pcg'    : iterative resolution of initial problem (Preconditioned conjugate gradients method)
                           % 'cgs'    : iterative resolution of initial problem (Conjugate gradients squared method)
                           % 'gmres'  : iterative resolution of initial problem (Generalized minimum residual method)
tol_initial_solver = 1e-15; % tolerance for iterative (pcg or cgs or gmres) solver of initial problem
maxiter_initial_solver = 100; % maximum number of iterations for iterative (pcg or cgs or gmres) solver of initial problem
precond_initial_solver = []; % (symmetric positive definite if pcg solver) preconditioner for solving initial problem
                             % []      : no preconditioner
                             % 'inv'   : Inverse of operator
                             % 'ichol' : Incomplete Cholesky factorization (if change_of_variable_initial_problem ~= 0)
                             % 'ilu'   : Sparse incomplete LU factorization (if change_of_variable_initial_problem ~= 0)
restart_initial_solver = []; % restarts gmres solver every restart inner iterations for initial problem

% Global-local iterative algorithm
display = true; % display error and stagnation indicators at each step
maxiter = 20; % maximum number of iterations
tol = 1e-13; % prescibed tolerance
optimal_rho = 1; % compute optimal relaxation parameter rho_opt and set rho = rho_opt
optimal_rho_approximation = 0; % compute approximation of optimal relaxation parameter rho_opt_approx and set rho = rho_opt_approx (unless optimal_rho = 1)
rho = 1; % relaxation parameter (unless optimal_rho_approximation = 1 or optimal_rho = 1)

% Global step
overlapping_domains = 1; % reformulation of global problem with overlapping domains
global_increment = 1;    % reformulation of global problem on increments (instead of current iterates)
global_solver = 'direct'; % 'direct' : direct resolution of global problem
                          % 'pcg'    : iterative resolution of global problem (Preconditioned conjugate gradients method)
                          % 'cgs'    : iterative resolution of global problem (Conjugate gradients squared method)
                          % 'gmres'  : iterative resolution of global problem (Generalized minimum residual method)
tol_global_solver = 1e-15; % tolerance for iterative (pcg or cgs or gmres) solver of global problem
maxiter_global_solver = 100; % maximum number of iterations for iterative (pcg or cgs or gmres) solver of global problem
precond_global_solver = []; % (symmetric positive definite if pcg solver) preconditioner for solving global problem
                            % []      : no preconditioner
                            % 'inv'   : Inverse of operator
                            % 'ichol' : Incomplete Cholesky factorization
                            % 'ilu'   : Sparse incomplete LU factorization
restart_global_solver = []; % restarts gmres solver every restart inner iterations for global problem
global_perturbation = 0; % introduction of a perturbation into global solution U through a truncated svd
tol_global_perturbation = 1e-3; % prescribed tolerance in truncated svd for approximation of global solution U
verification_global_prolongation = 0; % (only if overlapping_domains = 1) verification that prolongation of global solution U in the fictitious patch is uniquely defined and belongs to a suitable subspace

% Local step
change_of_variable_local_problem = 0; % reformulation of local problems by introducing a change of variable w=w_tilde+z
local_increment = 1;    % reformulation of local problems on increments (instead of current iterates)
local_solver = 'direct'; % 'direct' : direct resolution of local problems
                         % 'pcg'   : iterative resolution of local problems (Preconditioned conjugate gradients method) (if change_of_variable_local_problem ~= 0)
                         % 'cgs'   : iterative resolution of local problems (Conjugate gradients squared method)
                         % 'gmres' : iterative resolution of local problems (Generalized minimum residual method)
tol_local_solver = 1e-15; % tolerance for iterative (pcg or cgs or gmres) solver of local problems
maxiter_local_solver = 100; % maximum number of iterations for iterative (pcg or cgs or gmres) solver of local problems
precond_local_solver = []; % (symmetric positive definite if pcg solver) preconditioner for solving local problems
                           % []      : no preconditioner
                           % 'inv'   : Inverse of operator
                           % 'ichol' : Incomplete Cholesky factorization (if change_of_variable_local_problem ~= 0)
                           % 'ilu'   : Sparse incomplete LU factorization (if change_of_variable_local_problem ~= 0)
restart_local_solver = []; % restarts gmres solver every restart inner iterations for local problems
local_perturbation = 0; % introduction of a perturbation into local solutions w through a truncated svd
tol_local_perturbation = 1e-3; % prescribed tolerance in truncated svd for approximation of local solutions w

%% Domain and mesh definition

% Domain D
D = DOMAIN(dim,P(1),P(2));
% Patch D_patch{k}
D_patch = cell(1,nbpatch);
for k=1:nbpatch
    D_patch{k} = DOMAIN(dim,P_patch{k}(1),P_patch{k}(2));
end

% Coarse mesh S of domain D
if regular_mesh
    S = mesh(D,nbelem(1),nbelem(2));
    S = convertelem(S,'TRI3');
else
    S = gmsh(D,cl);
end
ddlprimal = DDL(DDLSCAL('u'));
ddldual = DDL(DDLSCAL('q'));
S = final(S,ddlprimal,ddldual);
% S = addcl(S,P(1),'u',0);
S = addcl(S,[],'u',0);

% Fine meshes S_patch{k} and S_patch_0{k} of patch D_patch{k}
S_patch = cell(1,nbpatch);
S_patch_0 = cell(1,nbpatch);
for k=1:nbpatch
    S_patch{k} = mesh(D_patch{k},nbelem_patch{k}(1),nbelem_patch{k}(2));
    S_patch{k} = convertelem(S_patch{k},'TRI3');
    S_patch{k} = final(S_patch{k},ddlprimal,ddldual);
    S_patch_0{k} = addcl(S_patch{k},[],'u',0);
end

% Partition of elements of mesh S
% 'partition' = 0 : group of elements of mesh S outside all patches D_patch{k} (or all associated meshes S_patch{k}) (not in a strict sense)
% 'partition' = k : group of elements of mesh S inside patch D_patch{k} (or associated mesh S_patch{k}) (in a strict sense)
S = setparamgroupelem(S,'partition',0);
for k=1:nbpatch
    [~,~,numelem] = intersect(S,D_patch{k},'strict',1);
    [S,newgroupelem] = separateelemwithnum(S,numelem);
    S = setparamgroupelem(S,'partition',k,newgroupelem);
end

figure('Name','Partition of elements of mesh T_H^Omega of domain Omega')
% set(gcf,'Name','Partition of elements of mesh T_H^Omega of domain Omega')
clf
plotparamelem(S,'partition');
set(gca,'FontSize',fontsize)
title('Partition of elements of mesh T_H^{\Omega} of domain \Omega')

% Coarse mesh S_out associated to 'partition' = 0
numgroupelem = getnumgroupelemwithparam(S,'partition',0);
S_out = keepgroupelem(S,numgroupelem);
S_out = removenodewithoutelem(S_out);
S_out = keepeleminnode(S_out);
S_out = final(S_out,ddlprimal,ddldual);
S_out = transfercl(S,S_out);

% Coarse meshes S_in{k} and S_in_0{k} associated to 'partition' = k
S_in = cell(1,nbpatch);
S_in_0 = cell(1,nbpatch);
for k=1:nbpatch
    numgroupelem = getnumgroupelemwithparam(S,'partition',k);
    S_in{k} = keepgroupelem(S,numgroupelem);
    S_in{k} = removenodewithoutelem(S_in{k});
    S_in{k} = keepeleminnode(S_in{k});
    S_in{k} = final(S_in{k},ddlprimal,ddldual);
    S_in_0{k} = addcl(S_in{k},[],'u',0);
    S_in{k} = transfercl(S,S_in{k});
end

figure('Name','Meshes')
% set(gcf,'Name','Meshes')
clf
plot(S_out,'facecolor',getfacecolor(1),'edgecolor','k');
leg = {'T_H^{\Omega \\ \Lambda} of domain \Omega \\ \Lambda'};
for k=1:nbpatch
    % plot(S_in{k},'facecolor',getfacecolor(k+1),'edgecolor','k');
    % leg = [leg, {['T_H^{\Lambda_{' num2str(k) '}} of patch \Lambda_{' num2str(k) '}']}];
    plot(S_patch{k},'facecolor',getfacecolor(k+1),'edgecolor','k');
    leg = [leg, {['T_h^{\Lambda_{' num2str(k) '}} of patch \Lambda_{' num2str(k) '}']}];
end
legend(leg{:})
set(gca,'FontSize',fontsize)
title(['Meshes T_H^{\Omega \\ \Lambda} of domain \Omega \\ \Lambda'...%', T_H^{\Lambda_k} of fictitious patches \Lambda_k'...
    ' and T_h^{\Lambda_k} of patches \Lambda_k'])

% Boundary mesh B of mesh S
B = create_boundary(S);
B = final(B,ddlprimal,ddldual);
B = transfercl(S,B);

% Boundary mesh parts B_out_1, B_out_2 of mesh S_out
D_out_1 = getedge(D,2);
B_out_1 = intersect(B,D_out_1);
B_out_1 = final(B_out_1,ddlprimal,ddldual);

D_out_2 = getedge(D,4);
B_out_2 = intersect(B,D_out_2);
B_out_2 = final(B_out_2,ddlprimal,ddldual);

% Boundary mesh B_in{k} of mesh S_in{k}
B_in = cell(1,nbpatch);
for k=1:nbpatch
    B_in{k} = create_boundary(S_in{k});
    B_in{k} = final(B_in{k},ddlprimal,ddldual);
    B_in{k} = transfercl(S_in{k},B_in{k});
end

% Boundary mesh B_patch{k} of mesh S_patch{k}
B_patch = cell(1,nbpatch);
for k=1:nbpatch
    B_patch{k} = create_boundary(S_patch{k});
    B_patch{k} = final(B_patch{k},ddlprimal,ddldual);
    B_patch{k} = transfercl(S_patch{k},B_patch{k});
end

%% Projection operator

% Projection operator from mesh S to mesh S_out
P_from_S_to_S_out = calc_P(S,S_out);
P_from_S_to_S_out = freevector(S,P_from_S_to_S_out,2);
P_from_S_to_S_out = freevector(S_out,P_from_S_to_S_out,1);

% Projection operator from mesh S to boundary mesh part B_out_1
P_from_S_to_B_out_1 = calc_P(S,B_out_1);
P_from_S_to_B_out_1 = freevector(S,P_from_S_to_B_out_1,2);

% Projection operator from mesh S to boundary mesh part B_out_2
P_from_S_to_B_out_2 = calc_P(S,B_out_2);
P_from_S_to_B_out_2 = freevector(S,P_from_S_to_B_out_2,2);

% Projection operator from mesh S_out to boundary mesh part B_out_1
P_from_S_out_to_B_out_1 = calc_P(S_out,B_out_1);
P_from_S_out_to_B_out_1 = freevector(S_out,P_from_S_out_to_B_out_1,2);

% Projection operator from mesh S_out to boundary mesh part B_out_2
P_from_S_out_to_B_out_2 = calc_P(S_out,B_out_2);
P_from_S_out_to_B_out_2 = freevector(S_out,P_from_S_out_to_B_out_2,2);

% Projection operator from mesh S_out to boundary mesh B_patch{k}
P_from_S_out_to_B_patch = cell(1,nbpatch);
numnode_coupling = cell(1,nbpatch);
I = speye(S_out.nbddl,S_out.nbddl);
for k=1:nbpatch
    if full_projection
        numnode_coupling{k} = 1:S_out.nbnode;
        P_from_S_out_to_B_patch{k} = eval_sol(S_out,I,B_patch{k},'u');
        P_from_S_out_to_B_patch{k} = sparse(squeeze(P_from_S_out_to_B_patch{k})');
    else
        P_from_S_out_to_B_patch{k} = sparse(B_patch{k}.nbddl,S_out.nbddl);
        [~,numnode_coupling{k},~] = intersect(S_out,D_patch{k},'strict',0);
        numddl_coupling = findddl(S_out,'all',numnode_coupling{k});
        P_from_S_out_to_B_patch_coupling = eval_sol(S_out,I(:,numddl_coupling),B_patch{k},'u');
        P_from_S_out_to_B_patch_coupling = sparse(squeeze(P_from_S_out_to_B_patch_coupling)');
        P_from_S_out_to_B_patch{k}(:,numddl_coupling) = P_from_S_out_to_B_patch_coupling;
    end
    P_from_S_out_to_B_patch{k} = freevector(S_out,P_from_S_out_to_B_patch{k},2);
end

figure('Name','Coupling nodes for projection operator from mesh T_H^Omega\Lambda of domain Omega\Lambda to boundary meshes T_h^Gamma_k of interfaces Gamma_k')
% set(gcf,'Name','Coupling nodes for projection operator from mesh T_H^Omega\Lambda of domain Omega\Lambda to boundary meshes T_h^Gamma_k of interfaces Gamma_k')
clf
plot(S_out,'facecolor',getfacecolor(1),'edgecolor','k');
leg = {'T_H^{\Omega \\ \Lambda} of domain \Omega \\ \Lambda'};
node_coupling = [];
for k=1:nbpatch
    plot(S_patch{k},'facecolor',getfacecolor(k+1),'edgecolor','k');
    leg = [leg, {['T_h^{\Lambda_{' num2str(k) '}} of patch \Lambda_{' num2str(k) '}']}];
    node_coupling = union(node_coupling,numnode_coupling{k});
end
plot(S_out.node(node_coupling),'o','color','b');
leg = [leg, {'i \in T_H^{\Omega \\ \Lambda}'}];
legend(leg{:})
set(gca,'FontSize',fontsize)
title('Coupling nodes for projection operator from mesh T_H^{\Omega \\ \Lambda} of domain \Omega \\ \Lambda to boundary meshes T_h^{\Gamma_k} of interfaces \Gamma_k')

% Projection operator from mesh S to boundary mesh B_patch{k}
P_from_S_to_B_patch = cell(1,nbpatch);
for k=1:nbpatch
    P_from_S_to_B_patch{k} = P_from_S_out_to_B_patch{k}*P_from_S_to_S_out;
end

% Projection operator from mesh S_patch{k} to boundary mesh B_patch{k}
P_from_S_patch_to_B_patch = cell(1,nbpatch);
for k=1:nbpatch
    P_from_S_patch_to_B_patch{k} = calc_P(S_patch{k},B_patch{k});
    P_from_S_patch_to_B_patch{k} = freevector(S_patch{k},P_from_S_patch_to_B_patch{k},2);
end

% Projection operator from mesh S_out to mesh S_patch{k}
P_from_S_out_to_S_patch = cell(1,nbpatch);
for k=1:nbpatch
    P_from_S_out_to_S_patch{k} = P_from_S_patch_to_B_patch{k}'*P_from_S_out_to_B_patch{k};
end

% Projection operator from mesh S to mesh S_patch{k}
P_from_S_to_S_patch = cell(1,nbpatch);
for k=1:nbpatch
    P_from_S_to_S_patch{k} = P_from_S_out_to_S_patch{k}*P_from_S_to_S_out;
end

% Projection operator from mesh S to mesh S_in{k}
P_from_S_to_S_in = cell(1,nbpatch);
for k=1:nbpatch
    P_from_S_to_S_in{k} = calc_P(S,S_in{k});
    P_from_S_to_S_in{k} = freevector(S,P_from_S_to_S_in{k},2);
    P_from_S_to_S_in{k} = freevector(S_in{k},P_from_S_to_S_in{k},1);
end

% Projection operator from mesh S_in{k} to boundary mesh B_in{k}
P_from_S_in_to_B_in = cell(1,nbpatch);
for k=1:nbpatch
    P_from_S_in_to_B_in{k} = calc_P(S_in{k},B_in{k});
    P_from_S_in_to_B_in{k} = freevector(S_in{k},P_from_S_in_to_B_in{k},2);
end

%% Mass matrix

% Bilinear form a
% a(u,v) = int( u.v )
a = BILINFORM(0,0,1);
a = setfree(a,0);

% Mass matrix M_B_patch{k} associated to boundary mesh B_patch{k}
M_B_patch = cell(1,nbpatch);
for k=1:nbpatch
    M_B_patch{k} = calc_matrix(a,B_patch{k});
    % M_B_patch{k} = a{B_patch{k}}(:,:);
end

% Mass matrix M_B_in{k} associated to boundary mesh B_in{k}
M_B_in = cell(1,nbpatch);
for k=1:nbpatch
    M_B_in{k} = calc_matrix(a,B_in{k});
    % M_B_in{k} = a{B_in{k}}(:,:);
end

%% Bilinear forms and linear forms associated to initial problem

% Bilinear form a
% a(u,v) = int( K.grad(u).grad(v) )
K = 1;
% a = BILINFORM(1,1,K);
a = DIFFUSIONFORM(K);
a = setfree(a,1);

% Linear form l
f = 1;
% l_out_1 = LINFORM(0,f);
% l_out_1 = setfree(l_out_1,1);
% l_out_2 = LINFORM(0,-f);
% l_out_2 = setfree(l_out_2,1);
% l = LINFORM(0,0);
% l(v) = int( f.v )
l = LINFORM(0,f);
l = setfree(l,1);

% Diffusion coefficients K_patch{k} evaluated on mesh S_patch{k} and K_in{k} evaluated on mesh S_in{k}
K_patch = cell(1,nbpatch);
K_in = cell(1,nbpatch);
for k=1:nbpatch
    % k_patch{k} = alpha*exp( -Amp*||x-C{k}||_2^2/L^2 ) if ||x-C{k}||_Inf < L
    %            = 0                                    if ||x-C{k}||_Inf >= L
    % K_patch{k}(x) = K +        k_patch{k}(x)
    % K_in{k}(x)    = K + beta * k_patch{k}(x)
    % alpha = 10;
    % C = getcenter(D_patch{k});
    % Amp = 2;
    % L = norm(getsize(D_patch{k}),Inf)/4;
    % k_patch = @(x) (distance(x,C,Inf)<L) * alpha * exp(-Amp*distance(x,C,2).^2/L^2);
    % P_patch = POINT(S_patch{k}.node);
    % K_patch{k} = K + k_patch(P_patch);
    % K_patch{k} = squeeze(K_patch{k});
    % 
    % beta = 0;
    % P_in = POINT(S_in{k}.node);
    % K_in{k} = K + beta * k_patch(P_in);
    % K_in{k} = squeeze(K_in{k});
    
    % k_patch{k} = 1 if ||x-C{k}||_Inf < L
    %            = 0 if ||x-C{k}||_Inf >= L
    % K_patch{k}(x) = 1 + k_patch{k}(x)
    % K_in{k}(x)    = 1
    C = getcenter(D_patch{k});
    L = norm(getsize(D_patch{k}),Inf)/4;
    k_patch = @(x) distance(x,C,Inf)<L;
    P_patch = POINT(S_patch{k}.node);
    K_patch{k} = ones(getnbnode(S_patch{k}),1) + double(squeeze(k_patch(P_patch)));
    
    K_in{k} = 1;
end

% Bilinear form a_patch{k} associated to diffusion coefficient kappa_patch{k}
% a_patch{k}(u,v) = int( kappa_patch{k}.grad(u).grad(v) )
a_patch = cell(1,nbpatch);
for k=1:nbpatch
    % a_patch{k} = BILINFORM(1,1,K_patch{k}); % uniform value
    % a_patch{k} = BILINFORM(1,1,K_patch{k},0); % nodal values
    a_patch{k} = DIFFUSIONFORM(K_patch{k});
    a_patch{k} = setfree(a_patch{k},1);
end

% Bilinear form l_patch{k}
% l_patch{k}(v) = l(v) = int( f.v )
l_patch = cell(1,nbpatch);
for k=1:nbpatch
    l_patch{k} = l;
end

% Bilinear form a_in{k} associated to fictitious diffusion coefficient kappa_in{k}
% a_in{k}(u,v) = int( kappa_in{k}.grad(u).grad(v) )
a_in = cell(1,nbpatch);
for k=1:nbpatch
    % a_in{k} = BILINFORM(1,1,K_in{k}); % uniform value
    % a_in{k} = BILINFORM(1,1,K_in{k},0); % nodal values
    a_in{k} = DIFFUSIONFORM(K_in{k});
    a_in{k} = setfree(a_in{k},1);
end

% Display diffusion coefficients K, K_patch{k}, K_in{k}
% figure('Name','Diffusion coefficients K')
% % set(gcf,'Name','Diffusion coefficients K')
% clf
% if overlapping_domains
%     nbrows = 3;
% else
%     nbrows = 2;
% end
% subplot(nbrows,nbpatch,1:nbpatch)
% plot(FENODEFIELD(K),removebc(S_out));%,'surface');
% % colormap('default')
% colorbar
% set(gca,'FontSize',fontsize)
% title('Diffusion coefficient K over domain \Omega \\ \Lambda')
% for k=1:nbpatch
%     subplot(nbrows,nbpatch,nbpatch+k)
%     plot(FENODEFIELD(K_patch{k}),S_patch{k});%,'surface');
%     % colormap('default')
%     colorbar
%     set(gca,'FontSize',fontsize)
%     title(['Diffusion coefficient K_{' num2str(k) '} over patch \Lambda_{' num2str(k) '}'])
%     if overlapping_domains
%         subplot(nbrows,nbpatch,2*nbpatch+k)
%         plot(FENODEFIELD(K_in{k}),S_in{k});%,'surface');
%         % colormap('default')
%         colorbar
%         set(gca,'FontSize',fontsize)
%         title(['Fictituous diffusion coefficient K_{' num2str(k) '} over fictitious patch \Lambda_{' num2str(k) '}'])
%     end
% end

%% Stiffness matrix and sollicitation vector

% Stiffness matrix A_out and sollicitation vector b_out associated to mesh S_out
% Corresponding prolongation A_out_tilde and b_out_tilde over mesh S for reformulation with overlapping domains
A_out = calc_matrix(a,S_out);
A_out_tilde = P_from_S_to_S_out'*A_out*P_from_S_to_S_out;
b_out = calc_vector(l,S_out);
b_out_tilde = P_from_S_to_S_out'*b_out;
% b_out_1 = calc_vector(l_out_1,B_out_1);
% b_out_2 = calc_vector(l_out_2,B_out_2);
% b_out = b_out + P_from_S_out_to_B_out_1'*b_out_1 + P_from_S_out_to_B_out_2'*b_out_2;
% b_out_tilde = b_out_tilde + P_from_S_to_B_out_1'*b_out_1 + P_from_S_to_B_out_2'*b_out_2;

% Stiffness matrix A_in{k} and sollicitation vector b_in{k} associated to mesh S_in{k}
% Corresponding prolongation A_in_tilde{k} and b_in_tilde{k} over mesh S for reformulation with overlapping domains
A_in = cell(1,nbpatch);
A_in_tilde = cell(1,nbpatch);
b_in = cell(1,nbpatch);
b_in_tilde = cell(1,nbpatch);
for k=1:nbpatch
    A_in{k} = calc_matrix(a_in{k},S_in{k});
    A_in_tilde{k} = P_from_S_to_S_in{k}'*A_in{k}*P_from_S_to_S_in{k};
    b_in{k} = calc_vector(l,S_in{k});
    b_in_tilde{k} = P_from_S_to_S_in{k}'*b_in{k};
end

% Stiffness matrix A_patch{k} and sollicitation vector b_patch{k} associated to mesh S_patch{k}
A_patch = cell(1,nbpatch);
b_patch = cell(1,nbpatch);
for k=1:nbpatch
    A_patch{k} = calc_matrix(a_patch{k},S_patch{k});
    b_patch{k} = calc_vector(l_patch{k},S_patch{k});
end

% Total matrix A_tilde associated to initial problem for reformulation with overlapping domains
A_tilde = A_out_tilde;
for k=1:nbpatch
    A_tilde = A_tilde + A_in_tilde{k};
end

% Total sollicitation vector b_tilde associated to initial problem for reformulation with overlapping domains
b_tilde = b_out_tilde;
for k=1:nbpatch
    b_tilde = b_tilde + b_in_tilde{k};
end

%% Direct resolution of initial problem based on domain decomposition
filename = 'reference_solution.mat';
fullfilename = fullfile(PathName,filename);
if solve_reference
    fprintf('\n --------------------------------------------------------------');
    fprintf('\n ------------ Direct resolution of initial problem ------------')
    fprintf('\n --------------------------------------------------------------\n');
    tic
    
    if ~change_of_variable_initial_problem
        fprintf('\nResolution of initial problem without change of variable and using a %s solver\n',initial_solver);
    else
        fprintf('\nResolution of initial problem with change of variable and using a %s solver\n',initial_solver);
    end
    
    if ~change_of_variable_initial_problem
        % Without change of variable and with Lagrange multipliers
        % Solution (U_ref,w_ref{k},lambda_ref{k})
        
        % Degrees of freedom
        nb_dof_U = getnbddlfree(S_out);
        nb_dof_w = cell(1,nbpatch);
        nb_dof_lambda = cell(1,nbpatch);
        for k=1:nbpatch
            nb_dof_w{k} = getnbddlfree(S_patch{k});
            nb_dof_lambda{k} = getnbddlfree(B_patch{k});
        end
        
        nb_block = 1+2*nbpatch;
        % Total matrix A associated to initial problem
        A = cell(nb_block);
        A{1,1} = A_out;
        for k=1:nbpatch
            for j=1:nbpatch
                A{2*k,2*j} = sparse(nb_dof_w{k},nb_dof_w{j});
                A{2*k+1,2*j+1} = sparse(nb_dof_lambda{k},nb_dof_lambda{j});
                A{2*k,2*j+1} = sparse(nb_dof_w{k},nb_dof_lambda{j});
                A{2*k+1,2*j} = sparse(nb_dof_lambda{k},nb_dof_w{j});
            end
        end
        for k=1:nbpatch
            A{2*k,2*k} = A_patch{k};
            A{1,2*k} = sparse(nb_dof_U,nb_dof_w{k});
            A{1,2*k+1} = P_from_S_out_to_B_patch{k}'*M_B_patch{k};
            A{2*k,2*k+1} = -P_from_S_patch_to_B_patch{k}'*M_B_patch{k};
            A{2*k,1} = A{1,2*k}'; % A{2*k,1} = sparse(nb_dof_w{k},nb_dof_U);
            A{2*k+1,1} = A{1,2*k+1}'; % A{2*k+1,1} = M_B_patch{k}*P_from_S_out_to_B_patch{k};
            A{2*k+1,2*k} = A{2*k,2*k+1}'; % A{2*k+1,2*k} = -M_B_patch{k}*P_from_S_patch_to_B_patch{k};
        end
        A = cell2mat(A);
        
        % Total sollicitation vector b associated to initial problem
        b = cell(nb_block,1);
        b{1} = b_out;
        for k=1:nbpatch
            b{2*k} = b_patch{k};
            b{2*k+1} = zeros(nb_dof_lambda{k},1);
        end
        b = cell2mat(b);
        
        % Total reconstructed solution vector sol_ref=[U_ref;w_ref{k};lambda_ref{k}] associated to initial problem
        if strcmp(initial_solver,'direct')
            sol_ref = solve(A,b);
        else
            if strcmp(precond_initial_solver,'inv')
                M1_initial_solver = inv(A);
                M2_initial_solver = [];
            elseif strcmp(precond_initial_solver,'ichol')
                L_chol = ichol(A);
                M1_initial_solver = L_chol;
                M2_initial_solver = L_chol';
            elseif strcmp(precond_initial_solver,'ilu')
                [L_lu,U_lu] = ilu(A);
                M1_initial_solver = L_lu;
                M2_initial_solver = U_lu;
            else
                M1_initial_solver = [];
                M2_initial_solver = [];
            end
            if strcmp(initial_solver,'pcg')
                sol_ref = pcg(A,b,tol_initial_solver,maxiter_initial_solver,M1_initial_solver,M2_initial_solver);
            elseif strcmp(initial_solver,'cgs')
                sol_ref = cgs(A,b,tol_initial_solver,maxiter_initial_solver,M1_initial_solver,M2_initial_solver);
            elseif strcmp(initial_solver,'gmres')
                sol_ref = gmres(A,b,restart_initial_solver,tol_initial_solver,maxiter_initial_solver,M1_initial_solver,M2_initial_solver);
            else
                error('Wrong solver for initial problem')
            end
        end
        rep = nb_dof_U;
        for k=1:nbpatch
            rep = [rep,nb_dof_w{k},nb_dof_lambda{k}];
        end
        sol_ref = mat2cell(sol_ref,rep);
        U_ref = sol_ref{1};
        w_ref = cell(1,nbpatch);
        lambda_ref = cell(1,nbpatch);
        for k=1:nbpatch
            w_ref{k} = sol_ref{2*k};
            lambda_ref{k} = sol_ref{2*k+1};
        end
    else
        % With change of variable and without Lagrange multipliers
        % Solution (U_ref,z_ref{k})
        
        % Degrees of freedom
        nb_dof_U = getnbddlfree(S_out);
        nb_dof_z = cell(1,nbpatch);
        for k=1:nbpatch
            nb_dof_z{k} = getnbddlfree(S_patch_0{k});
        end
        
        nb_block = 1+nbpatch;
        % Total matrix A associated to initial problem
        A = cell(nb_block);
        A{1,1} = A_out;
        for k=1:nbpatch
            for j=1:nbpatch
                A{k+1,j+1} = sparse(nb_dof_z{k},nb_dof_z{j});
            end
        end
        for k=1:nbpatch
            A{1,1} = A{1,1} + P_from_S_out_to_S_patch{k}'*A_patch{k}*P_from_S_out_to_S_patch{k};
            A{k+1,k+1} = freematrix(S_patch_0{k},A_patch{k});
            A{1,k+1} = P_from_S_out_to_S_patch{k}'*freevector(S_patch_0{k},A_patch{k},2);
            A{k+1,1} = freevector(S_patch_0{k},A_patch{k},1)*P_from_S_out_to_S_patch{k};
        end
        A = cell2mat(A);
        
        % Total sollicitation vector b associated to initial problem
        b = cell(nb_block,1);
        b{1} = b_out;
        for k=1:nbpatch
            b{1} = b{1} + P_from_S_out_to_S_patch{k}'*b_patch{k};
            b{k+1} = freevector(S_patch_0{k},b_patch{k});
        end
        b = cell2mat(b);
        
        % Total reconstructed solution vector sol_ref=[U_ref;z_ref{k}] associated to initial problem
        if strcmp(initial_solver,'direct')
            sol_ref = solve(A,b);
        else
            if strcmp(precond_initial_solver,'inv')
                M1_initial_solver = inv(A);
                M2_initial_solver = [];
            elseif strcmp(precond_initial_solver,'ichol')
                L_chol = ichol(A);
                M1_initial_solver = L_chol;
                M2_initial_solver = L_chol';
            elseif strcmp(precond_initial_solver,'ilu')
                [L_lu,U_lu] = ilu(A);
                M1_initial_solver = L_lu;
                M2_initial_solver = U_lu;
            else
                M1_initial_solver = [];
                M2_initial_solver = [];
            end
            if strcmp(initial_solver,'pcg')
                sol_ref = pcg(A,b,tol_initial_solver,maxiter_initial_solver,M1_initial_solver,M2_initial_solver);
            elseif strcmp(initial_solver,'cgs')
                sol_ref = cgs(A,b,tol_initial_solver,maxiter_initial_solver,M1_initial_solver,M2_initial_solver);
            elseif strcmp(initial_solver,'gmres')
                sol_ref = gmres(A,b,restart_initial_solver,tol_initial_solver,maxiter_initial_solver,M1_initial_solver,M2_initial_solver);
            else
                error('Wrong solver for initial problem')
            end
        end
        rep = nb_dof_U;
        for k=1:nbpatch
            rep = [rep,nb_dof_z{k}];
        end
        sol_ref = mat2cell(sol_ref,rep);
        U_ref = sol_ref{1};
        z_ref = cell(1,nbpatch);
        w_ref = cell(1,nbpatch);
        lambda_ref = cell(1,nbpatch);
        for k=1:nbpatch
            z_ref{k} = sol_ref{1+k};
            w_ref{k} = P_from_S_out_to_S_patch{k}*U_ref + unfreevector(S_patch_0{k},z_ref{k});
            
            % A = P_from_S_patch_to_B_patch{k}'*M_B_patch{k};
            % b = A_patch{k}*w_ref{k} - b_patch{k};
            A = M_B_patch{k};
            b = P_from_S_patch_to_B_patch{k}*(A_patch{k}*w_ref{k} - b_patch{k});
            lambda_ref{k} = solve(A,b);
        end
    end
    fprintf('\nElapsed time = %f s\n',toc);
    
    if save_reference
        % Save reference solution (U_ref,w_ref{k},lambda_ref{k})
        save(fullfilename,'U_ref','w_ref','lambda_ref');
    end
elseif load_reference
    if ~exist(filename,'file')
        error(['File ' filename ' does not exist in folder' PathName]);
    else
        % Load reference solution (U_ref,w_ref{k},lambda_ref{k})
        load(fullfilename,'U_ref','w_ref','lambda_ref');
    end
end

%% Display reference solution u_ref=(U_ref,w_ref)
figure('Name','Reference solution u_ref=(U_ref,w_ref)')
% set(gcf,'Name','Reference solution u_ref=(U_ref,w_ref)')
clf

subplot(1+nbpatch,2,1)
plot_sol(S_out,U_ref,'edgecolor','none');
for k=1:nbpatch
    plot_sol(S_patch{k},w_ref{k},'edgecolor','none');
    plot(D_patch{k},'facecolor','none','edgecolor','k');
    % gtext(['\Lambda_{' num2str(k) '}'],'FontSize',fontsize);
end
% colormap('default')
colorbar
ax=axis;
cax=caxis;
set(gca,'FontSize',fontsize)
title('Reference solution u^{ref}=(U^{ref},w^{ref}) over domain \Omega')

subplot(1+nbpatch,2,2)
plot_sol(S_out,U_ref,'edgecolor','none');
% colormap('default')
colorbar
axis(ax)
% caxis(cax)
set(gca,'FontSize',fontsize)
title('Global solution U^{ref} over mesh T_H^{\Omega \\ \Lambda} of domain \Omega \\ \Lambda')

for k=1:nbpatch
    subplot(1+nbpatch,2,2*k+1)
    plot_sol(S_patch{k},w_ref{k},'edgecolor','none');
    % colormap('default')
    colorbar
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',fontsize)
    title(['Local solution w^{ref}^{' num2str(k) '} over mesh T_h^{\Lambda_{' num2str(k) '}} of patch \Lambda_{' num2str(k) '}'])
    
    subplot(1+nbpatch,2,2*(k+1))
    plot_sol(B_patch{k},lambda_ref{k},'edgecolor','interp');
    % colormap('default')
    colorbar
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',fontsize)
    title(['Lagrange multiplier \lambda^{ref}^{' num2str(k) '} over boundary mesh T_h^{\Gamma_{' num2str(k) '}} of interface \Gamma_{' num2str(k) '}'])
end

%% Initial (resp. Reformulated) global-local iterative algorithm based on domain decomposition without (resp. with) overlapping domains

fprintf('\n ----------------------------------------------------------');
fprintf('\n ------------ Global-local iterative algorithm ------------')
fprintf('\n ----------------------------------------------------------\n');
tic
if ~overlapping_domains
    fprintf('\nInitial algorithm without overlapping domains\n');
else
    fprintf('\nReformulated algorithm with overlapping domains\n');
end
fprintf('\nResolution of global problems using a %s solver\n',global_solver);
if ~change_of_variable_local_problem
    fprintf('\nResolution of local problems without change of variable and using a %s solver\n',local_solver);
else
    fprintf('\nResolution of local problems with change of variable and using a %s solver\n',local_solver);
end
% Solution (U,w{k},lambda{k})

% Degrees of freedom
if overlapping_domains
    nb_dof_U = getnbddlfree(S);
else
    nb_dof_U = getnbddlfree(S_out);
end
nb_dof_w = cell(1,nbpatch);
nb_dof_lambda = cell(1,nbpatch);
for k=1:nbpatch
    nb_dof_w{k} = getnbddlfree(S_patch{k});
    nb_dof_lambda{k} = getnbddlfree(B_patch{k});
end
if change_of_variable_local_problem
    nb_dof_z = cell(1,nbpatch);
    for k=1:nbpatch
        nb_dof_z{k} = getnbddlfree(S_patch_0{k});
    end
end

% Optimal relaxation parameter rho_opt and its approximation rho_opt_approx
% Largest (resp. lowest) eigenvalue lambda_max (resp. lambda_min) of iteration operator A(V)=V-Phi(Psi(V))   without overlapping domains
%                                                                                       A(V)=V-Phi(Psi(V),V) with overlapping domains
if optimal_rho || optimal_rho_approximation
    
    A = @(V) calc_iteration_operator_old(V,nbpatch,S_patch_0,A_patch,A_out,A_in_tilde,A_tilde,M_B_patch,P_from_S_patch_to_B_patch,P_from_S_to_B_patch,P_from_S_out_to_B_patch,P_from_S_to_S_patch,P_from_S_out_to_S_patch,nb_dof_U,nb_dof_w,nb_dof_lambda,change_of_variable_local_problem,overlapping_domains);
    
    % Options for function eigs
    % opts.issym = false;
    % opts.isreal = true;
    % opts.tol = 1e-12;
    % opts.maxit = 300;
    % opts.v0 = rand(nb_dof_U,1);
    % opts.disp = 0;
    % opts.cholB = false;
    
    lambda_max = eigs(A,nb_dof_U,1,'lm');%,opts);
    if optimal_rho
        lambda_min = eigs(A,nb_dof_U,1,'sm');%,opts);
    end
    
    if optimal_rho
        rho = 2/(lambda_min + lambda_max);
        fprintf('\nOptimal relaxation parameter : rho = %f\n',rho);
    elseif optimal_rho_approximation
        rho = 1/lambda_max;
        fprintf('\nApproximation of optimal relaxation parameter : rho = %f\n',rho);
    end
else
    fprintf('\nRelaxation parameter : rho = %f\n',rho);
end

% Initialization
U_init = zeros(nb_dof_U,1);
w_init = cell(1,nbpatch);
lambda_init = cell(1,nbpatch);
if change_of_variable_local_problem
    z_init = cell(1,nbpatch);
end
for k=1:nbpatch
    w_init{k} = zeros(nb_dof_w{k},1);
    lambda_init{k} = zeros(nb_dof_lambda{k},1);
    if change_of_variable_local_problem
        z_init{k} = w_init{k};
        if overlapping_domains
            z_init{k} = z_init{k} - P_from_S_to_S_patch{k}*U_init;
        else
            z_init{k} = z_init{k} - P_from_S_out_to_S_patch{k}*U_init;
        end
        z_init{k} = freevector(S_patch_0{k},z_init{k});
    end
end

if overlapping_domains
    error_indicator_U_init = norm(P_from_S_to_S_out*U_init - U_ref)/norm(U_ref);
else
    error_indicator_U_init = norm(U_init - U_ref)/norm(U_ref);
end

error_indicator_w_init = cell(1,nbpatch);
error_indicator_lambda_init = cell(1,nbpatch);
for k=1:nbpatch
    error_indicator_w_init{k} = norm(w_init{k} - w_ref{k})/norm(w_ref{k});
    error_indicator_lambda_init{k} = norm(lambda_init{k} - lambda_ref{k})/norm(lambda_ref{k});
end

if display
    fprintf('\nInitialization : residual error = %.3e w.r.t. U',error_indicator_U_init);
    for k=1:nbpatch
        fprintf('\n                 residual error = %.3e w.r.t. w for patch #%u',error_indicator_w_init{k},k);
        fprintf('\n                 residual error = %.3e w.r.t. lambda for patch #%u\n',error_indicator_lambda_init{k},k);
    end
end

U = U_init;
w = w_init;
lambda = lambda_init;
if change_of_variable_local_problem
    z = z_init;
end

error_indicator_U = zeros(1,maxiter);
stagnation_indicator_U = zeros(1,maxiter);
error_indicator_w = cell(1,nbpatch);
stagnation_indicator_w = cell(1,nbpatch);
error_indicator_lambda = cell(1,nbpatch);
stagnation_indicator_lambda = cell(1,nbpatch);
for k=1:nbpatch
    error_indicator_w{k} = zeros(1,maxiter);
    stagnation_indicator_w{k} = zeros(1,maxiter);
    error_indicator_lambda{k} = zeros(1,maxiter);
    stagnation_indicator_lambda{k} = zeros(1,maxiter);
end

% Iteration step
tic;
for iter=1:maxiter
    
    % Update
    U_old = U;
    w_old = w;
    lambda_old = lambda;
    if change_of_variable_local_problem
        z_old = z;
    end
    
    % Global problem
    % Global solution U
    if overlapping_domains
        A = A_tilde;
        b = b_out_tilde;
        for k=1:nbpatch
            b = b - P_from_S_to_B_patch{k}'*M_B_patch{k}*lambda_old{k} + A_in_tilde{k}*U_old;
        end
    else
        A = A_out;
        b = b_out;
        for k=1:nbpatch
            b = b - P_from_S_out_to_B_patch{k}'*M_B_patch{k}*lambda_old{k};
        end
    end
    
    if ~global_increment
        if strcmp(global_solver,'direct')
            U = solve(A,b);
        else
            if strcmp(precond_global_solver,'inv')
                M1_global_solver = inv(A);
                M2_global_solver = [];
            elseif strcmp(precond_global_solver,'ichol')
                L_chol = ichol(A);
                M1_global_solver = L_chol;
                M2_global_solver = L_chol';
            elseif strcmp(precond_global_solver,'ilu')
                [L_lu,U_lu] = ilu(A);
                M1_global_solver = L_lu;
                M2_global_solver = U_lu;
            else
                M1_global_solver = [];
                M2_global_solver = [];
            end
            if strcmp(global_solver,'pcg')
                U = pcg(A,b,tol_global_solver,maxiter_global_solver,M1_global_solver,M2_global_solver);
            elseif strcmp(global_solver,'cgs')
                U = cgs(A,b,tol_global_solver,maxiter_global_solver,M1_global_solver,M2_global_solver);
            elseif strcmp(global_solver,'gmres')
                U = gmres(A,b,restart_global_solver,tol_global_solver,maxiter_global_solver,M1_global_solver,M2_global_solver);
            else
                error('Wrong solver for global problem')
            end
        end
    else
        b = b - A*U_old;
        if strcmp(global_solver,'direct')
            dU = solve(A,b);
        else
            if strcmp(precond_global_solver,'inv')
                M1_global_solver = inv(A);
                M2_global_solver = [];
            elseif strcmp(precond_global_solver,'ichol')
                L_chol = ichol(A);
                M1_global_solver = L_chol;
                M2_global_solver = L_chol';
            elseif strcmp(precond_global_solver,'ilu')
                [L_lu,U_lu] = ilu(A);
                M1_global_solver = L_lu;
                M2_global_solver = U_lu;
            else
                M1_global_solver = [];
                M2_global_solver = [];
            end
            if strcmp(global_solver,'pcg')
                dU = pcg(A,b,tol_global_solver,maxiter_global_solver,M1_global_solver,M2_global_solver);
            elseif strcmp(global_solver,'cgs')
                dU = cgs(A,b,tol_global_solver,maxiter_global_solver,M1_global_solver,M2_global_solver);
            elseif strcmp(global_solver,'gmres')
                dU = gmres(A,b,restart_global_solver,tol_global_solver,maxiter_global_solver,M1_global_solver,M2_global_solver);
            else
                error('Wrong solver for global problem')
            end
        end
        U = U_old + dU;
    end
    
    U = rho*U + (1-rho)*U_old;
    
    stagnation_indicator_U(iter) = norm(U - U_old)/norm(U);
    
    % Verification step : prolongation of global solution U in the fictitious patch is uniquely defined and belongs to a suitable subspace
    if overlapping_domains && verification_global_prolongation
        for k=1:nbpatch
            U_in_init = P_from_S_to_S_in{k}*U_init;
            
            P_from_S_in_to_S_in_0 = calc_P(S_in{k},S_in_0{k});
            P_from_S_in_to_S_in_0 = freevector(S_in{k},P_from_S_in_to_S_in_0,2);
            P_from_S_in_to_S_in_0 = freevector(S_in_0{k},P_from_S_in_to_S_in_0,1);
            
            A1 = P_from_S_in_to_S_in_0'*freevector(S_in_0{k},A_in{k},1);
            A2 = P_from_S_in_to_B_in{k}'*M_B_in{k}*P_from_S_in_to_B_in{k};
            A = A1 + A2;
            
            b1 = P_from_S_in_to_S_in_0'*freevector(S_in_0{k},A_in{k},1)*U_in_init;
            b2 = P_from_S_in_to_B_in{k}'*M_B_in{k}*P_from_S_to_B_in{k}*U;
            b = b1 + b2;
            
            U_in = solve(A,b);
            relative_error_prolongation_U = norm(P_from_S_to_S_in{k}*U - U_in)/norm(U_in);
            if relative_error_prolongation_U > tol
                fprintf('\nProlongation of U in fictitious patch #%u at iteration #%3.d : residual error = %.3e\n',k,iter,relative_error_prolongation_U);
            end
        end
    end
    
    if global_perturbation
        if overlapping_domains
            U = unfreevector(S,U);
        else
            U = unfreevector(S_out,U);
        end
        U = reshape(U,sqrt(length(U)),sqrt(length(U)));
        U = svdtruncate(U,tol_global_perturbation);
        U = U(:);
        if overlapping_domains
            U = freevector(S,U);
        else
            U = freevector(S_out,U);
        end
    end
    
    if overlapping_domains
        error_indicator_U(iter) = norm(P_from_S_to_S_out*U - U_ref)/norm(U_ref);
    else
        error_indicator_U(iter) = norm(U - U_ref)/norm(U_ref);
    end
    
    % Local problems
    % Local solutions (w{k},lambda{k}) without change of variable
    %                             z{k} with change of variable
    for k=1:nbpatch
        if ~change_of_variable_local_problem
            % Without change of variable and with Lagrange multipliers
            % Local solution (w{k},lambda{k})
            A11 = A_patch{k};
            A22 = sparse(nb_dof_lambda{k},nb_dof_lambda{k});
            A12 = -P_from_S_patch_to_B_patch{k}'*M_B_patch{k};
            A21 = M_B_patch{k}*P_from_S_patch_to_B_patch{k};
            A = [A11,A12;
                A21,A22];
            
            b1 = b_patch{k};
            if overlapping_domains
                b2 = M_B_patch{k}*P_from_S_to_B_patch{k}*U;
            else
                b2 = M_B_patch{k}*P_from_S_out_to_B_patch{k}*U;
            end
            b = [b1;b2];
            
            if ~local_increment
                if strcmp(local_solver,'direct')
                    local_sol = solve(A,b);
                else
                    if strcmp(precond_local_solver,'inv')
                        M1_local_solver = inv(A);
                        M2_local_solver = [];
                    elseif strcmp(precond_local_solver,'ichol')
                        L_chol = ichol(A);
                        M1_local_solver = L_chol;
                        M2_local_solver = L_chol';
                    elseif strcmp(precond_local_solver,'ilu')
                        [L_lu,U_lu] = ilu(A);
                        M1_local_solver = L_lu;
                        M2_local_solver = U_lu;
                    else
                        M1_local_solver = [];
                        M2_local_solver = [];
                    end
                    if strcmp(local_solver,'pcg')
                        local_sol = pcg(A,b,tol_local_solver,maxiter_local_solver,M1_local_solver,M2_local_solver);
                    elseif strcmp(local_solver,'cgs')
                        local_sol = cgs(A,b,tol_local_solver,maxiter_local_solver,M1_local_solver,M2_local_solver);
                    elseif strcmp(local_solver,'gmres')
                        local_sol = gmres(A,b,restart_local_solver,tol_local_solver,maxiter_local_solver,M1_local_solver,M2_local_solver);
                    else
                        error('Wrong solver for local problem')
                    end
                end
            else
                local_sol_old = [w_old{k};lambda_old{k}];
                b = b - A*local_sol_old;
                if strcmp(local_solver,'direct')
                    dlocal_sol = solve(A,b);
                else
                    if strcmp(precond_local_solver,'inv')
                        M1_local_solver = inv(A);
                        M2_local_solver = [];
                    elseif strcmp(precond_local_solver,'ichol')
                        L_chol = ichol(A);
                        M1_local_solver = L_chol;
                        M2_local_solver = L_chol';
                    elseif strcmp(precond_local_solver,'ilu')
                        [L_lu,U_lu] = ilu(A);
                        M1_local_solver = L_lu;
                        M2_local_solver = U_lu;
                    else
                        M1_local_solver = [];
                        M2_local_solver = [];
                    end
                    if strcmp(local_solver,'pcg')
                        dlocal_sol = pcg(A,b,tol_local_solver,maxiter_local_solver,M1_local_solver,M2_local_solver);
                    elseif strcmp(local_solver,'cgs')
                        dlocal_sol = cgs(A,b,tol_local_solver,maxiter_local_solver,M1_local_solver,M2_local_solver);
                    elseif strcmp(local_solver,'gmres')
                        dlocal_sol = gmres(A,b,restart_local_solver,tol_local_solver,maxiter_local_solver,M1_local_solver,M2_local_solver);
                    else
                        error('Wrong solver for local problem')
                    end
                end
                local_sol = local_sol_old + dlocal_sol;
            end
            
            w{k} = full(local_sol(1:nb_dof_w{k}));
            lambda{k} = full(local_sol(nb_dof_w{k}+1:end));
        else
            % With change of variable and without Lagrange multipliers
            % Local solution z{k}
            A = freematrix(S_patch_0{k},A_patch{k});
            
            b = freevector(S_patch_0{k},b_patch{k});
            if overlapping_domains
                b = b - freevector(S_patch_0{k},A_patch{k},1)*P_from_S_to_S_patch{k}*U;
            else
                b = b - freevector(S_patch_0{k},A_patch{k},1)*P_from_S_out_to_S_patch{k}*U;
            end
            
            if ~local_increment
                if strcmp(local_solver,'direct')
                    z{k} = solve(A,b);
                else
                    if strcmp(precond_local_solver,'inv')
                        M1_local_solver = inv(A);
                        M2_local_solver = [];
                    elseif strcmp(precond_local_solver,'ichol')
                        L_chol = ichol(A);
                        M1_local_solver = L_chol;
                        M2_local_solver = L_chol';
                    elseif strcmp(precond_local_solver,'ilu')
                        [L_lu,U_lu] = ilu(A);
                        M1_local_solver = L_lu;
                        M2_local_solver = U_lu;
                    else
                        M1_local_solver = [];
                        M2_local_solver = [];
                    end
                    if strcmp(local_solver,'pcg')
                        z{k} = pcg(A,b,tol_local_solver,maxiter_local_solver,M1_local_solver,M2_local_solver);
                    elseif strcmp(local_solver,'cgs')
                        z{k} = cgs(A,b,tol_local_solver,maxiter_local_solver,M1_local_solver,M2_local_solver);
                    elseif strcmp(local_solver,'gmres')
                        z{k} = gmres(A,b,restart_local_solver,tol_local_solver,maxiter_local_solver,M1_local_solver,M2_local_solver);
                    else
                        error('Wrong solver for local problem')
                    end
                end
            else
                b = b - A*z_old{k};
                if strcmp(local_solver,'direct')
                    dz = solve(A,b);
                else
                    if strcmp(precond_local_solver,'inv')
                        M1_local_solver = inv(A);
                        M2_local_solver = [];
                    elseif strcmp(precond_local_solver,'ichol')
                        L_chol = ichol(A);
                        M1_local_solver = L_chol;
                        M2_local_solver = L_chol';
                    elseif strcmp(precond_local_solver,'ilu')
                        [L_lu,U_lu] = ilu(A);
                        M1_local_solver = L_lu;
                        M2_local_solver = U_lu;
                    else
                        M1_local_solver = [];
                        M2_local_solver = [];
                    end
                    if strcmp(local_solver,'pcg')
                        dz = pcg(A,b,tol_local_solver,maxiter_local_solver,M1_local_solver,M2_local_solver);
                    elseif strcmp(local_solver,'cgs')
                        dz = cgs(A,b,tol_local_solver,maxiter_local_solver,M1_local_solver,M2_local_solver);
                    elseif strcmp(local_solver,'gmres')
                        dz = gmres(A,b,restart_local_solver,tol_local_solver,maxiter_local_solver,M1_local_solver,M2_local_solver);
                    else
                        error('Wrong solver for local problem')
                    end
                end
                z{k} = z_old{k} + dz;
            end
            
            w{k} = unfreevector(S_patch_0{k},z{k});
            if overlapping_domains
                w{k} = w{k} + P_from_S_to_S_patch{k}*U;
            else
                w{k} = w{k} + P_from_S_out_to_S_patch{k}*U;
            end
            
            % A = P_from_S_patch_to_B_patch{k}'*M_B_patch{k};
            % b = A_patch{k}*w{k} - b_patch{k};
            A = M_B_patch{k};
            b = P_from_S_patch_to_B_patch{k}*(A_patch{k}*w{k} - b_patch{k});
            lambda{k} = solve(A,b);
        end
        
        stagnation_indicator_w{k}(iter) = norm(w{k} - w_old{k})/norm(w{k});
        stagnation_indicator_lambda{k}(iter) = norm(lambda{k} - lambda_old{k})/norm(lambda{k});
        
        if local_perturbation
            w{k} = unfreevector(S_patch{k},w{k});
            w{k} = reshape(w{k},sqrt(length(w{k})),sqrt(length(w{k})));
            w{k} = svdtruncate(w{k},tol_local_perturbation);
            w{k} = w{k}(:);
            % A = P_from_S_patch_to_B_patch{k}'*M_B_patch{k};
            % b = A_patch{k}*w{k} - b_patch{k};
            A = M_B_patch{k};
            b = P_from_S_patch_to_B_patch{k}*(A_patch{k}*w{k} - b_patch{k});
            lambda{k} = solve(A,b);
        end
        
        error_indicator_w{k}(iter) = norm(w{k} - w_ref{k})/norm(w_ref{k});
        error_indicator_lambda{k}(iter) = norm(lambda{k} - lambda_ref{k})/norm(lambda_ref{k});
    end
    
    if error_indicator_U(iter) < tol
        break
    else
        if display
            fprintf('\nIteration #%3.d : stagnation = %.3e, residual error = %.3e w.r.t. U\n',iter,stagnation_indicator_U(iter),error_indicator_U(iter));
            for k=1:nbpatch
                fprintf('\n                 stagnation = %.3e, residual error = %.3e w.r.t. w for patch #%u',stagnation_indicator_w{k}(iter),error_indicator_w{k}(iter),k);
                fprintf('\n                 stagnation = %.3e, residual error = %.3e w.r.t. lambda for patch #%u\n',stagnation_indicator_lambda{k}(iter),error_indicator_lambda{k}(iter),k);
            end
        end
    end
end

if error_indicator_U(iter) < tol
    fprintf('\nAlgorithm converged at iteration #%d with residual error = %.3e w.r.t. U\n',iter,error_indicator_U(iter));
else
    fprintf('\nAlgorithm stopped at iteration #%d with residual error = %.3e w.r.t. U\n',iter,error_indicator_U(iter));
end
fprintf('\nElapsed time = %f s\n',toc);

%% Display reconstructed solution u=(U,w) at final iteration
figure('Name',['Reconstructed solution u=(U,w) at final iteration #' num2str(iter)])
% set(gcf,'Name',['Reconstructed solution u=(U,w) at final iteration #' num2str(iter)])
clf

subplot(1+nbpatch,2,1)
if overlapping_domains
    plot_sol(S_out,P_from_S_to_S_out*U,'edgecolor','none');
else
    plot_sol(S_out,U,'edgecolor','none');
end
for k=1:nbpatch
    plot_sol(S_patch{k},w{k},'edgecolor','none');
    plot(D_patch{k},'facecolor','none','edgecolor','k');
    % gtext(['\Lambda_{' num2str(k) '}'],'FontSize',fontsize);
end
% colormap('default')
colorbar
ax=axis;
cax=caxis;
set(gca,'FontSize',fontsize)
title(['Reconstructed solution u=(U,w) at final iteration #' num2str(iter) ' over domain \Omega'])

subplot(1+nbpatch,2,2)
if overlapping_domains
    plot_sol(S,U,'edgecolor','none');
else
    plot_sol(S_out,U,'edgecolor','none');
end
% colormap('default')
colorbar
axis(ax)
% caxis(cax)
set(gca,'FontSize',fontsize)
if overlapping_domains
    title(['Prolongation of global solution U at final iteration #' num2str(iter) ' over mesh T_H^{\Omega} of domain \Omega'])
else
    title(['Global solution U at final iteration #' num2str(iter) ' over mesh T_H{\Omega \\ \Lambda} of domain \Omega \\ \Lambda'])
end

for k=1:nbpatch
    subplot(1+nbpatch,2,2*k+1)
    plot_sol(S_patch{k},w{k},'edgecolor','none');
    % colormap('default')
    colorbar
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',fontsize)
    title(['Local solution w_{' num2str(k) '} at final iteration #' num2str(iter) ' over mesh T_h^{\Lambda_{' num2str(k) '}} of patch \Lambda_{' num2str(k) '}'])
    
    subplot(1+nbpatch,2,2*(k+1))
    plot_sol(B_patch{k},lambda{k},'edgecolor','interp');
    % colormap('default')
    colorbar
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',fontsize)
    title(['Lagrange multiplier \lambda_{' num2str(k) '} at final iteration #' num2str(iter) ' over boundary mesh T_h^{\Gamma_{' num2str(k) '}} of interface \Gamma_{' num2str(k) '}'])
end

% Display evolution of error indicators and stagnation indicators w.r.t. number of iterations
figure('Name','Evolution of error indicators epsilon(U;U_ref) w.r.t number of iterations') % , epsilon(w;w_ref) and epsilon(lambda;lambda_ref) w.r.t number of iterations')
% set(gcf,'Name','Evolution of error indicators epsilon(U;U_ref) w.r.t number of iterations') %, epsilon(w;w_ref) and epsilon(lambda;lambda_ref) w.r.t number of iterations')
clf
iter=1:iter;
semilogy([0,iter],[error_indicator_U_init,error_indicator_U(iter)],'-k')
hold on
if overlapping_domains
    leg = {'\epsilon_{\Omega}(U;U^{ref})'};
else
    leg = {'\epsilon_{\Omega \\ \Lambda}(U;U^{ref})'};
end
semilogy(iter,stagnation_indicator_U(iter),'--k')
if overlapping_domains
    leg = [leg, {'\epsilon_{\Omega}(U^n;U^{n-1})'}];
else
    leg = [leg, {'\epsilon_{\Omega \\ \Lambda}(U^n;U^{n-1})'}];
end
for k=1:nbpatch
    semilogy([0,iter],[error_indicator_w_init{k},error_indicator_w{k}(iter)],'Color',getfacecolor(2*k))
    leg = [leg, {['\epsilon_{\Lambda_{' num2str(k) '}}(w_{' num2str(k) '};w^{ref}_{' num2str(k) '})']}];
    semilogy(iter,stagnation_indicator_w{k}(iter),'--','Color',getfacecolor(2*k))
    leg = [leg, {['\epsilon_{\Lambda_{' num2str(k) '}}(w^n_{' num2str(k) '};w^{n-1}_{' num2str(k) '})']}];
    semilogy([0,iter],[error_indicator_lambda_init{k},error_indicator_lambda{k}(iter)],'-','Color',getfacecolor(2*k+1))
    leg = [leg, {['\epsilon_{\Gamma_{' num2str(k) '}}(\lambda_{' num2str(k) '};\lambda^{ref}_{' num2str(k) '})']}];
    semilogy(iter,stagnation_indicator_lambda{k}(iter),'--','Color',getfacecolor(2*k+1))
    leg = [leg, {['\epsilon_{\Gamma_{' num2str(k) '}}(\lambda^n_{' num2str(k) '};\lambda^{n-1}_{' num2str(k) '})']}];
end
hold off
grid on
set(gca,'FontSize',fontsize)
xlabel('Number of iterations')
ylabel('Error indicator')
legend(leg{:})
if overlapping_domains
    title('Evolution of error indicators \epsilon_{\Omega}(U;U^{ref}), \epsilon_{\Lambda}(w;w^{ref}) and \epsilon_{\Gamma}(\lambda;\lambda^{ref}) w.r.t number of iterations')
else
    title('Evolution of error indicators \epsilon_{\Omega \\ \Lambda}(U;U^{ref}), \epsilon_{\Lambda}(w;w^{ref}) and \epsilon_{\Gamma}(\lambda;\lambda^{ref}) w.r.t number of iterations')
end
