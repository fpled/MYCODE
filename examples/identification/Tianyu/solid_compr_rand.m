
clear all
close all
clc
format long
myparallel('start');

% 2D Hypothesis
% PLOption = 'PlaneStrain';
PLOption = 'PlaneStress';

% Boundary Condition Type
% LDOption = 'Dirichlet'; % Imposed displacement
LDOption = 'Neumann'; % Traction force density

% Loading Degree
DGOption = 'CST';
% DGOption = 'LIN';
% DGOption = 'QUA';

% DATA
% Element
m=4;
n=m*2;
%Domain
lx=1; ly=1;
nelx = 3; nely = 3;
dx=lx/nelx; dy=ly/nely;
% Coordinates
x=0:dx:lx;
y=0:dy:ly;
[Nx,Ny]=meshgrid(x,y);
% Quadrilateral mesh plot
quads=zeros(nelx*nely,m);
quads_c1=1:nelx*(nely+1)-1;
quads_c2=2+nelx:(nelx+1)*(nely+1)-1;
quads_c3=2+nelx+1:(nelx+1)*(nely+1);
quads_c4=2:nelx*(nely+1);
ind=find(mod(quads_c1,nely+1)==0);
quads_c1(ind)=[]; quads_c2(ind)=[];
quads_c3(ind)=[]; quads_c4(ind)=[];
quads(:,1)=quads_c1; quads(:,2)=quads_c2;
quads(:,3)=quads_c3; quads(:,4)=quads_c4;
d=quads(:,[1 2 3 4 1])';
% Degrees of freedom
ndof = 2*((nelx+1)*(nely+1));
alldofs = (1:2*(nely+1)*(nelx+1));
loaddofs = (2*(nely+1):2*(nely+1):2*(nely+1)*(nelx+1));
fixeddofs = union(1:2*(nely+1):2*(nely+1)*nelx+1,2:2*(nely+1):2*(nely+1)*nelx+2);
freedofs = setdiff(alldofs,union(fixeddofs,loaddofs));

% ASSEMBLE STIFFNESS MATRIX
% Indexing vectors
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
nodenrs = flipud(nodenrs);
edof=flipud(2*nodenrs(2:end,1:end-1)-1);
edofVec = reshape(edof,nelx*nely,1);
edofMat = repmat(edofVec,1,n)+repmat([0 1 2*nely+[2 3 4 5] +2 +3],nelx*nely,1);
iK = kron(edofMat,ones(n,1))';
jK = kron(edofMat,ones(1,n))';
% Material properties in the different elements
E = 1;
nu = 0.3;
kapa = E/(3*(1-2*nu));
eE = E*ones(nelx,nely);   eEVec=eE(:).'; 
emu = E/(2*(1+nu))*ones(nelx,nely);   emuVec=emu(:).';                                            % 2nd lame paramter = shear modulus (GPa - KN/mm2)
elambda = (kapa-2/3*E/(2*(1+nu)))*ones(nelx,nely);   elambdaVec=elambda(:).';       % 1st lame paramter (GPa - KN/mm2) 

% Imposed loading conditions
u = zeros(2*(nely+1)*(nelx+1),1);
switch LDOption
    case 'Dirichlet'
        [u] = DirichletSurfLoad(Nx,loaddofs,lx,u,DGOption);
        listdofs=freedofs;
        F=zeros(ndof,1);
        F(listdofs)=-(K(freedofs,loaddofs)*u(loaddofs));
        Fd=F(listdofs);
    case 'Neumann'
        [F] = NeumannSurfLoad(Nx,Ny,lx,ly,dx,ndof,DGOption);
        listdofs=setdiff(alldofs,fixeddofs);
        Fd=F(listdofs);
end

% Corresponding stiffness matrix entries
switch PLOption
    case 'PlaneStrain'
        [keLambda, keMu] = kePlaneStrain(dx/2,dy/2);
    case 'PlaneStress'
        [keE] = kePlaneStress(dx/2,dy/2,nu);
end

delta=0.31;
if delta < 0
    errordlg('valeur de delta < 0');
    return
elseif delta > sqrt((n+1)/(n+5))
    errordlg('valeur de delta > ((n+1)(n+5)^-1)^0.5')
    return
end
meanG = 0;
sigmaG = delta*(n+1)^(-0.5);   % standard_deviation

NomSimPower=1:0.2:5;
% NomSimPower=1;
NomSim=round(10.^NomSimPower);
uSol=zeros((nelx+1)*(nely+1)*2,max(NomSim));
% rng('default');    % r¨¦initialisation des NS premiers MC simulations
tic
parfor ger=1:max(NomSim)
    % M¨¦thode 4 : cas g¨¦n¨¦rale (lambda not an integer)
    L_hd = sigmaG*randn(n*(n-1)/2,1) + meanG;
    
    L_rand = triu(ones(n,n),1);
    L_rand(L_rand==1)=L_hd;
    
    A=(n+1)/(2*delta^2)+(1-(1:n)')/2;
    B=1;
    V_j = gamrnd(A,B);
    L_d = sigmaG*(2*V_j).^0.5;
    
    L_rand_d = diag(L_d);
    L_rand = L_rand + L_rand_d;
    
    G_rand=L_rand'*L_rand;
    
    switch PLOption
        case 'PlaneStrain'
            keLambda_rand = keLambda'*G_rand*keLambda;
            keMu_rand = keMu'*G_rand*keMu;
            % sK = keLambda*(kapa-2/3*E/(2*(1+nu)))+keMu*E/(2*(1+nu));
            sK = keLambda_rand(:)*elambdaVec + keMu_rand(:)*emuVec;
        case 'PlaneStress'
            keE_rand = keE'*G_rand*keE;
            sK = keE_rand(:)*eEVec;
    end

    K = sparse(iKVec, jKVec, sK(:), ndof, ndof);
    
    % SOLVE LINEAR SYSTEM
    % Staggered solution procedure
    v = K(listdofs,listdofs)\Fd;
    uSol(listdofs,ger)=v;
end
switch LDOption
    case 'Dirichlet'
        uSol(loaddofs,:) = repmat(u(loaddofs),1,size(uSol,2));
end
toc


%Post-traitement
u_moy=zeros((nelx+1)*(nely+1)*2,length(NomSim));
u_m2=zeros((nelx+1)*(nely+1)*2,length(NomSim));
u_var=zeros((nelx+1)*(nely+1)*2,length(NomSim));
u_moy_diff=zeros((nelx+1)*(nely+1)*2,length(NomSim)-1);
u_var_diff=zeros((nelx+1)*(nely+1)*2,length(NomSim)-1);
for sim=1:length(NomSim)
    NS=NomSim(sim);
    
    %     rng('default');    % r¨¦initialisation des NS premiers MC simulations
    
    u_moy(:,sim)=mean(uSol(:,1:NS),2);
    u_m2(:,sim)=mean(uSol(:,1:NS).^2,2);
    
    if sim>1
        u_moy_diff(:,sim-1)=u_moy(:,sim)-u_moy(:,sim-1);
    end
end
figure
est_u_moy=mean(u_moy_diff);
% plot(1:length(est_u_moy),est_u_moy)
semilogx(NomSim(2:end),est_u_moy)
tx = xlabel('Nombre Simulation (log)'); ty = ylabel('moyenne de E(u)');
set(tx,'Fontsize',10); set(ty,'Fontsize',10);
tt=title('Etude Convergence : esperance') ;
set(tt,'FontName','Bell MT','Fontsize',12); set(gca,'FontName','Times','FontSize',10) ;

for sim=1:length(NomSim)
    u_var(:,sim)=u_m2(:,sim)-u_moy(:,length(NomSim)).^2;
    
    if sim>1
        u_var_diff(:,sim-1)=u_var(:,sim)-u_var(:,sim-1);
    end
end
figure
est_u_var=mean(u_var_diff);
% plot(1:length(est_u_var),est_u_var)
semilogx(NomSim(2:end),est_u_var)
tx = xlabel('Nombre Simulation (log)'); ty = ylabel('moyenne de var(u)');
set(tx,'Fontsize',10); set(ty,'Fontsize',10);
tt=title('Etude Convergence : variance') ;
set(tt,'FontName','Bell MT','Fontsize',12); set(gca,'FontName','Times','FontSize',10) ;

delta_u=zeros(length(NomSim),1);
for sim=1:length(NomSim)
    NS=NomSim(sim);
    
    uumoy2=zeros(NS,1);
    for k=1:NS
        uumoy2(k,1) = (norm(uSol(:,k) - u_moy(:,length(NomSim)),'fro'))^2;
    end
    EM_uumoy2 = mean(uumoy2);
    delta_u(sim,1)=(EM_uumoy2/(norm(u_moy(:,length(NomSim)),'fro'))^2)^0.5;
end
figure
semilogx(NomSim,delta_u)
tx = xlabel('Nombre Simulation (log)'); ty = ylabel('\delta_u');
set(tx,'Fontsize',10); set(ty,'Fontsize',10);
tt=title('Etude Convergence : coefficient de variation') ;
set(tt,'FontName','Bell MT','Fontsize',12); set(gca,'FontName','Times','FontSize',10) ;

% save test
% save ('test.mat','u_diff_moy');

% % Final displacement meshing plot
% figure
% hold on
% plot(Nx(d), Ny(d),'b');
% axis equal; axis tight; axis off;
% ampl = lx/max(abs(u))/5;
% Nx_sol=Nx(:)+u(1:2:2*length(Nx(:)))*ampl;
% Ny_sol=Ny(:)+u(2:2:2*length(Nx(:)))*ampl;
% plot(Nx_sol(d), Ny_sol(d),'r');
% tt=title(['Final Displacement Profile']);
% set(tt,'FontName','Bell MT','Fontsize',12);
%
% % print('-r300','-dpng',[directory,'Final Displacement Profile'])

myparallel('stop');