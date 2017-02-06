
clear all
close all
clc
format long

% 2D Hypothesis
PLOption = 'PlaneStrain';
% PLOption = 'PlaneStress';

% Boundary Condition Type
% LDOption = 'Dirichlet'; % Imposed displacement
LDOption = 'Neumann'; % Traction force density

% Loading Degree
% DGOption = 'CST';
% DGOption = 'LIN';
DGOption = 'QUA';

tic
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
eE = E*ones(nelx,nely);
emu = E/(2*(1+nu))*ones(nelx,nely); % 2nd lame paramter = shear modulus (GPa - KN/mm2)
elambda = (kapa-2/3*E/(2*(1+nu)))*ones(nelx,nely); % 1st lame paramter (GPa - KN/mm2)
% Corresponding stiffness matrix entries
switch PLOption
    case 'PlaneStrain'
        [keLambda, keMu] = kePlaneStrain(dx/2,dy/2);
        % sK = keLambda*(kapa-2/3*E/(2*(1+nu)))+keMu*E/(2*(1+nu));
        sK = keLambda(:)*elambda(:)' + keMu(:)*emu(:)';
    case 'PlaneStress'
        [keE] = kePlaneStress(dx/2,dy/2,nu);
        sK = keE(:)*eE(:)';
end
% Assemblage
K = sparse(iK(:), jK(:), sK(:), ndof, ndof);

% Imposed loading conditions
u = zeros(2*(nely+1)*(nelx+1),1);
switch LDOption
    case 'Dirichlet'
        [u] = DirichletSurfLoad(Nx,loaddofs,lx,u,DGOption);
        listdofs=freedofs;
        F=zeros(ndof,1);
        F(listdofs)=-(K(freedofs,loaddofs)*u(loaddofs));
    case 'Neumann'
        [F] = NeumannSurfLoad(Nx,Ny,lx,ly,dx,ndof,DGOption);
        listdofs=setdiff(alldofs,fixeddofs);
end

% SOLVE LINEAR SYSTEM
% Staggered solution procedure
u(listdofs) = K(listdofs,listdofs)\F(listdofs);
toc

% save test
% save ('test.mat','u_diff_moy');

% % Final displacement meshing plot
figure
hold on
plot(Nx(d), Ny(d),'b');
ampl = lx/max(abs(u))/5;
Nx_sol=Nx(:)+u(1:2:2*length(Nx(:)))*ampl;
Ny_sol=Ny(:)+u(2:2:2*length(Nx(:)))*ampl;
plot(Nx_sol(d), Ny_sol(d),'r');
axis equal; axis off;axis tight; 
tt=title(['Final Displacement Profile']);
set(tt,'FontName','Bell MT','Fontsize',12);
% 
% % print('-r300','-dpng',[directory,'Final Displacement Profile'])

