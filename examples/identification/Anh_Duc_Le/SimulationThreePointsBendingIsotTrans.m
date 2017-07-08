clear all
close all
format long

%% Construction le maillage
% L=input('La longueur de la poutre L= ');
% h=input('La haute de la poutre h= ');
% nx=input('Nombre de noeud suivant x (impair): ');
% ny=input('Nombre de noeud suivant y: ');
L=200e-3;
l=1.5e-3;
h=25e-3;
b=25e-3;
nx=301;
ny=31;
k=0;
hx=L/(nx-1);
hy=h/(ny-1);
for i=1:nx
    for j=1:ny
        k=k+1;       
        Nx(k)=(i-1)*hx;
        Ny(k)=(j-1)*hy;
    end
end

List1=find(round((Nx(:) - l)* 1e5)/1e5 <= 0);
nn1=length(List1)/ny;
if nn1>1
    if l-hx*(nn1-1)<=hx/2
        Nx((nn1-1)*ny+1:nn1*ny)=l;
    elseif l-hx*(nn1-1)>hx/2
        Nx(nn1*ny+1:(nn1+1)*ny)=l;
    end
else
    Nx(ny+1:2*ny)=l;
end
    
List2=find(round((Nx(:) - (L-l))* 1e5)/1e5 >= 0);

nn2=length(List2)/ny;
if nn2>1
    if l-hx*(nn2-1)<=hx/2 
        Nx((nx-nn2)*ny+1:(nx-nn2+1)*ny)=L-l;
    elseif l-hx*(nn1-1)>hx/2
        Nx((nx-nn2-1)*ny+1:(nx-nn2)*ny)=L-l;
    end
else
    Nx((nx-2)*ny+1:(nx-1)*ny)=L-l;
end
    

List3=find(round((Nx(:) - L/2)* 1e5)/1e5 >= -l/2 & round((Nx(:) - L/2)* 1e5)/1e5...
    <= l/2);
nn3=length(List3)/ny;
if nn3>1
    if (l-hx*(nn3-1))/2<=hx/2 
        Nx(List3(1):List3(ny))=L/2-l/2;
        Nx(List3((nn3-1)*ny+1):List3(nn3*ny))=L/2+l/2;
    elseif (l-hx*(nn3-1))/2>hx/2
        Nx(List3(1)-ny:List3(1)-1)=L/2-l/2;
        Nx(List3(end)+1:List3(end)+ny)=L/2+l/2;
    end
else
     Nx(List3(1)-ny:List3(1)-1)=L/2-l/2;
     Nx(List3(end)+1:List3(end)+ny)=L/2+l/2;
end

hold on
plot(Nx,Ny,'k.')
axis equal

Connect=delaunay(Nx,Ny);
hold on
triplot(Connect,Nx,Ny,'b-');%plot the triangular mesh

%) Chercher les noeuds d'appuis

ListX11=find(round((Nx(:) - l)* 1e5)/1e5 <= 0 & round((Ny(:) - 0)* 1e5)/1e5 == 0);
ListX12=find(round((Nx(:) - (L-l))* 1e5)/1e5 >= 0  & round((Ny(:) - 0)* 1e5)/1e5 == 0);
ListX1=[ListX11; ListX12];


%) Chercher le noeud appliqué la force

ListX2=find(round((Nx(:) - L/2)* 1e5)/1e5 >= -l/2 & round((Nx(:) - L/2)* 1e5)/1e5...
    <= l/2 & round((Ny(:) - h)* 1e5)/1e5 == 0);

hold on
plot(Nx(ListX11),Ny(ListX11),'r^')
plot(Nx(ListX12),Ny(ListX12),'ro')
plot(Nx(ListX2),Ny(ListX2),'gv')

%% Define les properiétés élastique du matériau

E_l=11.550;%GPa
E_t=0.500;%GPa
nu_lt=0.4;
G_lt=0.550;%GPa 

%% Construct the assembly matrix

for i=1:size(Connect,1)
    
    MatAssemble(i,1)=Connect(i,1)*2-1;
    MatAssemble(i,2)=Connect(i,1)*2;
    MatAssemble(i,3)=Connect(i,2)*2-1;
    MatAssemble(i,4)=Connect(i,2)*2;
    MatAssemble(i,5)=Connect(i,3)*2-1;
    MatAssemble(i,6)=Connect(i,3)*2;
end

%% Construction la matrice K_glob

N=length(Nx)*2;
K_glob=sparse(N,N);

for i=1:size(Connect,1)
    n1=Connect(i,1);
    n2=Connect(i,2);
    n3=Connect(i,3);
    x1 = Nx(n1);
    x2 = Nx(n2);
    x3 = Nx(n3);
    y1 = Ny(n1);
    y2 = Ny(n2);
    y3 = Ny(n3);
    X=(Nx(n1)+Nx(n2)+Nx(n3))/3;
    Y=(Ny(n1)+Ny(n2)+Ny(n3))/3;
        
    % Calculer la section élementaire Ae
    
    ae = ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2;
    Ae = sqrt(ae*ae);
    
    % Calculer la matrice Be 
    
    A = [1 x1 y1;1 x2 y2;1 x3 y3];
    N1 = A \ [1 0 0]';
    N2 = A \ [0 1 0]';
    N3 = A \ [0 0 1]';
    Be = [N1(2) 0 N2(2) 0 N3(2) 0;0 N1(3) 0 N2(3) 0 N3(3);N1(3) N1(2) N2(3) N2(2) N3(3) N3(2)];
    
    % construction la matrice C (dans le cas isotrope transverse)
    
    MatS=[1/E_l -nu_lt/E_l 0
            -nu_lt/E_l 1/E_t 0
            0 0 1/G_lt];
    MatC=MatS^(-1);

    % % Construction la matrice Ke
    
   Ke=Be'*MatC*Be*Ae; 
  
   % Construction la matrice assemblage
  
   vecAssemble = MatAssemble(i,:);
    
  for ii=1:size(Ke,1)
        for jj=1:size(Ke,2)
            ligne=vecAssemble(ii);
            colonne=vecAssemble(jj);
            K_glob(ligne,colonne)=K_glob(ligne,colonne)+Ke(ii,jj);
        end
  end           
                   
end

%% Conditions aux limites

% F=input('La valeur de la force: ');
F=(-600/b)*1e-9;
vecF=zeros(length(Nx)*2,1);
vecF(ListX2*2)=F/(length(ListX2));
Uvec=zeros(length(Nx)*2,1);
Liste=1:length(Nx)*2;
Noeud_bloques = [(ListX11*2-1);(ListX11*2); (ListX12*2)];
%Noeud_bloques = [ListX1*2];
Uvec(Noeud_bloques)=0;
Liste = setdiff(Liste,Noeud_bloques);

%% Résoudre le problème

Uvec(Liste)= K_glob(Liste,Liste)\vecF(Liste);


%% Post-treatment
% plot deformed configuration

Nx2=Nx+20*Uvec(1:2:length(Uvec))';
Ny2=Ny+20*Uvec(2:2:length(Uvec))';

figure
triplot(Connect,Nx,Ny,'b')
hold on
triplot(Connect,Nx2,Ny2,'r')
axis equal
title('\bfDéformation forme de la poutre','FontSize',17, 'FontName','Times New Roman')





