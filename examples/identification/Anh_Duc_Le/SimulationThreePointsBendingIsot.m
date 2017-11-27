% clc
clearvars
close all

%% Construction le maillage
% L=input('La longueur de la poutre L= ');
% h=input('La haute de la poutre h= ');
% nx=input('Nombre de noeud suivant x (impair): ');
% ny=input('Nombre de noeud suivant y: ');
L=200e-3;
l=1.5e-3;
h=25e-3;
b=25e-3;
nx=201;
ny=21;
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


%) Chercher le noeud appliqu� la force

ListX2=find(round((Nx(:) - L/2)* 1e5)/1e5 >= -l/2 & round((Nx(:) - L/2)* 1e5)/1e5...
    <= l/2 & round((Ny(:) - h)* 1e5)/1e5 == 0);

hold on
plot(Nx(ListX11),Ny(ListX11),'r^')
plot(Nx(ListX12),Ny(ListX12),'ro')
plot(Nx(ListX2),Ny(ListX2),'gv')

%% Define les properi�t�s �lastique du mat�riau

E=11.550; %GPa
nu=0.4; 
G=E/(2*(1+nu));%GPa

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
% Contrainte_plan=1;

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
        
    % Calculer la section �lementaire Ae 
    
    ae = ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2;
    Ae = sqrt(ae*ae);
    
    % Calculer la matrice Be 
    
    A = [1 x1 y1;1 x2 y2;1 x3 y3];
    N1 = A \ [1 0 0]';
    N2 = A \ [0 1 0]';
    N3 = A \ [0 0 1]';
    Be = [N1(2) 0 N2(2) 0 N3(2) 0;0 N1(3) 0 N2(3) 0 N3(3);N1(3) N1(2) N2(3) N2(2) N3(3) N3(2)];
    
    % construction la matrice C (dans le cas isotrope)
    
    MatS=[1/E -nu/E 0
            -nu/E 1/E 0
            0 0 1/G];
    MatC=MatS^(-1);
   
    % Construction la matrice Ke
    
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
Noeud_bloques = [(ListX11*2-1); (ListX11*2); (ListX12*2)];
%Noeud_bloques = [ListX1*2];
Uvec(Noeud_bloques)=0;
Liste = setdiff(Liste,Noeud_bloques);

%% R�soudre le probl�me

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
title('\bfD�formation forme de la poutre','FontSize',17, 'FontName','latex')


%% Compte U_ef et U_theo
List1=find(round((Nx(:) - L/2)* 1e5)/1e5 == 0);
Uvec1=Uvec(List1*2);
U_max=mean(Uvec1)
I=b*h^3/12;
U_max_theo=-(600e-9*L^3)/(48*E*I)

% calculer la contrainte et la d�formation
% Chercher la quantit� d'�l�ments communs � chaque noeud
[a1,b1]=hist(Connect(:,1),unique(Connect));
[a2,b2]=hist(Connect(:,2),unique(Connect));
[a3,b3]=hist(Connect(:,3),unique(Connect));
num_element_per_node = (a1+a2+a3)';
% Trouver le nombre d'�lements communs � chaque noeud 
max = max(num_element_per_node);
node_to_element = zeros(length(Nx),max);
for i = 1 : length(Nx)
    num = num_element_per_node(i);
    node_to_element(i,1:num) = find(Connect(:,1)==i|Connect(:,2)==i|Connect(:,3)==i);
end
%calculer la contrainte et la d�formation
Mat_sigma = zeros(length(Nx),3);
Mat_epsilon = zeros(length(Nx),3);
for i = 1 : length(Nx)
    sigma = [0 0 0];
    epsilon = [0 0 0];
    V = 0;
    num_element = num_element_per_node(i);
    element_node = node_to_element(i,1:num_element);
    for j = 1 : num_element
        e = element_node(j);% nombre d'�lements au noeud i
        %3 noeudes d'�lements e
        n1 = Connect(e,1);
        n2 = Connect(e,2);
        n3 = Connect(e,3);
        TriNodes=[n1 n2 n3];
      
        x1 = Nx(n1);
        x2 = Nx(n2);
        x3 = Nx(n3);
        y1 = Ny(n1);
        y2 = Ny(n2);
        y3 = Ny(n3);
        X=(Nx(n1)+Nx(n2)+Nx(n3))/3;
        Y=(Ny(n1)+Ny(n2)+Ny(n3))/3;
        
    % calculer la section �lementaire Ae 
    
        ae = ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2;
        Ae = sqrt(ae*ae);
    
    % Calculer la matrice Be
    
        A = [1 x1 y1;1 x2 y2;1 x3 y3];
        N1 = A \ [1 0 0]';
        N2 = A \ [0 1 0]';
        N3 = A \ [0 0 1]';
        Be = [N1(2) 0 N2(2) 0 N3(2) 0;0 N1(3) 0 N2(3) 0 N3(3);N1(3) N1(2) N2(3) N2(2) N3(3) N3(2)];
    
        % Vecteur assemblage
        assemblage=MatAssemble(e,:);
        % Calculer la matrice d'�lasticite C_e
        MatS=[1/E -nu/E 0
            -nu/E 1/E 0
            0 0 1/G];
        MatC=MatS^(-1);
        % Calculer le champ de d�placement Ui
        U_i = Uvec(assemblage);
        sigma_terme = MatC*Be*U_i*Ae;
        epsilon_terme = Be*U_i*Ae;
        sigma = sigma + sigma_terme';
        epsilon = epsilon + epsilon_terme';
        V=V+Ae;
    end
    % Calculer la contrainte et la d�formation
    % Enregistrer des donn�es 
    sigma_moyen=(1e9/V)*sigma;
    epsilon_moyen=(1e9/V)*epsilon;
    Mat_sigma(i,:) = sigma_moyen(:);
    Mat_epsilon(i,:) = epsilon_moyen(:);
end
sigma_11 = Mat_sigma(:,1);
sigma_22 = Mat_sigma(:,2);
sigma_12 = Mat_sigma(:,3);
epsilon_11 = Mat_epsilon(:,1);
epsilon_22 = Mat_epsilon(:,2);
epsilon_12 = Mat_epsilon(:,3)/2;

Tracer2Dscalaire(sigma_11,Nx,Ny);
title('\sigma_{xx}');
Tracer2Dscalaire(sigma_22,Nx,Ny);
title('\sigma_{yy}');
Tracer2Dscalaire(sigma_12,Nx,Ny);
title('\sigma_{xy}');
Tracer2Dscalaire(epsilon_11,Nx,Ny);
title('\epsilon_{xx}');
Tracer2Dscalaire(epsilon_22,Nx,Ny);
title('\epsilon_{yy}');
Tracer2Dscalaire(epsilon_12,Nx,Ny);
title('\epsilon_{xy}');
%% Solution analytique
P=-600;
I=b*h^3/12;
for i=1:length(Nx)
    if Nx(i)<= (L/2)
        sigma2_11(i)=(P/(2*I))*Nx(i)*(Ny(i)-h/2);
    else
        sigma2_11(i)=(P/(2*I))*(L-Nx(i))*(Ny(i)-h/2);
    end
end
fontsize = 16;
linewidth = 1;
markersize = 36;
interpreter = 'latex';

Tracer2Dscalaire(sigma2_11,Nx,Ny);
set(gca,'FontSize',fontsize)
title('$\sigma_{xx} theorique$','Interpreter',interpreter,'FontWeight','bold');

Tracer2Dscalaire(sigma_11,Nx,Ny);
set(gca,'FontSize',fontsize)
title('$\sigma_{xx} numerique$','Interpreter',interpreter,'FontWeight','bold');

err=abs((sigma_11'-sigma2_11)./sigma2_11);
err(find(err==Inf))=0;
err(isnan(err)==1)=0;
Tracer2Dscalaire(err,Nx,Ny);
set(gca,'FontSize',fontsize)
title('$erreur \quad relative$','Interpreter',interpreter,'FontWeight','bold');