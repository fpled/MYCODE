function U = ThreePointsBendingIsot(x)
% function U = ThreePointsBendingIsot(x)

%% Generate mesh
L = 200e-3;
l = 1.5e-3;
h = 25e-3;
b = 25e-3;
nx = 51;
ny = 5;
k = 0;
hx = L/(nx-1);
hy = h/(ny-1);
for i=1:nx
    for j=1:ny
        k = k+1;
        Nx(k) = (i-1)*hx;
        Ny(k) = (j-1)*hy;
    end
end

List1 = find(round((Nx(:) - l)* 1e5)/1e5 <= 0);
nn1 = length(List1)/ny;
if nn1>1
    if l-hx*(nn1-1)<=hx/2
        Nx((nn1-1)*ny+1:nn1*ny) = l;
    elseif l-hx*(nn1-1)>hx/2
        Nx(nn1*ny+1:(nn1+1)*ny) = l;
    end
else
    Nx(ny+1:2*ny) = l;
end

List2 = find(round((Nx(:) - (L-l))* 1e5)/1e5 >= 0);
nn2 = length(List2)/ny;
if nn2>1
    if l-hx*(nn2-1)<=hx/2
        Nx((nx-nn2)*ny+1:(nx-nn2+1)*ny) = L-l;
    elseif l-hx*(nn1-1)>hx/2
        Nx((nx-nn2-1)*ny+1:(nx-nn2)*ny) = L-l;
    end
else
    Nx((nx-2)*ny+1:(nx-1)*ny) = L-l;
end

List3 = find(round((Nx(:) - L/2)* 1e5)/1e5 >= -l/2 & round((Nx(:) - L/2)* 1e5)/1e5...
    <= l/2);
nn3 = length(List3)/ny;
if nn3>1
    if (l-hx*(nn3-1))/2<=hx/2
        Nx(List3(1):List3(ny)) = L/2-l/2;
        Nx(List3((nn3-1)*ny+1):List3(nn3*ny)) = L/2+l/2;
    elseif (l-hx*(nn3-1))/2>hx/2
        Nx(List3(1)-ny:List3(1)-1) = L/2-l/2;
        Nx(List3(end)+1:List3(end)+ny) = L/2+l/2;
    end
else
    Nx(List3(1)-ny:List3(1)-1) = L/2-l/2;
    Nx(List3(end)+1:List3(end)+ny) = L/2+l/2;
end

Connect = delaunay(Nx,Ny);

% Nodes along Dirichlet boundaries
ListX11 = find(round((Nx(:) - l)* 1e5)/1e5 <= 0 & round((Ny(:) - 0)* 1e5)/1e5 == 0);
ListX12 = find(round((Nx(:) - (L-l))* 1e5)/1e5 >= 0  & round((Ny(:) - 0)* 1e5)/1e5 == 0);
ListX1 = [ListX11; ListX12];

% Nodes along Neumann boundaries
ListX2 = find(round((Nx(:) - L/2)* 1e5)/1e5 >= -l/2 & round((Nx(:) - L/2)* 1e5)/1e5...
    <= l/2 & round((Ny(:) - h)* 1e5)/1e5 == 0);

%% Material properties
E = x(1);
nu = x(2);
G = E/(2*(1+nu));

%% Assemble matrix
for i=1:size(Connect,1)
    MatAssemble(i,1) = Connect(i,1)*2-1;
    MatAssemble(i,2) = Connect(i,1)*2;
    MatAssemble(i,3) = Connect(i,2)*2-1;
    MatAssemble(i,4) = Connect(i,2)*2;
    MatAssemble(i,5) = Connect(i,3)*2-1;
    MatAssemble(i,6) = Connect(i,3)*2;
end

%% Construct rigidity matrix
N = length(Nx)*2;
K = sparse(N,N);
% Plane_stress = 1;

for i=1:size(Connect,1)
    n1 = Connect(i,1);
    n2 = Connect(i,2);
    n3 = Connect(i,3);
    x1 = Nx(n1);
    x2 = Nx(n2);
    x3 = Nx(n3);
    y1 = Ny(n1);
    y2 = Ny(n2);
    y3 = Ny(n3);
    X = (Nx(n1)+Nx(n2)+Nx(n3))/3;
    Y = (Ny(n1)+Ny(n2)+Ny(n3))/3;
    
    % Compute the area of triangle Ae
    ae = ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2;
    Ae = sqrt(ae*ae);
    
    % Compute matrix Be
    A = [1 x1 y1;1 x2 y2;1 x3 y3];
    N1 = A \ [1 0 0]';
    N2 = A \ [0 1 0]';
    N3 = A \ [0 0 1]';
    Be = [N1(2) 0 N2(2) 0 N3(2) 0;0 N1(3) 0 N2(3) 0 N3(3);N1(3) N1(2) N2(3) N2(2) N3(3) N3(2)];
    
    % Construct matrix S
    S = [1/E -nu/E 0
        -nu/E 1/E 0
        0 0 1/G];
    
    % Construct elementary matrices Ke
    Ke = Be'*(S\Be)*Ae;
    
    % Assemble global matrix K
    vecAssemble = MatAssemble(i,:);
    
    for ii=1:size(Ke,1)
        for jj=1:size(Ke,2)
            row = vecAssemble(ii);
            column = vecAssemble(jj);
            K(row,column) = K(row,column)+Ke(ii,jj);
        end
    end
    
end

%% Boundary conditions
% F = input('La valeur de la force: ');
F = (-600/b)*1e-9;
vecF = zeros(length(Nx)*2,1);
vecF(ListX2*2) = F/(length(ListX2));
U = zeros(length(Nx)*2,1);
List = 1:length(Nx)*2;
constrained_nodes = [(ListX11*2-1); (ListX11*2); (ListX12*2)];
% constrained_nodes = [ListX1*2];
U(constrained_nodes) = 0;
List = setdiff(List,constrained_nodes);

%% Solve
U(List)= K(List,List)\vecF(List);

end
