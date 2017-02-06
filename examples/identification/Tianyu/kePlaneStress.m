function [keE] = kePlaneStress(a,b,nu)

% Constitutive matrix contributions (Plane Stress)
D=1/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
% Two Gauss points in both directions
point=1/sqrt(3);
Gauss_point(1,1)=-point; Gauss_point(1,2)=-point;
Gauss_point(2,1)=point;  Gauss_point(2,2)=-point;
Gauss_point(3,1)=point;  Gauss_point(3,2)=point;
Gauss_point(4,1)=-point; Gauss_point(4,2)=point;
w=1;
% Initialize
keE = zeros(8,8); 
% node=[0,0;dx,0;dx,dy;0,dy];

for i=1:4
    % Integration point
    s=Gauss_point(i,1);t=Gauss_point(i,2);
    % Differentiated shape functions
    dNs = 1/4*[-(1-t)  (1-t) (1+t) -(1+t)];
    dNt = 1/4*[-(1-s) -(1+s) (1+s)  (1-s)];
    % Jacobian
    J = [dNs; dNt]*[-a a a -a; -b -b b b]';
    detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1);
    invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
    % Weight factor at this point
    weight = w*detJ;
    % Strain-displacement matrix
    dN_x_dN_y=invJ*[dNs; dNt];
    B1=[dN_x_dN_y(1,1),0;0,dN_x_dN_y(2,1);dN_x_dN_y(2,1),dN_x_dN_y(1,1)];
    B2=[dN_x_dN_y(1,2),0;0,dN_x_dN_y(2,2);dN_x_dN_y(2,2),dN_x_dN_y(1,2)];
    B3=[dN_x_dN_y(1,3),0;0,dN_x_dN_y(2,3);dN_x_dN_y(2,3),dN_x_dN_y(1,3)];
    B4=[dN_x_dN_y(1,4),0;0,dN_x_dN_y(2,4);dN_x_dN_y(2,4),dN_x_dN_y(1,4)];
    B=[B1,B2,B3,B4];
    % Element matrice
    keE=keE + weight*(B' * D * B);
end