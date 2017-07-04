function u_FEMU_in = ThreePointsBendingFEMU(ET,GL,x,Mesh,coorx,coory,u_exp_edge)
% function u_FEMU_in = ThreePointsBendingFEMU(ET,GL,x,Mesh,coorx,coory,u_exp_edge)

EL = x(1);
nuL = x(2);

S = [ 1/EL   -nuL/EL 0
     -nuL/EL 1/ET    0
      0      0     1/GL ];

K=zeros(2*Mesh.NNodeTot);
for i=1:Mesh.Nelem
    x123 = coorx( Mesh.TRI(i,:) );
    y123 = coory( Mesh.TRI(i,:) );
    
    A = [1 x123(1) y123(1)
         1 x123(2) y123(2)
         1 x123(3) y123(3)]\[1 0 0
                             0 1 0
                             0 0 1];
    
    B = [A(2,1) 0      A(2,2) 0      A(2,3) 0
         0      A(3,1) 0      A(3,2) 0      A(3,3)
         A(3,1) A(2,1) A(3,2) A(2,2) A(3,3) A(2,3)];
    
    Ae = 1/2*( (x123(2)-x123(1))*(y123(3)-y123(1))...
              -(y123(2)-y123(1))*(x123(3)-x123(1)) );
    
    Ke = B'*(S\(B*Ae));
    
    cnt = zeros(1,6);
    cnt(1:2:end)=2*Mesh.TRI(i,:)-1;
    cnt(2:2:end)=2*Mesh.TRI(i,:);
    
    K(cnt,cnt) = K(cnt,cnt) + Ke;
    
end

K1=K(1:2*(Mesh.NNodeTot-Mesh.NNodeEdge),...
     1:2*(Mesh.NNodeTot-Mesh.NNodeEdge));
K3=K(1:2*(Mesh.NNodeTot-Mesh.NNodeEdge),...
     2*(Mesh.NNodeTot-Mesh.NNodeEdge)+1:end);

u_FEMU_in = -K1\(K3*u_exp_edge);

end
