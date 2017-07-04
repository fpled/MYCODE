function u_ana = ThreePointsBendingAna(x,coorx,coory,F,k,Iz,Mesh,h)
% function u_ana = ThreePointsBendingAna(x,coorx,coory,F,k,Iz,Mesh,h)

ET = x(1);
GL = x(2);
Phi = x(3);
U0 = x(4);
V0 = x(5);

ux_ana = zeros(Mesh.NNodeTot,1);
uy_ana = zeros(Mesh.NNodeTot,1);

for i=1:Mesh.NNodeTot
    xi = coorx(i);
    yi = coory(i);
    
    ux_ana(i) = -F(k)*xi^2*yi/(4*ET*Iz)+F(k)*yi^3/(12*GL*Iz)+Phi*yi+U0;
    uy_ana(i) = F(k)*xi^3/(12*ET*Iz)-F(k)*h^2*xi/(16*GL*Iz)-Phi*xi+V0;
end

u_ana = [ux_ana; uy_ana];

end