function u_ana = solveThreePointBendingAna(param,coord,F,Iz,h)
% function u_ana = solveThreePointBendingAna(param,coord,F,Iz,h)

ET = param(1);
GL = param(2);
Phi = param(3);
U0 = param(4);
V0 = param(5);

ux_ana = @(x) -F*(x(:,1).^2).*x(:,2)./(4*ET*Iz)+F*(x(:,2).^3)./(12*GL*Iz)+Phi*x(:,2)+U0;
uy_ana = @(x) F*(x(:,1).^3)./(12*ET*Iz)-F*h^2*x(:,1)./(16*GL*Iz)-Phi*x(:,1)+V0;

u_ana = [ux_ana(coord) uy_ana(coord)]';
u_ana = u_ana(:);

end