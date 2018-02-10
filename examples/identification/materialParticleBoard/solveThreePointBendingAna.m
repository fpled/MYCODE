function u = solveThreePointBendingAna(param,coord,F,Iz,h)
% function u = solveThreePointBendingAna(param,coord,F,Iz,h)

ET = param(1);
GL = param(2);
Phi = param(3);
U0 = param(4);
V0 = param(5);

ux = @(x) -F*(x(:,1).^2).*x(:,2)./(4*ET*Iz)+F*(x(:,2).^3)./(12*GL*Iz)+Phi*x(:,2)+U0;
uy = @(x) F*(x(:,1).^3)./(12*ET*Iz)-F*h^2*x(:,1)./(16*GL*Iz)-Phi*x(:,1)+V0;

u = [ux(coord) uy(coord)]';
u = u(:);

end