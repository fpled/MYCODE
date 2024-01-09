function u = solveThreePointBendingAna(param,coord,F,Iz,h)
% function u = solveThreePointBendingAna(param,coord,F,Iz,h)

ET = param(1); % transverse Young modulus [MPa]
GL = param(2); % longitudinal shear modulus [MPa]
R0 = param(3); % rigid body rotation around z direction [rad]
U0 = param(4); % rigid body displacement along x direction [mm]
V0 = param(5); % rigid body displacement along y direction [mm]

ux = @(x) -F*(x(:,1).^2).*x(:,2)./(4*ET*Iz)+F*(x(:,2).^3)./(12*GL*Iz)+R0*x(:,2)+U0; % [mm]
uy = @(x) F*(x(:,1).^3)./(12*ET*Iz)-F*h^2*x(:,1)./(16*GL*Iz)-R0*x(:,1)+V0; % [mm]

u = [ux(coord) uy(coord)]';
u = u(:); % [mm]

end