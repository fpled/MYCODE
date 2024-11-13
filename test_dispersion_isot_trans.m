% Axial Transmission technique
% Mean transversely isotropic elasticity matrix C in Kelvin-Mandel notation
% Coefficient of variation delta_C of C as a fonction of dispersion parameter delta_S
% [Desceliers, Soize, Grimal, Talmant, Naili, 2009, JASA]
% [Desceliers, Soize, Naili, Haiat, 2012, MSSP] (containing an error: delta_S = 0.1029 and not delta_C)
EL = 17.717e9;
NUL = 0.3816;
GL = 4.7950e9;
ET = 9.8254e9;
NUT = 0.4495;
GT = ET/(2*(1+NUT));
kT = EL*ET/(2*(1-NUT)*EL-4*NUL^2*ET);

C11 = EL+4*NUL^2*kT; % C11 = EL^2*(1-NUT)/(EL-EL*NUT-2*ET*NUL^2);
C22 = kT+GT; % C22 = ET*(EL-ET*NUL^2)/((1+NUT)*(EL-EL*NUT-2*ET*NUL^2));
C12 = 2*kT*NUL; % C12 = ET*EL*NUL/(EL-EL*NUT-2*ET*NUL^2);
C23 = kT-GT; % C23 = ET*(EL*NUT+ET*NUL^2)/((1+NUT)*(EL-EL*NUT-2*ET*NUL^2));
C44 = 2*GT;
C55 = 2*GL;

C=[C11 C12 C12 0   0   0
   C12 C22 C23 0   0   0
   C12 C23 C22 0   0   0
   0   0   0   C44 0   0
   0   0   0   0   C55 0
   0   0   0   0   0   C55];
% Cl = 1e-6*C;
Cl = 0;
delta_S = 0.1029;
delta_C = delta_S/sqrt(7)*sqrt(1+trace(C-Cl)^2/trace((C-Cl)^2));