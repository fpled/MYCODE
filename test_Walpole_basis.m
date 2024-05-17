%% Tensorial basis for isotropic elasticity matrix
syms C1 C2
E1 = [1/3*ones(3) zeros(3); zeros(3) zeros(3)];

% in Voigt-Mandel notation
E2 = [eye(3) zeros(3); zeros(3) eye(3)/2] - E1;
C = 3*C1*E1 + 2*C2*E2;
C
eig(C)
det(C)

% in Kelvin-Mandel notation
E2 = eye(6) - E1;
C = 3*C1*E1 + 2*C2*E2;
C
eig(C)
det(C)

%% Tensorial basis for transversely isotropic elasticity matrix 
% with rotational symmetry axis n = (0,0,1)
syms C1 C2 C3 C4 C5
n = [0 0 1]';
P = n*n';
Q = eye(3) - P;

% in Voigt notation
I = [eye(3) zeros(3); zeros(3) eye(3)/2];
% E1 = [P zeros(3); zeros(3) zeros(3)];
E1 = [[0 0 0; 0 0 0; 0 0 1] zeros(3); zeros(3) zeros(3)];
E2 = [[1/2 1/2 0; 1/2 1/2 0; 0 0 0] zeros(3); zeros(3) zeros(3)];
E3 = [[0 0 1/sqrt(2); 0 0 1/sqrt(2); 1/sqrt(2) 1/sqrt(2) 0] zeros(3); zeros(3) zeros(3)];
E4 = [[1/2 -1/2 0; -1/2 1/2 0; 0 0 0] zeros(3); zeros(3) [0 0 0; 0 0 0; 0 0 1/2]];
E5 = I - E1 - E2 - E4;
% E5 = [zeros(3) zeros(3); zeros(3) [1/2 0 0; 0 1/2 0; 0 0 0]];
C = C1*E1 + C2*E2 + C3*E3 + C4*E4 + C5*E5;
C
eig(C)
simplify(det(C))

% in Kelvin-Mandel notation
I = eye(6);
E4 = [[1/2 -1/2 0; -1/2 1/2 0; 0 0 0] zeros(3); zeros(3) [0 0 0; 0 0 0; 0 0 1]];
E5 = I - E1 - E2 - E4;
% E5 = [zeros(3) zeros(3); zeros(3) [1 0 0; 0 1 0; 0 0 0]];
C = C1*E1 + C2*E2 + C3*E3 + C4*E4 + C5*E5;
C
eig(C)
simplify(det(C))

% in compact symbolic representation
C123 = [C1 C3; C3 C2];
C123
eig(C123)
simplify(det(C123))

%% [Guilleminot, Soize, 2012, IJNME]
% transversely isotropic mean elasticity matrix 
% with rotational symmetry axis n = (0,0,1)
% in Kelvin-Mandel notation
mC = [10.0735 0.5497 2.9745 0 0 0;
    0.5497 10.0735 2.9745 0 0 0;
    2.9745 2.9745 182.6657 0 0 0;
    0 0 0 14 0 0;
    0 0 0 0 14 0;
    0 0 0 0 0 9.5238];

mC1 = trace(mC'*E1)/norm(E1,'fro')^2; % mC1 = trace(mC'*E1)/trace(E1'*E1);
mC2 = trace(mC'*E2)/norm(E2,'fro')^2; % mC2 = trace(mC'*E2)/trace(E2'*E2);
mC3 = trace(mC'*E3)/norm(E3,'fro')^2; % mC3 = trace(mC'*E3)/trace(E3'*E3);
mC4 = trace(mC'*E4)/norm(E4,'fro')^2; % mC4 = trace(mC'*E4)/trace(E4'*E4);
mC5 = trace(mC'*E5)/norm(E5,'fro')^2; % mC5 = trace(mC'*E5)/trace(E5'*E5);

mC_data = [mC1 mC2 mC3 mC4 mC5];

%% Relations between coefficients (C1,...,C5) and engineering elastic moduli (ET,ET,GL,NUL,NUT) 
% for transversely isotropic elasticity matrix with rotational symmetry axis n = (1,0,0)
n = [1 0 0]';

syms C1 C2 C3 C4 C5
kT = C2/2;
NUL = (C3./kT)/(2*sqrt(2));
% NUL = C3./C2/sqrt(2);
EL = C1 - 4*(NUL.^2).*kT;
% EL = C1-(C3.^2)./C2;
GT = C4/2;
GL = C5/2;
ET = 4./(1./kT+1./GT+4*(NUL.^2)./EL); ET = simplify(ET);
% ET = 2*(C1.*C2-C3.^2).*C4./(C1.*C2-C3.^2+C1.*C4);
% ET = 2./(C1./(C1.*C2-C3.^2)+1./C4);
NUT = (ET./GT)/2-1; NUT = simplify(NUT);
% NUT = (C1.*C2-C3.^2-C1.*C4)./(C1.*C2-C3.^2+C1.*C4);
% NUT = (C1.*(C2-C4)-C3.^2)./(C1.*(C2+C4)-C3.^2);

syms EL ET GL NUL NUT
GT = ET./(2*(1+NUT));
kT = (EL.*ET)./(2*(1-NUT).*EL-4*ET.*(NUL).^2);
C1 = EL + 4*(NUL.^2).*kT;
C2 = 2*kT;
C3 = 2*sqrt(2)*kT.*NUL;
C4 = 2*GT;
C5 = 2*GL;
