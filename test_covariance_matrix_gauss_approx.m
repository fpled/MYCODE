% clc
% clearvars
clear all

%% Tensorial basis for transversely isotropic elasticity matrix 
% with rotational symmetry axis n = (1,0,0) in Kelvin-Mandel notation
syms c [5 1] matrix
syms la [6 1] matrix
% syms c1 c2 c3 c4 c5
% syms la1 la2 la3 la4 la5 la6
% syms 'c%d' [5 1] real
% syms 'la%d' [6 1] real
c = symmatrix2sym(c)
la = symmatrix2sym(la)
assume(c(1),"positive")
assume(c(2),"positive")
assumeAlso(c(1)*c(2)-c(3)^2,"positive")
assume(c(4),"positive")
assume(c(5),"positive")
assumptions_c = assumptions(c)

assume(la(1),"positive")
assume(la(2),"positive")
assumeAlso(2*sqrt(la(1)*la(2))-la(3),"positive")
assume(la(4),"positive")
assume(la(5),"positive")
% assume(la(6)<1/2)
assume(la(6)<0)
assumptions_lambda = assumptions(la)

n = [1 0 0]';
P = n*n';
Q = eye(3) - P;

% Kelvin-Mandel notation
I = eye(6);
% E1 = [P zeros(3); zeros(3) zeros(3)];
E1 = [[1 0 0; 0 0 0; 0 0 0] zeros(3); zeros(3) zeros(3)];
E2 = [[0 0 0; 0 1/2 1/2; 0 1/2 1/2] zeros(3); zeros(3) zeros(3)];
E3 = [[0 1/sqrt(2) 1/sqrt(2); 1/sqrt(2) 0 0; 1/sqrt(2) 0 0] zeros(3); zeros(3) zeros(3)];
E4 = [[0 0 0; 0 1/2 -1/2; 0 -1/2 1/2] zeros(3); zeros(3) [1 0 0; 0 0 0; 0 0 0]];
E5 = I - E1 - E2 - E4;
% E5 = [zeros(3) zeros(3); zeros(3) [0 0 0; 0 1 0; 0 0 1]];

C = simplify(c(1)*E1 + c(2)*E2 + c(3)*E3 + c(4)*E4 + c(5)*E5)
valC = eig(C)
detC = simplify(det(C))
phiC = log(detC)
PhiC = sum(la(1:5).*c) + la(6)*phiC;
% PhiC = la(1)*c(1) + la(2)*c(2) + la(3)*c(3) + la(4)*c(4) + la(5)*c(5) + la(6)*log((c(1)*c(2)-c(3)^2)*c(4)^2*c(5)^2);

g = diff(PhiC,symmatrix(c.')); g = symmatrix2sym(g)
g = gradient(PhiC,c)

h = diff(PhiC,symmatrix(c),symmatrix(c.')); h = simplify(symmatrix2sym(h))
h = simplify(hessian(PhiC,c))

c123 = c(1:3)
la123 = la(1:3)
C123 = [c(1) c(3); c(3) c(2)]
valC123 = eig(C123)
detC123 = det(C123)
phiC123 = log(detC123)
PhiC123 = sum(la(1:3).*c123) + la(6)*phiC123
% PhiC123 = la(1)*c(1) + la(2)*c(2) + la(3)*c(3) + la(6)*log(c(1)*c(2)-c(3)^2);

g123 = diff(PhiC123,symmatrix(c123.')); g123 = symmatrix2sym(g123)
g123 = gradient(PhiC123,c123)

H123 = diff(PhiC123,symmatrix(c123),symmatrix(c123.')); H123 = simplify(symmatrix2sym(H123))
H123 = simplify(hessian(PhiC123,c123))

Sigma123 = inv(H123)

% Sigma123 = -1/la(6)*[c(1)^2,    c(3)^2,    c(1)*c(3);
%                      c(3)^2,    c(2)^2,    c(2)*c(3);
%                      c(1)*c(3), c(2)*c(3), (c(1)*c(2)+c(3)^2)/2]
% H123 = simplify(inv(Sigma123))

assumptions_Sigma123 = assumptions(Sigma123)

valSigma123 = simplify(eig(Sigma123))
L123 = chol(Sigma123,'nocheck')
