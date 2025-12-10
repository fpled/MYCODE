% clc
% clearvars
clear all

%% Tensorial basis for isotropic elasticity matrix
% syms c1 c2 positive
syms 'c%d' [2 1] positive
assumptions_c = assumptions([c1 c2])

E1 = [1/3*ones(3) zeros(3); zeros(3) zeros(3)];

% Voigt-Mandel notation
E2 = [eye(3) zeros(3); zeros(3) eye(3)/2] - E1;
C = 3*c1*E1 + 2*c2*E2
valC = eig(C)
detC = det(C)

% Kelvin-Mandel notation
E2 = eye(6) - E1;
C = 3*c1*E1 + 2*c2*E2
valC = eig(C)
detC = det(C)

%% Tensorial basis for transversely isotropic elasticity matrix 
% with rotational symmetry axis n = (0,0,1)
% syms c1 c2 c3 c4 c5 positive
syms 'c%d' [5 1] positive
assume(c3,"clear")
assumeAlso(c1*c2-c3^2,"positive")
assumptions_C = assumptions([c1 c2 c3 c4 c5])

n = [0 0 1]';
P = n*n';
Q = eye(3) - P;

% Voigt notation
I = [eye(3) zeros(3); zeros(3) eye(3)/2];
% E1 = [P zeros(3); zeros(3) zeros(3)];
E1 = [[0 0 0; 0 0 0; 0 0 1] zeros(3); zeros(3) zeros(3)];
E2 = [[1/2 1/2 0; 1/2 1/2 0; 0 0 0] zeros(3); zeros(3) zeros(3)];
E3 = [[0 0 1/sqrt(2); 0 0 1/sqrt(2); 1/sqrt(2) 1/sqrt(2) 0] zeros(3); zeros(3) zeros(3)];
E4 = [[1/2 -1/2 0; -1/2 1/2 0; 0 0 0] zeros(3); zeros(3) [0 0 0; 0 0 0; 0 0 1/2]];
E5 = I - E1 - E2 - E4;
% E5 = [zeros(3) zeros(3); zeros(3) [1/2 0 0; 0 1/2 0; 0 0 0]];
C = c1*E1 + c2*E2 + c3*E3 + c4*E4 + c5*E5
valC = eig(C)
detC = simplify(det(C))
phiC = log(detC)

% Kelvin-Mandel notation
I = eye(6);
E4 = [[1/2 -1/2 0; -1/2 1/2 0; 0 0 0] zeros(3); zeros(3) [0 0 0; 0 0 0; 0 0 1]];
E5 = I - E1 - E2 - E4;
% E5 = [zeros(3) zeros(3); zeros(3) [1 0 0; 0 1 0; 0 0 0]];
C = c1*E1 + c2*E2 + c3*E3 + c4*E4 + c5*E5
valC = eig(C)
detC = simplify(det(C))
phiC = log(detC)

% compact symbolic representation
c123 = [c1; c2; c3]
C123 = [c1 c3; c3 c2]
valC123 = eig(C123)
detC123 = simplify(det(C123))
phiC123 = log(detC123)

%% [Guilleminot, Soize, 2012, IJNME]
% transversely isotropic mean elasticity matrix 
% with rotational symmetry axis n = (0,0,1)
% Kelvin-Mandel notation
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

mC_data = [mC1 mC2 mC3 mC4 mC5]'

%% Relations between coefficients (c1,...,c5) and engineering elastic moduli (EL,ET,GL,NUL,NUT) 
% for transversely isotropic elasticity matrix with rotational symmetry axis n = (1,0,0)

kT = c2/2;
NUL = (c3./kT)/(2*sqrt(2));
% NUL = c3./c2/sqrt(2);
EL = c1 - 4*(NUL.^2).*kT;
% EL = c1-(c3.^2)./c2;
GT = c4/2;
GL = c5/2;
ET = 4./(1./kT+1./GT+4*(NUL.^2)./EL); ET = simplify(ET);
% ET = 2*(c1.*c2-c3.^2).*c4./(c1.*c2-c3.^2+c1.*c4);
% ET = 2./(c1./(c1.*c2-c3.^2)+1./c4);
NUT = (ET./GT)/2-1; NUT = simplify(NUT);
% NUT = (c1.*c2-c3.^2-c1.*c4)./(c1.*c2-c3.^2+c1.*c4);
% NUT = (c1.*(c2-c4)-c3.^2)./(c1.*(c2+c4)-c3.^2);

syms EL ET GL positive
syms NUL NUT real
assume(NUL>-1 & NUL<0.5)
assume(NUT>-1 & NUT<0.5)
assumptions_moduli = assumptions([EL ET GL NUL NUT])

GT = ET./(2*(1+NUT));
kT = (EL.*ET)./(2*(1-NUT).*EL-4*ET.*(NUL).^2);
c1 = EL + 4*(NUL.^2).*kT;
c2 = 2*kT;
c3 = 2*sqrt(2)*kT.*NUL;
c4 = 2*GT;
c5 = 2*GL;
