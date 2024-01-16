%% Tensorial basis for isotropic elasticity matrix
syms C1 C2
E1 = [1/3*ones(3) zeros(3); zeros(3) zeros(3)];

% in Voigt-Mandel notation
E2 = [eye(3) zeros(3); zeros(3) eye(3)/2] - E1;
C = 3*C1*E1 + 2*C2*E2;
simplify(det(C))

% in Kelvin-Mandel notation
E2 = eye(6) - E1;
C = 3*C1*E1 + 2*C2*E2;
simplify(det(C))

%% Tensorial basis for transversely isotropic elasticity matrix 
syms C1 C2 C3 C4 C5
n = [0 0 1];
p = n'*n;
q = eye(3) - p;

% in Voigt notation
E1 = [[0 0 0; 0 0 0; 0 0 1] zeros(3); zeros(3) zeros(3)];
E2 = [[1/2 1/2 0; 1/2 1/2 0; 0 0 0] zeros(3); zeros(3) zeros(3)];
E3 = [[0 0 0; 0 0 0; 1/sqrt(2) 1/sqrt(2) 0] zeros(3); zeros(3) zeros(3)];
E4 = [[0 0 1/sqrt(2); 0 0 1/sqrt(2); 0 0 0] zeros(3); zeros(3) zeros(3)];
E5 = [[1/2 -1/2 0; -1/2 1/2 0; 0 0 0] zeros(3); zeros(3) [0 0 0; 0 0 0; 0 0 1/2]];
E6 = [eye(3) zeros(3); zeros(3) eye(3)/2] - E1 - E2 - E5;
% E6 = [zeros(3) zeros(3); zeros(3) [1/2 0 0; 0 1/2 0; 0 0 0]];
C = C1*E1 + C2*E2 + C3*(E3+E4) + C4*E5 + C5*E6;
simplify(det(C))

% in Kelvin-Mandel notation
E5 = [[1/2 -1/2 0; -1/2 1/2 0; 0 0 0] zeros(3); zeros(3) [0 0 0; 0 0 0; 0 0 1]];
E6 = eye(6) - E1 - E2 - E5;
% E6 = [zeros(3) zeros(3); zeros(3) [1 0 0; 0 1 0; 0 0 0]];
C = C1*E1 + C2*E2 + C3*(E3+E4) + C4*E5 + C5*E6;
simplify(det(C))

%% [Guilleminot, Soize, 2012, IJNME]
mC = [10.0735 0.5497 2.9745 0 0 0;
    0.5497 10.0735 2.9745 0 0 0;
    2.9745 2.9745 182.6657 0 0 0;
    0 0 0 14 0 0;
    0 0 0 0 14 0;
    0 0 0 0 0 9.5238];

mC1 = trace(mC'*E1)/norm(E1,'fro')^2;
mC2 = trace(mC'*E2)/norm(E2,'fro')^2;
mC3 = trace(mC'*E3)/norm(E3,'fro')^2;
% mC3 = trace(mC'*E4)/norm(E4,'fro')^2;
mC4 = trace(mC'*E5)/norm(E5,'fro')^2;
mC5 = trace(mC'*E6)/norm(E6,'fro')^2;

mC_data = [mC1 mC2 mC3 mC4 mC5];

