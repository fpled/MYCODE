% clc
clearvars

% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
% Elasticity matrix in reference material coordinate system [Pa]
% Orthotropic symmetry
Cref = 1e9*[65 20 0;
            20 260 0;
            0 0 30];
% Isotropic symmetry
% lambda = 121.15e9;
% mu = 80.77e9;
% Cref = [lambda+2*mu lambda, 0;
%         lambda lambda+2*mu 0;
%         0 0 mu];

ang = 30; % clockwise material orientation angle around z-axis [deg]
theta = deg2rad(ang); % clockwise material orientation angle around z-axis [rad]
c = cos(theta);
s = sin(theta);
% Transition matrix for elasticity matrix from reference material coordinate system to global coordinate system
P = [c^2 s^2 -c*s;
     s^2 c^2 c*s;
     2*c*s -2*c*s c^2-s^2];

% 2D rotation matrix around z-axis by angle -theta
% Transition matrix from global coordinate system to reference material coordinate system
% vGlob = R * vMat
% R = rotz(-ang); R = R(1:2,1:2);
c = cos(-theta);
s = sin(-theta);
R = [c -s;
    s c];
n1 = VECTOR(R(:,1));
n2 = VECTOR(R(:,2));
syscoord = CARTESIAN2D(n1,n2); % reference material coordinate system
Pmat = calc_Pmat(syscoord);
norm(P-Pmat)
% Elasticity matrix in global coordinate system [Pa]
C = P'*Cref*P;


