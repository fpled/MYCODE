function [F,J] = funlsqnonlinPF(d,A,b,varargin)
% function [F,J] = funlsqnonlinPF(d,A,b,varargin)

%% Least-squares minimization of squared 2-norm of the residual norm(R*d-c)^2
% with gradient 2*R'*(R*d-c) = 2*(A*d-b), and Hessian 2*R'*R = A

% Cholesky factorization for least-squares problem equivalent to quadratic programming problem
A = (A+A')/2;
R = chol(A); % A = R'*R
c = R'\b; % R'*c = b, c'*R = b', c' = b'/R

% Compute objective function
F = R*d - c;

% Least-squares minimization problem
% norm(R*d-c)^2 = (R*d)'*(R*d) - 2*c'*(R*d) + c'*c 
%               = d'*(R'*R)*d - 2*(R'*c)'*d + (b'/R)*(R'\b)
%               = d'*A*d - 2*(b'*d) + b'*(A\b)
% equivalent to convex quadratic minimization problem
% 1/2*(d'*A*d) - b'*d
% with gradient A*d-b, and Hessian A

if nargout > 1 % Jacobian required
    % Compute Jacobian of objective function
    J = R;
end

%% Mininimize squared 2-norm of the residual norm(A*d-b)^2
% with gradient A'*(A*d-b), and Hessian A'*A

% % Compute objective function
% F = A*d - b;
% 
% if nargout > 1 % Jacobian required
%     % Compute Jacobian of objective function
%     J = A;
% end