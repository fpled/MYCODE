function [F,J] = funlsqnonlinPF(d,A,b,varargin)
% function [F,J] = funlsqnonlinPF(d,A,b,varargin)

% Compute objective function
F = A*d - b;

if nargout > 1 % jacobian required
    % Compute jacobian of objective function
    J = A;
end