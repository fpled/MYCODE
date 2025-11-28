function [f, g, H] = funoptimPF(d,A,b,varargin)
% function [f, g, h] = funoptimPF(d,A,b,varargin)

% Compute objective function
f = 1/2*d'*A*d - d'*b;

if nargout > 1 % gradient required
    % Compute gradient of objective function
    g = A*d - b;
    
    if nargout > 2 % Hessian required
        % Compute Hessian of objective function
        H = A;
    end
end