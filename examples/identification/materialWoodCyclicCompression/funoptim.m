function f = funoptim(x,delta_exp,N,varargin)
% function f = funoptim(x,delta_exp,N,varargin)

delta = solveCyclicCompression(x,N);
f = norm(delta - delta_exp)^2;

end