function f = funlsqnonlin(x,delta_exp,N,varargin)
% function f = funlsqnonlin(x,delta_exp,N,varargin)

delta = solveCyclicCompression(x,N);
f = delta - delta_exp;

end