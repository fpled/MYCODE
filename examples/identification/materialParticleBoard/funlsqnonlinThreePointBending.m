function f = funlsqnonlinThreePointBending(x,fun,u_exp,varargin)
% function f = funlsqnonlinThreePointBending(x,fun,u_exp,varargin)

fun = fcnchk(fun);
u = fun(x,varargin{:});
f = u - u_exp;

end