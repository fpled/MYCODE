function f = funoptimThreePointBending(x,fun,u_exp,varargin)
% function f = funoptimThreePointBending(x,fun,u_exp,varargin)

fun = fcnchk(fun);
u = fun(x,varargin{:});
f = norm(u - u_exp)^2;

end