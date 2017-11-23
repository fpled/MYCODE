function f = funoptimNum(param,u_exp_in,S)
% function f = funoptimNum(param,u_exp_in,S)

u_in = solveTractionIsotTrans(param,S);
f = norm(u_in - u_exp_in)^2;

end