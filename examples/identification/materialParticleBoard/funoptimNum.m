function f = funoptimNum(param,u_exp_in,S,h)
% function f = funoptimNum(param,u_exp_in,S,h)

u_num_in = solveThreePointBendingNum(param,S,h);
f = norm(u_num_in - u_exp_in)^2;

end