function f = funlsqnonlinNum(param,u_exp_in,S)
% function f = funlsqnonlinNum(param,u_exp_in,S)

u_num_in = solveThreePointBendingNum(param,S);
f = u_num_in - u_exp_in;

end