function f = funlsqnonlinNum(param,u_exp_in,S)
% function f = funlsqnonlinNum(param,u_exp_in,S)

u_in = solveTractionIsotTrans(param,S);
f = u_in - u_exp_in;

end