function f = funlsqnonlinIsot(x,U_exp,noise)
% function f = funlsqnonlinIsot(x,U_exp,noise)

U_exp = U_exp + noise;
U = solveThreePointBendingIsot(x);
f = U-U_exp;

end