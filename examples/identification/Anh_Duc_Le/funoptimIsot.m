function f = funoptimIsot(x,U_exp,noise)
% function f = funoptimIsot(x,U_exp,noise)

U_exp = U_exp + noise;
U = solveThreePointBendingIsot(x);
f = norm(U-U_exp)^2;

end