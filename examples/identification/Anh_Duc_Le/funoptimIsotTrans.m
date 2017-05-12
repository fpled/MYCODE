function f = funoptimIsotTrans(x,U_exp,noise)
% function f = funoptimIsotTrans(x,U_exp,noise)

U_exp = U_exp + noise;
U = ThreePointsBendingIsotTrans(x);
f = norm(U-U_exp)^2;

end