function f = funlsqnonlinIsotTrans(x,U_exp,noise)
% function f = funlsqnonlinIsotTrans(x,U_exp,noise)

U_exp = U_exp + noise;
U = ThreePointsBendingIsotTrans(x);
f = U-U_exp;

end