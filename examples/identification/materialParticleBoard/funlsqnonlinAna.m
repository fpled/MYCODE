function f = funlsqnonlinAna(param,u_exp,coord,F,Iz,h)
% function f = funlsqnonlinAna(param,u_exp,coord,F,Iz,h)

u_ana = solveThreePointBendingAna(param,coord,F,Iz,h);
f = u_ana - u_exp;

end