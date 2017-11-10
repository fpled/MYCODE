function f = funoptimAna(param,u_exp,coord,F,Iz,h)
% function f = funoptimAna(param,u_exp,coord,F,Iz,h)

u_ana = solveThreePointBendingAna(param,coord,F,Iz,h);
f = norm(u_ana - u_exp)^2;

end