function f = funlsqnonlinAna(x,u_exp,coorx,coory,F,k,Iz,Mesh,h)
% function f = funlsqnonlinAna(x,U_exp,coorx,coory,F,k,Iz,Mesh,h)

u_ana = ThreePointsBendingAna(x,coorx,coory,F,k,Iz,Mesh,h);
f = u_ana-u_exp;

end