function f = funoptimAna(x,u_exp,coorx,coory,F,k,Iz,Mesh,h)
% function f = funoptimAna(x,u_exp,coorx,coory,F,k,Iz,Mesh,h)

u_ana = ThreePointsBendingAna(x,coorx,coory,F,k,Iz,Mesh,h);
f = norm(u_ana-u_exp)^2;

end