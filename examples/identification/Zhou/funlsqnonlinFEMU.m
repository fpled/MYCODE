function f = funlsqnonlinFEMU(ET,GL,x,u_exp_in,Mesh,coorx,coory,u_exp_edge)
% function f = funlsqnonlinFEMU(ET,GL,x,u_exp,Mesh,coorx,coory,ub)

u_FEMU_in = ThreePointsBendingFEMU(ET,GL,x,Mesh,coorx,coory,u_exp_edge);
f = u_FEMU_in-u_exp_in;

end