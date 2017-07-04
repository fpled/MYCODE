function f = funoptimFEMU(ET,GL,x,U_exp_in,Mesh,coorx,coory,u_exp_edge)
% function f = funoptimFEMU(ET,GL,x,U_exp_in,Mesh,coorx,coory,u_exp_edge)

u_FEMU_in = ThreePointsBendingFEMU(ET,GL,x,Mesh,coorx,coory,u_exp_edge);
f = norm(u_FEMU_in-U_exp_in)^2;

end