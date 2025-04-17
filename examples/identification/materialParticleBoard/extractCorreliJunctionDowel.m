function [u_exp_a,u_exp_b,coord_a,coord_b,cl_a,cl_b] = extractCorreliJunctionDowel(Job_a,Job_b,Mesh_a,Mesh_b,U_a,U_b,h)
% function [u_exp_a,u_exp_b,coord_a,coord_b,cl_a,cl_b] = extractCorreliJunctionDowel(Job_a,Job_b,Mesh_a,Mesh_b,U_a,U_b,h)

X_a = real(Mesh_a.Znode);
Y_a = imag(Mesh_a.Znode);
X_b = real(Mesh_b.Znode);
Y_b = imag(Mesh_b.Znode);
Coordx_a = (X_a+Job_a.ROI(1)-1);
Coordy_a = (Y_a+Job_a.ROI(2)-1);
Coordx_b = (X_b+Job_b.ROI(1)-1);
Coordy_b = (Y_b+Job_b.ROI(2)-1);
scaleFactorx_a = 2*h/(max(Coordx_a)-min(Coordx_a));
scaleFactory_a = h/(max(Coordy_a)-min(Coordy_a));
scaleFactorx_b = h/(max(Coordx_b)-min(Coordx_b));
scaleFactory_b = 3*h/(max(Coordy_b)-min(Coordy_b));
scaleFactor = mean([scaleFactorx_a scaleFactory_a scaleFactorx_b scaleFactory_b]);

coordx_a = Coordy_a*scaleFactor;
coordy_a = -Coordx_a*scaleFactor;
coordx_b = Coordy_b*scaleFactor;
coordy_b = -Coordx_b*scaleFactor;

mean_coordx_a = (max(coordx_a)+min(coordx_a))/2;
mean_coordy_b = (max(coordy_b)+min(coordy_b))/2;

coordx_a = coordx_a-mean_coordx_a;
coordy_b = coordy_b-mean_coordy_b;
coordy_a = coordy_a-max(coordy_a)+min(coordy_b);
coordx_b = coordx_b-mean_coordx_a; % PhD thesis Zhou Chen
% coordx_b = coordx_b-min(coordx_b)+min(coordx_a);

coord_a = [coordx_a coordy_a];
coord_b = [coordx_b coordy_b];

Ux_a = U_a(1:2:end);
Uy_a = U_a(2:2:end);
Ux_b = U_b(1:2:end);
Uy_b = U_b(2:2:end);
% ux_exp_a = Uy_a*scaleFactor;
% uy_exp_a = -Ux_a*scaleFactor;
% u_exp_a = [ux_exp_a uy_exp_a]';
% ux_exp_b = Uy_b*scaleFactor;
% uy_exp_b = -Ux_b*scaleFactor;
% u_exp_b = [ux_exp_b uy_exp_b]';
u_exp_a = [Uy_a -Ux_a]'*scaleFactor;
u_exp_a = u_exp_a(:);
u_exp_b = [Uy_b -Ux_b]'*scaleFactor;
u_exp_b = u_exp_b(:);

cl_a = Mesh_a.CharLength*scaleFactor;
cl_b = Mesh_b.CharLength*scaleFactor;

end

