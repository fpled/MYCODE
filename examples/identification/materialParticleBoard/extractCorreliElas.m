function [u_exp,coord] = extractCorreliElas(Job,Mesh,U,h,d)
% function [u_exp,coord] = extractCorreliElas(Job,Mesh,U,h,d)

X = real(Mesh.Znode);
Y = imag(Mesh.Znode);
Coordx = (X+Job.ROI(1)-1);
Coordy = (Y+Job.ROI(2)-1);
scaleFactor = h/(max(Coordx)-min(Coordx));
coordx = Coordy*scaleFactor;
coordy = -Coordx*scaleFactor;

coordx = coordx-min(coordx)+d;
coordy = coordy-(max(coordy)+min(coordy))/2;
coord = [coordx coordy];

Ux = U(1:2:end);
Uy = U(2:2:end);
% ux_exp = Uy*scaleFactor;
% uy_exp = -Ux*scaleFactor;
% u_exp = [ux_exp uy_exp]';
u_exp = [Uy -Ux]'*scaleFactor;
u_exp = u_exp(:);

end

