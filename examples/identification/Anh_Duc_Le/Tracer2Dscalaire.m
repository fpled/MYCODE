function []=Tracer2Dscalaire(Vecteur,Nx,Ny)

figure

X = linspace(min(Nx), max(Nx), 200);
Y = linspace(min(Ny), max(Ny), 20);
[X Y] = meshgrid(X, Y);
Z = griddata(Nx,Ny,Vecteur,X,Y);
[~,h]=contourf(X,Y,Z,'LineStyle','none');
axis equal
box on
colorbar
