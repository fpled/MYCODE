function [u] = DirichletSurfLoad(Nx,loaddofs,lx,u,DGOption)
% for a parabolic displacement compression test

switch DGOption
    case 'CST'
        u_max = -1e-2;
        u(loaddofs) = u_max;
    case 'LIN'
        u_l = -1e-2;
        u_r = 0;
        Nx_loaddofs=Nx(loaddofs/2);
        u_loadnode=(u_r-u_l)*Nx_loaddofs/lx+u_l;
        u(loaddofs) = u_loadnode;
    case 'QUA'
        u_max = -1e-2;
        Nx_loaddofs=Nx(loaddofs/2);
        u_loadnode=(4*abs(u_max)/lx^2)*(Nx_loaddofs-lx/2).^2+u_max;
        u(loaddofs) = u_loadnode;
end