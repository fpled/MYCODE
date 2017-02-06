function [F] = NeumannSurfLoad(Nx,Ny,lx,ly,dx,ndof,DGOption)
% linear force applied on the top surface
% two points de Gauss in 1D

x = Nx(Ny==ly);
numNode=find(Ny==ly);

ss=[-1/sqrt(3), 1/sqrt(3)]; 
ww=[1,1];

Nss1=[-1/2*ss(1)+1/2,0,1/2*ss(1)+1/2,0;0,-1/2*ss(1)+1/2,0,1/2*ss(1)+1/2];
Nss2=[-1/2*ss(2)+1/2,0,1/2*ss(2)+1/2,0;0,-1/2*ss(2)+1/2,0,1/2*ss(2)+1/2];

J = sqrt(1/4*dx^2);

Fe = zeros(length(x),1);

for i=1:(length(x)-1)
    switch DGOption
        case 'CST'
            f_max=-1;
            fss1=[0;f_max];
            fss2=[0;f_max];
        case 'LIN'
            f_l=1/lx*x(i)-1;
            f_r=1*x(i+1)-1;
            fss1=[0;((f_r-(f_l+f_r)/2)*ss(1)+(f_l+f_r)/2)];
            fss2=[0;((f_r-(f_l+f_r)/2)*ss(2)+(f_l+f_r)/2)];
        case 'QUA'
            f_max=-1;
            x_ss1=x(i)+(-1/2*x(i)+1/2*x(i+1))*abs(-1-ss(1));
            x_ss2=x(i)+(-1/2*x(i)+1/2*x(i+1))*abs(-1-ss(2));
            fss1=[0 ; (4*abs(f_max)/lx^2)*(x_ss1-lx/2).^2+f_max];
            fss2=[0 ; (4*abs(f_max)/lx^2)*(x_ss2-lx/2).^2+f_max];
    end

    fe=Nss1'*fss1*J*ww(1)+Nss2'*fss2*J*ww(2);
    fe=fe(fe~=0);
    Fe(i:i+1)=Fe(i:i+1)+fe;
end

ddlfe=[numNode*2];
F = sparse(ddlfe, 1, Fe(:),ndof,1);