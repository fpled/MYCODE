%% Comparison between the analytical displacement and the experimental data %%
%%--------------------------------------------------------------------------%%

clc
clear all
close all

samp_num = 'C9';
k = 10; % image number

F = applied_load(samp_num);
[b,h,dis_supp,Iz] = dim_sample(samp_num);

if k<10
    im = ['0' num2str(k)];
else
    im = num2str(k);
end

FileName = [samp_num '_00-' im '-Mesh'];
PathName = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','Zhou','results_data');
disp(['File = ',FileName]);
FullFileName = fullfile(PathName,FileName);
load(FullFileName);
X = real(Mesh.Znode);
Y = imag(Mesh.Znode);
ROI = Job.ROI;
Ux = U(1:2:end);
Uy = U(2:2:end);
Coorx = (X+ROI(1)-1);
Coory = (Y+ROI(2)-1);

% Change unit in mm and the coordinate system due to the use of Correli RT3
scale_factor = h/( max(Coorx) - min(Coorx) );
coorx = (Coory-min(Coory))*scale_factor+dis_supp;
coory = -(Coorx-1/2*(min(Coorx)+max(Coorx)))*scale_factor;
ux_exp = Uy*scale_factor;
uy_exp = -Ux*scale_factor;

loadPath = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materPropPartiBoards');
loadData = 'result_ET_GL.mat';
load([loadPath,filesep,loadData]);
eval( ['E=ET_' samp_num '(' im ')*1000;'] ); %Mpa
eval( ['G=GL_' samp_num '(' im ');'] ); %Mpa
eval( ['Phi=Phi_' samp_num '(' im ');'] );
eval( ['U0=U0_' samp_num '(' im ');'] );
eval( ['V0=V0_' samp_num '(' im ');'] );

ux_ana = @(x,y) F(k)*x^2*y/(4*E*Iz)+F(k)*y^3/(12*G*Iz)+Phi*y+U0;
uy_ana = @(x,y) F(k)*x^3/(12*E*Iz)-F(k)*h^2*x/(16*G*Iz)-Phi*x+V0;
Nelem = size(Mesh.TRI,1);

figure
subplot(2,1,1)
title('U_{x} experimental')
for i=1:Nelem
    ni=Mesh.TRI(i,1); nj=Mesh.TRI(i,2); nk=Mesh.TRI(i,3);
    xi=coorx(ni); yi=coory(ni);
    xj=coorx(nj); yj=coory(nj);
    xk=coorx(nk); yk=coory(nk);
    ui=ux_exp(ni);uj=ux_exp(nj);uk=ux_exp(nk);
    hold on
    fill([xi;xj;xk],[yi;yj;yk],[ui;uj;uk],'LineStyle','none');
end
axis equal
Unitx=' (mm)';
xlabel(['$x$ ',Unitx],'Interpreter','latex')
ylabel(['$y$ ',Unitx],'Interpreter','latex')
colorbar
kk=caxis;
subplot(2,1,2)
title('U_{x} analytical')
for i=1:Nelem
    ni=Mesh.TRI(i,1); nj=Mesh.TRI(i,2); nk=Mesh.TRI(i,3);
    xi=coorx(ni); yi=coory(ni);
    xj=coorx(nj); yj=coory(nj);
    xk=coorx(nk); yk=coory(nk);
    ui=ux_ana(xi,yi);
    uj=ux_ana(xj,yj);
    uk=ux_ana(xk,yk);
    hold on
    fill([xi;xj;xk],[yi;yj;yk],[ui;uj;uk],'LineStyle','none');
end
axis equal
Unitx=' (mm)';
xlabel(['$x$ ',Unitx],'Interpreter','latex')
ylabel(['$y$ ',Unitx],'Interpreter','latex')
caxis(kk)
colorbar

figure
subplot(2,1,1)
title('U_{y} experimental')
for i=1:Nelem
    ni=Mesh.TRI(i,1); nj=Mesh.TRI(i,2); nk=Mesh.TRI(i,3);
    xi=coorx(ni); yi=coory(ni);
    xj=coorx(nj); yj=coory(nj);
    xk=coorx(nk); yk=coory(nk);
    ui=uy_exp(ni);uj=uy_exp(nj);uk=uy_exp(nk);
    hold on
    fill([xi;xj;xk],[yi;yj;yk],[ui;uj;uk],'LineStyle','none');
end
axis equal
Unitx=' (mm)';
xlabel(['$x$ ',Unitx],'Interpreter','latex')
ylabel(['$y$ ',Unitx],'Interpreter','latex')
colorbar
subplot(2,1,2)
title('U_{y} analytical')
for i=1:Nelem
    ni=Mesh.TRI(i,1); nj=Mesh.TRI(i,2); nk=Mesh.TRI(i,3);
    xi=coorx(ni); yi=coory(ni);
    xj=coorx(nj); yj=coory(nj);
    xk=coorx(nk); yk=coory(nk);
    ui=uy_ana(xi,yi);
    uj=uy_ana(xj,yj);
    uk=uy_ana(xk,yk);
    hold on
    fill([xi;xj;xk],[yi;yj;yk],[ui;uj;uk],'LineStyle','none');
end
axis equal
Unitx=' (mm)';
xlabel(['$x$ ',Unitx],'Interpreter','latex')
ylabel(['$y$ ',Unitx],'Interpreter','latex')
colorbar
