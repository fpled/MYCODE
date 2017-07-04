%% Comparison between the numerical displacement and the experimental data %%
%%-------------------------------------------------------------------------%%

clc
clear all
close all

samp_num = 'C9';
k = 10; % image number

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

elemtype = 'TRI3';
% option = 'DEFO'; % plane strain
option = 'CONT'; % plane stress
S = MODEL('PLAN');
node = NODE([coorx(:),coory(:)],1:numel(coorx));
S = addnode(S,node);
elem = Mesh.TRI;
S = addelem(S,elemtype,elem,'option',option);
nodes_edgerd  = getnumber(getnode(create_boundary(S)));
nodes_inter = setdiff(1:Mesh.NNodeTot,nodes_edgerd);
N_node_egde = length(nodes_edgerd);
if N_node_egde ~= Mesh.NNodeEdge
    warning('Wrong Boundary Nodes')
elseif all(nodes_edgerd' ~= Mesh.NNodeTot-Mesh.NNodeEdge+1:Mesh.NNodeTot)
    warning('Wrong Global Matrix')
end
u_FEMU_edge=zeros(2*N_node_egde,1);
u_FEMU_edge(1:2:end)=ux_exp(nodes_edgerd);
u_FEMU_edge(2:2:end)=uy_exp(nodes_edgerd);

loadPath = fullfile(getfemobjectoptions('path'),'MYCODE','results',...
    'identification','materPropPartiBoards');
loadData1 = 'result_ET_GL.mat';
loadData2 = 'result_EL_nuL.mat';
load([loadPath,filesep,loadData1])
load([loadPath,filesep,loadData2])

eval( ['ET=ET_' samp_num '(' im ')*10^3;'] ); % MPa
eval( ['GL=GL_' samp_num '(' im ');'] ); % MPa
eval( ['EL=EL_f1_' samp_num '(' im ');'] ); % MPa
eval( ['nuL=nuL_f1_' samp_num '(' im ');'] );
S = [ 1/EL  -nuL/EL 0
     -nuL/EL 1/ET   0
      0      0      1/GL ];

K=zeros(2*Mesh.NNodeTot);
for i=1:Mesh.Nelem
    x123 = coorx( Mesh.TRI(i,:) );
    y123 = coory( Mesh.TRI(i,:) );
    
    A = [1 x123(1) y123(1)
        1 x123(2) y123(2)
        1 x123(3) y123(3)]\[1 0 0
        0 1 0
        0 0 1];
    
    B = [A(2,1) 0      A(2,2) 0      A(2,3) 0
        0      A(3,1) 0      A(3,2) 0      A(3,3)
        A(3,1) A(2,1) A(3,2) A(2,2) A(3,3) A(2,3)];
    
    Ae = 1/2*( (x123(2)-x123(1))*(y123(3)-y123(1))...
        -(y123(2)-y123(1))*(x123(3)-x123(1)) );
    
    Ke = B'*(S\(B*Ae));
    
    cnt = zeros(1,6);
    cnt(1:2:end)=2*Mesh.TRI(i,:)-1;
    cnt(2:2:end)=2*Mesh.TRI(i,:);
    
    K(cnt,cnt) = K(cnt,cnt) + Ke;
    
end
K1=K(1:2*(Mesh.NNodeTot-N_node_egde),1:2*(Mesh.NNodeTot-N_node_egde));
K3=K(1:2*(Mesh.NNodeTot-N_node_egde),2*(Mesh.NNodeTot-N_node_egde)+1:end);

u_FEMU_in = -K1\(K3*u_FEMU_edge);

u_FEMU = [u_FEMU_in
          u_FEMU_edge];
ux_FEMU = u_FEMU(1:2:end);
uy_FEMU = u_FEMU(2:2:end);

Nelem=size(Mesh.TRI,1);
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
title('U_{x} numerical')
for i=1:Nelem
    ni=Mesh.TRI(i,1); nj=Mesh.TRI(i,2); nk=Mesh.TRI(i,3);
    xi=coorx(ni); yi=coory(ni);
    xj=coorx(nj); yj=coory(nj);
    xk=coorx(nk); yk=coory(nk);
    ui=ux_FEMU(ni);
    uj=ux_FEMU(nj);
    uk=ux_FEMU(nk);
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
title('U_{y} numerical')
for i=1:Nelem
    ni=Mesh.TRI(i,1); nj=Mesh.TRI(i,2); nk=Mesh.TRI(i,3);
    xi=coorx(ni); yi=coory(ni);
    xj=coorx(nj); yj=coory(nj);
    xk=coorx(nk); yk=coory(nk);
    ui=uy_FEMU(ni);
    uj=uy_FEMU(nj);
    uk=uy_FEMU(nk);
    hold on
    fill([xi;xj;xk],[yi;yj;yk],[ui;uj;uk],'LineStyle','none');
end
axis equal
Unitx=' (mm)';
xlabel(['$x$ ',Unitx],'Interpreter','latex')
ylabel(['$y$ ',Unitx],'Interpreter','latex')
colorbar
