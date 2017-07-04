clc
clear all
close all

samp_num = 'C9';

F = applied_load(samp_num);
[b,h,dis_supp,Iz] = dim_sample(samp_num);

pathResult = fullfile(getfemobjectoptions('path'),'MYCODE','results',...
    'identification','materPropPartiBoards');
if ~exist(pathResult,'dir')
    mkdir(pathResult);
end

loadPath = fullfile(getfemobjectoptions('path'),'MYCODE','results',...
    'identification','materPropPartiBoards');
loadData = 'result_ET_GL.mat';
load([loadPath,filesep,loadData])

t1=tic;

% for k=1:length(F)
for k=4
    
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
    ecscale_factor = h/( max(Coorx) - min(Coorx) );
    coorx = (Coory-min(Coory))*ecscale_factor+dis_supp;
    coory = -(Coorx-1/2*(min(Coorx)+max(Coorx)))*ecscale_factor;
    ux = Uy*ecscale_factor;
    uy = -Ux*ecscale_factor;
    
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
    u_exp_edge=zeros(2*N_node_egde,1);
    u_exp_edge(1:2:end)=ux(nodes_edgerd);
    u_exp_edge(2:2:end)=uy(nodes_edgerd);
    u_exp_in = zeros(2*Mesh.NNodeIn,1);
    u_exp_in(1:2:end) = ux(nodes_inter);
    u_exp_in(2:2:end) = uy(nodes_inter);
    
%     set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',16,'DefaultTextInterpreter','tex');
%     figure('name','Imposed experimental displacement')
%     Scal=10;
%     Unitx=' (mm.)';
%     UnitU=' (mm.)';
%     plot(coorx(nodes_edgerd),coory(nodes_edgerd),'r*')
%     axis equal
%     xlabel(['$x$ ',Unitx],'Interpreter','latex')
%     ylabel(['$y$ ',Unitx],'Interpreter','latex')
%     hold on
%     triplot(Mesh.TRI,coorx,coory,'k');
%     legend('U_{exp}');
    
    eval( ['ET=ET_' samp_num '(' im ')*10^3;'] ); % MPa
    eval( ['GL=GL_' samp_num '(' im ');'] );      % MPa
    EL_series=linspace(100,700,10);
    nuL_series=linspace(0.01,0.1,10);
    
    err=zeros(length(EL_series),length(nuL_series));
    
    for m=1:length(EL_series)
        EL = EL_series(m);
        for n=1:length(nuL_series)
            nuL = nuL_series(n);

            S = [ 1/EL -nuL/EL 0
                 -nuL/EL 1/ET  0
                  0     0     1/GL ];
            
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
            
            u_FEMU_in = -K1\(K3*u_exp_edge);
            
            err(m,n) = norm(u_exp_in-u_FEMU_in);
              
        end
    end
    
    err = err./norm(u_exp_in);
    [errmin,I] = min(err);
    [errmin,c] = min(errmin);
    r = I(c);
    
    figure
    surfc(EL_series,nuL_series,err,'EdgeColor','none');
    colorbar
    hold on
    scatter3(EL_series(r),nuL_series(c),errmin,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    % set(gca,'ZScale','log')
    xlabel('$E^L$ (MPa)','Interpreter','latex')
    ylabel('$\nu^L$','Interpreter','latex')
    zlabel('$\varepsilon$','Interpreter','latex')
%     mysaveas(pathResult,['error3D_EL_nuL_' samp_num '_im' im],'fig');
    
    figure
    contourf(EL_series,nuL_series,err,30);
    colorbar
    hold on
    scatter(EL_series(r),nuL_series(c),'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    % set(gca,'ZScale','log')
    xlabel('$E^L$ (MPa)','Interpreter','latex')
    ylabel('$\nu^L$','Interpreter','latex')
%     mysaveas(pathResult,['error2D_EL_nuL_' samp_num '_im' im],'fig');
    
    toc(t1)
end


