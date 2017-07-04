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
    ech=h/( max(Coorx) - min(Coorx) );
    coorx = (Coory-min(Coory))*ech+dis_supp;
    coory = -(Coorx-1/2*(min(Coorx)+max(Coorx)))*ech;
    ux = Uy*ech;
    uy = -Ux*ech;
    u_exp = zeros(2*Mesh.NNodeTot,1);
    u_exp(1:2:end) = ux;
    u_exp(2:2:end) = uy;
    
    eval( ['Phi=Phi_' samp_num '(' im ');'] );
    eval( ['U0=U0_' samp_num '(' im ');'] );
    eval( ['V0=V0_' samp_num '(' im ');'] );
    ET_series = linspace(3000,5000,2000); % Mpa
    GL_series = linspace(100,200,100); % Mpa
    err = zeros(length(ET_series),length(GL_series));
    ux_ana = zeros(Mesh.NNodeTot,1);
    uy_ana = zeros(Mesh.NNodeTot,1);
    u_ana = zeros(2*Mesh.NNodeTot,1);
    
    for m=1:length(ET_series)
        ET = ET_series(m);
        for n=1:length(GL_series)
            GL = GL_series(n);
            
            for i = 1:Mesh.NNodeTot
                xi = coorx(i);
                yi = coory(i);
                ux_ana(i) = -F(k)*xi^2*yi/(4*ET*Iz)+F(k)*yi^3/(12*GL*Iz)+Phi*yi+U0;
                uy_ana(i) = F(k)*xi^3/(12*ET*Iz)-F(k)*h^2*xi/(16*GL*Iz)-Phi*xi+V0;
            end
            u_ana(1:2:end) = ux_ana;
            u_ana(2:2:end) = uy_ana;
            err(m,n) = norm(u_exp-u_ana)/norm(u_exp);
            
        end
    end
    [errmin,I] = min(err);
    [errmin,c] = min(errmin);
    r = I(c);
    
    figure
    surfc(GL_series,ET_series,err,'EdgeColor','none');
    colorbar
    hold on
    scatter3(GL_series(c),ET_series(r),errmin,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    % set(gca,'ZScale','log')
    xlabel('$G^L$ (MPa)','Interpreter','latex')
    ylabel('$E^T$ (MPa)','Interpreter','latex')
    zlabel('$\varepsilon$','Interpreter','latex')
%     mysaveas(pathResult,['error3D_ET_GL_' samp_num '_im' im],'fig');
    
    figure
    contourf(GL_series,ET_series,err,30);
    colorbar
    hold on
    scatter(GL_series(c),ET_series(r),'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    % set(gca,'ZScale','log')
    xlabel('$G^L$ (MPa)','Interpreter','latex')
    ylabel('$E^T$ (MPa)','Interpreter','latex')
%     mysaveas(pathResult,['error2D_ET_GL_' samp_num '_im' im],'fig');
    toc(t1)
    
    
end



