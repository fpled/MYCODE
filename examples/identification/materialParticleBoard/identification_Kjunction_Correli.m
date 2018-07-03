%% Identification of junction stiffness %%
%%--------------------------------------%%
% Call after Digital Image Correlation (DIC) RT3
% Coordinate system
% CORRELI:    right     down
% MATLAB:     right     up

clc
clearvars
close all

%% Input data
displaySolution = true;

filename = 'data_angle_Correli.mat';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
if ~exist(pathname,'dir')
    mkdir(pathname);
end
pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','materialParticleBoard','resultsDIC');
fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc2'};

%% Identification

numScrew = 4;
numDowel = 2;

d = 67.5; % mm
b = 113; % mm
% Same Dimensions for the two assembly junctions
% Moment per unit length: ml = F*d/b
% Junction stiffness per unit length: k = ml/var_angle

F_screw = cell(numScrew,1);
ml_screw = cell(numScrew,1);
var_angle_screw = cell(numScrew,1);
k_screw = cell(numScrew,1);
for j=1:numScrew
% for j=1
    
    numSample = ['S' num2str(j)];
    F = appliedLoad(numSample);
    numImages = length(F);
    var_angle = zeros(numImages,1);
    
    filenameDICa = ['S' num2str(j) 'a'];
    load(fullfile(pathnameDIC,filenameDICa));
    U_a = U;
    V_a = V;
    coorx_a = coorx;
    coory_a = coory;
    x_box_a = x_boite;
    y_box_a = y_boite;
    scaleFactor1_a = 45/(max(x_box_a)-min(x_box_a));
    scaleFactor2_a = 15/(max(y_box_a)-min(y_box_a));
    
    clear U V coorx coory x_boite y_boite
    filenameDICb = ['S' num2str(j) 'b'];
    load(fullfile(pathnameDIC,filenameDICb));
    U_b = U;
    V_b = V;
    coorx_b = coorx;
    coory_b = coory;
    x_box_b = x_boite;
    y_box_b = y_boite;
    scaleFactor1_b = 15/(max(x_box_b)-min(x_box_b));
    scaleFactor2_b = 30/(max(y_box_b)-min(y_box_b));
    
    scaleFactor = mean([scaleFactor1_a scaleFactor2_a...
        scaleFactor1_b scaleFactor2_b]);
    
    % Change of coordinate system from Correli to MATLAB
    step = (coorx(2)-coorx(1))*scaleFactor;
    coordx_a = coory_a*scaleFactor;
    coordy_a = -coorx_a*scaleFactor;
    coordx_b = coory_b*scaleFactor;
    coordy_b = -coorx_b*scaleFactor;
    
    [xi_a,yi_a] = meshgrid(coordx_a,coordy_a);
    xi_a = xi_a - min(y_box_a)*scaleFactor;
    yi_a = yi_a + max(x_box_a)*scaleFactor;
    [xi_b,yi_b] = meshgrid(coordx_b,coordy_b);
    xi_b = xi_b - min(y_box_a)*scaleFactor;
    yi_b = yi_b + max(x_box_a)*scaleFactor;
    
    TRI_a = delaunay(xi_a(:),yi_a(:));
    TRI_b = delaunay(xi_b(:),yi_b(:));
    
    for k=1:numImages
    % for k=3
        
        u_a = V_a(:,:,k)*scaleFactor;
        v_a = -U_a(:,:,k)*scaleFactor;
        u_b = V_b(:,:,k)*scaleFactor;
        v_b = -U_b(:,:,k)*scaleFactor;
        points_a = find( yi_a(:)>min(min(yi_b))-step &...
            yi_a(:)<max(max(yi_b)+step)...
            & xi_a(:) == max(max(xi_a)) );
        points_b = find( xi_b(:) == min(min(xi_b)) );
        
        % primary line1 and line2
        L1x0 = xi_a(points_a);
        L1y0 = yi_a(points_a);
        [L1y0_sort,index] = sort(L1y0);
        L1x0_sort = L1x0(index);
        fit10 = polyfit(L1y0_sort,L1x0_sort,1);
        val10 = polyval(fit10,L1y0_sort);
        
        L2x0 = xi_b(points_b);
        L2y0 = yi_b(points_b);
        [L2y0_sort,index] = sort(L2y0);
        L2x0_sort = L2x0(index);
        fit20 = polyfit(L2y0_sort,L2x0_sort,1);
        val20 = polyval(fit20,L2y0_sort);
        
        % deformed line1 and line2
        L1x = xi_a(points_a)+u_a(points_a);
        L1y = yi_a(points_a)+v_a(points_a);
        [L1y_sort,index] = sort(L1y);
        L1x_sort = L1x(index);
        fit1 = polyfit(L1y_sort,L1x_sort,1);
        val1 = polyval(fit1,L1y_sort);
        
        L2x = xi_b(points_b)+u_b(points_b);
        L2y = yi_b(points_b)+v_b(points_b);
        [L2y_sort,index] = sort(L2y);
        L2x_sort = L2x(index);
        fit2 = polyfit(L2y_sort,L2x_sort,1);
        val2 = polyval(fit2,L2y_sort);
        
        %-----------------------------
        % primary and deformed angle of junction
        %-----------------------------
        
        t = [val10(end)-val10(1) L1y0_sort(end)-L1y0_sort(1)];
        s = [val20(end)-val20(1) L2y0_sort(end)-L2y0_sort(1)];
        delta0 = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        t = [val1(end)-val1(1) L1y_sort(end)-L1y_sort(1)];
        s = [val2(end)-val2(1) L2y_sort(end)-L2y_sort(1)];
        delta = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        var_angle(k) = abs(delta0-delta);
        
        %-----------------------------
        % Reference and deformed mesh
        %-----------------------------
        set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',16,...
            'DefaultTextInterpreter','tex');
        Scal = 1;
        Unitx = '(mm)';
        UnitU = '(mm)';
        
%         figure('name','best fit line of initial mesh')
%         fig1 = triplot(TRI_a,xi_a,yi_a,'r');
%         hold on
%         fig2 = triplot(TRI_b,xi_b,yi_b,'r');
%         fig3 = plot(L1x0,L1y0,'k*',val10,L1y0_sort,'k');
%         fig4 = plot(L2x0,L2y0,'k*',val20,L2y0_sort,'k');
%         axis equal
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         hold off
        
        L1x = xi_a(points_a)+Scal*u_a(points_a);
        L1y = yi_a(points_a)+Scal*v_a(points_a);
        [L1y_sort,index] = sort(L1y);
        L1x_sort = L1x(index);
        fit1 = polyfit(L1y_sort,L1x_sort,1);
        val1 = polyval(fit1,L1y_sort);
        
        L2x = xi_b(points_b)+Scal*u_b(points_b);
        L2y = yi_b(points_b)+Scal*v_b(points_b);
        [L2y_sort,index] = sort(L2y);
        L2x_sort = L2x(index);
        fit2 = polyfit(L2y_sort,L2x_sort,1);
        val2 = polyval(fit2,L2y_sort);
        
%         figure('name','best fit line of deformed mesh')
%         fig1 = triplot(TRI_a,xi_a(:)+Scal*u_a(:),...
%             yi_a(:)+Scal*v_a(:),'k');
%         hold on
%         fig2 = triplot(TRI_b,xi_b(:)+Scal*u_b(:),...
%             yi_b(:)+Scal*v_b(:),'k');
%         fig3 = plot(xi_a(points_a)+Scal*u_a(points_a),...
%             yi_a(points_a)+Scal*v_a(points_a),'r*',...
%             val1,L1y_sort,'r');
%         fig4 = plot(xi_b(points_b)+Scal*u_b(points_b),...
%             yi_b(points_b)+Scal*v_b(points_b),'r*',...
%             val2,L2y_sort,'r');
%         axis equal
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         hold off
        
%         figure('name','Reference and deformed mesh')
%         fig1 = triplot(TRI_a,xi_a(:),yi_a(:),'r');
%         hold on
%         fig2 = triplot(TRI_b,xi_b(:),yi_b(:),'r');
%         fig3 = triplot(TRI_a,xi_a(:)+Scal*u_a(:),...
%             yi_a(:)+Scal*v_a(:),'k');
%         fig4 = triplot(TRI_b,xi_b(:)+Scal*u_b(:),...
%             yi_b(:)+Scal*v_b(:),'k');
%         axis equal
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         hold off
        
    end
    
    F_screw{j} = F;
    ml_screw{j} = F*d/b;
    var_angle_screw{j} = var_angle;
    k_screw{j} = (F*d)./(deg2rad(var_angle)'*b);
end

F_dowel = cell(numDowel,1);
ml_dowel = cell(numDowel,1);
var_angle_dowel = cell(numDowel,1);
k_dowel = cell(numDowel,1);
for j=1:numDowel
% for j=1
    
    numSample = ['D' num2str(j)];
    F = appliedLoad(numSample);
    numImages = length(F);
    var_angle = zeros(numImages,1);
    
    filenameDICa = ['D' num2str(j) 'a'];
    load(fullfile(pathnameDIC,filenameDICa));
    U_a = U;
    V_a = V;
    coorx_a = coorx;
    coory_a = coory;
    x_box_a = x_boite;
    y_box_a = y_boite;
    scaleFactor1_a = 45/(max(x_box_a)-min(x_box_a));
    scaleFactor2_a = 15/(max(y_box_a)-min(y_box_a));
    
    clear U V coorx coory x_boite y_boite
    filenameDICb = ['D' num2str(j) 'b'];
    load(fullfile(pathnameDIC,filenameDICb));
    U_b = U;
    V_b = V;
    coorx_b = coorx;
    coory_b = coory;
    x_box_b = x_boite;
    y_box_b = y_boite;
    scaleFactor1_b = 15/(max(x_box_b)-min(x_box_b));
    scaleFactor2_b = 30/(max(y_box_b)-min(y_box_b));
    
    scaleFactor = mean([scaleFactor1_a scaleFactor2_a...
        scaleFactor1_b scaleFactor2_b]);
    
    % Change of coordinate system from Correli to MATLAB
    step = (coorx(2)-coorx(1))*scaleFactor;
    coordx_a = coory_a*scaleFactor;
    coordy_a = -coorx_a*scaleFactor;
    coordx_b = coory_b*scaleFactor;
    coordy_b = -coorx_b*scaleFactor;
    
    [xi_a,yi_a] = meshgrid(coordx_a,coordy_a);
    xi_a = xi_a - min(y_box_a)*scaleFactor;
    yi_a = yi_a + max(x_box_a)*scaleFactor;
    [xi_b,yi_b] = meshgrid(coordx_b,coordy_b);
    xi_b = xi_b - min(y_box_a)*scaleFactor;
    yi_b = yi_b + max(x_box_a)*scaleFactor;
    
    TRI_a = delaunay(xi_a(:),yi_a(:));
    TRI_b = delaunay(xi_b(:),yi_b(:));
    
    for k=1:numImages
    % for k=3
        
        u_a = V_a(:,:,k)*scaleFactor;
        v_a = -U_a(:,:,k)*scaleFactor;
        u_b = V_b(:,:,k)*scaleFactor;
        v_b = -U_b(:,:,k)*scaleFactor;
        points_a = find( yi_a(:) == max(max(yi_a)) );
        points_b = find( yi_b(:) == min(min(yi_b)) &...
            xi_b(:)<max(max(xi_a))+step );
        
        % primary line1 and line2
        L1x0 = xi_a(points_a);
        L1y0 = yi_a(points_a);
        [L1x0_sort,index] = sort(L1x0);
        L1y0_sort = L1y0(index);
        fit10 = polyfit(L1x0_sort,L1y0_sort,1);
        val10 = polyval(fit10,L1x0_sort);
        
        L2x0 = xi_b(points_b);
        L2y0 = yi_b(points_b);
        [L2x0_sort,index] = sort(L2x0);
        L2y0_sort = L2y0(index);
        fit20 = polyfit(L2x0_sort,L2y0_sort,1);
        val20 = polyval(fit20,L2x0_sort);
        
        % deformed line1 and line2
        L1x = xi_a(points_a)+u_a(points_a);
        L1y = yi_a(points_a)+v_a(points_a);
        [L1x_sort,index] = sort(L1x);
        L1y_sort = L1y(index);
        fit1 = polyfit(L1x_sort,L1y_sort,1);
        val1 = polyval(fit1,L1x_sort);
        
        L2x = xi_b(points_b)+u_b(points_b);
        L2y = yi_b(points_b)+v_b(points_b);
        [L2x_sort,index] = sort(L2x);
        L2y_sort = L2y(index);
        fit2 = polyfit(L2x_sort,L2y_sort,1);
        val2 = polyval(fit2,L2x_sort);
        
        %-----------------------------
        % primary and deformed angle of junction
        %-----------------------------
        
        t = [L1x0_sort(end)-L1x0_sort(1) val10(end)-val10(1)];
        s = [L2x0_sort(end)-L2x0_sort(1) val20(end)-val20(1)];
        delta0 = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        t = [L1x_sort(end)-L1x_sort(1) val1(end)-val1(1)];
        s = [L2x_sort(end)-L2x_sort(1) val2(end)-val2(1)];
        delta = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        var_angle(k) = abs(delta0-delta);
        
        %-----------------------------
        % Reference and deformed mesh
        %-----------------------------
        set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',16,...
            'DefaultTextInterpreter','tex');
        Scal = 1;
        Unitx = '(mm)';
        UnitU = '(mm)';
        
%         figure('name','best fit line of initial mesh')
%         fig1 = triplot(TRI_a,xi_a,yi_a,'r');
%         hold on
%         fig2 = triplot(TRI_b,xi_b,yi_b,'r');
%         fig3 = plot(L1x0,L1y0,'k*',val10,L1y0_sort,'k');
%         fig4 = plot(L2x0,L2y0,'k*',val20,L2y0_sort,'k');
%         axis equal
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         hold off
        
        L1x = xi_a(points_a)+Scal*u_a(points_a);
        L1y = yi_a(points_a)+Scal*v_a(points_a);
        [L1x_sort,index] = sort(L1x);
        L1y_sort = L1y(index);
        fit1 = polyfit(L1x_sort,L1y_sort,1);
        val1 = polyval(fit1,L1x_sort);
        
        L2x = xi_b(points_b)+Scal*u_b(points_b);
        L2y = yi_b(points_b)+Scal*v_b(points_b);
        [L2x_sort,index] = sort(L2x);
        L2y_sort = L2y(index);
        fit2 = polyfit(L2x_sort,L2y_sort,1);
        val2 = polyval(fit2,L2x_sort);
        
%         figure('name','best fit line of deformed mesh')
%         fig1 = triplot(TRI_a,xi_a(:)+Scal*u_a(:),...
%             yi_a(:)+Scal*v_a(:),'k');
%         hold on
%         fig2 = triplot(TRI_b,xi_b(:)+Scal*u_b(:),...
%             yi_b(:)+Scal*v_b(:),'k');
%         fig3 = plot(xi_a(points_a)+Scal*u_a(points_a),...
%             yi_a(points_a)+Scal*v_a(points_a),'r*',...
%             L1x_sort,val1,'r');
%         fig4 = plot(xi_b(points_b)+Scal*u_b(points_b),...
%             yi_b(points_b)+Scal*v_b(points_b),'r*',...
%             L2x_sort,val2,'r');
%         axis equal
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         hold off
        
%         figure('name','Reference and deformed mesh')
%         fig1 = triplot(TRI_a,xi_a(:),yi_a(:),'r');
%         hold on
%         fig2 = triplot(TRI_b,xi_b(:),yi_b(:),'r');
%         fig3 = triplot(TRI_a,xi_a(:)+Scal*u_a(:),...
%             yi_a(:)+Scal*v_a(:),'k');
%         fig4 = triplot(TRI_b,xi_b(:)+Scal*u_b(:),...
%             yi_b(:)+Scal*v_b(:),'k');
%         axis equal
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         hold off
        
    end
    
    F_dowel{j} = F;
    ml_dowel{j} = F*d/b;
    var_angle_dowel{j} = var_angle;
    k_dowel{j} = (F*d)./(deg2rad(var_angle)'*b);
end

%% Save variables
save(fullfile(pathname,filename),'numScrew','numDowel',...
    'F_screw','ml_screw','var_angle_screw','k_screw',...
    'F_dowel','ml_dowel','var_angle_dowel','k_dowel');

%% Plot data
if displaySolution
    figure('name','Screw junction: Moment per unit length w.r.t. Variation of angle')
    for j = 1:numScrew
        plot(ml_screw{j},var_angle_screw{j},'+-')
        hold on
        grid on
        box on
        xlabel('Moment per unit length (N.mm/mm)','Interpreter',interpreter);
        ylabel('Variation of angle ($^{\circ}$)','Interpreter',interpreter);
    end
    
    figure('name','Screw junction: Applied load w.r.t. Junction stiffness per unit length')
    for j = 1:numScrew
        plot(F_screw{j},k_screw{j},'+-')
        hold on
        grid on
        box on
        xlabel('Applied load (N)','Interpreter',interpreter);
        ylabel('Junction stiffness per unit length (N/rad)','Interpreter',interpreter);    
    end
    
    figure('name','Dowel junction: Moment per unit length w.r.t. Variation of angle')
    for j = 1:numDowel
        plot(ml_dowel{j},var_angle_dowel{j},'+-')
        hold on
        grid on
        box on
        xlabel('Moment per unit length (N.mm/mm)','Interpreter',interpreter);
        ylabel('Variation of angle ($^{\circ}$)','Interpreter',interpreter);
    end
    
    figure('name','Dowel junction: Applied load w.r.t. Junction stiffness per unit length')
    for j = 1:numDowel
        plot(F_dowel{j},k_dowel{j},'+-')
        hold on
        grid on
        box on
        xlabel('Applied load (N)','Interpreter',interpreter);
        ylabel('Junction stiffness per unit length (N/rad)','Interpreter',interpreter);
    end
end
