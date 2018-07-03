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

filename = 'data_angle.mat';
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
    numSamplea = ['S' num2str(j) 'a'];
    numSampleb = ['S' num2str(j) 'b'];
    F = appliedLoad(numSample);
    numImages = length(F);
    var_angle = zeros(numImages,1);
    
    for k=1:numImages
    % for k=3
        
        numImage = num2str(k,'%02d');
        filenameDICa = [numSamplea '_00-' numImage '-Mesh'];
        filenameDICb = [numSampleb '_00-' numImage '-Mesh'];
        
        load(fullfile(pathnameDIC,filenameDICa));
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx_a_screw = (X+Job.ROI(1)-1);
        Coordy_a_screw = (Y+Job.ROI(2)-1);
        TRI_a_screw = Mesh.TRI;
        scaleFactor_1a = 45/(max(Coordx_a_screw)-min(Coordx_a_screw));
        scaleFactor_2a = 15/(max(Coordy_a_screw)-min(Coordy_a_screw));
        Ux_a_screw = U(1:2:end);
        Uy_a_screw = U(2:2:end);
        
        clear Mesh Job U
        load(fullfile(pathnameDIC,filenameDICb));
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx_b_screw = (X+Job.ROI(1)-1);
        Coordy_b_screw = (Y+Job.ROI(2)-1);
        TRI_b_screw = Mesh.TRI;
        scaleFactor_1b = 15/(max(Coordx_b_screw)-min(Coordx_b_screw));
        scaleFactor_2b = 30/(max(Coordy_b_screw)-min(Coordy_b_screw));
        Ux_b_screw = U(1:2:end);
        Uy_b_screw = U(2:2:end);
        
        scaleFactor = mean([scaleFactor_1a scaleFactor_2a...
            scaleFactor_1b scaleFactor_2b]);
        
        u_exp_a_screw = [Uy_a_screw -Ux_a_screw]'*scaleFactor;
        u_exp_a_screw = u_exp_a_screw(:);
        u_exp_b_screw = [Uy_b_screw -Ux_b_screw]'*scaleFactor;
        u_exp_b_screw = u_exp_b_screw(:);
        
        coordx_a_screw = Coordy_a_screw*scaleFactor;
        coordy_a_screw = -Coordx_a_screw*scaleFactor;
        min_coordx_a_screw = min(coordx_a_screw);
        min_coordy_a_screw = min(coordy_a_screw);
        coordx_a_screw(:) = coordx_a_screw(:) - min_coordx_a_screw;
        coordy_a_screw(:) = coordy_a_screw(:) - min_coordy_a_screw;
        coord_a_screw = [coordx_a_screw coordy_a_screw];
        
        coordx_b_screw = Coordy_b_screw*scaleFactor;
        coordy_b_screw = -Coordx_b_screw*scaleFactor;
        coordx_b_screw(:) = coordx_b_screw(:) - min_coordx_a_screw;
        coordy_b_screw(:) = coordy_b_screw(:) - min_coordy_a_screw;
        coord_b_screw = [coordx_b_screw coordy_b_screw];
                  
%         node_a_screw = NODE(coord_a_screw,1:size(coord_a_screw,1));
%         node_b_screw = NODE(coord_b_screw,1:size(coord_b_screw,1));
%         elem_a_screw = TRI_a_screw;
%         elem_b_screw = TRI_b_screw;
%         elemtype = 'TRI3';
%         S_a_screw = MODEL('PLAN');
%         S_b_screw = MODEL('PLAN');
%         S_a_screw = addnode(S_a_screw,node_a_screw);
%         S_b_screw = addnode(S_b_screw,node_b_screw);
%         S_a_screw = addelem(S_a_screw,elemtype,elem_a_screw);
%         S_b_screw = addelem(S_b_screw,elemtype,elem_b_screw);
%         S_a_screw = final(S_a_screw);
%         S_b_screw = final(S_b_screw);
%         numnode_bound_a_screw = getnumber(getnode(create_boundary(S_a_screw)));
%         numnode_bound_b_screw = getnumber(getnode(create_boundary(S_b_screw)));
%         figure
%         plot(create_boundary(S_a_screw));
%         hold on
%         plot(create_boundary(S_b_screw));
%         plot(coordx_a_screw(numnode_bound_a_screw),coordy_a_screw(numnode_bound_a_screw),'r*')
%         plot(coordx_b_screw(numnode_bound_b_screw),coordy_b_screw(numnode_bound_b_screw),'k*')
%         hold off

        points_a_screw = find(coordx_a_screw>max(coordx_a_screw)-Mesh.CharLength*scaleFactor/3 &...
            coordy_a_screw>min(coordy_b_screw));
        points_b_screw = find(coordx_b_screw<min(coordx_b_screw)+Mesh.CharLength*scaleFactor/3 &...
            coordy_b_screw>min(coordy_b_screw));
        
        % primary line1 and line2
        L1x0 = coordx_a_screw(points_a_screw);
        L1y0 = coordy_a_screw(points_a_screw);
        [L1y0_sort,index] = sort(L1y0);
        L1x0_sort = L1x0(index);
        fit10 = polyfit(L1y0_sort,L1x0_sort,1);
        val10 = polyval(fit10,L1y0_sort);
        
        L2x0 = coordx_b_screw(points_b_screw);
        L2y0 = coordy_b_screw(points_b_screw);
        [L2y0_sort,index] = sort(L2y0);
        L2x0_sort = L2x0(index);
        fit20 = polyfit(L2y0_sort,L2x0_sort,1);
        val20 = polyval(fit20,L2y0_sort);
        
        % deformed line1 and line2
        L1x = coordx_a_screw(points_a_screw)+u_exp_a_screw(2*points_a_screw-1);
        L1y = coordy_a_screw(points_a_screw)+u_exp_a_screw(2*points_a_screw);
        [L1y_sort,index] = sort(L1y);
        L1x_sort = L1x(index);
        fit1 = polyfit(L1y_sort,L1x_sort,1);
        val1 = polyval(fit1,L1y_sort);
        
        L2x = coordx_b_screw(points_b_screw)+u_exp_b_screw(2*points_b_screw-1);
        L2y = coordy_b_screw(points_b_screw)+u_exp_b_screw(2*points_b_screw);
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
%         fig1 = triplot(TRI_a_screw,coordx_a_screw,coordy_a_screw,'r');
%         hold on
%         fig2 = triplot(TRI_b_screw,coordx_b_screw,coordy_b_screw,'r');
%         fig3 = plot(L1x0,L1y0,'k*',val10,L1y0_sort,'k');
%         fig4 = plot(L2x0,L2y0,'k*',val20,L2y0_sort,'k');
%         axis equal
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         hold off
        
        L1x = coordx_a_screw(points_a_screw)+Scal*u_exp_a_screw(2*points_a_screw-1);
        L1y = coordy_a_screw(points_a_screw)+Scal*u_exp_a_screw(2*points_a_screw);
        [L1y_sort,index] = sort(L1y);
        L1x_sort = L1x(index);
        fit1 = polyfit(L1y_sort,L1x_sort,1);
        val1 = polyval(fit1,L1y_sort);
        
        L2x = coordx_b_screw(points_b_screw)+Scal*u_exp_b_screw(2*points_b_screw-1);
        L2y = coordy_b_screw(points_b_screw)+Scal*u_exp_b_screw(2*points_b_screw);
        [L2y_sort,index] = sort(L2y);
        L2x_sort = L2x(index);
        fit2 = polyfit(L2y_sort,L2x_sort,1);
        val2 = polyval(fit2,L2y_sort);
        
%         figure('name','best fit line of deformed mesh')
%         fig1 = triplot(TRI_a_screw,coordx_a_screw+Scal*u_exp_a_screw(1:2:end),...
%             coordy_a_screw+Scal*u_exp_a_screw(2:2:end),'k');
%         hold on
%         fig2 = triplot(TRI_b_screw,coordx_b_screw+Scal*u_exp_b_screw(1:2:end),...
%             coordy_b_screw+Scal*u_exp_b_screw(2:2:end),'k');
%         fig3 = plot(coordx_a_screw(points_a_screw)+Scal*u_exp_a_screw(2*points_a_screw-1),...
%             coordy_a_screw(points_a_screw)+Scal*u_exp_a_screw(2*points_a_screw),'r*',...
%             val1,L1y_sort,'r');
%         fig4 = plot(coordx_b_screw(points_b_screw)+Scal*u_exp_b_screw(2*points_b_screw-1),...
%             coordy_b_screw(points_b_screw)+Scal*u_exp_b_screw(2*points_b_screw),'r*',...
%             val2,L2y_sort,'r');
%         axis equal
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         hold off
        
%         figure('name','Reference and deformed mesh')
%         fig1 = triplot(TRI_a_screw,coordx_a_screw,coordy_a_screw,'r');
%         hold on
%         fig2 = triplot(TRI_b_screw,coordx_b_screw,coordy_b_screw,'r');
%         fig3 = triplot(TRI_a_screw,coordx_a_screw+Scal*u_exp_a_screw(1:2:end),...
%             coordy_a_screw+Scal*u_exp_a_screw(2:2:end),'k');
%         fig4 = triplot(TRI_b_screw,coordx_b_screw+Scal*u_exp_b_screw(1:2:end),...
%             coordy_b_screw+Scal*u_exp_b_screw(2:2:end),'k');
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
    numSamplea = ['D' num2str(j) 'a'];
    numSampleb = ['D' num2str(j) 'b'];
    F = appliedLoad(numSample);
    numImages = length(F);
    var_angle = zeros(numImages,1);
    
    for k=1:numImages
    % for k=3
        
        numImage = num2str(k,'%02d');
        filenameDICa = [numSamplea '_00-' numImage '-Mesh'];
        filenameDICb = [numSampleb '_00-' numImage '-Mesh'];
        
        clear Mesh Job U
        load(fullfile(pathnameDIC,filenameDICa));
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx_a_dowel = (X+Job.ROI(1)-1);
        Coordy_a_dowel = (Y+Job.ROI(2)-1);
        TRI_a_dowel = Mesh.TRI;
        scaleFactor_1a = 30/(max(Coordx_a_dowel)-min(Coordx_a_dowel));
        scaleFactor_2a = 15/(max(Coordy_a_dowel)-min(Coordy_a_dowel));
        Ux_a_dowel = U(1:2:end);
        Uy_a_dowel = U(2:2:end);
        
        clear Mesh Job U
        load(fullfile(pathnameDIC,filenameDICb));
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx_b_dowel = (X+Job.ROI(1)-1);
        Coordy_b_dowel = (Y+Job.ROI(2)-1);
        TRI_b_dowel = Mesh.TRI;
        scaleFactor_1b = 15/(max(Coordx_a_dowel)-min(Coordx_a_dowel));
        scaleFactor_2b = 45/(max(Coordy_a_dowel)-min(Coordy_a_dowel)); 
        Ux_b_dowel = U(1:2:end);
        Uy_b_dowel = U(2:2:end);
        
        scaleFactor = mean([scaleFactor_1a scaleFactor_2a...
            scaleFactor_1b scaleFactor_2b]);
        
        u_exp_a_dowel = [Uy_a_dowel -Ux_a_dowel]'*scaleFactor;
        u_exp_a_dowel = u_exp_a_dowel(:);
        u_exp_b_dowel = [Uy_b_dowel -Ux_b_dowel]'*scaleFactor;
        u_exp_b_dowel = u_exp_b_dowel(:);
        
        coordx_a_dowel = Coordy_a_dowel*scaleFactor;
        coordy_a_dowel = -Coordx_a_dowel*scaleFactor;
        min_coordx_a_dowel = min(coordx_a_dowel);
        min_coordy_a_dowel = min(coordy_a_dowel);
        coordx_a_dowel(:) = coordx_a_dowel(:) - min_coordx_a_dowel;
        coordy_a_dowel(:) = coordy_a_dowel(:) - min_coordy_a_dowel;
        coord_a_dowel = [coordx_a_dowel coordy_a_dowel];
        
        coordx_b_dowel = Coordy_b_dowel*scaleFactor;
        coordy_b_dowel = -Coordx_b_dowel*scaleFactor;
        coordx_b_dowel(:) = coordx_b_dowel(:) - min_coordx_a_dowel;
        coordy_b_dowel(:) = coordy_b_dowel(:) - min_coordy_a_dowel;
        coord_b_dowel = [coordx_b_dowel coordy_b_dowel];
        
%         node_a_dowel = NODE(coord_a_dowel,1:size(coord_a_dowel,1));
%         node_b_dowel = NODE(coord_b_dowel,1:size(coord_b_dowel,1));
%         elem_a_dowel = TRI_a_dowel;
%         elem_b_dowel = TRI_b_dowel;
%         elemtype = 'TRI3';
%         S_a_dowel = MODEL('PLAN');
%         S_b_dowel = MODEL('PLAN');
%         S_a_dowel = addnode(S_a_dowel,node_a_dowel);
%         S_b_dowel = addnode(S_b_dowel,node_b_dowel);
%         S_a_dowel = addelem(S_a_dowel,elemtype,elem_a_dowel);
%         S_b_dowel = addelem(S_b_dowel,elemtype,elem_b_dowel);
%         S_a_dowel = final(S_a_dowel);
%         S_b_dowel = final(S_b_dowel);
%         numnode_bound_a_dowel = getnumber(getnode(create_boundary(S_a_dowel)));
%         numnode_bound_b_dowel = getnumber(getnode(create_boundary(S_b_dowel)));
%         figure
%         plot(create_boundary(S_a_dowel));
%         hold on
%         plot(create_boundary(S_b_dowel));
%         plot(coordx_a_dowel(numnode_bound_a_dowel),coordy_a_dowel(numnode_bound_a_dowel),'r*')
%         plot(coordx_b_dowel(numnode_bound_b_dowel),coordy_b_dowel(numnode_bound_b_dowel),'k*')
%         hold off
        
        points_a_dowel = find(coordy_a_dowel>max(coordy_a_dowel)-Mesh.CharLength*scaleFactor/3);
        points_b_dowel = find(coordy_b_dowel<min(coordy_b_dowel)+Mesh.CharLength*scaleFactor/3 &...
            coordx_b_dowel<max(coordx_a_dowel));
        
        % primary line1 and line2
        L1x0 = coordx_a_dowel(points_a_dowel);
        L1y0 = coordy_a_dowel(points_a_dowel);
        [L1x0_sort,index] = sort(L1x0);
        L1y0_sort = L1y0(index);
        fit10 = polyfit(L1x0_sort,L1y0_sort,1);
        val10 = polyval(fit10,L1x0_sort);
        
        L2x0 = coordx_b_dowel(points_b_dowel);
        L2y0 = coordy_b_dowel(points_b_dowel);
        [L2x0_sort,index] = sort(L2x0);
        L2y0_sort = L2y0(index);
        fit20 = polyfit(L2x0_sort,L2y0_sort,1);
        val20 = polyval(fit20,L2x0_sort);
        
        % deformed line1 and lne2
        L1x = coordx_a_dowel(points_a_dowel)+u_exp_a_dowel(2*points_a_dowel-1);
        L1y = coordy_a_dowel(points_a_dowel)+u_exp_a_dowel(2*points_a_dowel);
        [L1x_sort,index] = sort(L1x);
        L1y_sort = L1y(index);
        fit1 = polyfit(L1x_sort,L1y_sort,1);
        val1 = polyval(fit1,L1x_sort);
        
        L2x = coordx_b_dowel(points_b_dowel)+u_exp_b_dowel(2*points_b_dowel-1);
        L2y = coordy_b_dowel(points_b_dowel)+u_exp_b_dowel(2*points_b_dowel);
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
%         fig1 = triplot(TRI_a_dowel,coordx_a_dowel,coordy_a_dowel,'r');
%         hold on
%         fig2 = triplot(TRI_b_dowel,coordx_b_dowel,coordy_b_dowel,'r');
%         fig3 = plot(L1x0,L1y0,'k*',L1x0_sort,val10,'k');
%         fig4 = plot(L2x0,L2y0,'k*',L2x0_sort,val20,'k');
%         axis equal
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         hold off
        
        L1x = coordx_a_dowel(points_a_dowel)+Scal*u_exp_a_dowel(2*points_a_dowel-1);
        L1y = coordy_a_dowel(points_a_dowel)+Scal*u_exp_a_dowel(2*points_a_dowel);
        [L1x_sort,index] = sort(L1x);
        L1y_sort = L1y(index);
        fit1 = polyfit(L1x_sort,L1y_sort,1);
        val1 = polyval(fit1,L1x_sort);
        
        L2x = coordx_b_dowel(points_b_dowel)+Scal*u_exp_b_dowel(2*points_b_dowel-1);
        L2y = coordy_b_dowel(points_b_dowel)+Scal*u_exp_b_dowel(2*points_b_dowel);
        [L2x_sort,index] = sort(L2x);
        L2y_sort = L2y(index);
        fit2 = polyfit(L2x_sort,L2y_sort,1);
        val2 = polyval(fit2,L2x_sort);
        
%         figure('name','best fit line of deformed mesh')
%         fig1 = triplot(TRI_a_dowel,coordx_a_dowel+Scal*u_exp_a_dowel(1:2:end),...
%             coordy_a_dowel+Scal*u_exp_a_dowel(2:2:end),'k');
%         hold on
%         fig2 = triplot(TRI_b_dowel,coordx_b_dowel+Scal*u_exp_b_dowel(1:2:end),...
%             coordy_b_dowel+Scal*u_exp_b_dowel(2:2:end),'k');
%         fig3 = plot(coordx_a_dowel(points_a_dowel)+Scal*u_exp_a_dowel(2*points_a_dowel-1),...
%             coordy_a_dowel(points_a_dowel)+Scal*u_exp_a_dowel(2*points_a_dowel),'r*',...
%             L1x_sort,val1,'r');
%         fig4 = plot(coordx_b_dowel(points_b_dowel)+Scal*u_exp_b_dowel(2*points_b_dowel-1),...
%             coordy_b_dowel(points_b_dowel)+Scal*u_exp_b_dowel(2*points_b_dowel),'r*',...
%             L2x_sort,val2,'r');
%         axis equal
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         hold off
        
%         figure('name','Reference and deformed mesh')
%         fig1 = triplot(TRI_a_dowel,coordx_a_dowel,coordy_a_dowel,'r');
%         hold on
%         fig2 = triplot(TRI_b_dowel,coordx_b_dowel,coordy_b_dowel,'r');
%         fig3 = triplot(TRI_a_dowel,coordx_a_dowel+Scal*u_exp_a_dowel(1:2:end),...
%             coordy_a_dowel+Scal*u_exp_a_dowel(2:2:end),'k');
%         fig4 = triplot(TRI_b_dowel,coordx_b_dowel+Scal*u_exp_b_dowel(1:2:end),...
%             coordy_b_dowel+Scal*u_exp_b_dowel(2:2:end),'k');
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
