%% Identification of bending stiffness for screw and dowel junctions %%
%%-------------------------------------------------------------------%%
% Call after Digital Image Correlation (DIC) RT3
% Coordinate system
% CORRELI:    right     down
% MATLAB:     right     up
% the cell of 50kN is too large for the Dowel junction which is fragile

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

filenameS = 'data_KS.mat';
filenameD = 'data_KD.mat';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
if ~exist(pathname,'dir')
    mkdir(pathname);
end
pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','materialParticleBoard','resultsDIC');

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};

%% Identification
if solveProblem
numScrewTotal = 16;
sampleScrewDeleted = [9 11 16]; % 3 samples to be deleted due to manufacturing defects
sampleNumScrew = setdiff(1:numScrewTotal,sampleScrewDeleted);
numScrew = length(sampleNumScrew);
numDowel = 8;

d = 67.5e-3; % [m]
b = 113e-3; % [m]
% Same Dimensions for the two assembly junctions
% Bending moment per unit length: ml = F*d/b
% Bending stiffness per unit length: k = ml/angle

FS = cell(numScrew,1);
mlS = cell(numScrew,1);
angleS = cell(numScrew,1);
kS = cell(numScrew,1);
mean_KS_data = zeros(numScrew,1);

for i=1:numScrew
% for i=6
    
    j = sampleNumScrew(i);
    
    numSample = ['S' num2str(j)];
    numSamplea = ['S' num2str(j) 'a'];
    numSampleb = ['S' num2str(j) 'b'];
    
    time = tic;
    
    F = appliedLoad(numSample);
    ml = F*d/b;
    numImages = length(F);
    angle = zeros(numImages,1);
    
    for k=1:numImages
    % for k=3
        
        numImage = num2str(k,'%02d');
        filenameDICa = [numSamplea '_00-' numImage '-Mesh'];
        filenameDICb = [numSampleb '_00-' numImage '-Mesh'];
        
        clear Mesh Job U
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
        
        angle(k) = abs(delta0-delta);
        
        %-----------------------------
        % Reference and deformed mesh
        %-----------------------------
        Scal = 1;
        Unitx = '[mm]';
        UnitU = '[mm]';
        
%         figure('name','best fit line of initial mesh')
%         triplot(TRI_a_screw,coordx_a_screw,coordy_a_screw,'r');
%         hold on
%         triplot(TRI_b_screw,coordx_b_screw,coordy_b_screw,'r');
%         plot(L1x0,L1y0,'k*',val10,L1y0_sort,'k');
%         plot(L2x0,L2y0,'k*',val20,L2y0_sort,'k');
%         hold off
%         axis equal
%         set(gca,'FontSize',fontsize)
%         xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
%         ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
%         mysaveas(pathname,['meshes_' numSample '_image_' numImage],formats);
        
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
%         triplot(TRI_a_screw,coordx_a_screw+Scal*u_exp_a_screw(1:2:end),...
%             coordy_a_screw+Scal*u_exp_a_screw(2:2:end),'k');
%         hold on
%         triplot(TRI_b_screw,coordx_b_screw+Scal*u_exp_b_screw(1:2:end),...
%             coordy_b_screw+Scal*u_exp_b_screw(2:2:end),'k');
%         plot(coordx_a_screw(points_a_screw)+Scal*u_exp_a_screw(2*points_a_screw-1),...
%             coordy_a_screw(points_a_screw)+Scal*u_exp_a_screw(2*points_a_screw),'r*',...
%             val1,L1y_sort,'r');
%         plot(coordx_b_screw(points_b_screw)+Scal*u_exp_b_screw(2*points_b_screw-1),...
%             coordy_b_screw(points_b_screw)+Scal*u_exp_b_screw(2*points_b_screw),'r*',...
%             val2,L2y_sort,'r');
%         hold off
%         axis equal
%         set(gca,'FontSize',fontsize)
%         xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
%         ylabel(['$z$ ',Unitx],'Interpreter',interpreter
%         mysaveas(pathname,['meshes_' numSample '_image_' numImage],formats);
        
        %-----------------------------
        % Reference and deformed mesh
        %-----------------------------
%         figure('name','Reference and deformed mesh')
%         triplot(TRI_a_screw,coordx_a_screw,coordy_a_screw,'r');
%         hold on
%         triplot(TRI_b_screw,coordx_b_screw,coordy_b_screw,'r');
%         triplot(TRI_a_screw,coordx_a_screw+Scal*u_exp_a_screw(1:2:end),...
%             coordy_a_screw+Scal*u_exp_a_screw(2:2:end),'k');
%         triplot(TRI_b_screw,coordx_b_screw+Scal*u_exp_b_screw(1:2:end),...
%             coordy_b_screw+Scal*u_exp_b_screw(2:2:end),'k');
%         hold off
%         axis equal
%         set(gca,'FontSize',fontsize)
%         xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
%         ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
%         mysaveas(pathname,['meshes_' numSample '_image_' numImage],formats);
        
    end
    
    %% Outputs
    fprintf('\n')
    disp('+-----------------------+')
    fprintf('| Sample S%2d            |\n',j)
    disp('+-------+---------------+-----------+------------------+')
    disp('| Load  | Moment p.u.l. |   Angle   | Stiffness p.u.l. |')
    disp('| F [N] |  ml [N.m/m]   | theta [°] |    k [N/rad]     |')
    disp('+-------+---------------+-----------+------------------+')
    for k=1:numImages
        fprintf('|  %3d  | %13.4f | %9.4f | %16.4e |\n',F(k),ml(k),angle(k),ml(k)./deg2rad(angle(k)))
    end
    disp('+-------+---------------+-----------+------------------+')
    
    toc(time)
    
    FS{i} = F; % [N]
    mlS{i} = ml; % [N.m/m]
    angleS{i} = angle; % [deg]
    kS{i} = ml./deg2rad(angle)'; % [N/rad]
    mean_KS_data(i) = mean( ml(2:end-1)./deg2rad(angle(2:end-1))' ); % [N/rad]
end

FD = cell(numDowel,1);
mlD = cell(numDowel,1);
angleD = cell(numDowel,1);
kD = cell(numDowel,1);
mean_KD_data = zeros(numDowel,1);

for j=1:numDowel
% for j=2
    
    numSample = ['D' num2str(j)];
    numSamplea = ['D' num2str(j) 'a'];
    numSampleb = ['D' num2str(j) 'b'];
    
    time = tic;
    
    F = appliedLoad(numSample);
    ml = F*d/b;
    numImages = length(F);
    angle = zeros(numImages,1);
    
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
        
        % deformed line1 and line2
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
        
        angle(k) = abs(delta0-delta);
        
        %-----------------------------
        % Reference and deformed mesh
        %-----------------------------
        Scal = 1;
        Unitx = '[mm]';
        UnitU = '[mm]';
        
%         figure('name','best fit line of initial mesh')
%         triplot(TRI_a_dowel,coordx_a_dowel,coordy_a_dowel,'r');
%         hold on
%         triplot(TRI_b_dowel,coordx_b_dowel,coordy_b_dowel,'r');
%         plot(L1x0,L1y0,'k*',L1x0_sort,val10,'k');
%         plot(L2x0,L2y0,'k*',L2x0_sort,val20,'k');
%         hold off
%         axis equal
%         set(gca,'FontSize',fontsize)
%         xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
%         ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
%         mysaveas(pathname,['meshes_' numSample '_image_' numImage],formats);
        
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
%         triplot(TRI_a_dowel,coordx_a_dowel+Scal*u_exp_a_dowel(1:2:end),...
%             coordy_a_dowel+Scal*u_exp_a_dowel(2:2:end),'k');
%         hold on
%         triplot(TRI_b_dowel,coordx_b_dowel+Scal*u_exp_b_dowel(1:2:end),...
%             coordy_b_dowel+Scal*u_exp_b_dowel(2:2:end),'k');
%         plot(coordx_a_dowel(points_a_dowel)+Scal*u_exp_a_dowel(2*points_a_dowel-1),...
%             coordy_a_dowel(points_a_dowel)+Scal*u_exp_a_dowel(2*points_a_dowel),'r*',...
%             L1x_sort,val1,'r');
%         plot(coordx_b_dowel(points_b_dowel)+Scal*u_exp_b_dowel(2*points_b_dowel-1),...
%             coordy_b_dowel(points_b_dowel)+Scal*u_exp_b_dowel(2*points_b_dowel),'r*',...
%             L2x_sort,val2,'r');
%         hold off
%         axis equal
%         set(gca,'FontSize',fontsize)
%         xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
%         ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
%         mysaveas(pathname,['meshes_' numSample '_image_' numImage],formats);
        
        %-----------------------------
        % Reference and deformed mesh
        %-----------------------------
%         figure('name','Reference and deformed mesh')
%         triplot(TRI_a_dowel,coordx_a_dowel,coordy_a_dowel,'r');
%         hold on
%         triplot(TRI_b_dowel,coordx_b_dowel,coordy_b_dowel,'r');
%         triplot(TRI_a_dowel,coordx_a_dowel+Scal*u_exp_a_dowel(1:2:end),...
%             coordy_a_dowel+Scal*u_exp_a_dowel(2:2:end),'k');
%         triplot(TRI_b_dowel,coordx_b_dowel+Scal*u_exp_b_dowel(1:2:end),...
%             coordy_b_dowel+Scal*u_exp_b_dowel(2:2:end),'k');
%         hold off
%         axis equal
%         set(gca,'FontSize',fontsize)
%         xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
%         ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
%         mysaveas(pathname,['meshes_' numSample '_image_' numImage],formats);
        
    end
    
    %% Outputs
    fprintf('\n')
    disp('+-----------------------+')
    fprintf('| Sample D%2d            |\n',j)
    disp('+-------+---------------+-----------+------------------+')
    disp('| Load  | Moment p.u.l. |   Angle   | Stiffness p.u.l. |')
    disp('| F [N] |  ml [N.m/m]   | theta [°] |    k [N/rad]     |')
    disp('+-------+---------------+-----------+------------------+')
    for k=1:numImages
        fprintf('|  %3d  | %13.4f | %9.4f | %16.4e |\n',F(k),ml(k),angle(k),ml(k)./deg2rad(angle(k)))
    end
    disp('+-------+---------------+-----------+------------------+')
    
    toc(time)
    
    FD{j} = F; % [N]
    mlD{j} = ml; % [N.m/m]
    angleD{j} = angle; % [deg]
    kD{j} = ml./deg2rad(angle)'; % [N/rad]
    mean_KD_data(j) = mean(kD{j}); % [N/rad]
end

%% Save variables
save(fullfile(pathname,filenameS),'numScrew','sampleNumScrew',...
    'FS','mlS','angleS','kS','mean_KS_data');
save(fullfile(pathname,filenameD),'numDowel',...
    'FD','mlD','angleD','kD','mean_KD_data');
else
%% Load variables
load(fullfile(pathname,filenameS),'numScrew','sampleNumScrew',...
    'FS','mlS','angleS','kS','mean_KS_data');
load(fullfile(pathname,filenameD),'numDowel',...
    'FD','mlD','angleD','kD','mean_KD_data');
end

%% Plot data
if displaySolution
    color = distinguishable_colors(max(numScrew,numDowel));
    
    for i=1:numScrew
        j = sampleNumScrew(i);
        numSample = ['S' num2str(j)];
        
        figure('name',['Screw junction ' numSample ': Variation of angle w.r.t. bending moment per unit length'])
        clf
        plot(mlS{i},angleS{i},'+b','LineWidth',linewidth)
        hold on
        plot(mlS{i},rad2deg(mlS{i}./mean_KS_data(i)),'-r','LineWidth',linewidth)
        hold off
        grid on
        set(gca,'FontSize',fontsize)
        xlabel('Bending moment per unit length [N.m/m]','Interpreter',interpreter);
        ylabel('Variation of angle [$^{\circ}$]','Interpreter',interpreter);
        %xlabel('Moment lin\''eique de flexion [N.m/m]','Interpreter',interpreter);
        %ylabel('Variation d''angle [$^{\circ}$]','Interpreter',interpreter);
        mysaveas(pathname,['data_angle_moment_' numSample],formats);
        mymatlab2tikz(pathname,['data_angle_moment_' numSample '.tex']);
        
        figure('name',['Screw junction ' numSample ': Bending stiffness per unit length w.r.t. applied load'])
        clf
        plot(FS{i},kS{i}*1e-3,'+b','LineWidth',linewidth)
        grid on
        set(gca,'FontSize',fontsize)
        xlabel('Applied load [N]','Interpreter',interpreter);
        ylabel('Bending stiffness per unit length [kN/rad]','Interpreter',interpreter);
        %xlabel('Chargement appliqu\''e [N]','Interpreter',interpreter);
        %ylabel('Rigidit\''e lin\''eique en flexion [kN/rad]','Interpreter',interpreter);
        mysaveas(pathname,['data_stiffness_load_' numSample],formats);
        mymatlab2tikz(pathname,['data_stiffness_load_' numSample '.tex']);
    end
    
    for j=1:numDowel
        numSample = ['D' num2str(j)];
        
        figure('name',['Dowel junction ' numSample ': Variation of angle w.r.t. bending moment per unit length'])
        clf
        plot(mlD{j},angleD{j},'+b','LineWidth',linewidth)
        hold on
        plot(mlD{j},rad2deg(mlD{j}./mean_KD_data(j)),'-r','LineWidth',linewidth)
        hold off
        grid on
        set(gca,'FontSize',fontsize)
        xlabel('Bending moment per unit length [N.m/m]','Interpreter',interpreter);
        ylabel('Variation of angle [$^{\circ}$]','Interpreter',interpreter);
        %xlabel('Moment lin\''eique de flexion [N.m/m]','Interpreter',interpreter);
        %ylabel('Variation d''angle [$^{\circ}$]','Interpreter',interpreter);
        mysaveas(pathname,['data_angle_moment_' numSample],formats);
        mymatlab2tikz(pathname,['data_angle_moment_' numSample '.tex']);
        
        figure('name',['Dowel junction ' numSample ': Bending stiffness per unit length w.r.t. applied load'])
        clf
        plot(FD{j},kD{j}*1e-3,'+b','LineWidth',linewidth)
        grid on
        set(gca,'FontSize',fontsize)
        xlabel('Applied load [N]','Interpreter',interpreter);
        ylabel('Bending stiffness per unit length [kN/rad]','Interpreter',interpreter);
        %xlabel('Chargement appliqu\''e [N]','Interpreter',interpreter);
        %ylabel('Rigidit\''e lin\''eique en flexion [kN/rad]','Interpreter',interpreter);
        mysaveas(pathname,['data_stiffness_load_' numSample],formats);
        mymatlab2tikz(pathname,['data_stiffness_load_' numSample '.tex']);
    end
    
    figure('name','Screw junctions: Variation of angle w.r.t. bending moment per unit length')
    clf
    leg = cell(numScrew,1);
    h = gobjects(numScrew,1);
    for i=1:numScrew
        plot(mlS{i},angleS{i},'+','Color',color(i,:),'LineWidth',linewidth)
        hold on
        h(i) = plot(mlS{i},rad2deg(mlS{i}./mean_KS_data(i)),'LineStyle','-','Color',color(i,:),'LineWidth',linewidth);
        leg{i} = ['Sample #' num2str(i)];
    end
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Bending moment per unit length [N.m/m]','Interpreter',interpreter);
    ylabel('Variation of angle [$^{\circ}$]','Interpreter',interpreter);
    %xlabel('Moment lin\''eique de flexion [N.m/m]','Interpreter',interpreter);
    %ylabel('Variation d''angle [$^{\circ}$]','Interpreter',interpreter);
    legend(h(:),leg{:},'Location','NorthEastOutside')
    mysaveas(pathname,'data_angle_moment_S',formats);
    mymatlab2tikz(pathname,'data_angle_moment_S.tex');
    
    figure('name','Screw junctions: Bending stiffness per unit length w.r.t. applied load')
    clf
    leg = cell(numScrew,1);
    for i=1:numScrew
        plot(FS{i},kS{i}*1e-3,'+','Color',color(i,:),'LineWidth',linewidth)
        hold on
        leg{i} = ['Sample #' num2str(i)];
    end
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Applied load [N]','Interpreter',interpreter);
    ylabel('Bending stiffness per unit length [kN/rad]','Interpreter',interpreter);
    %xlabel('Chargement appliqu\''e [N]','Interpreter',interpreter);
    %ylabel('Rigidit\''e lin\''eique en flexion [kN/rad]','Interpreter',interpreter);
    legend(leg{:},'Location','NorthEastOutside')
    mysaveas(pathname,'data_stiffness_load_S',formats);
    mymatlab2tikz(pathname,'data_stiffness_load_S.tex');
    
    figure('name','Dowel junctions: Variation of angle w.r.t. bending moment per unit length')
    clf
    leg = cell(numDowel,1);
    h = gobjects(numDowel,1);
    for j=1:numDowel
        plot(mlD{j},angleD{j},'+','Color',color(j,:),'LineWidth',linewidth)
        hold on
        h(j) = plot(mlD{j},rad2deg(mlD{j}./mean_KD_data(j)),'LineStyle','-','Color',color(j,:),'LineWidth',linewidth);
        leg{j} = ['Sample #' num2str(j)];
    end
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Bending moment per unit length [N.m/m]','Interpreter',interpreter);
    ylabel('Variation of angle [$^{\circ}$]','Interpreter',interpreter);
    %xlabel('Moment lin\''eique de flexion [N.m/m]','Interpreter',interpreter);
    %ylabel('Variation d''angle [$^{\circ}$]','Interpreter',interpreter);
    legend(h(:),leg{:},'Location','NorthEastOutside')
    mysaveas(pathname,'data_angle_moment_D',formats);
    mymatlab2tikz(pathname,'data_angle_moment_D.tex');
    
    figure('name','Dowel junctions: Bending stiffness per unit length w.r.t. applied load')
    clf
    leg = cell(numDowel,1);
    for j=1:numDowel
        plot(FD{j},kD{j}*1e-3,'+','Color',color(j,:),'LineWidth',linewidth)
        hold on
        leg{j} = ['Sample #' num2str(j)];
    end
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Applied load [N]','Interpreter',interpreter);
    ylabel('Bending stiffness per unit length [kN/rad]','Interpreter',interpreter);
    %xlabel('Chargement appliqu\''e [N]','Interpreter',interpreter);
    %ylabel('Rigidit\''e lin\''eique en flexion [kN/rad]','Interpreter',interpreter);
    legend(leg{:},'Location','NorthEastOutside')
    mysaveas(pathname,'data_stiffness_load_D',formats);
    mymatlab2tikz(pathname,'data_stiffness_load_D.tex');

    figure('Name','Screw junctions: Bending stiffness per unit length')
    clf
    bar(mean_KS_data*1e-3);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Bending stiffness per unit length $k^S$ [kN/rad]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Rigidit\''e lin\''eique en flexion $k^S$ [kN/rad]','Interpreter',interpreter);
    mysaveas(pathname,'data_KS',formats);
    mymatlab2tikz(pathname,'data_KS.tex');
    
    figure('Name','Dowel junctions: Bending stiffness per unit length')
    clf
    bar(mean_KD_data*1e-3);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Bending stiffness per unit length $k^D$ [kN/rad]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Rigidit\''e lin\''eique en flexion $k^D$ [kN/rad]','Interpreter',interpreter);
    mysaveas(pathname,'data_KD',formats);
    mymatlab2tikz(pathname,'data_KD.tex');
end
