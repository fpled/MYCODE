%% Identification of bending stiffness for screw and dowel junctions %%
%%-------------------------------------------------------------------%%
% Call after Digital Image Correlation (DIC) RT3
% Coordinate system
% CORRELI:    right     down
% MATLAB:     right     up
% The load cell capacity of 50kN is too large for the dowel junction which is fragile

% clc
clearvars
close all

%% Input data
solveProblem = true;
displayMesh = false;
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
renderer = 'OpenGL';

Scal = 1;
Unitx = '[mm]';
UnitU = '[mm]';

%% Identification
if solveProblem
% same geometric dimensions for both junction types
La = 142.5; % length of vertical column [mm]
Lb = 67.5; % length of horizontal beam [mm]
b = 113; % sample width [mm]
h = 15; % sample thickness [mm]

% Bending moment per unit length: ml = F*Lb/b
% Bending stiffness per unit length: k = ml/angle

%% Screw junctions
numScrewsTotal = 16;
% sampleScrewDeleted = [9 11 16]; % 3 samples to be deleted due to manufacturing defects
sampleScrewDeleted = []; % 3 samples to be deleted due to manufacturing defects
sampleNumScrew = setdiff(1:numScrewsTotal,sampleScrewDeleted);
numScrews = length(sampleNumScrew);

FS = cell(numScrews,1);
mlS = cell(numScrews,1);
angleS = cell(numScrews,1);
kS = cell(numScrews,1);
mean_KS_data = zeros(numScrews,1);

for i=1:numScrews
% for i=3
    
    j = sampleNumScrew(i);
    
    numSample = ['S' num2str(j)];
    numSamplea = ['S' num2str(j) 'a'];
    numSampleb = ['S' num2str(j) 'b'];
    
    time = tic;
    
    F = appliedLoad(numSample); % applied load [N]
    ml = F*Lb/b; % bending moment per unit length [N.m/m]
    
    numImages = length(F);
    angle = zeros(numImages,1);
    
    for k=1:numImages
    % for k=[4,numImages]
        
        numImage = num2str(k,'%02d');
        filenameDICa = [numSamplea '_00-' numImage '-Mesh'];
        filenameDICb = [numSampleb '_00-' numImage '-Mesh'];
        
        % clear Job Param Mesh U
        load(fullfile(pathnameDIC,filenameDICa));
        Job_a = Job;
        Mesh_a = Mesh;
        TRI_a = Mesh.TRI;
        U_a = U;
        
        % clear Job Param Mesh U
        load(fullfile(pathnameDIC,filenameDICb));
        Job_b = Job;
        Mesh_b = Mesh;
        TRI_b = Mesh.TRI;
        U_b = U;
        
        [u_exp_a,u_exp_b,coord_a,coord_b,cl_a,cl_b] = extractCorreliJunctionScrew(Job_a,Job_b,Mesh_a,Mesh_b,U_a,U_b,h); % [mm]
        coordx_a = coord_a(:,1);
        coordy_a = coord_a(:,2);
        coordx_b = coord_b(:,1);
        coordy_b = coord_b(:,2);
        
        %---------------------------------
        % Boundary nodes of reference mesh
        %---------------------------------
        if displayMesh && k==1
            node_a = NODE(coord_a,1:size(coord_a,1));
            node_b = NODE(coord_b,1:size(coord_b,1));
            elem_a = TRI_a;
            elem_b = TRI_b;
            elemtype = 'TRI3';
            S_a = MODEL('PLAN');
            S_b = MODEL('PLAN');
            S_a = addnode(S_a,node_a);
            S_b = addnode(S_b,node_b);
            S_a = addelem(S_a,elemtype,elem_a);
            S_b = addelem(S_b,elemtype,elem_b);
            S_a = final(S_a);
            S_b = final(S_b);
            numnode_bound_a = getnumber(getnode(create_boundary(S_a)));
            numnode_bound_b = getnumber(getnode(create_boundary(S_b)));
            
            % figure('name',['Sample ' numSample ' - Image ' numImage ': Boundary nodes of reference mesh'])
            figure('name',['Sample ' numSample ' - Image 00: Boundary nodes of reference mesh'])
            plot(create_boundary(S_a));
            hold on
            plot(create_boundary(S_b));
            plot(coordx_a(numnode_bound_a),coordy_a(numnode_bound_a),'k.')
            plot(coordx_b(numnode_bound_b),coordy_b(numnode_bound_b),'b.')
            hold off
            axis image
        end
        
        % points_a = find(coordx_a>max(coordx_a)-cl_a/3 &...
        %     coordy_a>min(coordy_b)); % PhD thesis Zhhou Chen
        % points_b = find(coordx_b<min(coordx_b)+cl_b/3 &...
        %     coordy_b>min(coordy_b)); % PhD thesis Zhhou Chen
        points_a = find(coordx_a>max(coordx_a)-cl_a/3 &...
            coordy_a>min(coordy_b));
        points_b = find(coordx_b<min(coordx_b)+cl_b/3);
        
        % initial line1 and line2
        L1x0 = coordx_a(points_a);
        L1y0 = coordy_a(points_a);
        [L1y0_sort,index] = sort(L1y0);
        L1x0_sort = L1x0(index);
        fit10 = polyfit(L1y0_sort,L1x0_sort,1);
        val10 = polyval(fit10,L1y0_sort);
        
        L2x0 = coordx_b(points_b);
        L2y0 = coordy_b(points_b);
        [L2y0_sort,index] = sort(L2y0);
        L2x0_sort = L2x0(index);
        fit20 = polyfit(L2y0_sort,L2x0_sort,1);
        val20 = polyval(fit20,L2y0_sort);
        
        % deformed line1 and line2
        L1x = coordx_a(points_a)+u_exp_a(2*points_a-1);
        L1y = coordy_a(points_a)+u_exp_a(2*points_a);
        [L1y_sort,index] = sort(L1y);
        L1x_sort = L1x(index);
        fit1 = polyfit(L1y_sort,L1x_sort,1);
        val1 = polyval(fit1,L1y_sort);
        
        L2x = coordx_b(points_b)+u_exp_b(2*points_b-1);
        L2y = coordy_b(points_b)+u_exp_b(2*points_b);
        [L2y_sort,index] = sort(L2y);
        L2x_sort = L2x(index);
        fit2 = polyfit(L2y_sort,L2x_sort,1);
        val2 = polyval(fit2,L2y_sort);
        
        %-----------------------------
        % initial and deformed angles of junction
        %-----------------------------
        t = [val10(end)-val10(1) L1y0_sort(end)-L1y0_sort(1)];
        s = [val20(end)-val20(1) L2y0_sort(end)-L2y0_sort(1)];
        delta0 = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        t = [val1(end)-val1(1) L1y_sort(end)-L1y_sort(1)];
        s = [val2(end)-val2(1) L2y_sort(end)-L2y_sort(1)];
        delta = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        angle(k) = abs(delta0-delta);
        
        %---------------------------------
        % Best fit line of reference mesh
        %---------------------------------
        if displayMesh && k==1
            % figure('name',['Sample ' numSample ' - Image ' numImage ': Best fit line of reference mesh'])
            figure('name',['Sample ' numSample ' - Image 00: Best fit line of reference mesh'])
            triplot(TRI_a,coordx_a,coordy_a,'k');
            hold on
            triplot(TRI_b,coordx_b,coordy_b,'k');
            plot(L1x0,L1y0,'b.',val10,L1y0_sort,'b','LineWidth',linewidth);
            plot(L2x0,L2y0,'b.',val20,L2y0_sort,'b','LineWidth',linewidth);
            hold off
            axis image
            grid on
            set(gca,'XLim',[min(coordx_a)-h/6,max(coordx_b)+h/6])
            set(gca,'YLim',[min(coordy_a)-h/6,max(coordy_b)+h/6])
            set(gca,'FontSize',fontsize)
            xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
            ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
            % mysaveas(pathname,['best_fit_line_mesh_init_' numSample '_' numImage],formats,renderer);
            mysaveas(pathname,['best_fit_line_mesh_init_' numSample '_00'],formats,renderer);
        end
        
        L1x = coordx_a(points_a)+Scal*u_exp_a(2*points_a-1);
        L1y = coordy_a(points_a)+Scal*u_exp_a(2*points_a);
        [L1y_sort,index] = sort(L1y);
        L1x_sort = L1x(index);
        fit1 = polyfit(L1y_sort,L1x_sort,1);
        val1 = polyval(fit1,L1y_sort);
        
        L2x = coordx_b(points_b)+Scal*u_exp_b(2*points_b-1);
        L2y = coordy_b(points_b)+Scal*u_exp_b(2*points_b);
        [L2y_sort,index] = sort(L2y);
        L2x_sort = L2x(index);
        fit2 = polyfit(L2y_sort,L2x_sort,1);
        val2 = polyval(fit2,L2y_sort);
        
        %--------------------------------
        % Best fit line of deformed mesh
        %--------------------------------
        if displayMesh
            figure('name',['Sample ' numSample ' - Image ' numImage ': Best fit line of deformed mesh'])
            triplot(TRI_a,coordx_a+Scal*u_exp_a(1:2:end),...
                coordy_a+Scal*u_exp_a(2:2:end),'r');
            hold on
            triplot(TRI_b,coordx_b+Scal*u_exp_b(1:2:end),...
                coordy_b+Scal*u_exp_b(2:2:end),'r');
            plot(coordx_a(points_a)+Scal*u_exp_a(2*points_a-1),...
                coordy_a(points_a)+Scal*u_exp_a(2*points_a),'b.',...
                val1,L1y_sort,'b','LineWidth',linewidth);
            plot(coordx_b(points_b)+Scal*u_exp_b(2*points_b-1),...
                coordy_b(points_b)+Scal*u_exp_b(2*points_b),'b.',...
                val2,L2y_sort,'b','LineWidth',linewidth);
            hold off
            axis image
            grid on
            set(gca,'XLim',[min(coordx_a)-h/6,max(coordx_b+Scal*u_exp_b(1:2:end))+h/6])
            set(gca,'YLim',[min(coordy_a)-h/6,max(coordy_b)+h/6])
            set(gca,'FontSize',fontsize)
            xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
            ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
            mysaveas(pathname,['best_fit_line_mesh_deformed_' numSample '_' numImage],formats,renderer);
        end
        
        %-------------------------------
        % Reference and deformed meshes
        %-------------------------------
        if displayMesh
            figure('name',['Sample ' numSample ' - Image ' numImage ': Reference and deformed meshes'])
            triplot(TRI_a,coordx_a,coordy_a,'k');
            hold on
            triplot(TRI_b,coordx_b,coordy_b,'k');
            triplot(TRI_a,coordx_a+Scal*u_exp_a(1:2:end),...
                coordy_a+Scal*u_exp_a(2:2:end),'r');
            triplot(TRI_b,coordx_b+Scal*u_exp_b(1:2:end),...
                coordy_b+Scal*u_exp_b(2:2:end),'r');
            hold off
            axis image
            grid on
            set(gca,'XLim',[min(coordx_a)-h/6,max(coordx_b+Scal*u_exp_b(1:2:end))+h/6])
            set(gca,'YLim',[min(coordy_a)-h/6,max(coordy_b)+h/6])
            set(gca,'FontSize',fontsize)
            xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
            ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
            mysaveas(pathname,['meshes_' numSample '_' numImage],formats,renderer);
        end
        
    end
    
    %% Outputs
    fprintf('\n')
    disp('+-----------------------+')
    fprintf('| Sample S%2d            |\n',j)
    disp('+-------+---------------+-----------+------------------+')
    disp('| Load  | Moment p.u.l. |   Angle   | Stiffness p.u.l. |')
    disp('| F [N] |  ml [N.m/m]   | theta [°] |    k [kN/rad]    |')
    disp('+-------+---------------+-----------+------------------+')
    for k=1:numImages
        fprintf('|  %3d  | %13.4f | %9.4f | %16.4f |\n',F(k),ml(k),angle(k),ml(k)./deg2rad(angle(k))*1e-3)
    end
    disp('+-------+---------------+-----------+------------------+')
    
    toc(time)
    
    FS{i} = F; % [N]
    mlS{i} = ml; % [N.m/m]
    angleS{i} = angle; % [deg]
    kS{i} = ml./deg2rad(angle)'; % [N/rad]
    mean_KS_data(i) = mean( ml(2:end-1)./deg2rad(angle(2:end-1))' ); % [N/rad]
end

%% Dowel junctions
numDowels = 8;

FD = cell(numDowels,1);
mlD = cell(numDowels,1);
angleD = cell(numDowels,1);
kD = cell(numDowels,1);
mean_KD_data = zeros(numDowels,1);

for j=1:numDowels
% for j=2
    
    numSample = ['D' num2str(j)];
    numSamplea = ['D' num2str(j) 'a'];
    numSampleb = ['D' num2str(j) 'b'];
    
    time = tic;
    
    F = appliedLoad(numSample); % applied load [N]
    ml = F*Lb/b; % bending moment per unit length [N.m/m]

    numImages = length(F);
    angle = zeros(numImages,1);
    
    for k=1:numImages
    % for k=[3,numImages]
        
        numImage = num2str(k,'%02d');
        filenameDICa = [numSamplea '_00-' numImage '-Mesh'];
        filenameDICb = [numSampleb '_00-' numImage '-Mesh'];
        
        % clear Job Param Mesh U
        load(fullfile(pathnameDIC,filenameDICa));
        [u_exp_a,coord_a,cl_a] = extractCorreliJunctionDowela(Job,Mesh,U,h); % [mm]
        coordx_a = coord_a(:,1);
        coordy_a = coord_a(:,2);
        TRI_a = Mesh.TRI;
        
        % clear Job Param Mesh U
        load(fullfile(pathnameDIC,filenameDICb));
        [u_exp_b,coord_b,cl_b] = extractCorreliJunctionDowelb(Job,Mesh,U,h); % [mm]
        coordx_b = coord_b(:,1);
        coordy_b = coord_b(:,2);
        TRI_b = Mesh.TRI;
        
        %---------------------------------
        % Boundary nodes of reference mesh
        %---------------------------------
        if displayMesh
            node_a = NODE(coord_a,1:size(coord_a,1));
            node_b = NODE(coord_b,1:size(coord_b,1));
            elem_a = TRI_a;
            elem_b = TRI_b;
            elemtype = 'TRI3';
            S_a = MODEL('PLAN');
            S_b = MODEL('PLAN');
            S_a = addnode(S_a,node_a);
            S_b = addnode(S_b,node_b);
            S_a = addelem(S_a,elemtype,elem_a);
            S_b = addelem(S_b,elemtype,elem_b);
            S_a = final(S_a);
            S_b = final(S_b);
            numnode_bound_a = getnumber(getnode(create_boundary(S_a)));
            numnode_bound_b = getnumber(getnode(create_boundary(S_b)));
            
            % figure('name',['Sample ' numSample ' - Image ' numImage ': Boundary nodes of reference mesh'])
            figure('name',['Sample ' numSample ' - Image 00: Boundary nodes of reference mesh'])
            plot(create_boundary(S_a));
            hold on
            plot(create_boundary(S_b));
            plot(coordx_a(numnode_bound_a),coordy_a(numnode_bound_a),'k.')
            plot(coordx_b(numnode_bound_b),coordy_b(numnode_bound_b),'b.')
            hold off
            axis image
        end
        
        points_a = find(coordy_a>=max(coordy_a)-cl_a/3);
        points_b = find(coordy_b<=min(coordy_b)+cl_b/3 &...
            coordx_b<=max(coordx_a));
        
        % initial line1 and line2
        L1x0 = coordx_a(points_a);
        L1y0 = coordy_a(points_a);
        [L1x0_sort,index] = sort(L1x0);
        L1y0_sort = L1y0(index);
        fit10 = polyfit(L1x0_sort,L1y0_sort,1);
        val10 = polyval(fit10,L1x0_sort);
        
        L2x0 = coordx_b(points_b);
        L2y0 = coordy_b(points_b);
        [L2x0_sort,index] = sort(L2x0);
        L2y0_sort = L2y0(index);
        fit20 = polyfit(L2x0_sort,L2y0_sort,1);
        val20 = polyval(fit20,L2x0_sort);
        
        % deformed line1 and line2
        L1x = coordx_a(points_a)+u_exp_a(2*points_a-1);
        L1y = coordy_a(points_a)+u_exp_a(2*points_a);
        [L1x_sort,index] = sort(L1x);
        L1y_sort = L1y(index);
        fit1 = polyfit(L1x_sort,L1y_sort,1);
        val1 = polyval(fit1,L1x_sort);
        
        L2x = coordx_b(points_b)+u_exp_b(2*points_b-1);
        L2y = coordy_b(points_b)+u_exp_b(2*points_b);
        [L2x_sort,index] = sort(L2x);
        L2y_sort = L2y(index);
        fit2 = polyfit(L2x_sort,L2y_sort,1);
        val2 = polyval(fit2,L2x_sort);
        
        %-----------------------------
        % initial and deformed angles of junction
        %-----------------------------
        t = [L1x0_sort(end)-L1x0_sort(1) val10(end)-val10(1)];
        s = [L2x0_sort(end)-L2x0_sort(1) val20(end)-val20(1)];
        delta0 = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        t = [L1x_sort(end)-L1x_sort(1) val1(end)-val1(1)];
        s = [L2x_sort(end)-L2x_sort(1) val2(end)-val2(1)];
        delta = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        angle(k) = abs(delta0-delta);
        
        %---------------------------------
        % Best fit line of reference mesh
        %---------------------------------
        if displayMesh && k==1
            % figure('name',['Sample ' numSample ' - Image ' numImage ': Best fit line of reference mesh'])
            figure('name',['Sample ' numSample ' - Image 00: Best fit line of reference mesh'])
            triplot(TRI_a,coordx_a,coordy_a,'k');
            hold on
            triplot(TRI_b,coordx_b,coordy_b,'k');
            plot(L1x0,L1y0,'b.',L1x0_sort,val10,'b','LineWidth',linewidth);
            plot(L2x0,L2y0,'b.',L2x0_sort,val20,'b','LineWidth',linewidth);
            hold off
            axis image
            grid on
            set(gca,'XLim',[min(coordx_a)-h/6,max(coordx_b)+h/6])
            set(gca,'YLim',[min(coordy_a)-h/6,max(coordy_b)+h/6])
            set(gca,'FontSize',fontsize)
            xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
            ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
            % mysaveas(pathname,['best_fit_line_mesh_init_' numSample '_' numImage],formats,renderer);
            mysaveas(pathname,['best_fit_line_mesh_init_' numSample '_00'],formats,renderer);
        end
        
        L1x = coordx_a(points_a)+Scal*u_exp_a(2*points_a-1);
        L1y = coordy_a(points_a)+Scal*u_exp_a(2*points_a);
        [L1x_sort,index] = sort(L1x);
        L1y_sort = L1y(index);
        fit1 = polyfit(L1x_sort,L1y_sort,1);
        val1 = polyval(fit1,L1x_sort);
        
        L2x = coordx_b(points_b)+Scal*u_exp_b(2*points_b-1);
        L2y = coordy_b(points_b)+Scal*u_exp_b(2*points_b);
        [L2x_sort,index] = sort(L2x);
        L2y_sort = L2y(index);
        fit2 = polyfit(L2x_sort,L2y_sort,1);
        val2 = polyval(fit2,L2x_sort);
        
        %--------------------------------
        % Best fit line of deformed mesh
        %--------------------------------
        if displayMesh
            figure('name',['Sample ' numSample ' - Image ' numImage ': Best fit line of deformed mesh'])
            triplot(TRI_a,coordx_a+Scal*u_exp_a(1:2:end),...
                coordy_a+Scal*u_exp_a(2:2:end),'r');
            hold on
            triplot(TRI_b,coordx_b+Scal*u_exp_b(1:2:end),...
                coordy_b+Scal*u_exp_b(2:2:end),'r');
            plot(coordx_a(points_a)+Scal*u_exp_a(2*points_a-1),...
                coordy_a(points_a)+Scal*u_exp_a(2*points_a),'b.',...
                L1x_sort,val1,'b','LineWidth',linewidth);
            plot(coordx_b(points_b)+Scal*u_exp_b(2*points_b-1),...
                coordy_b(points_b)+Scal*u_exp_b(2*points_b),'b.',...
                L2x_sort,val2,'b','LineWidth',linewidth);
            hold off
            axis image
            grid on
            set(gca,'XLim',[min(coordx_a)-h/6,max(coordx_b+Scal*u_exp_b(1:2:end))+h/6])
            set(gca,'YLim',[min(coordy_a)-h/6,max(coordy_b)+h/6])
            set(gca,'FontSize',fontsize)
            xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
            ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
            mysaveas(pathname,['best_fit_line_mesh_deformed_' numSample '_' numImage],formats,renderer);
        end
        
        %-------------------------------
        % Reference and deformed meshes
        %-------------------------------
        if displayMesh
            figure('name',['Sample ' numSample ' - Image ' numImage ': Reference and deformed meshes'])
            triplot(TRI_a,coordx_a,coordy_a,'k');
            hold on
            triplot(TRI_b,coordx_b,coordy_b,'k');
            triplot(TRI_a,coordx_a+Scal*u_exp_a(1:2:end),...
                coordy_a+Scal*u_exp_a(2:2:end),'r');
            triplot(TRI_b,coordx_b+Scal*u_exp_b(1:2:end),...
                coordy_b+Scal*u_exp_b(2:2:end),'r');
            hold off
            axis image
            grid on
            set(gca,'XLim',[min(coordx_a)-h/6,max(coordx_b+Scal*u_exp_b(1:2:end))+h/6])
            set(gca,'YLim',[min(coordy_a)-h/6,max(coordy_b)+h/6])
            set(gca,'FontSize',fontsize)
            xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
            ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
            mysaveas(pathname,['meshes_' numSample '_' numImage],formats,renderer);
        end
        
    end
    
    %% Outputs
    fprintf('\n')
    disp('+-----------------------+')
    fprintf('| Sample D%2d            |\n',j)
    disp('+-------+---------------+-----------+------------------+')
    disp('| Load  | Moment p.u.l. |   Angle   | Stiffness p.u.l. |')
    disp('| F [N] |  ml [N.m/m]   | theta [°] |    k [kN/rad]    |')
    disp('+-------+---------------+-----------+------------------+')
    for k=1:numImages
        fprintf('|  %3d  | %13.4f | %9.4f | %16.4f |\n',F(k),ml(k),angle(k),ml(k)./deg2rad(angle(k))*1e-3)
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
save(fullfile(pathname,filenameS),'numScrews','sampleNumScrew',...
    'FS','mlS','angleS','kS','mean_KS_data');
save(fullfile(pathname,filenameD),'numDowels',...
    'FD','mlD','angleD','kD','mean_KD_data');
else
%% Load variables
load(fullfile(pathname,filenameS),'numScrews','sampleNumScrew',...
    'FS','mlS','angleS','kS','mean_KS_data');
load(fullfile(pathname,filenameD),'numDowels',...
    'FD','mlD','angleD','kD','mean_KD_data');
end

%% Plot data
if displaySolution
    colors = distinguishable_colors(max(numScrews,numDowels));
    
    for i=1:numScrews
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
        ylabel('Variation of angle [deg]','Interpreter',interpreter);
        %xlabel('Moment lin\''eique de flexion [N.m/m]','Interpreter',interpreter);
        %ylabel('Variation d''angle [deg]','Interpreter',interpreter);
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
    
    for j=1:numDowels
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
        ylabel('Variation of angle [deg]','Interpreter',interpreter);
        %xlabel('Moment lin\''eique de flexion [N.m/m]','Interpreter',interpreter);
        %ylabel('Variation d''angle [deg]','Interpreter',interpreter);
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
    leg = cell(numScrews,1);
    h = gobjects(numScrews,1);
    for i=1:numScrews
        plot(mlS{i},angleS{i},'+','Color',colors(i,:),'LineWidth',linewidth)
        hold on
        h(i) = plot(mlS{i},rad2deg(mlS{i}./mean_KS_data(i)),'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth);
        leg{i} = ['Sample #' num2str(i)];
    end
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Bending moment per unit length [N.m/m]','Interpreter',interpreter);
    ylabel('Variation of angle [deg]','Interpreter',interpreter);
    %xlabel('Moment lin\''eique de flexion [N.m/m]','Interpreter',interpreter);
    %ylabel('Variation d''angle [deg]','Interpreter',interpreter);
    legend(h(:),leg{:},'Location','NorthEastOutside')
    mysaveas(pathname,'data_angle_moment_S',formats);
    mymatlab2tikz(pathname,'data_angle_moment_S.tex');
    
    figure('name','Dowel junctions: Variation of angle w.r.t. bending moment per unit length')
    clf
    leg = cell(numDowels,1);
    h = gobjects(numDowels,1);
    for j=1:numDowels
        plot(mlD{j},angleD{j},'+','Color',colors(j,:),'LineWidth',linewidth)
        hold on
        h(j) = plot(mlD{j},rad2deg(mlD{j}./mean_KD_data(j)),'LineStyle','-','Color',colors(j,:),'LineWidth',linewidth);
        leg{j} = ['Sample #' num2str(j)];
    end
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Bending moment per unit length [N.m/m]','Interpreter',interpreter);
    ylabel('Variation of angle [deg]','Interpreter',interpreter);
    %xlabel('Moment lin\''eique de flexion [N.m/m]','Interpreter',interpreter);
    %ylabel('Variation d''angle [deg]','Interpreter',interpreter);
    legend(h(:),leg{:},'Location','NorthEastOutside')
    mysaveas(pathname,'data_angle_moment_D',formats);
    mymatlab2tikz(pathname,'data_angle_moment_D.tex');
    
    figure('name','Screw junctions: Bending stiffness per unit length w.r.t. applied load')
    clf
    leg = cell(numScrews,1);
    for i=1:numScrews
        plot(FS{i},kS{i}*1e-3,'+','Color',colors(i,:),'LineWidth',linewidth)
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
    
    figure('name','Dowel junctions: Bending stiffness per unit length w.r.t. applied load')
    clf
    leg = cell(numDowels,1);
    for j=1:numDowels
        plot(FD{j},kD{j}*1e-3,'+','Color',colors(j,:),'LineWidth',linewidth)
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
    ylabel('Bending stiffness per unit length $k_S$ [kN/rad]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Rigidit\''e lin\''eique en flexion $k_S$ [kN/rad]','Interpreter',interpreter);
    mysaveas(pathname,'data_KS',formats);
    mymatlab2tikz(pathname,'data_KS.tex');
    
    figure('Name','Dowel junctions: Bending stiffness per unit length')
    clf
    bar(mean_KD_data*1e-3);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Bending stiffness per unit length $k_D$ [kN/rad]','Interpreter',interpreter);
    %xlabel('Num\''ero d''\''echantillon','Interpreter',interpreter);
    %ylabel('Rigidit\''e lin\''eique en flexion $k_D$ [kN/rad]','Interpreter',interpreter);
    mysaveas(pathname,'data_KD',formats);
    mymatlab2tikz(pathname,'data_KD.tex');
end
