%% Identification of bending stiffness for screw junction %%
%%--------------------------------------------------------%%
% Call after Digital Image Correlation (DIC) RT3
% Coordinate system
% CORRELI:    right     down
% MATLAB:     right     up

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = true;

filename = 'data_KS.mat';
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
numScrewTotal = 16;
sampleScrewDeleted = [9 11 16]; % 3 samples to be deleted due to manufacturing defects
sampleNumScrew = setdiff(1:numScrewTotal,sampleScrewDeleted);
numScrew = length(sampleNumScrew);

d = 67.5e-3; % [m]
b = 113e-3; % [m]
% Bending moment per unit length: ml = F*d/b
% Bending stiffness per unit length: k = ml/angle

FS = cell(numScrew,1);
mlS = cell(numScrew,1);
angleS = cell(numScrew,1);
kS = cell(numScrew,1);
mean_KS_data = zeros(numScrew,1);

for i=1:numScrew
% for i=4
    
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
    % for k=numImages
        
        numImage = num2str(k,'%02d');
        filenameDICa = [numSamplea '_00-' numImage '-Mesh'];
        filenameDICb = [numSampleb '_00-' numImage '-Mesh'];
        
        clear Mesh Job U
        load(fullfile(pathnameDIC,filenameDICa));
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx_a = (X+Job.ROI(1)-1);
        Coordy_a = (Y+Job.ROI(2)-1);
        TRI_a = Mesh.TRI;
        scaleFactor_1a = 45/(max(Coordx_a)-min(Coordx_a));
        scaleFactor_2a = 15/(max(Coordy_a)-min(Coordy_a));
        Ux_a = U(1:2:end);
        Uy_a = U(2:2:end);
        
        clear Mesh Job U
        load(fullfile(pathnameDIC,filenameDICb));
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx_b = (X+Job.ROI(1)-1);
        Coordy_b = (Y+Job.ROI(2)-1);
        TRI_b = Mesh.TRI;
        scaleFactor_1b = 15/(max(Coordx_b)-min(Coordx_b));
        scaleFactor_2b = 30/(max(Coordy_b)-min(Coordy_b));
        Ux_b = U(1:2:end);
        Uy_b = U(2:2:end);
        
        scaleFactor = mean([scaleFactor_1a scaleFactor_2a...
            scaleFactor_1b scaleFactor_2b]);
        
        u_exp_a = [Uy_a -Ux_a]'*scaleFactor;
        u_exp_a = u_exp_a(:);
        u_exp_b = [Uy_b -Ux_b]'*scaleFactor;
        u_exp_b = u_exp_b(:);
        
        coordx_a = Coordy_a*scaleFactor;
        coordy_a = -Coordx_a*scaleFactor;
        min_coordx_a = min(coordx_a);
        min_coordy_a = min(coordy_a);
        coordx_a(:) = coordx_a(:) - min_coordx_a;
        coordy_a(:) = coordy_a(:) - min_coordy_a;
        coord_a = [coordx_a coordy_a];
        
        coordx_b = Coordy_b*scaleFactor;
        coordy_b = -Coordx_b*scaleFactor;
        coordx_b(:) = coordx_b(:) - min_coordx_a;
        coordy_b(:) = coordy_b(:) - min_coordy_a;
        coord_b = [coordx_b coordy_b];
        
        % node_a = NODE(coord_a,1:size(coord_a,1));
        % node_b = NODE(coord_b,1:size(coord_b,1));
        % elem_a = TRI_a;
        % elem_b = TRI_b;
        % elemtype = 'TRI3';
        % S_a = MODEL('PLAN');
        % S_b = MODEL('PLAN');
        % S_a = addnode(S_a,node_a);
        % S_b = addnode(S_b,node_b);
        % S_a = addelem(S_a,elemtype,elem_a);
        % S_b = addelem(S_b,elemtype,elem_b);
        % S_a = final(S_a);
        % S_b = final(S_b);
        % numnode_bound_a = getnumber(getnode(create_boundary(S_a)));
        % numnode_bound_b = getnumber(getnode(create_boundary(S_b)));
        % figure
        % plot(create_boundary(S_a));
        % hold on
        % plot(create_boundary(S_b));
        % plot(coordx_a(numnode_bound_a),coordy_a(numnode_bound_a),'r*')
        % plot(coordx_b(numnode_bound_b),coordy_b(numnode_bound_b),'k*')
        % hold off
        
        points_a = find(coordx_a>max(coordx_a)-Mesh.CharLength*scaleFactor/3 &...
            coordy_a>min(coordy_b));
        points_b = find(coordx_b<min(coordx_b)+Mesh.CharLength*scaleFactor/3 &...
            coordy_b>min(coordy_b));
        
        % primary line1 and line2
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
        % primary and deformed angle of junction
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
        % figure('name',['Sample ' numSample ' - Image ' numImage ': Best fit line of reference mesh'])
        % triplot(TRI_a,coordx_a,coordy_a,'k');
        % hold on
        % triplot(TRI_b,coordx_b,coordy_b,'k');
        % plot(L1x0,L1y0,'b.',val10,L1y0_sort,'b','LineWidth',linewidth);
        % plot(L2x0,L2y0,'b.',val20,L2y0_sort,'b','LineWidth',linewidth);
        % hold off
        % axis equal
        % grid on
        % set(gca,'FontSize',fontsize)
        % xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
        % ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
        % mysaveas(pathname,['best_fit_line_mesh_init_' numSample '_image_' numImage],formats,renderer);
        
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
        % figure('name',['Sample ' numSample ' - Image ' numImage ': Best fit line of deformed mesh'])
        % triplot(TRI_a,coordx_a+Scal*u_exp_a(1:2:end),...
        %     coordy_a+Scal*u_exp_a(2:2:end),'r');
        % hold on
        % triplot(TRI_b,coordx_b+Scal*u_exp_b(1:2:end),...
        %     coordy_b+Scal*u_exp_b(2:2:end),'r');
        % plot(coordx_a(points_a)+Scal*u_exp_a(2*points_a-1),...
        %     coordy_a(points_a)+Scal*u_exp_a(2*points_a),'b.',...
        %     val1,L1y_sort,'b','LineWidth',linewidth);
        % plot(coordx_b(points_b)+Scal*u_exp_b(2*points_b-1),...
        %     coordy_b(points_b)+Scal*u_exp_b(2*points_b),'b.',...
        %     val2,L2y_sort,'b','LineWidth',linewidth);
        % hold off
        % axis equal
        % grid on
        % set(gca,'FontSize',fontsize)
        % xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
        % ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
        % mysaveas(pathname,['best_fit_line_mesh_deformed_' numSample '_image_' numImage],formats,renderer);
        
        %-------------------------------
        % Reference and deformed meshes
        %-------------------------------
        % figure('name',['Sample ' numSample ' - Image ' numImage ': Reference and deformed meshes'])
        % triplot(TRI_a,coordx_a,coordy_a,'k');
        % hold on
        % triplot(TRI_b,coordx_b,coordy_b,'k');
        % triplot(TRI_a,coordx_a+Scal*u_exp_a(1:2:end),...
        %     coordy_a+Scal*u_exp_a(2:2:end),'r');
        % triplot(TRI_b,coordx_b+Scal*u_exp_b(1:2:end),...
        %     coordy_b+Scal*u_exp_b(2:2:end),'r');
        % hold off
        % axis equal
        % grid on
        % set(gca,'FontSize',fontsize)
        % xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
        % ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
        % mysaveas(pathname,['meshes_' numSample '_image_' numImage],formats,renderer);
        
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

%% Save variables
save(fullfile(pathname,filename),'numScrew','sampleNumScrew',...
    'FS','mlS','angleS','kS','mean_KS_data');
else
%% Load variables
load(fullfile(pathname,filename),'numScrew','sampleNumScrew',...
    'FS','mlS','angleS','kS','mean_KS_data');
end

%% Plot data
if displaySolution
    color = distinguishable_colors(numScrew);
    
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
    ylabel('Variation of angle [deg]','Interpreter',interpreter);
    %xlabel('Moment lin\''eique de flexion [N.m/m]','Interpreter',interpreter);
    %ylabel('Variation d''angle [deg]','Interpreter',interpreter);
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
end
