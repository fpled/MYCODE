%% Identification of bending stiffness for dowel junction %%
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
displaySolution = false;

filename = 'data_KD.mat';
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

Scal = 1;
Unitx = '[mm]';
UnitU = '[mm]';

%% Identification
if solveProblem
numDowels = 8;

d = 67.5; % [mm]
b = 113; % [mm]
h = 15; % [mm]
% Bending moment per unit length: ml = F*d/b
% Bending stiffness per unit length: k = ml/angle

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
    
    F = appliedLoad(numSample);
    ml = F*d/b;
    numImages = length(F);
    angle = zeros(numImages,1);
    
    for k=1:numImages
    % for k=[3,numImages]
        
        numImage = num2str(k,'%02d');
        filenameDICa = [numSamplea '_00-' numImage '-Mesh'];
        filenameDICb = [numSampleb '_00-' numImage '-Mesh'];
        
        % clear Job Param Mesh U
        load(fullfile(pathnameDIC,filenameDICa));
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx_a = (X+Job.ROI(1)-1);
        Coordy_a = (Y+Job.ROI(2)-1);
        TRI_a = Mesh.TRI;
        scaleFactor_1a = 30/(max(Coordx_a)-min(Coordx_a));
        scaleFactor_2a = 15/(max(Coordy_a)-min(Coordy_a));
        Ux_a = U(1:2:end);
        Uy_a = U(2:2:end);
        
        % clear Job Param Mesh U
        load(fullfile(pathnameDIC,filenameDICb));
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx_b = (X+Job.ROI(1)-1);
        Coordy_b = (Y+Job.ROI(2)-1);
        TRI_b = Mesh.TRI;
        scaleFactor_1b = 15/(max(Coordx_a)-min(Coordx_a));
        scaleFactor_2b = 45/(max(Coordy_a)-min(Coordy_a));
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
        % plot(coordx_a(numnode_bound_a),coordy_a(numnode_bound_a),'k.')
        % plot(coordx_b(numnode_bound_b),coordy_b(numnode_bound_b),'b.')
        % hold off
        
        points_a = find(coordy_a>max(coordy_a)-Mesh.CharLength*scaleFactor/3);
        points_b = find(coordy_b<min(coordy_b)+Mesh.CharLength*scaleFactor/3 &...
            coordx_b<max(coordx_a));
        
        % primary line1 and line2
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
        % primary and deformed angle of junction
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
        % figure('name',['Sample ' numSample ' - Image ' numImage ': Best fit line of reference mesh'])
        % triplot(TRI_a,coordx_a,coordy_a,'k');
        % hold on
        % triplot(TRI_b,coordx_b,coordy_b,'k');
        % plot(L1x0,L1y0,'b.',L1x0_sort,val10,'b','LineWidth',linewidth);
        % plot(L2x0,L2y0,'b.',L2x0_sort,val20,'b','LineWidth',linewidth);
        % hold off
        % axis equal
        % grid on
        % set(gca,'FontSize',fontsize)
        % xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
        % ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
        % mysaveas(pathname,['best_fit_line_mesh_init_' numSample '_' numImage],formats);
        
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
        % figure('name',['Sample ' numSample ' - Image ' numImage ': Best fit line of deformed mesh'])
        % triplot(TRI_a,coordx_a+Scal*u_exp_a(1:2:end),...
        %     coordy_a+Scal*u_exp_a(2:2:end),'r');
        % hold on
        % triplot(TRI_b,coordx_b+Scal*u_exp_b(1:2:end),...
        %     coordy_b+Scal*u_exp_b(2:2:end),'r');
        % plot(coordx_a(points_a)+Scal*u_exp_a(2*points_a-1),...
        %     coordy_a(points_a)+Scal*u_exp_a(2*points_a),'b.',...
        %     L1x_sort,val1,'b','LineWidth',linewidth);
        % plot(coordx_b(points_b)+Scal*u_exp_b(2*points_b-1),...
        %     coordy_b(points_b)+Scal*u_exp_b(2*points_b),'b.',...
        %     L2x_sort,val2,'b','LineWidth',linewidth);
        % hold off
        % axis equal
        % grid on
        % set(gca,'FontSize',fontsize)
        % xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
        % ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
        % mysaveas(pathname,['best_fit_line_mesh_deformed_' numSample '_' numImage],formats);
        
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
        % mysaveas(pathname,['meshes_' numSample '_' numImage],formats);
        
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
save(fullfile(pathname,filename),'numDowels',...
    'FD','mlD','angleD','kD','mean_KD_data');
else
%% Load variables
load(fullfile(pathname,filename),'numDowels',...
    'FD','mlD','angleD','kD','mean_KD_data');
end

%% Statistics
fprintf('\nDowel junctions: Bending stiffness per unit length kD\n');
fprintf('mean(kD) = %g kN/rad\n',mean(mean_KD_data));
fprintf('var(kD)  = %g (kN/rad)^2\n',var(mean_KD_data));
fprintf('std(kD)  = %g kN/rad\n',std(mean_KD_data));
fprintf('cv(kD)   = %g\n',std(mean_KD_data)/mean(mean_KD_data));

%% Plot data
if displaySolution
    colors = distinguishable_colors(numDowels);
    
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
