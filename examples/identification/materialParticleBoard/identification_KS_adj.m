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
displayMesh = true;
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
% formats = {'fig','epsc'};
formats = {'fig','epsc','png'};

Scal = 1;
Unitx = '[mm]';
UnitU = '[mm]';

%% Identification
if solveProblem
% geometric dimensions
La = 142.5; % length of vertical column [mm]
Lb = 67.5; % length of horizontal beam [mm]
b = 113; % sample width [mm]
h = 15; % sample thickness [mm]

% Bending moment per unit length: ml = F*Lb/b
% Bending stiffness per unit length: k = ml/angle

numScrews = 16;

FS = cell(numScrews,1);
mlS = cell(numScrews,1);
angleS = cell(numScrews,1);
kS = cell(numScrews,1);
pS = zeros(numScrews,2);
mean_KS_data = zeros(numScrews,1);

for j=1:numScrews
    
    numSample = ['S' num2str(j)];
    numSamplea = ['S' num2str(j) 'a'];
    numSampleb = ['S' num2str(j) 'b'];
    
    time = tic;
    
    F = appliedLoad(numSample); % applied load [N]
    ml = F*Lb/b; % bending moment per unit length [N.m/m]
    
    numImages = length(F);
    angle = zeros(numImages,1);
    
    for k=1:numImages
        
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
            S_a = MODEL('PLANE');
            S_b = MODEL('PLANE');
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
        
        points_a = find(coordx_a>max(coordx_a)-cl_a/3 &...
            coordy_a>=min(coordy_b));
        points_b = find(coordx_b<min(coordx_b)+cl_b/3);
        
        % initial line a and line b
        x0a = coordx_a(points_a);
        y0a = coordy_a(points_a);
        [y0a_sort,index] = sort(y0a);
        x0a_sort = x0a(index);
        p0a = polyfit(y0a_sort,x0a_sort,1);
        v0a = polyval(p0a,y0a_sort);
        
        x0b = coordx_b(points_b);
        y0b = coordy_b(points_b);
        [y0b_sort,index] = sort(y0b);
        x0b_sort = x0b(index);
        p0b = polyfit(y0b_sort,x0b_sort,1);
        v0b = polyval(p0b,y0b_sort);
        
        % deformed line a and line b
        xa = coordx_a(points_a)+u_exp_a(2*points_a-1);
        ya = coordy_a(points_a)+u_exp_a(2*points_a);
        [ya_sort,index] = sort(ya);
        xa_sort = xa(index);
        pa = polyfit(ya_sort,xa_sort,1);
        va = polyval(pa,ya_sort);
        
        xb = coordx_b(points_b)+u_exp_b(2*points_b-1);
        yb = coordy_b(points_b)+u_exp_b(2*points_b);
        [yb_sort,index] = sort(yb);
        xb_sort = xb(index);
        pb = polyfit(yb_sort,xb_sort,1);
        vb = polyval(pb,yb_sort);
        
        %-----------------------------
        % initial and deformed angles of junction
        %-----------------------------
        t0 = [v0a(end)-v0a(1) y0a_sort(end)-y0a_sort(1)];
        s0 = [v0b(end)-v0b(1) y0b_sort(end)-y0b_sort(1)];
        delta0 = acosd(abs(t0*s0')/(norm(t0)*norm(s0)));
        
        t = [va(end)-va(1) ya_sort(end)-ya_sort(1)];
        s = [vb(end)-vb(1) yb_sort(end)-yb_sort(1)];
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
            plot(x0a,y0a,'b.',v0a,y0a_sort,'b','LineWidth',linewidth);
            plot(x0b,y0b,'b.',v0b,y0b_sort,'b','LineWidth',linewidth);
            hold off
            axis image
            grid on
            set(gca,'XLim',[min(coordx_a)-h/6,max(coordx_b)+h/6])
            set(gca,'YLim',[min(coordy_a)-h/6,max(coordy_b)+h/6])
            set(gca,'FontSize',fontsize)
            xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
            ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
            % mysaveas(pathname,['best_fit_line_mesh_init_' numSample '_' numImage],formats);
            mysaveas(pathname,['best_fit_line_mesh_init_' numSample '_00'],formats);
        end
        
        xa = coordx_a(points_a)+Scal*u_exp_a(2*points_a-1);
        ya = coordy_a(points_a)+Scal*u_exp_a(2*points_a);
        [ya_sort,index] = sort(ya);
        xa_sort = xa(index);
        pa = polyfit(ya_sort,xa_sort,1);
        va = polyval(pa,ya_sort);
        
        xb = coordx_b(points_b)+Scal*u_exp_b(2*points_b-1);
        yb = coordy_b(points_b)+Scal*u_exp_b(2*points_b);
        [yb_sort,index] = sort(yb);
        xb_sort = xb(index);
        pb = polyfit(yb_sort,xb_sort,1);
        vb = polyval(pb,yb_sort);
        
        %--------------------------------
        % Best fit line of deformed mesh
        %--------------------------------
        if displayMesh && k==min(6,numImages)
            figure('name',['Sample ' numSample ' - Image ' numImage ': Best fit line of deformed mesh'])
            triplot(TRI_a,coordx_a+Scal*u_exp_a(1:2:end),...
                coordy_a+Scal*u_exp_a(2:2:end),'r');
            hold on
            triplot(TRI_b,coordx_b+Scal*u_exp_b(1:2:end),...
                coordy_b+Scal*u_exp_b(2:2:end),'r');
            plot(xa,ya,'b.',va,ya_sort,'b','LineWidth',linewidth);
            plot(xb,yb,'b.',vb,yb_sort,'b','LineWidth',linewidth);
            hold off
            axis image
            grid on
            set(gca,'XLim',[min(coordx_a)-h/6,max(coordx_b+Scal*u_exp_b(1:2:end))+h/6])
            set(gca,'YLim',[min(coordy_a)-h/6,max(coordy_b)+h/6])
            set(gca,'FontSize',fontsize)
            xlabel(['$y$ ',Unitx],'Interpreter',interpreter)
            ylabel(['$z$ ',Unitx],'Interpreter',interpreter)
            mysaveas(pathname,['best_fit_line_mesh_deformed_' numSample '_' numImage],formats);
        end
        
        %-------------------------------
        % Reference and deformed meshes
        %-------------------------------
        if displayMesh && k==min(6,numImages)
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
            mysaveas(pathname,['meshes_' numSample '_' numImage],formats);
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
    
    % keptImages = 2:numImages-1;
    if j==1 || j==2 || j==4 || j==8 || j==9 || j==13 || j==14 || j==15 || j==16
        keptImages = 2:numImages-1;
    elseif j==3 || j==6 || j==12
        keptImages = 1:numImages-2;
    elseif j==11
        keptImages = 2:numImages-2;
    else
        keptImages = 1:numImages-1;
    end
    % if j==2 || j==3 || j==6 || j==8 || j==9 || j==10 || j==11 || j==12
    %     keptImages = 1:numImages-2;
    % else
    %     keptImages = 1:numImages-1;
    % end
    FS{j} = F(keptImages); % [N]
    mlS{j} = ml(keptImages); % [N.m/m]
    angleS{j} = angle(keptImages); % [deg]
    kS{j} = mlS{j}./deg2rad(angleS{j}); % [N/rad]
    mean_KS_data(j) = mean(kS{j}); % [N/rad]
    pS(j,:) = [mean_KS_data(j) 0];
    % pS(j,:) = polyfit(deg2rad(angleS{j}),mlS{j},1); % linear fit
    % mean_KS_data(j) = pS(j,1); % [N/rad]
end

%% Save variables
save(fullfile(pathname,filename),'numScrews',...
    'FS','mlS','angleS','kS','pS','mean_KS_data');
else
%% Load variables
load(fullfile(pathname,filename),'numScrews',...
    'FS','mlS','angleS','kS','pS','mean_KS_data');
end

%% Statistics
fprintf('\nScrew junctions: Bending stiffness per unit length kS\n');
fprintf('mean(kS) = %g kN/rad\n',mean(mean_KS_data)*1e-3);
fprintf('var(kS)  = %g (kN/rad)^2\n',var(mean_KS_data)*1e-6);
fprintf('std(kS)  = %g kN/rad\n',std(mean_KS_data)*1e-3);
fprintf('cv(kS)   = %g\n',std(mean_KS_data)/mean(mean_KS_data));

%% Plot data
if displaySolution
    colors = distinguishable_colors(numScrews);
    
    for j=1:numScrews
        numSample = ['S' num2str(j)];
        
        figure('name',['Screw junction ' numSample ': Bending moment per unit length w.r.t. variation of angle'])
        clf
        plot(angleS{j},mlS{j},'+b','LineWidth',linewidth)
        hold on
        plot(angleS{j},polyval(pS(j,:),deg2rad(angleS{j})),'-r','LineWidth',linewidth)
        hold off
        grid on
        set(gca,'FontSize',fontsize)
        set(gca,'XLim',[0,inf])
        set(gca,'YLim',[0,inf])
        xlabel('Variation of angle [deg]','Interpreter',interpreter);
        ylabel('Bending moment per unit length [N.m/m]','Interpreter',interpreter);
        %xlabel('Variation d''angle [deg]','Interpreter',interpreter);
        %ylabel('Moment lin\''eique de flexion [N.m/m]','Interpreter',interpreter);
        mysaveas(pathname,['data_moment_angle_' numSample],formats);
        mymatlab2tikz(pathname,['data_moment_angle_' numSample '.tex']);
        
        % figure('name',['Screw junction ' numSample ': Bending stiffness per unit length w.r.t. applied load'])
        % clf
        % plot(FS{j},kS{j}*1e-3,'+b','LineWidth',linewidth)
        % grid on
        % set(gca,'FontSize',fontsize)
        % xlabel('Applied load [N]','Interpreter',interpreter);
        % ylabel('Bending stiffness per unit length [kN/rad]','Interpreter',interpreter);
        % %xlabel('Chargement appliqu\''e [N]','Interpreter',interpreter);
        % %ylabel('Rigidit\''e lin\''eique en flexion [kN/rad]','Interpreter',interpreter);
        % mysaveas(pathname,['data_stiffness_load_' numSample],formats);
        % mymatlab2tikz(pathname,['data_stiffness_load_' numSample '.tex']);
    end
    
    figure('name','Screw junctions: Bending moment per unit length w.r.t. variation of angle')
    clf
    leg = cell(numScrews,1);
    h = gobjects(numScrews,1);
    for j=1:numScrews
        plot(angleS{j},mlS{j},'+','Color',colors(j,:),'LineWidth',linewidth)
        hold on
        h(j) = plot(angleS{j},polyval(pS(j,:),deg2rad(angleS{j})),'LineStyle','-','Color',colors(j,:),'LineWidth',linewidth);
        leg{j} = ['Sample #' num2str(j)];
    end
    grid on
    set(gca,'FontSize',fontsize)
    set(gca,'XLim',[0,inf])
    set(gca,'YLim',[0,inf])
    xlabel('Variation of angle [deg]','Interpreter',interpreter);
    ylabel('Bending moment per unit length [N.m/m]','Interpreter',interpreter);
    %xlabel('Variation d''angle [deg]','Interpreter',interpreter);
    %ylabel('Moment lin\''eique de flexion [N.m/m]','Interpreter',interpreter);
    legend(h(:),leg{:},'Location','NorthEastOutside')
    mysaveas(pathname,'data_moment_angle_S',formats);
    mymatlab2tikz(pathname,'data_moment_angle_S.tex');
    
    % figure('name','Screw junctions: Bending stiffness per unit length w.r.t. applied load')
    % clf
    % leg = cell(numScrews,1);
    % for j=1:numScrews
    %     plot(FS{j},kS{j}*1e-3,'+','Color',colors(j,:),'LineWidth',linewidth)
    %     hold on
    %     leg{j} = ['Sample #' num2str(j)];
    % end
    % grid on
    % set(gca,'FontSize',fontsize)
    % xlabel('Applied load [N]','Interpreter',interpreter);
    % ylabel('Bending stiffness per unit length [kN/rad]','Interpreter',interpreter);
    % %xlabel('Chargement appliqu\''e [N]','Interpreter',interpreter);
    % %ylabel('Rigidit\''e lin\''eique en flexion [kN/rad]','Interpreter',interpreter);
    % legend(leg{:},'Location','NorthEastOutside')
    % mysaveas(pathname,'data_stiffness_load_S',formats);
    % mymatlab2tikz(pathname,'data_stiffness_load_S.tex');
    
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
