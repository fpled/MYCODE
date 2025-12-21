%% Identification of bending stiffness for dowel junction %%
%%--------------------------------------------------------%%
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
displayMesh = true;
displaySolution = true;

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

numDowels = 8;

FD = cell(numDowels,1);
mlD = cell(numDowels,1);
angleD = cell(numDowels,1);
kD = cell(numDowels,1);
pD = zeros(numDowels,2);
mean_KD_data = zeros(numDowels,1);

for j=1:numDowels
    
    numSample = ['D' num2str(j)];
    numSamplea = ['D' num2str(j) 'a'];
    numSampleb = ['D' num2str(j) 'b'];
    
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
        
        [u_exp_a,u_exp_b,coord_a,coord_b,cl_a,cl_b] = extractCorreliJunctionDowel(Job_a,Job_b,Mesh_a,Mesh_b,U_a,U_b,h); % [mm]
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
        
        points_a = find(coordy_a>max(coordy_a)-cl_a/3);
        points_b = find(coordy_b<min(coordy_b)+cl_b/3 &...
            coordx_b<=max(coordx_a));
        
        % initial line a and line b
        x0a = coordx_a(points_a);
        y0a = coordy_a(points_a);
        [x0a_sort,index] = sort(x0a);
        y0a_sort = y0a(index);
        p0a = polyfit(x0a_sort,y0a_sort,1);
        v0a = polyval(p0a,x0a_sort);
        
        x0b = coordx_b(points_b);
        y0b = coordy_b(points_b);
        [x0b_sort,index] = sort(x0b);
        y0b_sort = y0b(index);
        p0b = polyfit(x0b_sort,y0b_sort,1);
        v0b = polyval(p0b,x0b_sort);
        
        % deformed line a and line b
        xa = coordx_a(points_a)+u_exp_a(2*points_a-1);
        ya = coordy_a(points_a)+u_exp_a(2*points_a);
        [xa_sort,index] = sort(xa);
        ya_sort = ya(index);
        pa = polyfit(xa_sort,ya_sort,1);
        va = polyval(pa,xa_sort);
        
        xb = coordx_b(points_b)+u_exp_b(2*points_b-1);
        yb = coordy_b(points_b)+u_exp_b(2*points_b);
        [xb_sort,index] = sort(xb);
        yb_sort = yb(index);
        pb = polyfit(xb_sort,yb_sort,1);
        vb = polyval(pb,xb_sort);
        
        %-----------------------------
        % initial and deformed angles of junction
        %-----------------------------
        t0 = [x0a_sort(end)-x0a_sort(1) v0a(end)-v0a(1)];
        s0 = [x0b_sort(end)-x0b_sort(1) v0b(end)-v0b(1)];
        delta0 = acosd(abs(t0*s0')/(norm(t0)*norm(s0)));
        
        t = [xa_sort(end)-xa_sort(1) va(end)-va(1)];
        s = [xb_sort(end)-xb_sort(1) vb(end)-vb(1)];
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
            plot(x0a,y0a,'b.',x0a_sort,v0a,'b','LineWidth',linewidth);
            plot(x0b,y0b,'b.',x0b_sort,v0b,'b','LineWidth',linewidth);
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
        [xa_sort,index] = sort(xa);
        ya_sort = ya(index);
        pa = polyfit(xa_sort,ya_sort,1);
        va = polyval(pa,xa_sort);
        
        xb = coordx_b(points_b)+Scal*u_exp_b(2*points_b-1);
        yb = coordy_b(points_b)+Scal*u_exp_b(2*points_b);
        [xb_sort,index] = sort(xb);
        yb_sort = yb(index);
        pb = polyfit(xb_sort,yb_sort,1);
        vb = polyval(pb,xb_sort);
        
        %--------------------------------
        % Best fit line of deformed mesh
        %--------------------------------
        if displayMesh && k==min(3,numImages)
            figure('name',['Sample ' numSample ' - Image ' numImage ': Best fit line of deformed mesh'])
            triplot(TRI_a,coordx_a+Scal*u_exp_a(1:2:end),...
                coordy_a+Scal*u_exp_a(2:2:end),'r');
            hold on
            triplot(TRI_b,coordx_b+Scal*u_exp_b(1:2:end),...
                coordy_b+Scal*u_exp_b(2:2:end),'r');
            plot(xa,ya,'b.',xa_sort,va,'b','LineWidth',linewidth);
            plot(xb,yb,'b.',xb_sort,vb,'b','LineWidth',linewidth);
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
        if displayMesh && k==min(3,numImages)
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
    fprintf('Sample D%d\n',j)
    disp('+-------+---------------+-----------+------------------+')
    disp('| Load  | Moment p.u.l. |   Angle   | Stiffness p.u.l. |')
    disp('| F [N] |  ml [N.m/m]   | theta [°] |    k [kN/rad]    |')
    disp('+-------+---------------+-----------+------------------+')
    for k=1:numImages
        fprintf('|  %3d  | %13.4f | %9.4f | %16.4f |\n',F(k),ml(k),angle(k),ml(k)./deg2rad(angle(k))*1e-3)
    end
    disp('+-------+---------------+-----------+------------------+')
    
    toc(time)
    
    if j==1 || j==3 || j==5 || j==6
        keptImages = 1:numImages-1;
    elseif j==2
        keptImages = 1:numImages-2;
    else
        keptImages = 1:numImages;
    end
    FD{j} = F(keptImages); % [N]
    mlD{j} = ml(keptImages); % [N.m/m]
    angleD{j} = angle(keptImages); % [deg]
    kD{j} = mlD{j}./deg2rad(angleD{j}); % [N/rad]
    mean_KD_data(j) = mean(kD{j}); % [N/rad]
    pD(j,:) = [mean_KD_data(j) 0];
    % pD(j,:) = polyfit(deg2rad(angleD{j}),mlD{j},1); % linear fit
    % mean_KD_data(j) = pD(j,1); % [N/rad]
end

%% Save variables
save(fullfile(pathname,filename),'numDowels',...
    'FD','mlD','angleD','kD','pD','mean_KD_data');
else
%% Load variables
load(fullfile(pathname,filename),'numDowels',...
    'FD','mlD','angleD','kD','pD','mean_KD_data');
end

%% Statistics
fprintf('\n');
fprintf('Dowel junctions: Bending stiffness per unit length kD\n');
fprintf('mean(kD) = %g kN/rad\n',mean(mean_KD_data)*1e-3);
fprintf('var(kD)  = %g (kN/rad)^2\n',var(mean_KD_data)*1e-6);
fprintf('std(kD)  = %g kN/rad\n',std(mean_KD_data)*1e-3);
fprintf('cv(kD)   = %g\n',std(mean_KD_data)/mean(mean_KD_data));

%% Plot data
if displaySolution
    colors = distinguishable_colors(numDowels);
    
    for j=1:numDowels
        numSample = ['D' num2str(j)];
        
        figure('name',['Dowel junction ' numSample ': bending moment per unit length w.r.t. variation of angle'])
        clf
        plot(angleD{j},mlD{j},'+b','LineWidth',linewidth)
        hold on
        plot(angleD{j},polyval(pD(j,:),deg2rad(angleD{j})),'-r','LineWidth',linewidth)
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
        
        % figure('name',['Dowel junction ' numSample ': Bending stiffness per unit length w.r.t. applied load'])
        % clf
        % plot(FD{j},kD{j}*1e-3,'+b','LineWidth',linewidth)
        % grid on
        % set(gca,'FontSize',fontsize)
        % xlabel('Applied load [N]','Interpreter',interpreter);
        % ylabel('Bending stiffness per unit length [kN/rad]','Interpreter',interpreter);
        % %xlabel('Chargement appliqu\''e [N]','Interpreter',interpreter);
        % %ylabel('Rigidit\''e lin\''eique en flexion [kN/rad]','Interpreter',interpreter);
        % mysaveas(pathname,['data_stiffness_load_' numSample],formats);
        % mymatlab2tikz(pathname,['data_stiffness_load_' numSample '.tex']);
    end
    
    figure('name','Dowel junctions: Bending moment per unit length w.r.t. variation of angle')
    clf
    leg = cell(numDowels,1);
    h = gobjects(numDowels,1);
    for j=1:numDowels
        plot(angleD{j},mlD{j},'+','Color',colors(j,:),'LineWidth',linewidth)
        hold on
        h(j) = plot(angleD{j},polyval(pD(j,:),deg2rad(angleD{j})),'LineStyle','-','Color',colors(j,:),'LineWidth',linewidth);
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
    mysaveas(pathname,'data_moment_angle_D',formats);
    mymatlab2tikz(pathname,'data_moment_angle_D.tex');
    
    % figure('name','Dowel junctions: Bending stiffness per unit length w.r.t. applied load')
    % clf
    % leg = cell(numDowels,1);
    % for j=1:numDowels
    %     plot(FD{j},kD{j}*1e-3,'+','Color',colors(j,:),'LineWidth',linewidth)
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
    % mysaveas(pathname,'data_stiffness_load_D',formats);
    % mymatlab2tikz(pathname,'data_stiffness_load_D.tex');
    
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
