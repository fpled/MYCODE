%% Plot error with respect to ET and GL %%
%%--------------------------------------%%

% clc
clearvars
close all

%% Input Data
filename = 'data_ET_GL.mat';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathname,filename));
pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','materialParticleBoard','resultsDIC');

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};

%% Plot relative error
% geometric dimensions
b = 50; % sample width [mm]
h = 15; % sample thickness [mm]
Iz = b*h^3/12; % planar second moment of area (or planar area moment of inertia) [mm^4]
d = 20; % distance between the support and the region of interest (ROI) [mm]

numSamples = 27;
for j=1:numSamples
    
    numSample = ['B' num2str(j)];
    F = appliedLoad(numSample); % applied load [N]
    
    numImages = length(F);
    % for k=1:numImages
    for k=[6,numImages]
        
        t = tic;
        
        numImage = num2str(k,'%02d');
        filenameDIC = [numSample '_00-' numImage '-Mesh'];
        disp(['File = ',filenameDIC]);
        load(fullfile(pathnameDIC,filenameDIC));
        
        [u_exp,coord] = extractCorreliElas(Job,Mesh,U,h,d); % [mm]
        
        ET = ET_data{j}(k); % [MPa]
        GL = GL_data{j}(k); % [MPa]
        R0 = R0_data{j}(k); % [rad]
        U0 = U0_data{j}(k); % [mm]
        V0 = V0_data{j}(k); % [mm]
        err = err_ana{j}(k);
        
        ET_series = linspace(ET*0.5,ET*1.5,1e2); % [MPa]
        GL_series = linspace(GL*0.5,GL*1.5,1e2); % [MPa]
        
        err_series = zeros(length(ET_series),length(GL_series));
        for m=1:length(ET_series)
            ET_m = ET_series(m);
            for n=1:length(GL_series)
                GL_n = GL_series(n);
                x = [ET_m GL_n R0 U0 V0];
                u_ana = solveThreePointBendingAna(x,coord,F(k),Iz,h);
                err_series(m,n) = norm(u_exp - u_ana);
            end
        end
        err_series = err_series./norm(u_exp);
        % [err_min,I] = min(err_series);
        % [err_min,c] = min(err_min);
        % r = I(c);
        % GL_min = GL_series(c);
        % ET_min = ET_series(r);
        
        figure('Name','Surface plot: Relative error with respect to ET and GL')
        clf
        surfc(GL_series,ET_series,err_series,'EdgeColor','none');
        colorbar
        set(gca,'ColorScale','log')
        %view(-37.5,30) % default view
        view(-30,30)
        hold on
        % scatter3(GL_min,ET_min,err_min,'MarkerEdgeColor','k','MarkerFaceColor','r');
        scatter3(GL,ET,err,'MarkerEdgeColor','k','MarkerFaceColor','r');
        hold off
        set(gca,'FontSize',fontsize)
        set(gca,'ZScale','log')
        xlabel('$G_L$ [MPa]','Interpreter',interpreter)
        ylabel('$E_T$ [MPa]','Interpreter',interpreter)
        zlabel('Relative error','Interpreter',interpreter)
        %zlabel('Erreur relative','Interpreter',interpreter)
        mysaveas(pathname,['error_ET_GL_' numSample '_' numImage '_3D'],formats);
        close(gcf)
        
        figure('Name','Contour plot: Relative error with respect to ET and GL')
        clf
        contourf(GL_series,ET_series,err_series,50);
        colorbar
        set(gca,'ColorScale','log')
        hold on
        % scatter(GL_min,ET_min,'MarkerEdgeColor','k','MarkerFaceColor','r');
        scatter(GL,ET,'MarkerEdgeColor','k','MarkerFaceColor','r');
        hold off
        set(gca,'FontSize',fontsize)
        xlabel('$G_L$ [MPa]','Interpreter',interpreter)
        ylabel('$E_T$ [MPa]','Interpreter',interpreter)
        mysaveas(pathname,['error_ET_GL_' numSample '_' numImage '_2D'],formats);
        close(gcf)
        
        toc(t)
    end
end
