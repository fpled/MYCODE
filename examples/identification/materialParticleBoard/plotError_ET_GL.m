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

%% Plot error
sample = 'B';
j = 14; % sample number
numSample = [sample num2str(j)];

F = appliedLoad(numSample);
[b,h,d,Iz] = dimSample(numSample);

numImages = length(F);
% for k=1:numImages
for k=6
    
    t = tic;
    
    numImage = num2str(k,'%02d');
    filenameDIC = [numSample '_00-' numImage '-Mesh'];
    disp(['File = ',filenameDIC]);
    load(fullfile(pathnameDIC,filenameDIC));
    
    [u_exp,coord] = extractCorreli(Job,Mesh,U,h,d); % [mm]
    
    ET = ET_data{j}(k); % [MPa]
    GL = GL_data{j}(k); % [MPa]
    R0 = R0_data{j}(k); % [rad]
    U0 = U0_data{j}(k); % [mm]
    V0 = V0_data{j}(k); % [mm]
    err = err_ana_data{j}(k);
    
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
    
    figure('Name','Surface plot: Error with respect to ET and GL')
    clf
    surfc(GL_series,ET_series,err_series,'EdgeColor','none');
    colorbar
    hold on
    % scatter3(GL_min,ET_min,err_min,'MarkerEdgeColor','k','MarkerFaceColor','r');
    scatter3(GL,ET,err,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    % set(gca,'ZScale','log')
    xlabel('$G^L$ [MPa]','Interpreter',interpreter)
    ylabel('$E^T$ [MPa]','Interpreter',interpreter)
    zlabel('Error','Interpreter',interpreter)
    %zlabel('Erreur','Interpreter',interpreter)
    mysaveas(pathname,['error_ET_GL_' numSample '_image_' numImage '_3D'],formats);
    
    figure('Name','Contour plot: Error with respect to ET and GL')
    clf
    contourf(GL_series,ET_series,err_series,30);
    colorbar
    hold on
    % scatter(GL_min,ET_min,'MarkerEdgeColor','k','MarkerFaceColor','r');
    scatter(GL,ET,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    % set(gca,'ZScale','log')
    xlabel('$G^L$ [MPa]','Interpreter',interpreter)
    ylabel('$E^T$ [MPa]','Interpreter',interpreter)
    mysaveas(pathname,['error_ET_GL_' numSample '_image_' numImage '_2D'],formats);
    
    toc(t)
end
