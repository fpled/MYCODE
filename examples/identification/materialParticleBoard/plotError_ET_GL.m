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

fontsize = 12;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

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
    
    [u_exp,coord] = extractCorreli(Job,Mesh,U,h,d);
    
    ET = ET_data{j}(k); % MPa
    GL = GL_data{j}(k); % MPa
    Phi = Phi_data{j}(k);
    U0 = U0_data{j}(k); % mm
    V0 = V0_data{j}(k); % mm
    ET_series = linspace(ET*0.5,ET*1.5,1e2); % MPa
    GL_series = linspace(GL*0.5,GL*1.5,1e2); % MPa
    
    err = zeros(length(ET_series),length(GL_series));
    for m=1:length(ET_series)
        ET = ET_series(m);
        for n=1:length(GL_series)
            GL = GL_series(n);
            x = [ET GL Phi U0 V0];
            u_ana = solveThreePointBendingAna(x,coord,F(k),Iz,h);
            err(m,n) = norm(u_exp - u_ana);
        end
    end
    err = err./norm(u_exp);
    [errmin,I] = min(err);
    [errmin,c] = min(errmin);
    r = I(c);
    
    figure
    surfc(GL_series,ET_series,err,'EdgeColor','none');
    colorbar
    hold on
    scatter3(GL_series(c),ET_series(r),errmin,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    % set(gca,'ZScale','log')
    xlabel('$G^L$ [MPa]','Interpreter',interpreter)
    ylabel('$E^T$ [MPa]','Interpreter',interpreter)
    zlabel('Error','Interpreter',interpreter)
    mysaveas(pathname,['error_ET_GL_' numSample '_image_' numImage '_3D'],formats,renderer);
    
    figure
    contourf(GL_series,ET_series,err,30);
    colorbar
    hold on
    scatter(GL_series(c),ET_series(r),'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    % set(gca,'ZScale','log')
    xlabel('$G^L$ [MPa]','Interpreter',interpreter)
    ylabel('$E^T$ [MPa]','Interpreter',interpreter)
    mysaveas(pathname,['error_ET_GL_' numSample '_image_' numImage '_2D'],formats,renderer);
    
    toc(t)
end
