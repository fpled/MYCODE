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

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Plot error
sampleNum = 'B8';

F = appliedLoad(sampleNum);
[b,h,d,Iz] = dimSample(sampleNum);

for k = 1:length(F)
    
    if k<10
        imageNum = ['0' num2str(k)];
    else
        imageNum = num2str(k);
    end
    
    t = tic;
    
    filenameDIC = [sampleNum '_00-' imageNum '-Mesh'];
    pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'examples','identification','materialParticleBoard','resultsDIC');
    disp(['File = ',filenameDIC]);
    load(fullfile(pathnameDIC,filenameDIC));
    
    [u_exp,coord] = extractCorreli(Job,Mesh,U,h,d);
    
    ET = eval(['ET_' sampleNum '(' imageNum ');'])*1e3; % MPa
    GL = eval(['GL_' sampleNum '(' imageNum ');']); % MPa
    Phi = eval(['Phi_' sampleNum '(' imageNum ');']);
    U0 = eval(['U0_' sampleNum '(' imageNum ');']);
    V0 = eval(['V0_' sampleNum '(' imageNum ');']);
    ET_series = linspace(ET*0.5,ET*1.5,1e2); % MPa
    GL_series = linspace(GL*0.5,GL*1.5,1e2); % MPa
    
    err = zeros(length(ET_series),length(GL_series));
    for m=1:length(ET_series)
        ET = ET_series(m);
        for n=1:length(GL_series)
            GL = GL_series(n);
            param = [ET GL Phi U0 V0];
            u_ana = solveThreePointBendingAna(param,coord,F(k),Iz,h);
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
    xlabel('$G^L$ (MPa)','Interpreter',interpreter)
    ylabel('$E^T$ (MPa)','Interpreter',interpreter)
    zlabel('$\varepsilon$','Interpreter',interpreter)
    mysaveas(pathname,['error_ET_GL_' sampleNum '_image_' imageNum '_3D'],formats,renderer);
    
    figure
    contourf(GL_series,ET_series,err,30);
    colorbar
    hold on
    scatter(GL_series(c),ET_series(r),'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    % set(gca,'ZScale','log')
    xlabel('$G^L$ (MPa)','Interpreter',interpreter)
    ylabel('$E^T$ (MPa)','Interpreter',interpreter)
    mysaveas(pathname,['error_ET_GL_' sampleNum '_image_' imageNum '_2D'],formats,renderer);
    
    toc(t)
end
