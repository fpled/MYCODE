%% Identification of ET and GL - Analytical beam model %%
%%-----------------------------------------------------%%
% Call after Digital Image Correlation (DIC) RT3
% Coordinate system
% CORRELI:    right     down
% MATLAB:     right     up

% clc
clearvars
close all

%% Input data
displaySolution = true;

filename = 'data_ET_GL.mat';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc2'};

%% Identification
sample = 'B';
for j = 1:27
    
    sampleNum = [sample num2str(j)];
    
    F = appliedLoad(sampleNum);
    [b,h,d,Iz] = dimSample(sampleNum);
    
    t = tic;
    
    ET = zeros(length(F),1);
    GL = zeros(length(F),1);
    Phi = zeros(length(F),1);
    U0 = zeros(length(F),1);
    V0 = zeros(length(F),1);
    err = zeros(length(F),1);
    
    for k = 1:length(F)
        
        if k<10
            imageNum = ['0' num2str(k)];
        else
            imageNum = num2str(k);
        end
        
        filenameDIC = [sampleNum '_00-' imageNum '-Mesh'];
        pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
            'examples','identification','materialParticleBoard','resultsDIC');
        load(fullfile(pathnameDIC,filenameDIC));
        
        [u_exp,coord] = extractCorreli(Job,Mesh,U,h,d);
        
        A = zeros(5);
        B = zeros(5,1);
        for i=1:Mesh.NNodeTot
            xi = coord(i,1);
            yi = coord(i,2);
            ui = u_exp(2*i-1);
            vi = u_exp(2*i);
            
            A_UE = -F(k)/(4*Iz)*yi*xi^2;
            A_UG = F(k)/(12*Iz)*yi^3;
            A_VE = F(k)/(12*Iz)*xi^3;
            A_VG = -F(k)*h^2/(16*Iz)*xi;
            
            A = A + [ A_UE^2+A_VE^2        A_UE*A_UG+A_VE*A_VG  A_UE*yi-A_VE*xi  A_UE  A_VE
                A_UE*A_UG+A_VE*A_VG  A_UG^2+A_VG^2        A_UG*yi-A_VG*xi  A_UG  A_VG
                A_UE*yi-A_VE*xi      A_UG*yi-A_VG*xi      yi^2+xi^2        yi    -xi
                A_UE                 A_UG                 yi               1     0
                A_VE                 A_VG                 -xi              0     1 ];
            
            B = B + [ A_UE*ui + A_VE*vi
                A_UG*ui + A_VG*vi
                yi*ui   - xi*vi
                ui
                vi ];
        end
        
        param = A\B;
        
        ET(k)= 1/param(1)*1e-3; % GPa
        GL(k) = 1/param(2); % MPa
        Phi(k) = param(3);
        U0(k) = param(4);
        V0(k) = param(5);
        
        x = [ET(k)*1e3 GL(k) Phi(k) U0(k) V0(k)];
        u_ana = solveThreePointBendingAna(x,coord,F(k),Iz,h);
        err(k) = norm(u_ana - u_exp);
        err(k) = err(k)./norm(u_exp);
    end
    toc(t)
    
    %% Outputs
    fprintf('\n')
    disp('+---------------+')
    fprintf('| Sample %s%2d    |\n',sample,j)
    disp('+---------------+---------------+-----------------+')
    disp('| Young modulus | Shear modulus |  Error between  |')
    disp('|    ET (GPa)   |    GL (MPa)   | U_ana and U_exp |')
    disp('+---------------+---------------+-----------------+')
    for k=1:length(F)
        fprintf('| %13.4f | %13.4f | %15.4e |\n',ET(k),GL(k),err(k))
    end
    disp('+---------------+---------------+-----------------+')
    
    eval(['ET_' sampleNum '= ET;']);
    eval(['GL_' sampleNum '= GL;']);
    eval(['Phi_' sampleNum '= Phi;']);
    eval(['U0_' sampleNum '= U0;']);
    eval(['V0_' sampleNum '= V0;']);
    eval(['err_ana_' sampleNum '= err;']);
    
    initIm = 8;
    eval(['ET_' sampleNum '_data = ET_' sampleNum '(initIm:end);']);
    eval(['GL_' sampleNum '_data = GL_' sampleNum '(initIm:end);']);
    eval(['mean_ET_' sampleNum '_data = mean(ET_' sampleNum '_data);']);
    eval(['mean_GL_' sampleNum '_data = mean(GL_' sampleNum '_data);']);
    eval(['std_ET_' sampleNum '_data = std(ET_' sampleNum '_data);']);
    eval(['std_GL_' sampleNum '_data = std(GL_' sampleNum '_data);']);
    
    %% Save variables
    if isempty(dir(fullfile(pathname,filename)))
        save(fullfile(pathname,filename),['ET_' sampleNum],['GL_' sampleNum],...
            ['Phi_' sampleNum],['U0_' sampleNum],['V0_' sampleNum],['err_ana_' sampleNum],...
            ['ET_' sampleNum '_data'],['GL_' sampleNum '_data'],...
            ['mean_ET_' sampleNum '_data'],['mean_GL_' sampleNum '_data']);
    else
        save(fullfile(pathname,filename),['ET_' sampleNum],['GL_' sampleNum],...
            ['Phi_' sampleNum],['U0_' sampleNum],['V0_' sampleNum],['err_ana_' sampleNum],...
            ['ET_' sampleNum '_data'],['GL_' sampleNum '_data'],...
            ['mean_ET_' sampleNum '_data'],['mean_GL_' sampleNum '_data'],'-append');
    end
    
    %% Plot data
    if displaySolution
        figure
        eval(['bar(ET_' sampleNum '_data);']);
        grid on
        set(gca,'FontSize',fontsize)
        legend(sampleNum,'Location','NorthEastOutside');
        xlabel('Image number','Interpreter',interpreter);
        ylabel('Young modulus $E^T$ (GPa)','Interpreter',interpreter);
        mysaveas(pathname,['data_ET_' sampleNum],formats);
        mymatlab2tikz(pathname,['data_ET_' sampleNum '.tex']);
        
        figure
        eval(['bar(GL_' sampleNum '_data);']);
        grid on
        set(gca,'FontSize',fontsize)
        legend(sampleNum,'Location','NorthEastOutside');
        xlabel('Image number','Interpreter',interpreter);
        ylabel('Shear modulus $G^L$ (MPa)','Interpreter',interpreter);
        mysaveas(pathname,['data_GL_' sampleNum],formats);
        mymatlab2tikz(pathname,['data_GL_' sampleNum '.tex']);
    end
end
