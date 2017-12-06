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
x0 = [1e3 1e2 0 0 0];
lb = [0 0 -Inf -Inf -Inf];
ub = [Inf Inf Inf Inf Inf];
tolX = 1e-14;
tolFun = 1e-14;
display = 'off';

optionslsqnonlin  = optimoptions('lsqnonlin','Display',display,'TolX',tolX,'TolFun',tolFun);
optionsfminsearch = optimset('Display',display,'TolX',tolX,'TolFun',tolFun);
optionsfminunc    = optimoptions('fminunc','Display',display,'TolX',tolX,'TolFun',tolFun,'Algorithm','quasi-newton');
optionsfmincon    = optimoptions('fmincon','Display',display,'TolX',tolX,'TolFun',tolFun);

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
        
        funlsqnonlin = @(x) funlsqnonlinAna(x,u_exp,coord,F(k),Iz,h);
        % funoptim = @(x) funoptimAna(x,u_exp,coord,F(imageList),Iz,h);
        
        [x,err(k),~,exitflag,output] = lsqnonlin(funlsqnonlin,x0,lb,ub,optionslsqnonlin);
        % [x,err(k),exitflag,output] = fminsearch(funoptim,x0,optionsfminsearch);
        % [x,err(k),exitflag,output] = fminunc(funoptim,x0,optionsfminunc);
        % [x,err(k),exitflag,output] = fmincon(funoptim,x0,[],[],[],[],lb,ub,[],optionsfmincon);
        
        ET(k) = x(1)*1e-3; % GPa
        GL(k) = x(2); % MPa
        Phi(k) = x(3);
        U0(k) = x(4);
        V0(k) = x(5);
        err(k) = sqrt(err(k))./norm(u_exp);
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
    
    initIm = 3;
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
            ['mean_ET_' sampleNum '_data'],['mean_GL_' sampleNum '_data'],...
            ['std_ET_' sampleNum '_data'],['std_GL_' sampleNum '_data']);
    else
        save(fullfile(pathname,filename),['ET_' sampleNum],['GL_' sampleNum],...
            ['Phi_' sampleNum],['U0_' sampleNum],['V0_' sampleNum],['err_ana_' sampleNum],...
            ['ET_' sampleNum '_data'],['GL_' sampleNum '_data'],...
            ['mean_ET_' sampleNum '_data'],['mean_GL_' sampleNum '_data'],...
            ['std_ET_' sampleNum '_data'],['std_GL_' sampleNum '_data'],'-append');
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
