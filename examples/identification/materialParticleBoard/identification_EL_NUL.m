%% Identification of EL and NUL - Finite Element Model Updating (FEMU) method %%
%%----------------------------------------------------------------------------%%
% Call after Digital Image Correlation (DIC) RT3 and Identification of ET and GL

% clc
clearvars
close all

%% Input data
displaySolution = true;

filenameAna = 'data_ET_GL.mat';
filenameNum = 'data_EL_NUL.mat';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathname,filenameAna));

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc2'};

%% Identification
x0 = [1e2 1e-3];
lb = [0 0];
ub = [Inf 0.5];
tolX = 1e-14;
tolFun = 1e-14;
display = 'off';

optionslsqnonlin  = optimoptions('lsqnonlin','Display',display,'TolX',tolX,'TolFun',tolFun);
optionsfminsearch = optimset('Display',display,'TolX',tolX,'TolFun',tolFun);
optionsfminunc    = optimoptions('fminunc','Display',display,'TolX',tolX,'TolFun',tolFun,'Algorithm','quasi-newton');
optionsfmincon    = optimoptions('fmincon','Display',display,'TolX',tolX,'TolFun',tolFun);

sample = 'C';
for j = 1:20
    
    sampleNum = [sample num2str(j)];
    
    F = appliedLoad(sampleNum);
    [b,h,d,Iz] = dimSample(sampleNum);
    
    t = tic;
    
    EL = zeros(length(F),1);
    NUL = zeros(length(F),1);
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
        
        node = NODE(coord,1:size(coord,1));
        elem = Mesh.TRI;
        elemtype = 'TRI3';
        % option = 'DEFO'; % plane strain
        option = 'CONT'; % plane stress
        S = MODEL('PLAN');
        S = addnode(S,node);
        S = addelem(S,elemtype,elem,'option',option);
        
        ET = eval(['ET_' sampleNum '(' imageNum ')'])*1e3; % MPa
        GL = eval(['GL_' sampleNum '(' imageNum ')']); % MPa
        mat = ELAS_ISOT_TRANS('EL',[],'ET',ET,'NUL',[],'GL',GL,'DIM3',h);
        mat = setnumber(mat,1);
        S = setmaterial(S,mat);
        S = final(S);
        I = create_boundary(S);
        I = final(I);
        P = calc_P(S,I);
        u_exp_b = P*u_exp;
        S = addcl(S,[],'U',u_exp_b);
        u_exp_in = freevector(S,u_exp);
        
        funlsqnonlin = @(x) funlsqnonlinNum(x,u_exp_in,S);
        % funoptim = @(x) funoptimNum(x,u_exp_in,S);
        
        [x,err(k),~,exitflag,output] = lsqnonlin(funlsqnonlin,x0,lb,ub,optionslsqnonlin);
        % [x,err(k),exitflag,output] = fminsearch(funoptim,x0,optionsfminsearch);
        % [x,err(k),exitflag,output] = fminunc(funoptim,x0,optionsfminunc);
        % [x,err(k),exitflag,output] = fmincon(funoptim,x0,[],[],[],[],lb,ub,[],optionsfmincon);
        
        EL(k) = x(1); % MPa
        NUL(k) = x(2);
        err(k) = sqrt(err(k))./norm(u_exp_in);
    end
    toc(t)
    
    %% Outputs
    fprintf('\n')
    disp('+---------------+')
    fprintf('| Sample %s%2d    |\n',sample,j)
    disp('+---------------+---------------+-----------------+')
    disp('| Young modulus | Poisson ratio |  Error between  |')
    disp('|    EL (MPa)   |      NUL      | U_num and U_exp |')
    disp('+---------------+---------------+-----------------+')
    for k=1:length(F)
        fprintf('| %13.4f | %13.4f | %15.4e |\n',EL(k),NUL(k),err(k))
    end
    disp('+---------------+---------------+-----------------+')
    
    eval(['EL_' sampleNum '= EL;']);
    eval(['NUL_' sampleNum '= NUL;']);
    eval(['err_num_' sampleNum '= err;']);
    
    imageInit = 8;
    eval(['EL_' sampleNum '_data = EL_' sampleNum '(imageInit:end);']);
    eval(['NUL_' sampleNum '_data = NUL_' sampleNum '(imageInit:end);']);
    eval(['mean_EL_' sampleNum '_data = mean(EL_' sampleNum '_data);']);
    eval(['mean_NUL_' sampleNum '_data = mean(NUL_' sampleNum  '_data);']);
    eval(['std_EL_' sampleNum '_data = std(EL_' sampleNum '_data);']);
    eval(['std_NUL_' sampleNum '_data = std(NUL_' sampleNum '_data);']);
    
    %% Save variables
    if isempty(dir(fullfile(pathname,filenameNum)))
        save(fullfile(pathname,filenameNum),['EL_' sampleNum],['NUL_' sampleNum],['err_num_' sampleNum],...
            ['EL_' sampleNum '_data'],['GL_' sampleNum '_data'],...
            ['mean_EL_' sampleNum '_data'],['mean_NUL_' sampleNum '_data']);
    else
        save(fullfile(pathname,filenameNum),['EL_' sampleNum],['NUL_' sampleNum],['err_num_' sampleNum],...
            ['EL_' sampleNum '_data'],['GL_' sampleNum '_data'],...
            ['mean_EL_' sampleNum '_data'],['mean_NUL_' sampleNum '_data'],'-append');
    end
    
    %% Plot data
    if displaySolution
        figure
        eval(['bar(EL_' sampleNum '_data);']);
        grid on
        set(gca,'FontSize',fontsize)
        legend(sampleNum,'Location','NorthEastOutside');
        xlabel('Image number','Interpreter',interpreter);
        ylabel('Young modulus $E^L$ (MPa)','Interpreter',interpreter);
        mysaveas(pathname,['data_EL_' sampleNum],formats);
        mymatlab2tikz(pathname,['data_EL_' sampleNum '.tex']);
        
        figure
        eval(['bar(NUL_' sampleNum '_data);']);
        grid on
        set(gca,'FontSize',fontsize)
        legend(sampleNum,'Location','NorthEastOutside');
        xlabel('Image number','Interpreter',interpreter);
        ylabel('Poisson ratio $\nu^L$','Interpreter',interpreter);
        mysaveas(pathname,['data_NUL_' sampleNum],formats);
        mymatlab2tikz(pathname,['data_NUL_' sampleNum '.tex']);
    end
end
