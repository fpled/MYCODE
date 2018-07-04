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
pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','materialParticleBoard','resultsDIC');

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};

%% Identification
EL0 = 1e2; % MPa
NUL0 = 1e-2;
x0 = [EL0 NUL0];
lb = [0 0];
ub = [Inf 0.5];
optimFun = 'lsqnonlin';
% optimFun = 'fminsearch';
% optimFun = 'fminunc';
% optimFun = 'fmincon';
display = 'off';
tolX = 1e-14;
tolFun = 1e-14;

switch optimFun
    case {'lsqnonlin','fminunc','fmincon'}
        options  = optimoptions(optimFun,'Display',display,'TolX',tolX,'TolFun',tolFun);
    case 'fminsearch'
        options = optimset('Display',display,'TolX',tolX,'TolFun',tolFun);
    otherwise
        error(['Wrong optimization function' optimFun])
end

sample = 'B';
numSamples = 27;
EL_data = cell(numSamples,1);
NUL_data = cell(numSamples,1);
err_num_data = cell(numSamples,1);

mean_EL_data = zeros(numSamples,1);
mean_NUL_data = zeros(numSamples,1);
std_EL_data = zeros(numSamples,1);
std_NUL_data = zeros(numSamples,1);

for j=1:numSamples
    
    numSample = [sample num2str(j)];
    F = appliedLoad(numSample);
    [b,h,d,Iz] = dimSample(numSample);
    
    t = tic;
    
    numImages = length(F);
    EL = zeros(numImages,1);
    NUL = zeros(numImages,1);
    err = zeros(numImages,1);
    
    for k=1:numImages
        
        numImage = num2str(k,'%02d');
        filenameDIC = [numSample '_00-' numImage '-Mesh'];
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
        
        ET = ET_data{j}(k); % MPa
        GL = GL_data{j}(k); % MPa
        mat = ELAS_ISOT_TRANS('AXISL',[0;1],'AXIST',[1;0],'EL',[],'ET',ET,'NUL',[],'GL',GL,'DIM3',h);
        mat = setnumber(mat,1);
        S = setmaterial(S,mat);
        S = final(S);
        I = create_boundary(S);
        I = final(I);
        P = calc_P(S,I);
        u_exp_b = P*u_exp;
        S = addcl(S,[],'U',u_exp_b);
        u_exp_in = freevector(S,u_exp);
        
        switch optimFun
            case 'lsqnonlin'
                fun = @(x) funlsqnonlin(x,@solveThreePointBendingNum,u_exp_in,S);
                [x,err(k),~,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
            case 'fminsearch'
                fun = @(x) funoptim(x,@solveThreePointBendingNum,u_exp_in,S);
                [x,err(k),exitflag,output] = fminsearch(fun,x0,options);
            case 'fminunc'
                fun = @(x) funoptim(x,@solveThreePointBendingNum,u_exp_in,S);
                [x,err(k),exitflag,output] = fminunc(fun,x0,options);
            case 'fmincon'
                fun = @(x) funoptim(x,@solveThreePointBendingNum,u_exp_in,S);
                [x,err(k),exitflag,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
        end
        
        EL(k) = x(1); % MPa
        NUL(k) = x(2);
        err(k) = sqrt(err(k))./norm(u_exp_in);
    end
    
    %% Outputs
    fprintf('\n')
    disp('+-----------------+')
    fprintf('| Sample %s%2d      |\n',sample,j)
    disp('+-----------------+-----------------+-----------------+')
    disp('| Young''s modulus | Poisson''s ratio |  Error between  |')
    disp('|     EL (MPa)    |       NUL       | U_num and U_exp |')
    disp('+-----------------+-----------------+-----------------+')
    for k=1:numImages
        fprintf('| %15.4f | %15.4f | %15.4e |\n',EL(k),NUL(k),err(k))
    end
    disp('+-----------------+-----------------+-----------------+')
    
    toc(t)
    
    EL_data{j} = EL;
    NUL_data{j} = NUL;
    err_num_data{j} = err;
    mean_EL_data(j) = mean(EL(initImage:end));
    mean_NUL_data(j) = mean(NUL(initImage:end));
    std_EL_data(j) = std(EL(initImage:end));
    std_NUL_data(j) = std(NUL(initImage:end));
    
end

%% Save variables
save(fullfile(pathname,filenameNum),'EL_data','NUL_data','err_num_data',...
    'mean_EL_data','mean_NUL_data','std_EL_data','std_NUL_data');

%% Plot data
if displaySolution
%     for j=1:numSamples
%         numSample = ['B' num2str(j)];
%         
%         figure
%         clf
%         bar(EL_data{j}(initImage:end));
%         grid on
%         set(gca,'FontSize',fontsize)
%         legend(numSample,'Location','NorthEastOutside');
%         xlabel('Image number','Interpreter',interpreter);
%         ylabel('Young''s modulus $E^L$ (MPa)','Interpreter',interpreter);
%         mysaveas(pathname,['data_EL_' numSample],formats);
%         mymatlab2tikz(pathname,['data_EL_' numSample '.tex']);
%         
%         figure
%         clf
%         bar(NUL_data{j}(initImage:end));
%         grid on
%         set(gca,'FontSize',fontsize)
%         legend(numSample,'Location','NorthEastOutside');
%         xlabel('Image number','Interpreter',interpreter);
%         ylabel('Poisson''s ratio $\nu^L$','Interpreter',interpreter);
%         mysaveas(pathname,['data_NUL_' numSample],formats);
%         mymatlab2tikz(pathname,['data_NUL_' numSample '.tex']);
%     end
    
    figure
    clf
    bar(mean_EL_data,'LineWidth',linewidth);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Young''s modulus $E^L$ (MPa)','Interpreter',interpreter);
    mysaveas(pathname,'data_EL',formats);
    mymatlab2tikz(pathname,'data_EL.tex');
    
    figure
    clf
    bar(mean_NUL_data,'LineWidth',linewidth);
    grid on
    set(gca,'FontSize',fontsize)
    xlabel('Sample number','Interpreter',interpreter);
    ylabel('Poisson''s ratio $\nu^L$','Interpreter',interpreter);
    mysaveas(pathname,'data_NUL',formats);
    mymatlab2tikz(pathname,'data_NUL.tex');
end
