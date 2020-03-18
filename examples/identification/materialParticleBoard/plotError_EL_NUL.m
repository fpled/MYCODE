%% Plot error with respect to EL and NUL %%
%%---------------------------------------%%

% clc
clearvars
close all

%% Input Data
filenameAna = 'data_ET_GL.mat';
filenameNum = 'data_EL_NUL.mat';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathname,filenameAna));
load(fullfile(pathname,filenameNum));
pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','materialParticleBoard','resultsDIC');

fontsize = 16;
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
    
    [u_exp,coord] = extractCorreli(Job,Mesh,U,h,d); % [mm]
    
    node = NODE(coord,1:size(coord,1));
    elem = Mesh.TRI;
    elemtype = 'TRI3';
    % option = 'DEFO'; % plane strain
    option = 'CONT'; % plane stress
    S = MODEL('PLAN');
    S = addnode(S,node);
    S = addelem(S,elemtype,elem,'option',option);
    
    ET = ET_data{j}(k); % [MPa]
    GL = GL_data{j}(k); % [MPa]
    EL = EL_data{j}(k); % [MPa]
    NUL = NUL_data{j}(k);
    mat = ELAS_ISOT_TRANS('AXISL',[0;1],'AXIST',[1;0],'EL',EL,'ET',ET,'NUL',NUL,'GL',GL,'DIM3',h);
    mat = setnumber(mat,1);
    S = setmaterial(S,mat);
    S = final(S);
    I = create_boundary(S);
    I = final(I);
    P = calc_P(S,I);
    u_exp_b = P*u_exp;
    S = addcl(S,[],'U',u_exp_b);
    u_exp_in = freevector(S,u_exp);
    
    % nodes_b = getnumber(getnode(create_boundary(S)));
    % u_exp_b = [ux_exp(nodes_b) uy_exp(nodes_b)]';
    % u_exp_b = u_exp_b(:);
    
    % nodes_in = setdiff(1:Mesh.NNodeTot,nodes_b);
    % u_exp_in = [ux_exp(nodes_in) uy_exp(nodes_in)]';
    % u_exp_in = u_exp_in(:);
    
    EL_series=linspace(EL*0.5,EL*1.5,1e2);
    NUL_series=linspace(NUL*0.5,NUL*1.5,1e2);
    
    err=zeros(length(EL_series),length(NUL_series));
    for m=1:length(EL_series)
        EL = EL_series(m);
        for n=1:length(NUL_series)
            NUL = NUL_series(n);
            x = [EL NUL];
            u_num_in = solveThreePointBendingNum(x,S);
            err(m,n) = norm(u_exp_in - u_num_in);
        end
    end
    err = err./norm(u_exp_in);
    [errmin,I] = min(err);
    [errmin,c] = min(errmin);
    r = I(c);
    
    figure
    surfc(NUL_series,EL_series,err,'EdgeColor','none');
    colorbar
    hold on
    scatter3(NUL_series(c),EL_series(r),errmin,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    % set(gca,'ZScale','log')
    xlabel('$\nu^L$','Interpreter',interpreter)
    ylabel('$E^L$ [MPa]','Interpreter',interpreter)
    %zlabel('Error','Interpreter',interpreter)
    zlabel('Erreur','Interpreter',interpreter)
    mysaveas(pathname,['error_EL_NUL_' numSample '_image_' numImage '_3D'],formats,renderer);
    
    figure
    contourf(NUL_series,EL_series,err,30);
    colorbar
    hold on
    scatter(NUL_series(c),EL_series(r),'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    % set(gca,'ZScale','log')
    xlabel('$\nu^L$','Interpreter',interpreter)
    ylabel('$E^L$ [MPa]','Interpreter',interpreter)
    mysaveas(pathname,['error_EL_NUL_' numSample '_image_' numImage '_2D'],formats,renderer);
    
    toc(t)
end
