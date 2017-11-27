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

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Plot error
sampleNum = 'C3';

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
    
    node = NODE(coord,1:size(coord,1));
    elem = Mesh.TRI;
    elemtype = 'TRI3';
    % option = 'DEFO'; % plane strain
    option = 'CONT'; % plane stress
    S = MODEL('PLAN');
    S = addnode(S,node);
    S = addelem(S,elemtype,elem,'option',option);
    
    ET = eval(['ET_' sampleNum '(' imageNum ');'])*1e3; % MPa
    GL = eval(['GL_' sampleNum '(' imageNum ');']); % MPa
    EL = eval(['EL_' sampleNum '(' imageNum ');']); % MPa
    NUL = eval(['NUL_' sampleNum '(' imageNum ');']);
    mat = ELAS_ISOT_TRANS('EL',EL,'ET',ET,'NUL',NUL,'GL',GL,'DIM3',h);
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
    
    EL_series=linspace(EL*0.5,EL*1.5,10);
    NUL_series=linspace(NUL*0.5,NUL*1.5,10);
    
    err=zeros(length(EL_series),length(NUL_series));
    for m=1:length(EL_series)
        EL = EL_series(m);
        for n=1:length(NUL_series)
            NUL = NUL_series(n);
            param = [EL NUL];
            u_num_in = solveThreePointBendingNum(param,S);
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
    ylabel('$E^L$ (MPa)','Interpreter',interpreter)
    zlabel('$\varepsilon$','Interpreter',interpreter)
    mysaveas(pathname,['error_EL_NUL_' sampleNum '_image_' imageNum '_3D'],formats,renderer);
    
    figure
    contourf(NUL_series,EL_series,err,30);
    colorbar
    hold on
    scatter(NUL_series(c),EL_series(r),'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    set(gca,'FontSize',fontsize)
    % set(gca,'ZScale','log')
    xlabel('$\nu^L$','Interpreter',interpreter)
    ylabel('$E^L$ (MPa)','Interpreter',interpreter)
    mysaveas(pathname,['error_EL_NUL_' sampleNum '_image_' imageNum '_2D'],formats,renderer);
    
    toc(t)
end
