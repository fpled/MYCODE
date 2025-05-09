%% Plot numerical solution and experimental data %%
%%-----------------------------------------------%%

% clc
clearvars
close all

%% Input data
filenameAna = 'data_ET_GL.mat';
filenameNum = 'data_EL_NUL.mat';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathname,filenameAna));
load(fullfile(pathname,filenameNum));

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'painters';

%% Compute analytical and numerical solutions
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
        numImage = num2str(k,'%02d');
        
        filenameDIC = [numSample '_00-' numImage '-Mesh'];
        pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
            'examples','identification','materialParticleBoard','resultsDIC');
        disp(['File = ',filenameDIC]);
        load(fullfile(pathnameDIC,filenameDIC));
        
        [u_exp,coord] = extractCorreliElas(Job,Mesh,U,h,d); % [mm]
        
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
        R0 = R0_data{j}(k); % [rad]
        U0 = U0_data{j}(k); % [mm]
        V0 = V0_data{j}(k); % [mm]
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
        
        x_ana = [ET GL R0 U0 V0];
        u_ana = solveThreePointBendingAna(x_ana,coord,F(k),Iz,h); % [mm]
        
        x_num = [EL NUL];
        [u_in,S] = solveThreePointBendingNum(x_num,S); % [mm]
        u_num = unfreevector(S,u_in);
        
        %% Display solutions
        ampl = 0;
        v_exp = calc_init_dirichlet(S);
        figure('Name','Imposed experimental displacement')
        clf
        hM = plot(S,'FaceColor','w','LineWidth',0.5);
        hold on
        [hD,legD] = vectorplot(S,'U',v_exp,ampl,'r','LineWidth',1);
        hold off
        set(gca,'FontSize',fontsize)
        hg = hggroup;
        set([hM(:),hD],'Parent',hg);
        axis image
        l = legend(hD,'{\boldmath$u$}$^{\mathrm{exp}}$','Location','NorthEastOutside');
        % l = legend(hD,legD,'Location','NorthEastOutside');
        set(l,'Interpreter',interpreter);
        mysaveas(pathname,['u_exp_imposed_' numSample '_' numImage],formats,renderer);
        
        ampl = 0;
        % ampl = getsize(S)/max(max(max(abs(u_ana)),max(abs(u_num))),max(abs(u_exp)))/5;
        ddl = [findddl(S,'UX'),findddl(S,'UY')];
        
        for i=1:2
            u_ana_i = u_ana(ddl(:,i));
            u_num_i = u_num(ddl(:,i));
            u_exp_i = u_exp(ddl(:,i));
            u_min = min(min(min(u_ana_i),min(u_num_i)),min(u_exp_i));
            u_max = max(max(max(u_ana_i),max(u_num_i)),max(u_exp_i));
            
            plotSolution(S,u_ana,'displ',i,'ampl',ampl,'FontSize',fontsize);
            set(gca,'CLim',[u_min,u_max])
            mysaveas(pathname,['u_' num2str(i) '_ana_' numSample '_' numImage],formats,renderer);
            
            plotSolution(S,u_num,'displ',i,'ampl',ampl,'FontSize',fontsize);
            set(gca,'CLim',[u_min,u_max])
            mysaveas(pathname,['u_' num2str(i) '_num_' numSample '_' numImage],formats,renderer);
            
            plotSolution(S,u_exp,'displ',i,'ampl',ampl,'FontSize',fontsize);
            set(gca,'CLim',[u_min,u_max])
            mysaveas(pathname,['u_' num2str(i) '_exp_' numSample '_' numImage],formats,renderer);
            
            figure('Name',['Solution u_' num2str(i)])
            tiledlayout(3,1)
            nexttile
            plot_sol(S,u_ana,'displ',i,'ampl',ampl);
            if i==1
                t = '$u^{\mathrm{ana}}$';
            elseif i==2
                t = '$v^{\mathrm{ana}}$';
            end
            title(t,'FontSize',fontsize,'Interpreter',interpreter)
            set(gca,'CLim',[u_min,u_max])
            set(gca,'FontSize',fontsize)
            nexttile
            plot_sol(S,u_num,'displ',i,'ampl',ampl);
            if i==1
                t = '$u^{\mathrm{num}}$';
            elseif i==2
                t = '$v^{\mathrm{num}}$';
            end
            title(t,'FontSize',fontsize,'Interpreter',interpreter)
            set(gca,'CLim',[u_min,u_max])
            set(gca,'FontSize',fontsize)
            nexttile
            plot_sol(S,u_exp,'displ',i,'ampl',ampl);
            if i==1
                t = '$u^{\mathrm{exp}}$';
            elseif i==2
                t = '$v^{\mathrm{exp}}$';
            end
            title(t,'FontSize',fontsize,'Interpreter',interpreter)
            set(gca,'CLim',[u_min,u_max])
            set(gca,'FontSize',fontsize)
            cb = colorbar;
            cb.Layout.Tile = 'east';
            mysaveas(pathname,['u_' num2str(i) '_' numSample '_' numImage],formats,renderer);
        end
        
        % for i=1:3
        %     plotSolution(S,u_ana,'epsilon',i,'ampl',ampl);
        %     mysaveas(pathname,['eps_' num2str(i) '_ana_' numSample '_' numImage],formats,renderer);
        %     plotSolution(S,u_num,'epsilon',i,'ampl',ampl);
        %     mysaveas(pathname,['eps_' num2str(i) '_num_' numSample '_' numImage],formats,renderer);
        %     plotSolution(S,u_exp,'epsilon',i,'ampl',ampl);
        %     mysaveas(pathname,['eps_' num2str(i) '_exp_' numSample '_' numImage],formats,renderer);
        %
        %     plotSolution(S,u_ana,'sigma',i,'ampl',ampl);
        %     mysaveas(pathname,['sig_' num2str(i) '_num_' numSample '_' numImage],formats,renderer);
        %     plotSolution(S,u_num,'sigma',i,'ampl',ampl);
        %     mysaveas(pathname,['sig_' num2str(i) '_num_' numSample '_' numImage],formats,renderer);
        %     plotSolution(S,u_exp,'sigma',i,'ampl',ampl);
        %     mysaveas(pathname,['sig_' num2str(i) '_exp_' numSample '_' numImage],formats,renderer);
        % end
        % 
        % plotSolution(S,u_ana,'epsilon','mises','ampl',ampl);
        % mysaveas(pathname,['eps_von_mises_ana_' numSample '_' numImage],formats,renderer);
        % plotSolution(S,u_num,'epsilon','mises','ampl',ampl);
        % mysaveas(pathname,['eps_von_mises_num_' numSample '_' numImage],formats,renderer);
        % plotSolution(S,u_exp,'epsilon','mises','ampl',ampl);
        % mysaveas(pathname,['eps_von_mises_exp_' numSample '_' numImage],formats,renderer);
        % 
        % plotSolution(S,u_ana,'sigma','mises','ampl',ampl);
        % mysaveas(pathname,['sig_von_mises_ana_' numSample '_' numImage],formats,renderer);
        % plotSolution(S,u_num,'sigma','mises','ampl',ampl);
        % mysaveas(pathname,['sig_von_mises_num_' numSample '_' numImage],formats,renderer);
        % plotSolution(S,u_exp,'sigma','mises','ampl',ampl);
        % mysaveas(pathname,['sig_von_mises_exp_' numSample '_' numImage],formats,renderer);
    end
    close all
end
