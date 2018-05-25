%% Identification of Kjunction %%
%%-----------------------------%%
% Call after Digital Image Correlation (DIC) RT3
% Coordinate system
% CORRELI:    right     down
% MATLAB:     right     up

clc
clear all
close all

%% Input data
displaySolution = true;

filename = 'data_angle.mat';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc2'};

%% Identification

NumScrew = 4;
NumDowel = 2;

F_screw = cell(NumScrew,1);
Var_angle_screw = cell(NumScrew,1);
for j = 1:NumScrew
    
    sampleNum = ['S' num2str(j)];
    sampleNuma = ['S' num2str(j) 'a'];
    sampleNumb = ['S' num2str(j) 'b'];
    
    F = appliedLoad(sampleNum);
    
    Var_angle = zeros(length(F),1);
    
    for k = 1:length(F)
        
        if k<10
            imageNum = ['0' num2str(k)];
        else
            imageNum = num2str(k);
        end
        
        filenameDICa = [sampleNuma '_00-' imageNum '-Mesh'];
        filenameDICb = [sampleNumb '_00-' imageNum '-Mesh'];
        pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
            'examples','identification','materialParticleBoard','resultsDIC');
        load(fullfile(pathnameDIC,filenameDICa));
        
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx = (X+Job.ROI(1)-1);
        Coordy = (Y+Job.ROI(2)-1);
        TRIa = Mesh.TRI;
        scaleFactor_1a = 45/(max(Coordx)-min(Coordx));
        scaleFactor_2a = 15/(max(Coordy)-min(Coordy));
        scaleFactor_a = mean([scaleFactor_1a scaleFactor_2a]);
        min_Coordy_a = min(Coordy);
        max_Coordx_a = max(Coordx);
        coordx_a = (Coordy-min(Coordy))*scaleFactor_a;
        coordy_a = (max(Coordx)-Coordx)*scaleFactor_a;
        coord_a = [coordx_a coordy_a];
        Ux_a = U(1:2:end);
        Uy_a = U(2:2:end);
        u_exp_a = [Uy_a -Ux_a]'*scaleFactor_a;
        u_exp_a = u_exp_a(:);
        
        % primary line1
        L1x0 = coordx_a;
        L1y0 = coordy_a;
        [L1y0_sort index] = sort(L1y0);
        L1x0_sort = L1x0(index);
        fit10 = polyfit(L1y0_sort,L1x0_sort,1);
        val10=polyval(fit10,L1y0_sort);
        
        % deformed line1
        L1x = coordx_a+Uy_a;
        L1y = coordy_a-Ux_a;
        [L1y_sort index] = sort(L1y);
        L1x_sort = L1x(index);
        fit1 = polyfit(L1y_sort,L1x_sort,1);
        val1=polyval(fit1,L1y_sort);
        
        clear Mesh Job
        load(fullfile(pathnameDIC,filenameDICb));
        
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx = (X+Job.ROI(1)-1);
        Coordy = (Y+Job.ROI(2)-1);
        TRIb = Mesh.TRI;
        scaleFactor_1b = 15/(max(Coordx)-min(Coordx));
        scaleFactor_2b = 30/(max(Coordy)-min(Coordy));
        scaleFactor_b = mean([scaleFactor_1b scaleFactor_2b]);
        scaleFactor = mean([scaleFactor_a scaleFactor_b]);
        min_Coordy_b = min(Coordy);
        max_Coordx_b = max(Coordx);
        coordx_b = (Coordy-min(Coordy))*scaleFactor_b+15;
        coordy_b = (max(Coordx)-Coordx)*scaleFactor_b+30;
        coord_b = [coordx_b coordy_b];
        Ux_b = U(1:2:end);
        Uy_b = U(2:2:end);
        u_exp_b = [Uy_b -Ux_b]'*scaleFactor_b;
        u_exp_b = u_exp_b(:);
        
        % primary line2
        L2x0 = coordx_b;
        L2y0 = coordy_b;
        [L2x0_sort index] = sort(L2x0);
        L2y0_sort = L2y0(index);
        fit20 = polyfit(L2x0_sort,L2y0_sort,1);
        val20=polyval(fit20,L2x0_sort);
        
        % deformed line2
        L2x = coordx_b+Uy_b;
        L2y = coordy_b-Ux_b;
        [L2x_sort index] = sort(L2x);
        L2y_sort = L2y(index);
        fit2 = polyfit(L2x_sort,L2y_sort,1);
        val2=polyval(fit2,L2x_sort);
        
        %-----------------------------
        % primary and deformed angle of junction
        %-----------------------------
        
        t = [val10(end)-val10(1) L1y0_sort(end)-L1y0_sort(1)];
        s = [L2x0_sort(end)-L2x0_sort(1) val20(end)-val20(1)];
        delta0 = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        t = [val1(end)-val1(1) L1y_sort(end)-L1y_sort(1)];
        s = [L2x_sort(end)-L2x_sort(1) val2(end)-val2(1)];
        delta = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        %-----------------------------
        % Reference and deformed mesh
        %-----------------------------
%         set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',16,...
%             'DefaultTextInterpreter','tex');
%         Scal=2;
%         Unitx = '(mm.)';
%         UnitU = '(mm.)';
%         
%         figure('name','best fit line')
%         plot(L1x,L1y,'r*',val1,L1y_sort,'b')
%         hold on
%         plot(L2x,L2y,'k*',L2x_sort,val2,'g')
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         axis equal
%         
%         figure('name','Reference and deformed mesh')
%         fig1 = triplot(TRIa,coordx_a,coordy_a,'r');
%         hold on
%         fig2 = triplot(TRIb,coordx_b,coordy_b,'r');
%         fig3 = triplot(TRIa,coordx_a+Scal*u_exp_a(1:2:end),...
%                        coordy_a+Scal*u_exp_a(2:2:end),'k');
%         fig4 = triplot(TRIb,coordx_b+Scal*u_exp_b(1:2:end),...
%                        coordy_b+Scal*u_exp_b(2:2:end),'k');
%         axis equal
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         hold off
        
        Var_angle(k) = abs(delta0-delta);
        
        
    end
    
    F_screw{j} = F;
    Var_angle_screw{j} = Var_angle;
end


F_dowel = cell(NumDowel,1);
Var_angle_dowel = cell(NumDowel,1);
for j = 1:NumDowel
    
    sampleNum = ['D' num2str(j)];
    sampleNuma = ['D' num2str(j) 'a'];
    sampleNumb = ['D' num2str(j) 'b'];
    
    F = appliedLoad(sampleNum);
    
    Var_angle = zeros(length(F),1);
    
    for k = 1:length(F)
        
        if k<10
            imageNum = ['0' num2str(k)];
        else
            imageNum = num2str(k);
        end
        
        filenameDICa = [sampleNuma '_00-' imageNum '-Mesh'];
        filenameDICb = [sampleNumb '_00-' imageNum '-Mesh'];
        pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
            'examples','identification','materialParticleBoard','resultsDIC');
        load(fullfile(pathnameDIC,filenameDICa));
        
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx = (X+Job.ROI(1)-1);
        Coordy = (Y+Job.ROI(2)-1);
        TRIa = Mesh.TRI;
        scaleFactor_1a = 30/(max(Coordx)-min(Coordx));
        scaleFactor_2a = 15/(max(Coordy)-min(Coordy));
        scaleFactor_a = mean([scaleFactor_1a scaleFactor_2a]);
        coordx_a = (Coordy-min(Coordy))*scaleFactor_a;
        coordy_a = (max(Coordx)-Coordx)*scaleFactor_a;
        coord_a = [coordx_a coordy_a];
        Ux_a = U(1:2:end);
        Uy_a = U(2:2:end);
        u_exp_a = [Uy_a -Ux_a]'*scaleFactor_a;
        u_exp_a = u_exp_a(:);
        
        % primary line1
        L1x0 = coordx_a;
        L1y0 = coordy_a;
        [L1y0_sort index] = sort(L1y0);
        L1x0_sort = L1x0(index);
        fit10 = polyfit(L1y0_sort,L1x0_sort,1);
        val10=polyval(fit10,L1y0_sort);
        
        % deformed line1
        L1x = coordx_a+Uy_a;
        L1y = coordy_a-Ux_a;
        [L1y_sort index] = sort(L1y);
        L1x_sort = L1x(index);
        fit1 = polyfit(L1y_sort,L1x_sort,1);
        val1=polyval(fit1,L1y_sort);
        
        clear Mesh Job
        load(fullfile(pathnameDIC,filenameDICb));
        
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx = (X+Job.ROI(1)-1);
        Coordy = (Y+Job.ROI(2)-1);
        TRIb = Mesh.TRI;
        scaleFactor_1b = 15/(max(Coordx)-min(Coordx));
        scaleFactor_2b = 45/(max(Coordy)-min(Coordy));
        scaleFactor_b = mean([scaleFactor_1b scaleFactor_2b]);
        coordx_b = (Coordy-min(Coordy))*scaleFactor_b;
        coordy_b = (max(Coordx)-Coordx)*scaleFactor_b+30;
        coord_b = [coordx_b coordy_b];
        Ux_b = U(1:2:end);
        Uy_b = U(2:2:end);
        u_exp_b = [Uy_b -Ux_b]'*scaleFactor_b;
        u_exp_b = u_exp_b(:);
        
        % primary line2
        L2x0 = coordx_b;
        L2y0 = coordy_b;
        [L2x0_sort index] = sort(L2x0);
        L2y0_sort = L2y0(index);
        fit20 = polyfit(L2x0_sort,L2y0_sort,1);
        val20=polyval(fit20,L2x0_sort);
        
        % deformed line2
        L2x = coordx_b+Uy_b;
        L2y = coordy_b-Ux_b;
        [L2x_sort index] = sort(L2x);
        L2y_sort = L2y(index);
        fit2 = polyfit(L2x_sort,L2y_sort,1);
        val2=polyval(fit2,L2x_sort);
        
        %-----------------------------
        % primary and deformed angle of junction
        %-----------------------------
        
        t = [val10(end)-val10(1) L1y0_sort(end)-L1y0_sort(1)];
        s = [L2x0_sort(end)-L2x0_sort(1) val20(end)-val20(1)];
        delta0 = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        t = [val1(end)-val1(1) L1y_sort(end)-L1y_sort(1)];
        s = [L2x_sort(end)-L2x_sort(1) val2(end)-val2(1)];
        delta = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        %-----------------------------
        % Reference and deformed mesh
        %-----------------------------
%         set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',16,...
%             'DefaultTextInterpreter','tex');
%         Scal=2;
%         Unitx = '(mm.)';
%         UnitU = '(mm.)';
%         
%         figure('name','best fit line')
%         plot(L1x,L1y,'r*',val1,L1y_sort,'b')
%         hold on
%         plot(L2x,L2y,'k*',L2x_sort,val2,'g')
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         axis equal
%         
%         figure('name','Reference and deformed mesh')
%         fig1 = triplot(TRIa,coordx_a,coordy_a,'r');
%         hold on
%         fig2 = triplot(TRIb,coordx_b,coordy_b,'r');
%         fig3 = triplot(TRIa,coordx_a+Scal*u_exp_a(1:2:end),...
%                        coordy_a+Scal*u_exp_a(2:2:end),'k');
%         fig4 = triplot(TRIb,coordx_b+Scal*u_exp_b(1:2:end),...
%                        coordy_b+Scal*u_exp_b(2:2:end),'k');
%         axis equal
%         xlabel(['$x$ ',Unitx],'Interpreter','latex')
%         ylabel(['$y$ ',Unitx],'Interpreter','latex')
%         hold off
        
        Var_angle(k) = abs(delta0-delta);
        
        
    end
    
    F_dowel{j} = F;
    Var_angle_dowel{j} = Var_angle;
end

%% Save variables
save(fullfile(pathname,filename),'NumScrew','NumDowel',...
    'F_screw','Var_angle_screw','F_dowel','Var_angle_dowel');


%% Plot data
if displaySolution
    
    
    figure('name','JoncVis: variation of angle VS applied load')
    for j = 1:NumScrew
        plot(F_screw{j},Var_angle_screw{j},'*-')
        hold on
        xlabel('Applied load (N)','Interpreter',interpreter);
        ylabel('Variation of angle','Interpreter',interpreter);
    end
    
    figure('name','JoncTou: variation of angle VS applied load')
    for j = 1:NumDowel
        plot(F_dowel{j},Var_angle_dowel{j},'*-')
        hold on
        xlabel('Applied load (N)','Interpreter',interpreter);
        ylabel('Variation of angle','Interpreter',interpreter);
    end
    
    
    
end
