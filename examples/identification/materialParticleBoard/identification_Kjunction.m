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

fontsize = 19;
interpreter = 'latex';
formats = {'fig','epsc2'};

%% Identification

NumScrew = 4;
NumDowel = 2;

F_screw = cell(NumScrew,1);
Var_angle_screw = cell(NumScrew,1);    
for j = 1:NumScrew
    
    sampleNum = ['S' num2str(j)];
    
    F = appliedLoad(sampleNum);
    
    Var_angle = zeros(length(F),1);
    
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
        
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx = (X+Job.ROI(1)-1);
        Coordy = (Y+Job.ROI(2)-1);
        scaleFactor1 = 45/(max(Coordx)-min(Coordx));
        scaleFactor2 = 45/(max(Coordy)-min(Coordy));
        scaleFactor = mean([scaleFactor1 scaleFactor2]);
        coordx = (Coordy-min(Coordy))*scaleFactor;
        coordy = (max(Coordx)-Coordx)*scaleFactor;
        coord = [coordx coordy];
        Ux = U(1:2:end);
        Uy = U(2:2:end);
        u_exp = [Uy -Ux]'*scaleFactor;
        u_exp = u_exp(:);
        
        P1 = find(coordx<15+scaleFactor*Mesh.CharLength/6&...
                  coordy<45+scaleFactor*Mesh.CharLength/6);
        
        P2 = find(coordx>15-scaleFactor*Mesh.CharLength/6&...
                  coordy>30-scaleFactor*Mesh.CharLength/6);
        
        %-----------------------------
        % primary and deformed angle of junction
        %-----------------------------
        % primary angle
        L1x0 = coordx(P1);
        L1y0 = coordy(P1);
        [L1y0_sort index] = sort(L1y0); 
        L1x0_sort = L1x0(index);
        fit10 = polyfit(L1y0_sort,L1x0_sort,1);
        val10=polyval(fit10,L1y0_sort); 
        
        L2x0 = coordx(P2);
        L2y0 = coordy(P2);
        [L2x0_sort index] = sort(L2x0); 
        L2y0_sort = L2y0(index);
        fit20 = polyfit(L2x0_sort,L2y0_sort,1);
        val20=polyval(fit20,L2x0_sort);
        
        t = [val10(end)-val10(1) L1y0_sort(end)-L1y0_sort(1)];
        s = [L2x0_sort(end)-L2x0_sort(1) val20(end)-val20(1)];
        delta0 = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        % deformed angle
        L1x = coordx(P1)+Uy(P1);
        L1y = coordy(P1)-Ux(P1);
        [L1y_sort index] = sort(L1y); 
        L1x_sort = L1x(index);
        fit1 = polyfit(L1y_sort,L1x_sort,1);
        val1=polyval(fit1,L1y_sort); 
        
        L2x = coordx(P2)+Uy(P2);
        L2y = coordy(P2)-Ux(P2);
        [L2x_sort index] = sort(L2x); 
        L2y_sort = L2y(index);
        fit2 = polyfit(L2x_sort,L2y_sort,1);
        val2=polyval(fit2,L2x_sort);
        
        t = [val1(end)-val1(1) L1y_sort(end)-L1y_sort(1)];
        s = [L2x_sort(end)-L2x_sort(1) val2(end)-val2(1)];
        delta = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        %-----------------------------
        % Reference and deformed mesh
        %-----------------------------
%             set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',16,...
%             'DefaultTextInterpreter','tex');
%             Scal=3;
%             Unitx = '(mm.)';
%             UnitU = '(mm.)';
%             
%             figure('name','best fit line')
%             plot(L1x,L1y,'r*',val1,L1y_sort,'r') 
%             hold on
%             plot(L2x,L2y,'k*',L2x_sort,val2,'k')
%             xlabel(['$x$ ',Unitx],'Interpreter','latex')
%             ylabel(['$y$ ',Unitx],'Interpreter','latex')
%             axis equal
%         
%             figure('name','Reference and deformed mesh')
%             fig1 = triplot(Mesh.TRI,coordx,coordy,'r');
%             hold on
%             fig2 = triplot(Mesh.TRI,coordx+Scal*u_exp(1:2:end),...
%                    coordy+Scal*u_exp(2:2:end),'k');
%             axis equal
%             xlabel(['$x$ ',Unitx],'Interpreter','latex')
%             ylabel(['$y$ ',Unitx],'Interpreter','latex')
%             hold off
        
            Var_angle(k) = abs(delta0-delta);
            
            
    end
   
    F_screw{j} = F;
    Var_angle_screw{j} = Var_angle;
end
 

F_dowel = cell(NumDowel,1);
Var_angle_dowel = cell(NumDowel,1);    
for j = 1:NumDowel
    
    sampleNum = ['D' num2str(j)];
    
    F = appliedLoad(sampleNum);
    
    Var_angle = zeros(length(F),1);
    
    for k = 1:length(F)
        
        if k<10
            imageNum = ['0' num2str(k)];
        else
            imageNum = num2str(k);
        end
        
        filenameDIC = [sampleNum '_00-' imageNum '-Mesh'];
        pathnameDIC = [cd '/resultsDIC'];
        load(fullfile(pathnameDIC,filenameDIC));
        
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        Coordx = (X+Job.ROI(1)-1);
        Coordy = (Y+Job.ROI(2)-1);
        scaleFactor1 = 45/(max(Coordx)-min(Coordx));
        scaleFactor2 = 45/(max(Coordy)-min(Coordy));
        scaleFactor = mean([scaleFactor1 scaleFactor2]);
        coordx = (Coordy-min(Coordy))*scaleFactor;
        coordy = (max(Coordx)-Coordx)*scaleFactor;
        coord = [coordx coordy];
        Ux = U(1:2:end);
        Uy = U(2:2:end);
        u_exp = [Uy -Ux]'*scaleFactor;
        u_exp = u_exp(:);
        
        P1 = find(coordx<15+scaleFactor*Mesh.CharLength/6&...
                  coordy<30+scaleFactor*Mesh.CharLength/6);
        
        P2 = find(coordx>-scaleFactor*Mesh.CharLength/6&...
                  coordy>30-scaleFactor*Mesh.CharLength/6);
        
        %-----------------------------
        % primary and deformed angle of junction
        %-----------------------------
        % primary angle
        L1x0 = coordx(P1);
        L1y0 = coordy(P1);
        [L1y0_sort index] = sort(L1y0); 
        L1x0_sort = L1x0(index);
        fit10 = polyfit(L1y0_sort,L1x0_sort,1);
        val10=polyval(fit10,L1y0_sort); 
        
        L2x0 = coordx(P2);
        L2y0 = coordy(P2);
        [L2x0_sort index] = sort(L2x0); 
        L2y0_sort = L2y0(index);
        fit20 = polyfit(L2x0_sort,L2y0_sort,1);
        val20=polyval(fit20,L2x0_sort);
        
        t = [val10(end)-val10(1) L1y0_sort(end)-L1y0_sort(1)];
        s = [L2x0_sort(end)-L2x0_sort(1) val20(end)-val20(1)];
        delta0 = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        % deformed angle
        L1x = coordx(P1)+Uy(P1);
        L1y = coordy(P1)-Ux(P1);
        [L1y_sort index] = sort(L1y); 
        L1x_sort = L1x(index);
        fit1 = polyfit(L1y_sort,L1x_sort,1);
        val1=polyval(fit1,L1y_sort); 
        
        L2x = coordx(P2)+Uy(P2);
        L2y = coordy(P2)-Ux(P2);
        [L2x_sort index] = sort(L2x); 
        L2y_sort = L2y(index);
        fit2 = polyfit(L2x_sort,L2y_sort,1);
        val2=polyval(fit2,L2x_sort);
        
        t = [val1(end)-val1(1) L1y_sort(end)-L1y_sort(1)];
        s = [L2x_sort(end)-L2x_sort(1) val2(end)-val2(1)];
        delta = acosd(abs(t*s')/(norm(t)*norm(s)));
        
        
        %-----------------------------
        % Reference and deformed mesh
        %-----------------------------
%             set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',16,...
%             'DefaultTextInterpreter','tex');
%             Scal=3;
%             Unitx = '(mm.)';
%             UnitU = '(mm.)';
%             
%             figure('name','best fit line')
%             plot(L1x,L1y,'r*',val1,L1y_sort,'r') 
%             hold on
%             plot(L2x,L2y,'k*',L2x_sort,val2,'k')
%             xlabel(['$x$ ',Unitx],'Interpreter','latex')
%             ylabel(['$y$ ',Unitx],'Interpreter','latex')
%             axis equal
%         
%             figure('name','Reference and deformed mesh')
%             fig1 = triplot(Mesh.TRI,coordx,coordy,'r');
%             hold on
%             fig2 = triplot(Mesh.TRI,coordx+Scal*u_exp(1:2:end),...
%                    coordy+Scal*u_exp(2:2:end),'k');
%             axis equal
%             xlabel(['$x$ ',Unitx],'Interpreter','latex')
%             ylabel(['$y$ ',Unitx],'Interpreter','latex')
%             hold off
        
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
   xlabel('Applied load (N)','Interpreter',interpreter,'fontsize',fontsize);
   ylabel('Variation of angle ($^\circ$)','Interpreter',interpreter,'fontsize',fontsize); 
end

figure('name','JoncTou: variation of angle VS applied load')
for j = 1:NumDowel
   plot(F_dowel{j},Var_angle_dowel{j},'*-')
   hold on   
   xlabel('Applied load (N)','Interpreter',interpreter,'fontsize',fontsize);
   ylabel('Variation of angle ($^\circ$)','Interpreter',interpreter,'fontsize',fontsize); 
end



end
