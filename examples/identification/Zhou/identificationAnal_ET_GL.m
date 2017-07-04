%% Identification of ET and GL             %%
% Call after Digital Image Correlation RT3 %%
% Coordinate system                        %%
% CORRELI:    right     down               %%
% MATLAB:     right     up                 %%
%%-----------------------------------------%%

clc
clear all
close all

resultData = 'result_ET_GL.mat';
pathResult = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materPropPartiBoards');
if ~exist(pathResult,'dir')
    mkdir(pathResult);
end

for sample_list = 1:20 
    
    sample = 'C';
    samp_num = [sample num2str(sample_list)];
    F = applied_load(samp_num);
    [b,h,dis_supp,Iz] = dim_sample(samp_num);
    
    err = zeros(length(F),1);
    ET = zeros(length(F),1);
    GL = zeros(length(F),1);
    Phi = zeros(length(F),1);
    U0 = zeros(length(F),1);
    V0 = zeros(length(F),1);
    for k=1:length(F)
        
        if k<10
            im = ['0' num2str(k)];
        else
            im = num2str(k);
        end
        
        fileName = [samp_num '_00-' im '-Mesh'];
        pathName = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'examples','identification','Zhou','results_data');
        fullfileName = fullfile(pathName,fileName);
        load(fullfileName);
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        ROI = Job.ROI;
        Ux = U(1:2:end);
        Uy = U(2:2:end);
        Coorx = (X+ROI(1)-1);
        Coory = (Y+ROI(2)-1);
        
        % Change unit in mm and the coordinate system due to the use of Correli RT3
        scale_factor=h/( max(Coorx) - min(Coorx) );
        coorx = (Coory-min(Coory))*scale_factor+dis_supp;
        coory = -(Coorx-1/2*(min(Coorx)+max(Coorx)))*scale_factor;
        ux = Uy*scale_factor;
        uy = -Ux*scale_factor;
        
        %-----------------------------
        % Reference and deformed mesh
        %-----------------------------
        %     set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',16,...
        %     'DefaultTextInterpreter','tex');
        %     Scal=10;
        %     Unitx = '(mm.)';
        %     UnitU = '(mm.)';
        %     figure('name','Reference and deformed mesh')
        %     fig1 = triplot(Mesh.TRI,coorx,coory,'r');
        %     hold on
        %     fig2 = triplot(Mesh.TRI,coorx+Scal*ux,coory+Scal*uy,'k');
        %     axis equal
        %     xlabel(['$x$ ',Unitx],'Interpreter','latex')
        %     ylabel(['$y$ ',Unitx],'Interpreter','latex')
        %     hold off
        
        A = zeros(5);
        B = zeros(5,1);
        for i=1:Mesh.NNodeTot
            xi = coorx(i);
            yi = coory(i);
            ui = ux(i);
            vi = uy(i);
            
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
        
        resu = A\B;
        
        ET(k)= 1/(resu(1)*1000); % Gpa
        GL(k) = 1/(resu(2)); % Mpa
        Phi(k) = resu(3);
        U0(k) = resu(4);
        V0(k) = resu(5);
        
        %% Error between U_ana and U_exp
        ux_ana = @(x,y) -F(k)*(x.^2).*y./(4*ET(k)*1000*Iz)+F(k)*(y.^3)./(12*GL(k)*Iz)+Phi(k)*y+U0(k);
        uy_ana = @(x,y) F(k)*(x.^3)./(12*ET(k)*1000*Iz)-F(k)*h^2*x./(16*GL(k)*Iz)-Phi(k)*x+V0(k);
        
        u_ana = [ux_ana(coorx,coory);uy_ana(coorx,coory)];
        u_exp = [ux;uy];
        
        err(k) = norm(u_ana - u_exp)/norm(u_exp);
    end
    
    disp('+-----------------+')
    fprintf('| Sample %s%2d      |\n',sample,sample_list)
    disp('+-----------------+-----------------+-----------------+')
    disp('| Young''s modulus |  Shear modulus  |  Error between  |')
    disp('|      ET (GPa)   |     GL (MPa)    | U_ana and U_exp |')
    disp('+-----------------+-----------------+-----------------+')
    for k=1:length(F)
        fprintf('| %15.4f | %15.4f | %15.4f |\n',ET(k),GL(k),err(k))
    end
    disp('+-----------------+-----------------+-----------------+')
    
    eval( ['ET_' samp_num '= zeros(1,length(ET));'] );
    eval( ['GL_' samp_num '= zeros(1,length(GL));'] );
    eval( ['Phi_' samp_num '= zeros(1,length(Phi));'] );
    eval( ['U0_' samp_num '= zeros(1,length(U0));'] );
    eval( ['V0_' samp_num '= zeros(1,length(V0));'] );
    eval( ['err_' samp_num '= zeros(1,length(err));'] );
    
    for i=1:length(F)
        eval( ['ET_' samp_num '(i)=' num2str(ET(i)) ';'] );
        eval( ['GL_' samp_num '(i)=' num2str(GL(i)) ';'] );
        eval( ['Phi_' samp_num '(i)=' num2str(Phi(i)) ';'] );
        eval( ['U0_' samp_num '(i)=' num2str(U0(i)) ';'] );
        eval( ['V0_' samp_num '(i)=' num2str(V0(i)) ';'] );
        eval( ['err_' samp_num '(i)=' num2str(err(i)) ';'] );
    end
    
    if isempty( dir([pathResult,filesep,resultData]) )
        save([pathResult,filesep,resultData],['ET_' samp_num],['GL_' samp_num],['Phi_' samp_num],...
            ['U0_' samp_num],['V0_' samp_num],['err_' samp_num]);
    else
        save([pathResult,filesep,resultData],['ET_' samp_num],['GL_' samp_num],['Phi_' samp_num],...
            ['U0_' samp_num],['V0_' samp_num],['err_' samp_num],'-append');
    end

end


