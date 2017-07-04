%% Identification of ET and GL                %%
% Call after Digital Image Correlation RT3    %%
% Nonlinear least-squares function: lsqnonlin %%
% Coordinate system                           %%
% CORRELI:    right     down                  %%
% MATLAB:     right     up                    %%
%%--------------------------------------------%%

clc
clear all
close all

loadData = 'result_ET_GL.mat';
loadPath = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materPropPartiBoards');
load([loadPath,filesep,loadData])

for sample_list = 1:20
    
    sample = 'C';
    samp_num = [sample num2str(sample_list)];
    
    F = applied_load(samp_num);
    [b,h,dis_supp,Iz] = dim_sample(samp_num);
    
    lb = [0 0 -Inf -Inf -Inf];
    ub = [Inf Inf Inf Inf Inf];
    tolX = 1e-14;
    tolFun = 1e-14;
    display = 'off';
    
    optionslsqnonlin  = optimoptions('lsqnonlin','Display',display,'TolX',tolX,'TolFun',tolFun);
    optionsfminsearch = optimset('Display',display,'TolX',tolX,'TolFun',tolFun);
    optionsfminunc    = optimoptions('fminunc','Display',display,'TolX',tolX,'TolFun',tolFun,'Algorithm','quasi-newton');
    optionsfmincon    = optimoptions('fmincon','Display',display,'TolX',tolX,'TolFun',tolFun);
    
    for k=1:length(F)
              
        if k<10
            im = ['0' num2str(k)];
        else
            im = num2str(k);
        end
        
        FileName = [samp_num '_00-' im '-Mesh'];
        PathName = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'examples','identification','Zhou','results_data');
        FullFileName = fullfile(PathName,FileName);
        load(FullFileName);
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        ROI = Job.ROI;
        Ux = U(1:2:end);
        Uy = U(2:2:end);
        Coorx = (X+ROI(1)-1);
        Coory = (Y+ROI(2)-1);
        
        % Change unit in mm and the coordinate system due to the use of Correli RT3
        ecscale_factor = h/( max(Coorx) - min(Coorx) );
        coorx = (Coory-min(Coory))*ecscale_factor+dis_supp;
        coory = -(Coorx-1/2*(min(Coorx)+max(Coorx)))*ecscale_factor;
        ux = Uy*ecscale_factor;
        uy = -Ux*ecscale_factor;
        u_exp = [ux;uy];
        
        eval( ['ET0=ET_' samp_num '(' im ')*1000/2;'] ); % Mpa
        eval( ['GL0=GL_' samp_num '(' im ')/2;'] ); % Mpa
        eval( ['Phi0=Phi_' samp_num '(' im ')/2;'] );
        eval( ['U00=U0_' samp_num '(' im ')/2;'] );
        eval( ['V00=V0_' samp_num '(' im ')/2;'] );
        
        x0 = [ET0 GL0 Phi0 U00 V00];
        
        fprintf('\n--------------------');
        fprintf('\nFile = %s',FileName);
        fprintf('\n--------------------');
        fprintf('\nAnalytical values');
        fprintf('\n-----------------\n');
        fprintf('ET_ana = %g MPa\n',eval(['ET_' samp_num '(' im ')*1e3']) );
        fprintf('GL_ana = %g MPa\n',eval(['GL_' samp_num '(' im ')']) );
        fprintf('err_ana = %g \n',eval(['err_' samp_num '(' im ')']) );
        
        funlsqnonlin = @(x) funlsqnonlinAna(x,u_exp,coorx,coory,F,k,Iz,Mesh,h);
        % funoptim = @(x) funoptimAna(x,u_exp,coorx,coory,F,k,Iz,Mesh,h);
        
        t1 = tic;
        [x1,err1,~,exitflag1,output1] = lsqnonlin(@(x) funlsqnonlin(x),x0,lb,ub,optionslsqnonlin);
        fprintf('');
        fprintf('\nlsqnonlin');
        fprintf('\n---------\n');
        fprintf('ET  = %g MPa\n',x1(1));
        fprintf('GL  = %g MPa\n',x1(2));
        fprintf('err = %g\n',sqrt(err1)/norm(u_exp));
        % fprintf('exitflag = %g\n',exitflag1);
        % disp(output1);
        toc(t1)
        
        %     t2 = tic;
        %     [x2,err2,exitflag2,output2] = fminsearch(@(x) funoptim(x),x0,optionsfminsearch);
        %     fprintf('\nfminsearch');
        %     fprintf('\n----------\n');
        %     fprintf('ET  = %g MPa\n',x2(1));
        %     fprintf('GL  = %g MPa\n',x2(2));
        %     fprintf('err = %g\n',sqrt(err2)/norm(u_exp));
        %     % fprintf('exitflag = %g\n',exitflag2);
        %     % disp(output2);
        %     toc(t2)
        
        %     t3 = tic;
        %     [x3,err3,exitflag3,output3] = fminunc(@(x) funoptim(x),x0,optionsfminunc);
        %     fprintf('\nfminunc');
        %     fprintf('\n-------\n');
        %     fprintf('ET  = %g MPa\n',x3(1));
        %     fprintf('GL  = %g MPa\n',x3(2));
        %     fprintf('err = %g\n',sqrt(err3)/norm(u_exp));
        %     % fprintf('exitflag = %g\n',exitflag3);
        %     % disp(output3);
        %     toc(t3)
        
        %     t4 = tic;
        %     [x4,err4,exitflag4,output4] = fmincon(@(x) funoptim(x),x0,[],[],[],[],lb,ub,[],optionsfmincon);
        %     fprintf('\nfmincon');
        %     fprintf('\n-------\n');
        %     fprintf('ET  = %g MPa\n',x4(1));
        %     fprintf('GL  = %g MPa\n',x4(2));
        %     fprintf('err = %g\n',sqrt(err4)/norm(u_exp));
        %     % fprintf('exitflag = %g\n',exitflag4);
        %     % disp(output4);
        %     toc(t4)
    end
end