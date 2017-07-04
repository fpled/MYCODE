%% Identification of EL and nuL              %%
% FEMU Method                                %%
% with ET and GL already identified          %%
% Nonlinear least-squares funtion: lsqnonlin %%
%%-------------------------------------------%%

clc
clear all
close all

%%
loadPath = fullfile(getfemobjectoptions('path'),'MYCODE','results',...
    'identification','materPropPartiBoards');
loadData = 'result_ET_GL.mat';
load([loadPath,filesep,loadData])

lb = [0 0];
ub = [Inf 0.5];
tolX = 1e-14;
tolFun = 1e-14;
display = 'off';

optionslsqnonlin  = optimoptions('lsqnonlin','Display',display,'TolX',tolX,'TolFun',tolFun);
optionsfminsearch = optimset('Display',display,'TolX',tolX,'TolFun',tolFun);
optionsfminunc    = optimoptions('fminunc','Display',display,'TolX',tolX,'TolFun',tolFun,'Algorithm','quasi-newton');
optionsfmincon    = optimoptions('fmincon','Display',display,'TolX',tolX,'TolFun',tolFun);

% for sample_list = 1:20
for sample_list = 3
    sample = 'C';
    num_sam = [sample num2str(sample_list)];
    F = applied_load(num_sam);
    [b,h,dis_supp,Iz] = dim_sample(num_sam);
    
    eval( ['EL_f1_' num_sam '= zeros(1,length(F));' ] );
    eval( ['nuL_f1_' num_sam '= zeros(1,length(F));' ] );
    eval( ['EL_f2_' num_sam '= zeros(1,length(F));' ] );
    eval( ['nuL_f2_' num_sam '= zeros(1,length(F));' ] );
    eval( ['EL_f3_' num_sam '= zeros(1,length(F));' ] );
    eval( ['nuL_f3_' num_sam '= zeros(1,length(F));' ] );
    eval( ['EL_f4_' num_sam '= zeros(1,length(F));' ] );
    eval( ['nuL_f4_' num_sam '= zeros(1,length(F));' ] );
    
    for k=1:length(F)
        % for k = 5
        
        if k<10
            im = ['0' num2str(k)];
        else
            im = num2str(k);
        end
        
        FileName = [num_sam '_00-' im '-Mesh'];
        PathName = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'examples','identification','Zhou','results_data');
        fprintf('\n**************************\n');
        disp(['File = ',FileName]);
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
        
        elemtype = 'TRI3';
        % option = 'DEFO'; % plane strain
        option = 'CONT'; % plane stress
        S = MODEL('PLAN');
        node = NODE([coorx(:),coory(:)],1:numel(coorx));
        S = addnode(S,node);
        elem = Mesh.TRI;
        S = addelem(S,elemtype,elem,'option',option);
        nodes_edgerd  = getnumber(getnode(create_boundary(S)));
        nodes_inter = setdiff(1:Mesh.NNodeTot,nodes_edgerd);
        u_exp_edge = zeros(2*Mesh.NNodeEdge,1);
        u_exp_in = zeros(2*Mesh.NNodeIn,1);
        u_exp_in(1:2:end) = ux(nodes_inter);
        u_exp_in(2:2:end) = uy(nodes_inter);
        u_exp_edge(1:2:end) = ux(nodes_edgerd);
        u_exp_edge(2:2:end) = uy(nodes_edgerd);
        N_node_egde = length(nodes_edgerd);
        if N_node_egde ~= Mesh.NNodeEdge
            warning('Wrong Bord Nodes')
        elseif all(nodes_edgerd' ~= Mesh.NNodeTot-Mesh.NNodeEdge+1:Mesh.NNodeTot)
            warning('Wrong Global Matrix')
        end
        
        eval( ['ET=ET_' num_sam '(' im ')*1000;'] ); %Mpa
        eval( ['GL=GL_' num_sam '(' im ');'] ); %Mpa
        
        EL0 = 100;
        nuL0 = 0.001;
        x0 = [EL0 nuL0];
        
        funlsqnonlin = @(x) funlsqnonlinFEMU(ET,GL,x,u_exp_in,Mesh,coorx,coory,u_exp_edge);
        % funoptim = @(x) funoptimFEMU(ET,GL,x,u_exp_in,Mesh,coorx,coory,u_exp_edge);
        
        t1 = tic;
        [x1,err1,~,exitflag1,output1] = lsqnonlin(@(x) funlsqnonlin(x),x0,lb,ub,optionslsqnonlin);
        fprintf('\nlsqnonlin');
        fprintf('\n---------\n');
        fprintf('EL  = %g MPa\n',x1(1));
        fprintf('nuL = %g\n',x1(2));
        fprintf('err = %g\n',sqrt(err1)/norm(u_exp_in));
        % fprintf('exitflag = %g\n',exitflag1);
        % disp(output1);
        toc(t1)
        eval( ['EL_f1_' num_sam '(' num2str(k) ')=' num2str(x1(1)) ';' ]);
        eval( ['nuL_f1_' num_sam '(' num2str(k) ')=' num2str(x1(2)) ';' ]);
        
        %         t2 = tic;
        %         [x2,err2,exitflag2,output2] = fminsearch(@(x) funoptim(x),x0,optionsfminsearch);
        %         fprintf('\nfminsearch');
        %         fprintf('\n----------\n');
        %         fprintf('EL  = %g MPa\n',x2(1));
        %         fprintf('nuL = %g\n',x2(2));
        %         fprintf('err = %g\n',sqrt(err2)/norm(u_exp_in));
        %         % fprintf('exitflag = %g\n',exitflag2);
        %         % disp(output2);
        %         toc(t2)
        %         eval( ['EL_f2_' num_sam '(' num2str(k) ')=' num2str(x2(1)) ';' ]);
        %         eval( ['nuL_f2_' num_sam '(' num2str(k) ')=' num2str(x2(2)) ';' ]);
        %
        %         t3 = tic;
        %         [x3,err3,exitflag3,output3] = fminunc(@(x) funoptim(x),x0,optionsfminunc);
        %         fprintf('\nfminunc');
        %         fprintf('\n-------\n');
        %         fprintf('EL  = %g MPa\n',x3(1));
        %         fprintf('nuL = %g\n',x3(2));
        %         fprintf('err = %g\n',sqrt(err3)/norm(u_exp_in));
        %         % fprintf('exitflag = %g\n',exitflag3);
        %         % disp(output3);
        %         toc(t3)
        %         eval( ['EL_f3_' num_sam '(' num2str(k) ')=' num2str(x3(1)) ';' ]);
        %         eval( ['nuL_f3_' num_sam '(' num2str(k) ')=' num2str(x3(2)) ';' ]);
        %
        %         t4 = tic;
        %         [x4,err4,exitflag4,output4] = fmincon(@(x) funoptim(x),x0,[],[],[],[],lb,ub,[],optionsfmincon);
        %         fprintf('\nfmincon');
        %         fprintf('\n-------\n');
        %         fprintf('EL  = %g MPa\n',x4(1));
        %         fprintf('nuL = %g\n',x4(2));
        %         fprintf('err = %g\n',sqrt(err4)/norm(u_exp_in));
        %         % fprintf('exitflag = %g\n',exitflag4);
        %         % disp(output4);
        %         toc(t4)
        %         eval( ['EL_f4_' num_sam '(' num2str(k) ')=' num2str(x4(1)) ';' ]);
        %         eval( ['nuL_f4_' num_sam '(' num2str(k) ')=' num2str(x4(2)) ';' ]);
        
    end
    
    dataResult = 'result_EL_nuL.mat';
    pathResult = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','identification','materPropPartiBoards');
    if ~exist(pathResult,'dir')
        mkdir(pathResult);
    end
    
    %     if isempty( dir([pathResult,filesep,dataResult]) )
    %         save([pathResult,filesep,dataResult],...
    %             ['EL_f1_' num_sam],['nuL_f1_' num_sam],...
    %             ['EL_f2_' num_sam],['nuL_f2_' num_sam],...
    %             ['EL_f3_' num_sam],['nuL_f3_' num_sam],...
    %             ['EL_f4_' num_sam],['nuL_f4_' num_sam]);
    %     else
    %         save([pathResult,filesep,dataResult],...
    %             ['EL_f1_' num_sam],['nuL_f1_' num_sam],...
    %             ['EL_f2_' num_sam],['nuL_f2_' num_sam],...
    %             ['EL_f3_' num_sam],['nuL_f3_' num_sam],...
    %             ['EL_f4_' num_sam],['nuL_f4_' num_sam],'-append');
    %     end
    
end