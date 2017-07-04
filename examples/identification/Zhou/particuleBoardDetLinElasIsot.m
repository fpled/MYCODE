%% Particule Board - Determinsitic isotropic linear elasticity problem %%
%%---------------------------------------------------------------------%%
% Call after Digital Image Correlation RT3

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

filename = 'particuleBoardDetLinElasIsot';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc2'};
renderer = 'OpenGL';


% for sample_list = 1:20
for sample_list = 20
    
    sample = 'C';
    sample_number = [sample num2str(sample_list)];
    F = applied_load(sample_number);
    [b,h,d,Iz] = dim_sample(sample_number);
    
    % for im=1:length(F)
    for im=15
        
        FileName = [sample_number '_00-' num2str(im) '-Mesh'];
        PathName = fullfile(getfemobjectoptions('path'),'MYCODE',...
            'examples','identification','Zhou','results_data');
        disp(['File = ',FileName]);
        FullFileName = fullfile(PathName,FileName);
        load(FullFileName);
        X = real(Mesh.Znode);
        Y = imag(Mesh.Znode);
        ROI = Job.ROI;
        Ux = U(1:2:end);
        Uy = U(2:2:end);
        Coordx = (X+ROI(1)-1);
        Coordy = (Y+ROI(2)-1);
		% Change unit in mm and the coordinate system due to the use of Correli RT3
        scale_factor = h/( max(Coordx) - min(Coordx) );
        coordx = (Coordy-min(Coordy))*scale_factor+d;
        coordy = -(Coordx-1/2*(min(Coordx)+max(Coordx)))*scale_factor;
        ux_exp = Uy*scale_factor;
        uy_exp = -Ux*scale_factor;
        u_exp = zeros(2*Mesh.NNodeTot,1);
        u_exp(1:2:end) = ux_exp;
        u_exp(2:2:end) = uy_exp;
        
        %% Problem
        if setProblem
            %% Meshes
            
            elemtype = 'TRI3';
            % option = 'DEFO'; % plane strain
            option = 'CONT'; % plane stress
            S = MODEL('PLAN');
            node = NODE([coordx(:),coordy(:)],1:numel(coordx));
            S = addnode(S,node);
            elem = Mesh.TRI;
            S = addelem(S,elemtype,elem,'option',option);
            
            nodes_edge = getnumber(getnode(create_boundary(S)));
            N_node_egde = length(nodes_edge);
            
            u_exp_edge = zeros(S.dim*N_node_egde,1);
            u_exp_edge(1:2:end) = ux_exp(nodes_edge);
            u_exp_edge(2:2:end) = uy_exp(nodes_edge);
            % ux_exp_edge = ux_exp(nodes_edge);
            % uy_exp_edge = uy_exp(nodes_edge);
            
            %% Materials
            % Poisson ratio
            NU = 0.34;
            % Thickness
            DIM3 = 1;
            % Density
            RHO = 1;
            % Young modulus
            E = 12e9;
            
            % Material
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
            mat = setnumber(mat,1);
            S = setmaterial(S,mat);
            
            %% Dirichlet boundary conditions
            
            S = final(S);
            S = addcl(S,[],'U',u_exp_edge);
            % S = addcl(S,[],'UX',ux_exp_edge);
            % S = addcl(S,[],'UY',uy_exp_edge);
            %% Stiffness matrices and sollicitation vectors
            [A,b] = calc_rigi(S,'noforce');
            b = -b;
            
            %% Save variables
            save(fullfile(pathname,'problem.mat'),'S','b');
        else
            load(fullfile(pathname,'problem.mat'),'S','b');
        end
        
        %% Solution
        if solveProblem
            t = tic;
            u = A\b;
            time = toc(t);
            
            e = calc_epsilon(S,u);
            s = calc_sigma(S,u);
            
            save(fullfile(pathname,'solution.mat'),'u','time','e','s');
        else
            load(fullfile(pathname,'solution.mat'),'u','time','e','s');
        end
        
        %% Outputs
        fprintf('\nSquare specimen\n');
        fprintf(['mesh     : ' elemtype ' elements\n']);
        fprintf('nb elements = %g\n',getnbelem(S));
        fprintf('nb nodes    = %g\n',getnbnode(S));
        fprintf('nb dofs     = %g\n',getnbddl(S));
        fprintf('elapsed time = %f s\n',time);
        fprintf('\n');
        
        %% Display
        if displaySolution
            %% Display domains, boundary conditions and meshes
            % [hD,legD] = plotBoundaryConditions(S,'legend',false);
            % legend(hD,legD)
            ampl = 0.5;
            v_exp = calc_init_dirichlet(S);
            figure('Name','Imposed experimental displacement')
            clf
            h = plot(S,'FaceColor','w','LineWidth',0.5);
            hold on
            [hD,legD] = vectorplot(S,'U',v_exp,ampl,'r','LineWidth',1);
            hold off
            hg = hggroup;
            set([h(:),hD],'Parent',hg);
            axis image
            legend(hD,'U_{exp}',-1);
            % legend(hD,legD)
            
            mysaveas(pathname,'boundary_conditions',formats,renderer);
            
            % plotModel(S,'legend',false);
            % mysaveas(pathname,'mesh',formats,renderer);
            
            plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
            mysaveas(pathname,'mesh',formats,renderer);
            
            ampl = getsize(S)/max(abs(u))/5;
            plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
            mysaveas(pathname,'mesh_deflected',formats,renderer);
            
            figure('Name','Meshes')
            clf
            plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
            plot(S+ampl*unfreevector(S,u),'Color','b','FaceColor','b','FaceAlpha',0.1);
            mysaveas(pathname,'meshes_deflected',formats,renderer);
            
            %% Display solution
            ampl = 0;
            % ampl = getsize(S)/max(abs(u))/5;
            
            for i=1:2
                plotSolution(S,u,'displ',i,'ampl',ampl);
                mysaveas(pathname,['u_' num2str(i)],formats,renderer);
            end
            
            for i=1:3
                plotSolution(S,u,'epsilon',i,'ampl',ampl);
                mysaveas(pathname,['eps_' num2str(i)],formats,renderer);
                
                plotSolution(S,u,'sigma',i,'ampl',ampl);
                mysaveas(pathname,['sig_' num2str(i)],formats,renderer);
            end
            
            plotSolution(S,u,'epsilon','mises','ampl',ampl);
            mysaveas(pathname,'eps_von_mises',formats,renderer);
            
            plotSolution(S,u,'sigma','mises','ampl',ampl);
            mysaveas(pathname,'sig_von_mises',formats,renderer);
            
            % u = unfreevector(S,u);
            %
            % figure('Name','Solution eps_xx')
            % clf
            % plot(e,S+ampl*u,'compo','EPXX')
            % colorbar
            % set(gca,'FontSize',16)
            % mysaveas(pathname,'eps_xx',formats,renderer);
            %
            % figure('Name','Solution eps_yy')
            % clf
            % plot(e,S+ampl*u,'compo','EPYY')
            % colorbar
            % set(gca,'FontSize',16)
            % mysaveas(pathname,'eps_yy',formats,renderer);
            %
            % figure('Name','Solution eps_xy')
            % clf
            % plot(e,S+ampl*u,'compo','EPXY')
            % colorbar
            % set(gca,'FontSize',16)
            % mysaveas(pathname,'eps_xy',formats,renderer);
            
            % figure('Name','Solution sig_xx')
            % clf
            % plot(s,S+ampl*u,'compo','SMXX')
            % colorbar
            % set(gca,'FontSize',16)
            % mysaveas(pathname,'sig_xx',formats,renderer);
            %
            % figure('Name','Solution sig_yy')
            % clf
            % plot(s,S+ampl*u,'compo','SMYY')
            % colorbar
            % set(gca,'FontSize',16)
            % mysaveas(pathname,'sig_yy',formats,renderer);
            %
            % figure('Name','Solution sig_xy')
            % clf
            % plot(s,S+ampl*u,'compo','SMXY')
            % colorbar
            % set(gca,'FontSize',16)
            % mysaveas(pathname,'sig_xy',formats,renderer);
            
        end
    end
end
