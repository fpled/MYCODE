%% Specimen under 3-point bending - Determinsitic linear elasticity problem %%
%%--------------------------------------------------------------------------%%
% Call after Digital Image Correlation RT3 and Identification of ET, GL, EL and NUL

% clc
clearvars
close all

%% Input data
setProblem = true;
solveProblem = true;
displaySolution = true;

filename = 'specimenThreePointBendingDetLinElas';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

filenameAna = 'data_ET_GL.mat';
filenameNum = 'data_EL_NUL.mat';
pathnameIdentification = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathnameIdentification,filenameAna));
load(fullfile(pathnameIdentification,filenameNum));

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc2'};
renderer = 'OpenGL';

sample = 'C';
% for j = 1:20
for j = 3
    
    sampleNum = [sample num2str(j)];
    
    F = appliedLoad(sampleNum);
    [b,h,d,Iz] = dimSample(sampleNum);
    
    % for k = 1:length(F)
    for k = 1
        
        if k<10
            imageNum = ['0' num2str(k)];
        else
            imageNum = num2str(k);
        end

        filenameDIC = [sampleNum '_00-' num2str(imageNum) '-Mesh'];
        pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
            'examples','identification','materialParticleBoard','resultsDIC');
        disp(['File = ',filenameDIC]);
        load(fullfile(pathnameDIC,filenameDIC));
        
        %% Problem
        if setProblem
            %% Meshes
            [u_exp,coord] = extractCorreli(Job,Mesh,U,h,d);
            node = NODE(coord,1:size(coord,1));
            elem = Mesh.TRI;
            elemtype = 'TRI3';
            % option = 'DEFO'; % plane strain
            option = 'CONT'; % plane stress
            S = MODEL('PLAN');
            S = addnode(S,node);
            S = addelem(S,elemtype,elem,'option',option);
            
            %% Materials
            % Thickness
            DIM3 = h;
            % Density
            RHO = 1;
            
            % Material Symmetry
            materialSym = 'isotTrans';
            
            switch lower(materialSym)
                case 'isot'
                    % Young modulus
                    E = eval(['ET_' sampleNum '(' imageNum ')*1e3']); % MPa
                    % Shear modulus
                    G = eval(['GL_' sampleNum '(' imageNum ')']); % MPa
                    % Poisson ratio
                    NU = E/G/2-1;
                    % Material
                    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',DIM3);
                case 'isottrans'
                    % Transverse Young modulus
                    ET = eval(['ET_' sampleNum '(' imageNum ')'])*1e3; % MPa
                    % Longitudinal shear modulus
                    GL = eval(['GL_' sampleNum '(' imageNum ')']); % MPa
                    % Longitudinal Young modulus
                    EL = eval(['EL_' sampleNum '(' imageNum ');']); % MPa
                    % Longitudinal Poisson ratio
                    NUL = eval(['NUL_' sampleNum '(' imageNum ');']);
                    % Transverse Poisson ratio
                    NUT = 0.25;
                    % Material
                    mat = ELAS_ISOT_TRANS('EL',EL,'ET',ET,'NUL',NUL,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',DIM3);
                otherwise
                    error('Wrong material symmetry !')
            end
            mat = setnumber(mat,1);
            S = setmaterial(S,mat);
            
            %% Dirichlet boundary conditions
            S = final(S);
            
            I = create_boundary(S);
            I = final(I);
            
            P = calc_P(S,I);
            u_exp_b = P*u_exp;
            S = addcl(S,[],'U',u_exp_b);
            
            %% Stiffness matrices and sollicitation vectors
            [A,b] = calc_rigi(S);
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
            % legend(hD,legD,'Location','northeastoutside')
            ampl = 0.5;
            v_exp = calc_init_dirichlet(S);
            figure('Name','Imposed experimental displacement')
            clf
            h = plot(S,'FaceColor','w','LineWidth',0.5);
            hold on
            [hD,legD] = vectorplot(S,'U',v_exp,ampl,'r','LineWidth',1);
            hold off
            set(gca,'FontSize',fontsize)
            hg = hggroup;
            set([h(:),hD],'Parent',hg);
            axis image
            l = legend(hD,'$U_{\mathrm{exp}}$','Location','northeastoutside');
            % l = legend(hD,legD,'Location','northeastoutside');
            set(l,'Interpreter',interpreter);
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
            % set(gca,'FontSize',fontsize)
            % mysaveas(pathname,'eps_xx',formats,renderer);
            %
            % figure('Name','Solution eps_yy')
            % clf
            % plot(e,S+ampl*u,'compo','EPYY')
            % colorbar
            % set(gca,'FontSize',fontsize)
            % mysaveas(pathname,'eps_yy',formats,renderer);
            %
            % figure('Name','Solution eps_xy')
            % clf
            % plot(e,S+ampl*u,'compo','EPXY')
            % colorbar
            % set(gca,'FontSize',fontsize)
            % mysaveas(pathname,'eps_xy',formats,renderer);
            
            % figure('Name','Solution sig_xx')
            % clf
            % plot(s,S+ampl*u,'compo','SMXX')
            % colorbar
            % set(gca,'FontSize',fontsize)
            % mysaveas(pathname,'sig_xx',formats,renderer);
            %
            % figure('Name','Solution sig_yy')
            % clf
            % plot(s,S+ampl*u,'compo','SMYY')
            % colorbar
            % set(gca,'FontSize',fontsize)
            % mysaveas(pathname,'sig_yy',formats,renderer);
            %
            % figure('Name','Solution sig_xy')
            % clf
            % plot(s,S+ampl*u,'compo','SMXY')
            % colorbar
            % set(gca,'FontSize',fontsize)
            % mysaveas(pathname,'sig_xy',formats,renderer);
            
        end
    end
end
