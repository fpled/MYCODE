% clc
clearvars
close all
% rng('default');

%% Input data
setProblem = false;
solveProblem = false;
displayModel = false;
displaySolution = false;
makeMovie = false;
saveParaview = false;

% test = true; % coarse mesh
test = false; % fine mesh

numWorkers = 100;
% numWorkers = 1; maxNumCompThreads(1); % mono-thread computation

% Deterministic model parameters
Dim = 2; % space dimension Dim = 2, 3
symmetry = 'Isot'; % 'Isot', 'MeanIsot', 'Anisot'. Material symmetry
ang = 45; % clockwise material orientation angle around z-axis for anisotopic material [deg]
loading = 'Shear'; % 'Tension' or 'Shear'
PFmodel = 'Miehe'; % 'Isot', 'Amor', 'Miehe', 'He'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 1; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
initialCrack = 'GeometricCrack'; % 'GeometricCrack', 'GeometricNotch', 'InitialPhaseField'
FEmesh = 'Optim'; % 'Unif' or 'Optim'

% Random model parameters
N = 500; % number of samples
% N = numWorkers;
coeff_gc = 1.0;
randMat = struct('delta',0.2,'lcorr',1e-4); % random material parameters model
% aGc = 0;
% bGc = 0;
gc = 2.7e3;
aGc = 0.6*gc;
bGc = 1.4*gc;
% aGc = [0.7,1.2]*gc;
% bGc = [0.8,1.3]*gc;
randPF = struct('aGc',aGc,'bGc',bGc,'lcorr',Inf); % random phase field parameters model

foldername = ['singleEdgeCrack' loading '_' num2str(Dim) 'D'];
filename = ['linElas' symmetry];
if strcmpi(symmetry,'anisot') % anisotropic material
    filename = [filename num2str(ang) 'deg'];
end
filename = [filename PFmodel PFsplit PFregularization PFsolver initialCrack...% 'MaxIter' num2str(maxIter) 'Tol' num2str(tolConv)
    'Mesh' FEmesh '_' num2str(N) 'samples'];
if any(randMat.delta)
    filename = [filename '_RandMat_Delta' num2str(randMat.delta,'_%g') '_Lcorr' num2str(randMat.lcorr,'_%g')];
end
if any(randPF.aGc) && any(randPF.bGc)
    gcbounds = [randPF.aGc(:),randPF.bGc(:)]';
    filename = [filename '_RandPF_Gc' num2str(gcbounds(:)','_%g') '_Lcorr' num2str(randPF.lcorr,'_%g')];
end
filename = [filename '_coeffgc' num2str(coeff_gc,'_%g')];

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','phasefieldSto',foldername,filename);
if test
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','phasefieldSto_test',foldername,filename);
end
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'epsc'};
renderer = 'OpenGL';

load(fullfile(pathname,'problem.mat'),'T');
load(fullfile(pathname,'solution.mat'),'N','ft',...
    'ft_mean','ft_std','ft_ci','probs','time',...
    'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
    'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw');

[t,~] = gettevol(T);
    
%% Display force-displacement curve
figure('Name','Force vs displacement')
clf
plot(t*1e3,ft_mean*((Dim==2)*1e-6+(Dim==3)*1e-3),'-b','LineWidth',linewidth)
hold on
ciplot(ft_ci(1,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),ft_ci(2,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),t*1e3,'b');
alpha(0.2)
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Displacement [mm]','Interpreter',interpreter)
ylabel('Force [kN]','Interpreter',interpreter)
legend('mean function',...
    ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
    'Location','NorthWest','Interpreter',interpreter)
mysaveas(pathname,'force_displacement',formats);
mymatlab2tikz(pathname,'force_displacement.tex');

colors = distinguishable_colors(N);
figure('Name','Forces vs displacement')
clf
for i=1:N
    plot(t*1e3,ft(i,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
    hold on
end
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Displacement [mm]'...,'Interpreter',interpreter
    )
ylabel('Force [kN]'...,'Interpreter',interpreter
    )
mysaveas(pathname,'forces_displacement',formats);
mymatlab2tikz(pathname,'forces_displacement.tex');

