%% Plot analytical solution and experimental data %%
%%------------------------------------------------%%

% clc
clearvars
close all

%% Input data
j = 14; % sample number
k = 6; % image number
numSample = ['B' num2str(j)];
numImage = num2str(k,'%02d'); 

F = appliedLoad(numSample);
[b,h,d,Iz] = dimSample(numSample);

filenameDIC = [numSample '_00-' numImage '-Mesh'];
pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','materialParticleBoard','resultsDIC');
disp(['File = ',filenameDIC]);
load(fullfile(pathnameDIC,filenameDIC));

filename = 'data_ET_GL.mat';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathname,filename));

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Compute analytical solution
[u_exp,coord] = extractCorreli(Job,Mesh,U,h,d); % [mm]

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
R0 = R0_data{j}(k); % [rad]
U0 = U0_data{j}(k); % [mm]
V0 = V0_data{j}(k); % [mm]
mat = ELAS_ISOT_TRANS('AXISL',[0;1],'AXIST',[1;0],'EL',[],'ET',ET,'NUL',[],'GL',GL,'DIM3',h);
mat = setnumber(mat,1);
S = setmaterial(S,mat);
S = final(S);

x = [ET GL R0 U0 V0];
u = solveThreePointBendingAna(x,coord,F(k),Iz,h); % [mm]

%% Display solutions
ampl = 0;
% ampl = getsize(S)/max(max(abs(u)),max(abs(u_exp)))/5;

for i=1:2
    plotSolution(S,u,'displ',i,'ampl',ampl);
    set(gcf,'position',[100,100,500,150])
    mysaveas(pathname,['u_' num2str(i) '_ana'],formats,renderer);
    plotSolution(S,u_exp,'displ',i,'ampl',ampl);
    set(gcf,'position',[100,100,500,150])
    mysaveas(pathname,['u_' num2str(i) '_exp'],formats,renderer);
end

% for i=1:3
%     plotSolution(S,u,'epsilon',i,'ampl',ampl);
%     mysaveas(pathname,['eps_' num2str(i) '_ana'],formats,renderer);
%     plotSolution(S,u_exp,'epsilon',i,'ampl',ampl);
%     mysaveas(pathname,['eps_' num2str(i) '_exp'],formats,renderer);
%     
%     plotSolution(S,u,'sigma',i,'ampl',ampl);
%     mysaveas(pathname,['sig_' num2str(i) '_ana'],formats,renderer);
%     plotSolution(S,u_exp,'sigma',i,'ampl',ampl);
%     mysaveas(pathname,['sig_' num2str(i) '_exp'],formats,renderer);
% end
% 
% plotSolution(S,u,'epsilon','mises','ampl',ampl);
% mysaveas(pathname,'eps_von_mises_ana',formats,renderer);
% plotSolution(S,u_exp,'epsilon','mises','ampl',ampl);
% mysaveas(pathname,'eps_von_mises_exp',formats,renderer);
% 
% plotSolution(S,u,'sigma','mises','ampl',ampl);
% mysaveas(pathname,'sig_von_mises_ana',formats,renderer);
% plotSolution(S,u_exp,'sigma','mises','ampl',ampl);
% mysaveas(pathname,'sig_von_mises_exp',formats,renderer);
