%% Plot analytical solution and experimental data %%
%%------------------------------------------------%%

% clc
clear all
close all

%% Input data
sampleNum = 'C3'; % sample number
imageNum = '01'; % image number

F = appliedLoad(sampleNum);
[b,h,d,Iz] = dimSample(sampleNum);

filenameDIC = [sampleNum '_00-' imageNum '-Mesh'];
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
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Compute analytical solution
X = real(Mesh.Znode);
Y = imag(Mesh.Znode);
Coordx = (X+Job.ROI(1)-1);
Coordy = (Y+Job.ROI(2)-1);
scaleFactor = h/(max(Coordx)-min(Coordx));
coordx = (Coordy-min(Coordy))*scaleFactor+d;
coordy = -(Coordx-1/2*(min(Coordx)+max(Coordx)))*scaleFactor;
coord = [coordx coordy];
Ux = U(1:2:end);
Uy = U(2:2:end);
u_exp = [Uy -Ux]'*scaleFactor;
u_exp = u_exp(:);

node = NODE([coordx,coordy],1:numel(coordx));
elem = Mesh.TRI;
elemtype = 'TRI3';
% option = 'DEFO'; % plane strain
option = 'CONT'; % plane stress
S = MODEL('PLAN');
S = addnode(S,node);
S = addelem(S,elemtype,elem,'option',option);

ET = eval(['ET_' sampleNum '(' imageNum ');'])*1e3; % MPa
GL = eval(['GL_' sampleNum '(' imageNum ');']); % MPa
Phi = eval(['Phi_' sampleNum '(' imageNum ');']);
U0 = eval(['U0_' sampleNum '(' imageNum ');']); % mm
V0 = eval(['V0_' sampleNum '(' imageNum ');']); % mm
mat = ELAS_ISOT_TRANS('EL',[],'ET',ET,'NUL',[],'GL',GL,'DIM3',h);
mat = setnumber(mat,1);
S = setmaterial(S,mat);
S = final(S);

x = [ET GL Phi U0 V0];
u_ana = solveThreePointBendingAna(x,coord,F(str2double(imageNum)),Iz,h);

%% Display solutions
ampl = 0;
% ampl = getsize(S)/max(max(abs(u_ana)),max(abs(u_exp)))/5;

for i=1:2
    plotSolution(S,u_ana,'displ',i,'ampl',ampl);
    mysaveas(pathname,['u_ana_' num2str(i)],formats,renderer);
    plotSolution(S,u_exp,'displ',i,'ampl',ampl);
    mysaveas(pathname,['u_exp_' num2str(i)],formats,renderer);
end

% for i=1:3
%     plotSolution(S,u_ana,'epsilon',i,'ampl',ampl);
%     mysaveas(pathname,['eps_ana_' num2str(i)],formats,renderer);
%     plotSolution(S,u_exp,'epsilon',i,'ampl',ampl);
%     mysaveas(pathname,['eps_exp_' num2str(i)],formats,renderer);
% end
% 
% plotSolution(S,u_ana,'epsilon','mises','ampl',ampl);
% mysaveas(pathname,'eps_von_mises_ana',formats,renderer);
% plotSolution(S,u_exp,'epsilon','mises','ampl',ampl);
% mysaveas(pathname,'eps_von_mises_exp',formats,renderer);
