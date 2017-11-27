%% Plot numerical solution and experimental data %%
%%-----------------------------------------------%%

% clc
clearvars
close all

%% Input data
sampleNum = 'C3'; % sample number
imageNum = '01'; % image number

[b,h,d,Iz] = dimSample(sampleNum);

filenameDIC = [sampleNum '_00-' imageNum '-Mesh'];
pathnameDIC = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','identification','materialParticleBoard','resultsDIC');
disp(['File = ',filenameDIC]);
load(fullfile(pathnameDIC,filenameDIC));

filenameAna = 'data_ET_GL.mat';
filenameNum = 'data_EL_NUL.mat';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','identification','materialParticleBoard');
load(fullfile(pathname,filenameAna));
load(fullfile(pathname,filenameNum));

fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc2'};
renderer = 'OpenGL';

%% Compute numerical solution
[u_exp,coord] = extractCorreli(Job,Mesh,U,h,d);

node = NODE(coord,1:size(coord,1));
elem = Mesh.TRI;
elemtype = 'TRI3';
% option = 'DEFO'; % plane strain
option = 'CONT'; % plane stress
S = MODEL('PLAN');
S = addnode(S,node);
S = addelem(S,elemtype,elem,'option',option);

ET = eval(['ET_' sampleNum '(' imageNum ');'])*1e3; % MPa
GL = eval(['GL_' sampleNum '(' imageNum ');']); % MPa
EL = eval(['EL_' sampleNum '(' imageNum ');']); % MPa
NUL = eval(['NUL_' sampleNum '(' imageNum ');']);
mat = ELAS_ISOT_TRANS('EL',EL,'ET',ET,'NUL',NUL,'GL',GL,'DIM3',h);
mat = setnumber(mat,1);
S = setmaterial(S,mat);
S = final(S);
I = create_boundary(S);
I = final(I);
P = calc_P(S,I);
u_exp_b = P*u_exp;
S = addcl(S,[],'U',u_exp_b);
u_exp_in = freevector(S,u_exp);

% [A,b] = calc_rigi(S);
% b = -b;
% u_num_in = A\b;
% u_num = unfreevector(S,u_num_in);

param = [EL NUL];
[u_num_in,S] = solveThreePointBendingNum(param,S);
u_num = unfreevector(S,u_num_in);

%% Display solutions
ampl = 0;
% ampl = getsize(S)/max(max(abs(u_num)),max(abs(u_exp)))/5;

for i=1:2
    plotSolution(S,u_num,'displ',i,'ampl',ampl);
    mysaveas(pathname,['u_num_' num2str(i)],formats,renderer);
    plotSolution(S,u_exp,'displ',i,'ampl',ampl);
    mysaveas(pathname,['u_exp_' num2str(i)],formats,renderer);
end

% for i=1:3
%     plotSolution(S,u_num,'epsilon',i,'ampl',ampl);
%     mysaveas(pathname,['eps_num_' num2str(i)],formats,renderer);
%     plotSolution(S,u_exp,'epsilon',i,'ampl',ampl);
%     mysaveas(pathname,['eps_exp_' num2str(i)],formats,renderer);
%     
%     plotSolution(S,u_num,'sigma',i,'ampl',ampl);
%     mysaveas(pathname,['sig_num_' num2str(i)],formats,renderer);
%     plotSolution(S,u_exp,'sigma',i,'ampl',ampl);
%     mysaveas(pathname,['sig_exp_' num2str(i)],formats,renderer);
% end
% 
% plotSolution(S,u_num,'epsilon','mises','ampl',ampl);
% mysaveas(pathname,'eps_von_mises_num',formats,renderer);
% plotSolution(S,u_exp,'epsilon','mises','ampl',ampl);
% mysaveas(pathname,'eps_von_mises_exp',formats,renderer);
% 
% plotSolution(S,u_num,'sigma','mises','ampl',ampl);
% mysaveas(pathname,'sig_von_mises_num',formats,renderer);
% plotSolution(S,u_exp,'sigma','mises','ampl',ampl);
% mysaveas(pathname,'sig_von_mises_exp',formats,renderer);
