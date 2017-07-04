clc
clear all
close all

samp_num = 'C9';

fontsize = 16;
linewidth = 1;
markersize = 36;
interpreter = 'latex';

dataResult = 'result_ET_GL_EL_nuL_mean_value.mat';
pathResult = fullfile(getfemobjectoptions('path'),'MYCODE','results',...
    'identification','materPropPartiBoards');
if ~exist(pathResult,'dir')
   mkdir(pathResult);
end

loadPath = fullfile(getfemobjectoptions('path'),'MYCODE','results',...
                    'identification','materPropPartiBoards');
loadData1 = 'result_ET_GL.mat';
loadData2 = 'result_EL_nuL.mat';
load([loadPath,filesep,loadData1])
load([loadPath,filesep,loadData2])

eval( ['ET_' samp_num '_rest = ET_' samp_num '(8:end);'] );
eval( ['GL_' samp_num '_rest = GL_' samp_num '(8:end);'] ); 
eval( ['EL_' samp_num '_rest = EL_f1_' samp_num '(8:end);'] );
eval( ['nuL_' samp_num '_rest = nuL_f1_' samp_num '(8:end);'] ); 

eval( ['mean_ET_' samp_num '_rest = mean(ET_' samp_num '_rest);'] );
eval( ['std_ET_' samp_num '_rest = std(ET_' samp_num '_rest);'] );
eval( ['mean_GL_' samp_num '_rest = mean(GL_' samp_num '_rest);'] );
eval( ['std_GL_' samp_num '_rest = std(GL_' samp_num '_rest);'] );
eval( ['mean_EL_' samp_num '_rest = mean(EL_' samp_num '_rest);'] );
eval( ['std_EL_' samp_num '_rest = std(EL_' samp_num '_rest);'] );
eval( ['mean_nuL_' samp_num '_rest = mean(nuL_' samp_num  '_rest);'] );
eval( ['std_nuL_' samp_num '_rest = std(nuL_' samp_num '_rest);'] );


figure
eval( ['bar(ET_' samp_num '_rest);'] );
grid on;
set(gca,'FontSize',fontsize)
legend(samp_num,-1);
xlabel({'Image number'},'Interpreter',interpreter);
ylabel({'Young''s modulus $E^T$ (GPa)'},'Interpreter',interpreter);
% mysaveas(pathResult,['ET_' samp_num],'fig');
% mymatlab2tikz(pathResult,['ET_' samp_num '.tex']);


figure
eval( ['bar(GL_' samp_num '_rest);'] );
grid on;
set(gca,'FontSize',fontsize)
legend(samp_num,-1);
xlabel({'Image number'},'Interpreter',interpreter);
ylabel({'Shear modulus $G^L$ (MPa)'},'Interpreter',interpreter);
% mysaveas(pathResult,['GL_' samp_num],'fig');
% mymatlab2tikz(pathResult,['GL_' samp_num '.tex']);


figure
eval( ['bar(EL_' samp_num '_rest);'] );
grid on;
set(gca,'FontSize',fontsize)
legend(samp_num,-1);
xlabel({'Image number'},'Interpreter',interpreter);
ylabel({'Young''s modulus $E^L$ (MPa)'},'Interpreter',interpreter);
% mysaveas(pathResult,['EL_' samp_num],'fig');
% mymatlab2tikz(pathResult,['EL_' samp_num '.tex']);


figure
eval( ['bar(nuL_' samp_num '_rest);'] );
grid on;
set(gca,'FontSize',fontsize)
legend(samp_num,-1);
xlabel({'Image number'},'Interpreter',interpreter);
ylabel({'Poisson''s ratio $\nu^L$'},'Interpreter',interpreter);
% mysaveas(pathResult,['nuL_' samp_num],'fig');
% mymatlab2tikz(pathResult,['nuL_' samp_num '.tex']);


% if isempty( dir([pathResult,filesep,dataResult]) )
% save([pathResult,filesep,dataResult],['mean_ET_' samp_num '_rest'],...
%     ['mean_GL_' samp_num '_rest'],['mean_EL_' samp_num '_rest'],...
%     ['mean_nuL_' samp_num '_rest']);
% else
% save([pathResult,filesep,dataResult],['mean_ET_' samp_num '_rest'],...
%     ['mean_GL_' samp_num '_rest'],['mean_EL_' samp_num '_rest'],...
%     ['mean_nuL_' samp_num '_rest'],'-append');   
% end
