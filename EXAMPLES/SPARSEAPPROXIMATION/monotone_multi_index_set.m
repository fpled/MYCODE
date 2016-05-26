%% Monotone multi-index set %%
%%--------------------------%%

% clc
clear all
close all

%% Input data
filename = 'sparse_approx_monotone_multi_index_set';
pathname = fullfile(getfemobjectoptions('path'),'MYCODE',filesep,'RESULTS',filesep,filename,filesep);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
% set(0,'DefaultFigureVisible','off'); % change the default figure properties of the MATLAB root object

%% Random variables
M = 2; % number of random variables
rv = RVUNIFORM(0,1);
RV = RANDVARS(repmat({rv},1,M));

%% Polynomial chaos basis
PC = POLYCHAOS(RV,0,'typebase',1);
indices = [0 0 0;
    1 0 1;
    0 1 1;
    2 0 2;
    1 1 2;
    0 2 2;
    3 0 3;
    2 1 3;
    0 3 3;
    4 0 4];
PC = setindices(PC,indices,'sort','update');

%% Display multi-index set
figure('Name','Multi-index set')
clf
plot_indices(PC,'dim',1:2,'max','marg','marg_red')
mysaveas(pathname,filename,'fig');
mymatlab2tikz(pathname,[filename '.tex']);
