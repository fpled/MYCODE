%% Phase-field fracture model - deterministic linear elasticity problem %%
%  Plate with a central circular hole under compression                 %%
%%----------------------------------------------------------------------%%
% [Noel, Pled, Chevalier, Wilquin, EFM, 2025] (anisotropic phase-field model of He et al.)

% clc
clearvars
close all
% myparallel('start');

%% Experimental data from [Noel, Pled, Chevalier, Wilquin, EFM, 2025]
% filenameData   = '_data_instron';
filenameTests  = '_tests';
filenameAngles = '_angles';
filenameElas   = '_elastic';
filenameGc     = '_gc';

pathnameExp = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'examples','phasefield','dataPlateWithHoleWoodNoel');
% filenameData    = fullfile(pathnameExp,filenameData);
filenameTests   = fullfile(pathnameExp,filenameTests);
filenameAngles = fullfile(pathnameExp,filenameAngles);
filenameElas   = fullfile(pathnameExp,filenameElas);
filenameGc     = fullfile(pathnameExp,filenameGc);

% optsData  = detectImportOptions(filenameData);
optsTests  = detectImportOptions(filenameTests);
optsAngles = detectImportOptions(filenameAngles);
optsElas   = detectImportOptions(filenameElas);
optsGc     = detectImportOptions(filenameGc);

% optsData.VariableNames           = {'Test','Temps','Deplacements','Forces'};
optsTests.SelectedVariableNames  = {'Test','Tree_ringAngle_deg_','CrackForce_kN_','CrackImage','CrackThreshold','Forces_kN_','Displacements_mm_','ForcesRedim_kN_','DisplacementsRedim_mm_'};
optsAngles.SelectedVariableNames = {'Test','Angles'};
optsElas.SelectedVariableNames   = {'EL','ET','GL','vL'};
optsGc.SelectedVariableNames     = {'threshold','inc0','inc1','f_exp_kN_','f_kN_','err','Gc_mJ_mm2_','l0_mm_'};

% T_Data   = readtable(filenameData,optsData); T_Data.Test = fillmissing(T_Data.Test, 'previous'); T_Data = rmmissing(T_Data);
T_Tests  = readtable(filenameTests,optsTests);
T_Angles = readtable(filenameAngles,optsAngles);
T_Elas   = readtable(filenameElas,optsElas);
T_Gc     = readtable(filenameGc,optsGc);

% T_DataGroup = groupcounts(T_Data,'Test');
% [G,ID] = findgroups(T_Data.Test);
% % 1x18 cell: each cell is a column vector
% t_exp_data  = splitapply(@(x){x}, T_Data.Temps, G); % [s]
% ut_exp_data_brut = splitapply(@(x){x}, T_Data.Deplacements, G); % [mm] (brut)
% ft_exp_data_brut = splitapply(@(x){x}, T_Data.Forces, G); % [kN] (brut)
% [fmax_exp_data_brut,idmax_exp_brut] = splitapply(@max, T_Data.Forces, G); % [kN] (brut)
% udmax_exp_brut = arrayfun(@(i) ut_exp_data_brut{i}(idmax_exp_brut(i)),1:numel(ID)).'; % [mm] (brut)

vars = {'Forces_kN_','Displacements_mm_','ForcesRedim_kN_','DisplacementsRedim_mm_'};
for k = 1:numel(vars)
    s = T_Tests.(vars{k});
    s = erase(erase(s,'['),']');   % remove brackets
    s = replace(s, newline, ' ');  % remove line breaks
    T_Tests.(vars{k}) = cellfun(@(t) sscanf(t,'%f'), s, 'UniformOutput', false);
end

test_data           = T_Tests.Test; % test names
angle_data          = T_Tests.Tree_ringAngle_deg_; % [deg]
fc_exp_data         = T_Tests.CrackForce_kN_; % [kN]
image_data          = T_Tests.CrackImage;
crackthreshold_data = T_Tests.CrackThreshold;
ft_exp_data_noredim = T_Tests.Forces_kN_; % [kN] (no redim) (same as ft_exp_data_brut for #2, #13 and #16: ft_exp_data_noredim{i}(1) = 0 for i=1:18, but ft_exp_data_brut{i}(1) = 0 only for i=[2,13,16])
ut_exp_data_noredim = T_Tests.Displacements_mm_; % [mm] (no redim) (same as ut_exp_data_brut except for #14: ut_exp_data_noredim{14} = ut_exp_data_brut{14} - 1e-4, since ut_exp_data_brut{14}(1) = 1e-4 mm instead of 0 mm)
ft_exp_data         = T_Tests.ForcesRedim_kN_; % [kN] (redim)
ut_exp_data         = T_Tests.DisplacementsRedim_mm_; % [mm] (redim)

numTests = numel(test_data); % number of tests

% Redimensionning
% f1 = 15; f2 = 25; % [kN]
% idf1_exp = arrayfun(@(i) find(ft_exp_data_noredim{i}>=f1,1),1:numTests).'; % [kN] (no redim)
% idf2_exp = arrayfun(@(i) find(ft_exp_data_noredim{i}>=f2,1),1:numTests).'; % [kN] (no redim)
% f1_exp = arrayfun(@(i) ft_exp_data_noredim{i}(idf1_exp(i)),1:numTests).'; % [kN] (no redim)
% f2_exp = arrayfun(@(i) ft_exp_data_noredim{i}(idf2_exp(i)),1:numTests).'; % [kN] (no redim)
% u1_exp = arrayfun(@(i) ut_exp_data_noredim{i}(idf1_exp(i)),1:numTests).'; % [mm] (no redim)
% u2_exp = arrayfun(@(i) ut_exp_data_noredim{i}(idf2_exp(i)),1:numTests).'; % [mm] (no redim)
% % p_exp  = arrayfun(@(i) polyfit([u1_exp(i),u2_exp(i)],[f1_exp(i),f2_exp(i)],1),1:numTests,'UniformOutput',false).';
% % p_exp  = arrayfun(@(u1,u2,f1,f2) polyfit([u1,u2],[f1,f2],1),u1_exp,u2_exp,f1_exp,f2_exp,'UniformOutput',false);
% % a_exp  = cellfun(@(p) p(1),p_exp);
% % b_exp  = cellfun(@(p) p(2),p_exp);
% a_exp = (f2_exp-f1_exp)./(u2_exp-u1_exp);
% b_exp = f1_exp - a_exp .* u1_exp;
% ft_exp_data_redim  = ft_exp_data_noredim; % [kN] (redim)
% ut_exp_data_redim  = ut_exp_data_noredim; % [mm] (redim)
% for i=1:numTests
%     ft_i = ft_exp_data_redim{i};
%     ut_i = ut_exp_data_redim{i};
%     idf1 = idf1_exp(i);
%     ft_i(1:idf1) = a_exp(i) * ut_exp_data_noredim{i}(1:idf1) + b_exp(i);
%     idsup0 = find(ft_i>=0,1);
%     ft_i = ft_i(idsup0:end);
%     ut_i = ut_i(idsup0:end);
%     ft_exp_data_redim{i} = ft_i - ft_i(1);
%     ut_exp_data_redim{i} = ut_i - ut_i(1);
% end
% diff_ft_exp_data = cellfun(@(ft_redim,ft) max(abs(ft_redim-ft)./ft),ft_exp_data_redim,ft_exp_data);
% diff_ut_exp_data = cellfun(@(ut_redim,ut) max(abs(ut_redim-ut)./ut),ut_exp_data_redim,ut_exp_data);

[fmax_exp_data_noredim,idmax_exp_noredim] = cellfun(@max,ft_exp_data_noredim); % [kN] (no redim)
udmax_exp_data_noredim = arrayfun(@(i) ut_exp_data_noredim{i}(idmax_exp_noredim(i)),1:numTests).'; % [mm] (no redim)

[fmax_exp_data,idmax_exp_data] = cellfun(@max,ft_exp_data); % [kN] (redim)
udmax_exp_data = arrayfun(@(i) ut_exp_data{i}(idmax_exp_data(i)),1:numTests).'; % [mm] (redim)

% idc_exp_data_noredim = arrayfun(@(i) find(ft_exp_data_noredim{i} == fc_exp_data(i),1),1:numTests,'UniformOutput',false).'; % [kN] (no redim)
idc_exp_data_noredim = arrayfun(@(i) find(abs(ft_exp_data_noredim{i}-fc_exp_data(i))<=eps(fc_exp_data(i)),1),1:numTests).'; % [kN] (no redim)
udc_exp_data_noredim = arrayfun(@(i) ut_exp_data_noredim{i}(idc_exp_data_noredim(i)),1:numTests).'; % [mm] (no redim)

idc_exp_data = arrayfun(@(i) find(ft_exp_data{i}>=fc_exp_data(i),1),1:numTests).'; % [kN] (redim)
% udc_exp_data = arrayfun(@(i) ut_exp_data{i}(idc_exp_data(i)),1:numTests).'; % [mm] (redim)
udc_exp_data = arrayfun(@(i) interp1(ft_exp_data{i}([idc_exp_data(i)-1,idc_exp_data(i)]),ut_exp_data{i}([idc_exp_data(i)-1,idc_exp_data(i)]),fc_exp_data(i),'linear'),1:numTests).'; % [mm] (redim)

sample_data = T_Angles.Test; % sample numbers from 0 to 17
% angle_data  = T_Angles.Angles; % [deg] (same as T_Tests.Tree_ringAngle_deg_)

EL_data  = T_Elas.EL; % [MPa]
ET_data  = T_Elas.ET; % [MPa]
GL_data  = T_Elas.GL; % [MPa]
NUL_data = T_Elas.vL;

threshold_data = T_Gc.threshold;
inc0_data      = T_Gc.inc0; % [mm]
inc1_data      = T_Gc.inc1; % [mm]
% fc_exp_data    = T_Gc.f_exp_kN_; % [kN] (same as T_Tests.CrackForce_kN_ within tolerance 1e-14)
fc_data        = T_Gc.f_kN_; % [kN]
% err_fc_data    = T_Gc.err; % relative error (abs(fc_data-fc_exp_data)./fc_exp_data)
gc_data        = T_Gc.Gc_mJ_mm2_; % [mJ/mm^2]=[kN/m]
l_data         = T_Gc.l0_mm_; % [mm]

samplesA = sample_data(abs(angle_data)<=37.5); % samplesA = [0,5,7:8,12:16]'; % samples from crossbeam A
samplesB = sample_data(abs(angle_data)>37.5);  % samplesB = [1:4,6,9:11,17]'; % samples from crossbeam B
samples  = sample_data; % samples = union(samplesA,samplesB); % samples = (0:17)'; % set of all samples
samplesRemoved = [2,8,13,17]'; % set of removed samples
samplesKept  = setdiff(samples,samplesRemoved); % samplesKept = [1,3:7,9:16]'; % set of kept samples
samplesKeptA = setdiff(samplesA,samplesRemoved); % samplesKeptA = [0,5,7,12,14:16]'; % selected samples from crossbeam A
samplesKeptB = setdiff(samplesB,samplesRemoved); % samplesKeptB = [1,3:4,6,9:11]'; % selected samples from crossbeam B
numSamples = numel(samples); % number of samples/tests
numSample = 4; % sample number

%% Statistics of elastic and fracture properties: mean value and biased standard deviation
filenameStat = fullfile(pathname,'stat_data.txt');
fidStat = fopen(filenameStat,'w');
fprintf(fidStat,['Samples : all = [' sprintf(' %d ',samples) '], kept = [' sprintf(' %d ',samplesKept) ']\n']);
fprintf(fidStat,['Group A : all = [' sprintf(' %d ',samplesA) '], kept = [' sprintf(' %d ',samplesKeptA) ']\n']);
fprintf(fidStat,['Group B : all = [' sprintf(' %d ',samplesB) '], kept = [' sprintf(' %d ',samplesKeptB) ']\n']);

fprintf(fidStat,'\n');
fprintf(fidStat,'Longitudinal Young modulus EL\n');
fprintf(fidStat,'All     : mean = %g GPa, std = %g MPa, cv = %g %%\n',mean(EL_data(samplesKept+1))*1e-3,std(EL_data(samplesKept+1),1),std(EL_data(samplesKept+1),1)/mean(EL_data(samplesKept+1))*1e2);
fprintf(fidStat,'Group A : mean = %g GPa, std = %g MPa, cv = %g %%\n',mean(EL_data(samplesKeptA+1))*1e-3,std(EL_data(samplesKeptA+1),1),std(EL_data(samplesKeptA+1),1)/mean(EL_data(samplesKeptA+1))*1e2);
fprintf(fidStat,'Group B : mean = %g GPa, std = %g MPa, cv = %g %%\n',mean(EL_data(samplesKeptB+1))*1e-3,std(EL_data(samplesKeptB+1),1),std(EL_data(samplesKeptB+1),1)/mean(EL_data(samplesKeptB+1))*1e2);

fprintf(fidStat,'\n');
fprintf(fidStat,'Transverse Young modulus ET\n');
fprintf(fidStat,'All     : mean = %g MPa, std = %g MPa, cv = %g %%\n',mean(ET_data(samplesKept+1)),std(ET_data(samplesKept+1),1),std(ET_data(samplesKept+1),1)/mean(ET_data(samplesKept+1))*1e2);
fprintf(fidStat,'Group A : mean = %g MPa, std = %g MPa, cv = %g %%\n',mean(ET_data(samplesKeptA+1)),std(ET_data(samplesKeptA+1),1),std(ET_data(samplesKeptA+1),1)/mean(ET_data(samplesKeptA+1))*1e2);
fprintf(fidStat,'Group B : mean = %g MPa, std = %g MPa, cv = %g %%\n',mean(ET_data(samplesKeptB+1)),std(ET_data(samplesKeptB+1),1),std(ET_data(samplesKeptB+1),1)/mean(ET_data(samplesKeptB+1))*1e2);

fprintf(fidStat,'\n');
fprintf(fidStat,'Longitudinal shear modulus GL\n');
fprintf(fidStat,'All     : mean = %g MPa, std = %g MPa, cv = %g %%\n',mean(GL_data(samplesKept+1)),std(GL_data(samplesKept+1),1),std(GL_data(samplesKept+1),1)/mean(GL_data(samplesKept+1))*1e2);
fprintf(fidStat,'Group A : mean = %g MPa, std = %g MPa, cv = %g %%\n',mean(GL_data(samplesKeptA+1)),std(GL_data(samplesKeptA+1),1),std(GL_data(samplesKeptA+1),1)/mean(GL_data(samplesKeptA+1))*1e2);
fprintf(fidStat,'Group B : mean = %g MPa, std = %g MPa, cv = %g %%\n',mean(GL_data(samplesKeptB+1)),std(GL_data(samplesKeptB+1),1),std(GL_data(samplesKeptB+1),1)/mean(GL_data(samplesKeptB+1))*1e2);

fprintf(fidStat,'\n');
fprintf(fidStat,'Longitudinal Poisson ratio NUL\n');
fprintf(fidStat,'All     : mean = %g, std = %g, cv = %g %%\n',mean(NUL_data(samplesKept+1)),std(NUL_data(samplesKept+1),1),std(NUL_data(samplesKept+1),1)/mean(NUL_data(samplesKept+1))*1e2);
fprintf(fidStat,'Group A : mean = %g, std = %g, cv = %g %%\n',mean(NUL_data(samplesKeptA+1)),std(NUL_data(samplesKeptA+1),1),std(NUL_data(samplesKeptA+1),1)/mean(NUL_data(samplesKeptA+1))*1e2);
fprintf(fidStat,'Group B : mean = %g, std = %g, cv = %g %%\n',mean(NUL_data(samplesKeptB+1)),std(NUL_data(samplesKeptB+1),1),std(NUL_data(samplesKeptB+1),1)/mean(NUL_data(samplesKeptB+1))*1e2);

fprintf(fidStat,'\n');
fprintf(fidStat,'Critical energy release rate Gc\n');
fprintf(fidStat,'All     : mean = %g N/m, std = %g N/m, cv = %g %%\n',mean(gc_data)*1e3,std(gc_data,1)*1e3,std(gc_data,1)/mean(gc_data)*1e2);
fprintf(fidStat,'Group A : mean = %g N/m, std = %g N/m, cv = %g %%\n',mean(gc_data(samplesA+1))*1e3,std(gc_data(samplesA+1),1)*1e3,std(gc_data(samplesA+1),1)/mean(gc_data(samplesA+1))*1e2);
fprintf(fidStat,'Group B : mean = %g N/m, std = %g N/m, cv = %g %%\n',mean(gc_data(samplesB+1))*1e3,std(gc_data(samplesB+1),1)*1e3,std(gc_data(samplesB+1),1)/mean(gc_data(samplesB+1))*1e2);

fprintf(fidStat,'\n');
fprintf(fidStat,'Experimental maximum compression force Fm_exp (no redim)\n');
fprintf(fidStat,'All     : mean = %g kN, std = %g kN, cv = %g %%\n',mean(fmax_exp_data_noredim),std(fmax_exp_data_noredim,1),std(fmax_exp_data_noredim,1)/mean(fmax_exp_data_noredim)*1e2);
fprintf(fidStat,'Group A : mean = %g kN, std = %g kN, cv = %g %%\n',mean(fmax_exp_data_noredim(samplesA+1)),std(fmax_exp_data_noredim(samplesA+1),1),std(fmax_exp_data_noredim(samplesA+1),1)/mean(fmax_exp_data_noredim(samplesA+1))*1e2);
fprintf(fidStat,'Group B : mean = %g kN, std = %g kN, cv = %g %%\n',mean(fmax_exp_data_noredim(samplesB+1)),std(fmax_exp_data_noredim(samplesB+1),1),std(fmax_exp_data_noredim(samplesB+1),1)/mean(fmax_exp_data_noredim(samplesB+1))*1e2);

fprintf(fidStat,'\n');
fprintf(fidStat,'Experimental maximum compression force Fm_exp (redim)\n');
fprintf(fidStat,'All     : mean = %g kN, std = %g kN, cv = %g %%\n',mean(fmax_exp_data),std(fmax_exp_data,1),std(fmax_exp_data,1)/mean(fmax_exp_data)*1e2);
fprintf(fidStat,'Group A : mean = %g kN, std = %g kN, cv = %g %%\n',mean(fmax_exp_data(samplesA+1)),std(fmax_exp_data(samplesA+1),1),std(fmax_exp_data(samplesA+1),1)/mean(fmax_exp_data(samplesA+1))*1e2);
fprintf(fidStat,'Group B : mean = %g kN, std = %g kN, cv = %g %%\n',mean(fmax_exp_data(samplesB+1)),std(fmax_exp_data(samplesB+1),1),std(fmax_exp_data(samplesB+1),1)/mean(fmax_exp_data(samplesB+1))*1e2);

fprintf(fidStat,'\n');
fprintf(fidStat,'Experimental crack initiation force Fc_exp (no redim)\n');
fprintf(fidStat,'All     : mean = %g kN, std = %g kN, cv = %g %%\n',mean(fc_exp_data),std(fc_exp_data,1),std(fc_exp_data,1)/mean(fc_exp_data)*1e2);
fprintf(fidStat,'Group A : mean = %g kN, std = %g kN, cv = %g %%\n',mean(fc_exp_data(samplesA+1)),std(fc_exp_data(samplesA+1),1),std(fc_exp_data(samplesA+1),1)/mean(fc_exp_data(samplesA+1))*1e2);
fprintf(fidStat,'Group B : mean = %g kN, std = %g kN, cv = %g %%\n',mean(fc_exp_data(samplesB+1)),std(fc_exp_data(samplesB+1),1),std(fc_exp_data(samplesB+1),1)/mean(fc_exp_data(samplesB+1))*1e2);

fprintf(fidStat,'\n');
fprintf(fidStat,'Numerical crack initiation force Fc\n');
fprintf(fidStat,'All     : mean = %g kN, std = %g kN, cv = %g %%\n',mean(fc_data),std(fc_data,1),std(fc_data,1)/mean(fc_data)*1e2);
fprintf(fidStat,'Group A : mean = %g kN, std = %g kN, cv = %g %%\n',mean(fc_data(samplesA+1)),std(fc_data(samplesA+1),1),std(fc_data(samplesA+1),1)/mean(fc_data(samplesA+1))*1e2);
fprintf(fidStat,'Group B : mean = %g kN, std = %g kN, cv = %g %%\n',mean(fc_data(samplesB+1)),std(fc_data(samplesB+1),1),std(fc_data(samplesB+1),1)/mean(fc_data(samplesB+1))*1e2);

fprintf(fidStat,'\n');
fprintf(fidStat,'Sample #%d\n',numSample);
fprintf(fidStat,'EL     = %g GPa\n',EL_data(numSample+1)*1e-3);
fprintf(fidStat,'ET     = %g MPa\n',ET_data(numSample+1));
fprintf(fidStat,'GL     = %g MPa\n',GL_data(numSample+1));
fprintf(fidStat,'NUL    = %g\n',NUL_data(numSample+1));
fprintf(fidStat,'Gc     = %g N/m\n',gc_data(numSample+1)*1e3);
fprintf(fidStat,'Fm_exp = %g kN (no redim)\n',fmax_exp_data_noredim(numSample+1));
fprintf(fidStat,'Fm_exp = %g kN (redim)\n',fmax_exp_data(numSample+1));
fprintf(fidStat,'Fc_exp = %g kN (no redim)\n',fc_exp_data(numSample+1));
fprintf(fidStat,'Fc     = %g kN\n',fc_data(numSample+1));
fclose(fidStat);
type(filenameStat) % fprintf('%s', fileread(filenameExp))

%% Input data
setProblem = true;
solveProblem = true;
displayModel = false;
displaySolution = false;
makeMovie = false;
saveParaview = false;

test = true; % coarse mesh
% test = false; % fine mesh

Dim = 2; % space dimension Dim = 2, 3
symmetry = 'Anisot'; % 'Isot' or 'Anisot'. Material symmetry
PFmodel = 'HeFreddi'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT1'; % 'AT1' or 'AT2'
PFsolver = 'HistoryFieldElem'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = Inf; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
critConv = 'Energy'; % 'Solution', 'Residual', 'Energy'
FEmesh = 'Optim'; % 'Unif' or 'Optim'

% PFmodels = {'Bourdin','Amor','Miehe','HeAmor','HeFreddi','Zhang'};
% PFsplits = {'Strain','Stress'};
% PFregularizations = {'AT1','AT2'};
% PFsolvers = {'HistoryFieldElem','BoundConstrainedOptim'};
% maxIters = [1,Inf];

% for iPFmodel=1:length(PFmodels)
% PFmodel = PFmodels{iPFmodel};
% for iPFsplit=1:length(PFsplits)
% PFsplit = PFsplits{iPFsplit};
% for iPFRegularization=1:length(PFregularizations)
% PFregularization = PFregularizations{iPFRegularization};
% for iPFsolver=1:length(PFsolvers)
% PFsolver = PFsolvers{iPFsolver};
% for imaxIter=1:length(maxIters)
% maxIter = maxIters(imaxIter);
% close all

suffix = '';

foldername = ['plateWithHoleWoodNoel_' num2str(Dim) 'D'];
testname = test_data{numSample+1};
filename = ['linElas' symmetry];
filename = [filename PFmodel PFsplit PFregularization PFsolver...
    'MaxIter' num2str(maxIter)];
if maxIter>1
    filename = [filename 'Tol' num2str(tolConv) num2str(critConv)];
end
filename = [filename 'Mesh' FEmesh suffix];

pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
    'results','phasefieldDet',foldername,testname,filename);
if test
    pathname = fullfile(getfemobjectoptions('path'),'MYCODE',...
        'results','phasefieldDet_test',foldername,testname,filename);
end
if ~exist(pathname,'dir')
    mkdir(pathname);
end

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    % [Noel, Pled, Chevalier, Wilquin, EFM, 2025]
    L = 45e-3; % length
    H = 2*L; % height
    e = 20e-3; % thickness
    r = 5e-3; % radius of the hole
    l = L/100; % regularization length (width of the smeared crack)
    if Dim==2
        D = DOMAIN(2,[0.0,0.0],[L,H]);
        C = CIRCLE(L/2,H/2,r);
    elseif Dim==3
        D = DOMAIN(3,[0.0,0.0,0.0],[L,H,e]);
        C = CIRCLE(L/2,H/2,0.0,r);
        C = CYLINDER(C,e); % CYLINDER(L/2,H/2,0.0,r,e);
    end
    
    switch lower(FEmesh)
        case 'unif' % uniform mesh
            cl = l/2;
            if test
                cl = l;
            end
            clD = cl; % characteristic length for domain
            clC = cl; % characteristic length for circular hole
            B = [];
        case 'optim' % optimized mesh
            % [Noel, Pled, Chevalier, Wilquin, EFM, 2025]
            clD = 3*l/2; % characteristic length for domain
            clC = l/2; % characteristic length for circular hole
            if test
                clC = l;
            end
            VIn = clC; VOut = clD;
            switch lower(PFmodel)
                case {'bourdin','amor','heamor'}
                    XMin = 0; XMax = L;
                    YMin = H/2-2*r; YMax = H/2+2*r;
                    % Thickness = H/2-2*r;
                otherwise
                    XMin = L/2-2*r; XMax = L/2+2*r;
                    YMin = 0; YMax = H;
                    % Thickness = L/2-2*r;
            end
            Thickness = 0;
            if Dim==2
                B = struct('VIn',VIn,'VOut',VOut,'XMin',XMin,'XMax',XMax,'YMin',YMin,'YMax',YMax,'Thickness',Thickness);
            elseif Dim==3
                ZMin = 0; ZMax = e;
                B = struct('VIn',VIn,'VOut',VOut,'XMin',XMin,'XMax',XMax,'YMin',YMin,'YMax',YMax,'ZMin',ZMin,'ZMax',ZMax,'Thickness',Thickness);
            end
        otherwise
            error('Wrong FE mesh')
    end
    if Dim==2
        S_phase = gmshDomainWithHole(D,C,clD,clC,fullfile(pathname,'gmsh_plate_with_hole_wood_Noel'),Dim,'Box',B);
        % S_phase = gmsh2femobject(2,fullfile(getfemobjectoptions('path'),'MYCODE','examples','phasefield','gmsh_plate_with_hole_wood_Noel_unif.msh'),2);
        % S_phase = gmsh2femobject(2,fullfile(getfemobjectoptions('path'),'MYCODE','examples','phasefield','gmsh_plate_with_hole_wood_Noel_unif_test.msh'),2);
        % S_phase = gmsh2femobject(2,fullfile(getfemobjectoptions('path'),'MYCODE','examples','phasefield','gmsh_plate_with_hole_wood_Noel_optim.msh'),2);
        % S_phase = gmsh2femobject(2,fullfile(getfemobjectoptions('path'),'MYCODE','examples','phasefield','gmsh_plate_with_hole_wood_Noel_optim_test.msh'),2);
    elseif Dim==3
        S_phase = gmshDomainWithHole(D,C,clD,clC,fullfile(pathname,'gmsh_plate_with_hole_wood_Noel'),Dim,'Box',B,'extrude');
    end
    S = S_phase;
    
    %% Phase-field problem
    %% Material
    % Critical energy release rate (or fracture toughness)
    gc = gc_data(numSample+1)*1e3; % [J/m^2]=[N/m]
    % Small artificial residual stiffness
    % k = 1e-12;
    k = 0;
    
    % Material
    [K,R,Qn] = setphasefieldparam(gc,l,PFregularization);
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'DIM3',e,'PFregularization',PFregularization);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BRight = LINE([L,0.0],[L,H]);
        BLeft = LINE([0.0,0.0],[0.0,H]);
    elseif Dim==3
        BRight = PLANE([L,0.0,0.0],[L,H,0.0],[L,0.0,e]);
        BLeft = PLANE([0.0,0.0,0.0],[0.0,H,0.0],[0.0,0.0,e]);
    end
    
    findddlboundary = @(S_phase) union(findddl(S_phase,'T',BRight),findddl(S_phase,'T',BLeft));
    
    S_phase = final(S_phase);
    
    %% Stiffness matrices and sollicitation vectors
    % a_phase = BILINFORM(1,1,K); % uniform values
    % % a_phase = DIFFUSIONFORM(K);
    % a_phase = setfree(a_phase,0);
    % K_phase = calc_matrix(a_phase,S_phase); % quadorder=0, nbgauss=1
    % % K_phase = calc_matrix(a_phase,S_phase,[],[],'quadorder',2);
    % b_phase = calc_nonhomogeneous_vector(S_phase,K_phase);
    % K_phase = freematrix(S_phase,K_phase);
    
    % r_phase = BILINFORM(0,0,R); % uniform values
    % R_phase = calc_matrix(r_phase,S_phase); % quadorder=2, nbgauss=3
    % A_phase = K_phase + R_phase;
    
    % l_phase = LINFORM(0,Qn); % uniform values
    % l_phase = setfree(l_phase,1);
    % b_phase = -b_phase + calc_vector(l_phase,S_phase);
    
    % [A_phase,b_phase] = calc_rigi(S_phase);
    % b_phase = -b_phase + bodyload(S_phase,[],'QN',Qn);
    
    %% Linear elastic displacement field problem
    %% Materials
    % Option
    % option = 'DEFO'; % plane strain [Romani, Bornert, Leguillon, Roy, Sab, 2015, EJMS], [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Luo, Chen, Wang, Li, 2022, CM]
    option = 'CONT'; % plane stress [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Noel, Pled, Chevalier, Wilquin, EFM, 2025]
    switch lower(symmetry)
        case 'isot' % isotropic material
            % Lame coefficients
            % Young modulus and Poisson ratio
            E = EL_data(numSample+1)*1e6; % [Pa]
            NU = 0.44;
        case 'anisot' % anisotropic material
            % Longitudinal Young modulus
            EL = EL_data(numSample+1)*1e6; % [Pa]
            % Transverse Young modulus
            ET = ET_data(numSample+1)*1e6; % [Pa]
            % Longitudinal shear modulus
            GL = GL_data(numSample+1)*1e6; % [Pa]
            % Longitudinal Poisson ratio
            NUL = NUL_data(numSample+1);
            % Transverse Poisson ratio
            NUT = 0.44;
        otherwise
            error('Wrong material symmetry class');
    end
    % Energetic degradation function
    g = @(d) (1-d).^2;
    % Density
    RHO = 1;
    
    % Material
    d = calc_init_dirichlet(S_phase);
    switch lower(symmetry)
        case 'isot' % isotropic material
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
        case 'anisot' % anisotropic material
            mat = ELAS_ISOT_TRANS('AXISL',[0;1;0],'AXIST',[1;0;0],'EL',EL,'ET',ET,'NUL',NUL,'NUT',NUT,'GL',GL,'RHO',RHO,'DIM3',e);
            mat = setnumber(mat,1);
            S = setoption(S,option);
            S = setmaterial(S,mat);
            Cmat = double(double(calc_opmat(S)));
            mat = ELAS_ANISOT('C',Cmat,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit);
        otherwise
            error('Wrong material symmetry class');
    end
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BU = LINE([0.0,H],[L,H]);
        BL = LINE([0.0,0.0],[L,0.0]);
    elseif Dim==3
        BU = PLANE([0.0,H,0.0],[L,H,0.0],[0.0,H,e]);
        BL = PLANE([0.0,0.0,0.0],[L,0.0,0.0],[0.0,0.0,e]);
    end
    PD = getvertices(D);
    P1 = POINT(PD{1});
    P2 = POINT(PD{2});
    
    addbc = @(S,ud) addbcPlateWithHole(S,ud,BU,BL,P1,P2);
    findddlforce = @(S) findddl(S,'UY',BU);
    
    S = final(S);
    
    ud = 0;
    S = addbc(S,ud);
    
    u = calc_init_dirichlet(S);
    mats = MATERIALS(S);
    for m=1:length(mats)
        mats{m} = setparam(mats{m},'u',u);
    end
    S = actualisematerials(S,mats);
    
    %% Stiffness matrices and sollicitation vectors
    % [A,b] = calc_rigi(S);
    % b = -b;
    
    %% Time scheme
    % [Noel, Pled, Chevalier, Wilquin, EFM, 2025]
    % du = 24e-5 mm during the first stage (until the phase-field reaches the threshold value)
    % du = 2e-5 mm during the last stage (as soon as the phase-field exceeds the threshold value)
    dt0 = 24e-8; % dt0 = inc0_data(numSample)*1e-3; [m]
    dt1 = 2e-8;  % dt1 = inc1_data(numSample)*1e-3; [m]
    tf = 0.31e-3; % tf = 0.31 mm
    dth = 0.2; % dth = threshold_data(numSample);
    if test
        % du = 48e-5 mm during the first stage (until the phase-field reaches the threshold value)
        % du = 6e-5 mm during the last stage (as soon as the phase-field exceeds the threshold value)
        dt0 = 48e-8;
        dt1 = 6e-8;
        tf = 0.31e-3; % tf = 0.31 mm
        dth = 0.2;
    end
    T = struct('dt0',dt0,'dt1',dt1,'tf',tf,'dth',dth);
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','addbc','findddlforce','findddlboundary');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','addbc','findddlforce','findddlboundary');
end

%% Solution
if solveProblem
    tTotal = tic;
    
    displayIter = true;
    displaySol  = false;
    
    switch lower(PFsolver)
        case {'historyfieldelem','historyfieldnode'}
            [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElasThreshold(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,...
                'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
                'displayiter',displayIter,'displaysol',displaySol);
        otherwise
            [dt,ut,ft,~,Edt,Eut,output] = solvePFDetLinElasThreshold(S_phase,S,T,PFsolver,addbc,findddlforce,findddlboundary,...
                'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
                'displayiter',displayIter,'displaysol',displaySol);
    end
    % switch lower(PFsolver)
    %     case {'historyfieldelem','historyfieldnode'}
    %         [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElasPlateWithHoleThreshold(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P1,P2,...
    %             'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
    %             'displayiter',displayIter,'displaysol',displaySol);
    %     otherwise
    %         [dt,ut,ft,~,Edt,Eut,output] = solvePFDetLinElasPlateWithHoleThreshold(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P1,P2,...
    %             'maxiter',maxIter,'tol',tolConv,'crit',critConv,...
    %             'displayiter',displayIter,'displaysol',displaySol);
    % end
    
    T = gettimemodel(dt);
    t = gettevol(T);
    dt_val = getvalue(dt);
    dmaxt = max(dt_val);
    % idc = find(dmaxt>=min(0.75,max(dmaxt)),1);
    idc = find(dmaxt>=min(1,max(dmaxt)),1);
    fc = ft(idc);
    udc = t(idc);
    [fmax,idmax] = max(ft,[],2);
    udmax = t(idmax);
    
    time = toc(tTotal);
    
    save(fullfile(pathname,'solution.mat'),'dt','ut','ft','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','T','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        save(fullfile(pathname,'solution.mat'),'Ht','-append');
    end
else
    load(fullfile(pathname,'solution.mat'),'dt','ut','ft','Edt','Eut','output','dmaxt','fmax','udmax','fc','udc','T','time');
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        load(fullfile(pathname,'solution.mat'),'Ht');
    end
end

ft_exp_noredim = ft_exp_data_noredim{numSample+1};
ut_exp_noredim  = ut_exp_data_noredim{numSample+1};

ft_exp = ft_exp_data{numSample+1};
ut_exp  = ut_exp_data{numSample+1};

fmax_exp_noredim = fmax_exp_data_noredim(numSample+1);
udmax_exp_noredim = udmax_exp_data_noredim(numSample+1);

udc_exp_noredim = udc_exp_data_noredim(numSample+1);

idmax_exp = idmax_exp_data(numSample+1);
fmax_exp  = fmax_exp_data(numSample+1);
udmax_exp = udmax_exp_data(numSample+1);

idc_exp = idc_exp_data(numSample+1);
fc_exp  = fc_exp_data(numSample+1);
udc_exp = udc_exp_data(numSample+1);

% Experimental elastic stiffness for redim data
fe = 15;
ide_exp = find(ft_exp>fe,1)-1;
fte_exp = ft_exp(1:ide_exp); % [kN] (redim)
ute_exp = ut_exp(1:ide_exp); % [mm] (redim)
k_exp = (fte_exp(end)-fte_exp(1))./(ute_exp(end)-ute_exp(1)); % [kN/mm] (redim)
% pe_exp = polyfit([ute_exp(1),ute_exp(end)],[fte_exp(1),fte_exp(end)],1);
% k_exp  = pe_exp(1); % [kN/mm] (redim)

% Numerical elastic stiffness
t = gettevol(T);
ide = find(ft*1e-3>fe,1)-1;
fte = ft(1:ide)*1e-3; % [kN]
ute = t(1:ide)*1e3; % [mm] (redim)
k_mat = (fte(end)-fte(1))./(ute(end)-ute(1)); % [kN/mm]
% pe_mat = polyfit([ute(1),ute(end)],[fte(1),fte(end)],1);
% k_mat  = pe_mat(1); % [kN/mm]

% Rescaled experimental displacements
k_setup = 1/(1/k_exp-1/k_mat);
ut_setup = ft_exp/k_setup;
t_exp = ut_exp - ut_setup;

udmax_setup = ut_setup(idmax_exp);
udmax_exp_rescale = udmax_exp - udmax_setup;

udc_setup = interp1(ft_exp([idc_exp-1,idc_exp]),ut_setup([idc_exp-1,idc_exp]),fc_exp,'linear'); % [mm] (redim)
udc_exp_rescale = udc_exp - udc_setup;
% udc_exp_rescale = interp1(ft_exp([idc_exp-1,idc_exp]),t_exp([idc_exp-1,idc_exp]),fc_exp,'linear'); % [mm] (redim)

%% Outputs
if solveProblem
    filenameResults = fullfile(pathname,'results.txt');
    fid = fopen(filenameResults,'w');
    fprintf(fid,'Plate with hole\n');
    fprintf(fid,'\n');
    fprintf(fid,'Sample #%d\n',numSample);
    fprintf(fid,'\n');
    fprintf(fid,'dim      = %d\n',Dim);
    fprintf(fid,'mat sym  = %s\n',symmetry);
    fprintf(fid,'PF model = %s\n',PFmodel);
    fprintf(fid,'PF split = %s\n',PFsplit);
    fprintf(fid,'PF regularization = %s\n',PFregularization);
    fprintf(fid,'PF solver = %s\n',PFsolver);
    fprintf(fid,'nb elements = %g\n',getnbelem(S));
    fprintf(fid,'nb nodes    = %g\n',getnbnode(S));
    fprintf(fid,'nb dofs     = %g\n',getnbddl(S));
    fprintf(fid,'nb time dofs = %g\n',getnbtimedof(T));
    fprintf(fid,'elapsed time = %f s\n',time);
    
    fprintf(fid,'\n');
    fprintf(fid,'fmax      = %g kN/m\n',fmax*1e-3);
    fprintf(fid,'fmax_exp  = %g kN/m\n',fmax_exp);
    fprintf(fid,'fc        = %g kN/m\n',fc*1e-3);
    fprintf(fid,'fc_exp    = %g kN/m\n',fc_exp);
    fprintf(fid,'udmax     = %g mm\n',udmax*1e3);
    fprintf(fid,'udmax_exp = %g mm\n',udmax_exp_rescale);
    fprintf(fid,'udc       = %g mm\n',udc*1e3);
    fprintf(fid,'udc_exp   = %g mm\n',udc_exp_rescale);
    fclose(fid);
    type(filenameResults) % fprintf('%s', fileread(filenameResults))
end

%% Display
if Dim==2
    facealpha = 0.1;
    facecolor = 'k';
    facecolordef = 'b';
elseif Dim==3
    facealpha = 1;
    facecolor = 'w';
    facecolordef = 'w';
end

if displayModel
    [t,rep] = gettevol(T);
    
    %% Display domains, boundary conditions and meshes
    plotDomain({D,C},'legend',false);
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 0.5;
    v = calc_init_dirichlet(S);
    [hN,legN] = vectorplot(S,'U',v,ampl,'r','LineWidth',1);
    % legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions_displacement',formats,renderer);
    
    [hD_phase,legD_phase] = plotBoundaryConditions(S_phase,'legend',false);
    % legend(hD_phase,legD_phase,'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions_damage',formats,renderer);
    
    % plotModel(S,'legend',false);
    % mysaveas(pathname,'mesh',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor',facecolor,'FaceAlpha',facealpha,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
    
    u = getmatrixatstep(ut,rep(end));
    ampl = getsize(S)/max(abs(u))/20;
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor',facecolordef,'FaceAlpha',facealpha,'legend',false);
    mysaveas(pathname,'mesh_deflected',formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor',facecolor,'FaceAlpha',facealpha);
    plot(S+ampl*unfreevector(S,u),'Color','b','FaceColor',facecolordef,'FaceAlpha',facealpha);
    mysaveas(pathname,'meshes_deflected',formats,renderer);
end

%% Display solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display experimental force-displacement curves (no redim)
    colors = distinguishable_colors(numSamples);
    leg = arrayfun(@(i) num2str(i),samples,'UniformOutput',false).';
    % leg = cell(numSamples,1);
    figure('Name','Experimental forces vs displacement (no redim)')
    clf
    for i=1:numSamples
        plot(ut_exp_data_noredim{i},ft_exp_data_noredim{i},'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
        % leg{i} = [num2str(samples(i))];
        % leg{i} = ['Sample #' num2str(samples(i))];
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    limsy = get(gca,'YLim');
    set(gca,'YLim',[-1,limsy(2)])
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    legend(leg{:},'Location','NorthEastOutside')
    mysaveas(pathname,'forces_displacement_noredim',formats);
    mymatlab2tikz(pathname,'forces_displacement_noredim.tex');
    
    %% Display experimental force-displacement curves (redim)
    % leg = cell(numSamples,1);
    figure('Name','Experimental forces vs displacements (redim)')
    clf
    for i=1:numSamples
        plot(ut_exp_data{i},ft_exp_data{i},'LineStyle','-','Color',colors(i,:),'LineWidth',linewidth)
        % leg{i} = [num2str(samples(i))];
        % leg{i} = ['Sample #' num2str(samples(i))];
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    limsy = get(gca,'YLim');
    set(gca,'YLim',[-1,limsy(2)])
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    legend(leg{:},'Location','NorthEastOutside')
    mysaveas(pathname,'forces_displacement_redim',formats);
    mymatlab2tikz(pathname,'forces_displacement_redim.tex');
    
    %% Display experimental force-displacement curve
    figure('Name','Experimental force vs displacement (no redim)')
    clf
    plot(ut_exp_noredim,ft_exp_noredim,'-k','LineWidth',linewidth)
    hold on
    scatter(udmax_exp_noredim,fmax_exp_noredim,'Marker','x','MarkerEdgeColor','k','LineWidth',linewidth)
    scatter(udc_exp_noredim,fc_exp,'Marker','x','MarkerEdgeColor','r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    limsy = get(gca,'YLim');
    set(gca,'YLim',[-1,limsy(2)])
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    legend('experimental',...
        ['$F_m^{\mathrm{exp}} = ' sprintf('%.02f',fmax_exp_noredim) '$ kN'],...
        ['$F_c^{\mathrm{exp}} = ' sprintf('%.02f',fc_exp) '$ kN'],...
        'Location','SouthEast','Interpreter',interpreter)
    mysaveas(pathname,'force_displacement_noredim',formats);
    mymatlab2tikz(pathname,'force_displacement_noredim.tex');
    
    %% Display experimental force-displacement curve
    figure('Name','Experimental force vs displacement (redim)')
    clf
    plot(ut_exp,ft_exp,'-k','LineWidth',linewidth)
    hold on
    scatter(udmax_exp,fmax_exp,'Marker','x','MarkerEdgeColor','k','LineWidth',linewidth)
    scatter(udc_exp,fc_exp,'Marker','x','MarkerEdgeColor','r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    limsy = get(gca,'YLim');
    set(gca,'YLim',[-1,limsy(2)])
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    legend('experimental',...
        ['$F_m^{\mathrm{exp}} = ' sprintf('%.02f',fmax_exp) '$ kN'],...
        ['$F_c^{\mathrm{exp}} = ' sprintf('%.02f',fc_exp) '$ kN'],...
        'Location','SouthEast','Interpreter',interpreter)
    mysaveas(pathname,'force_displacement_redim',formats);
    mymatlab2tikz(pathname,'force_displacement_redim.tex');

    %% Display numerical and experimental force-displacement curve
    figure('Name','Force vs displacement')
    clf
    plot(t*1e3,ft*1e-3,'-b','LineWidth',linewidth)
    hold on
    plot(t_exp,ft_exp,'--k','LineWidth',linewidth)
    scatter(udmax*1e3,fmax*1e-3,'Marker','+','MarkerEdgeColor','b','LineWidth',linewidth)
    scatter(udmax_exp_rescale,fmax_exp,'Marker','x','MarkerEdgeColor','k','LineWidth',linewidth)
    scatter(udc*1e3,fc*1e-3,'Marker','+','MarkerEdgeColor','r','LineWidth',linewidth)
    scatter(udc_exp_rescale,fc_exp,'Marker','x','MarkerEdgeColor','r','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    legend('numerical','experimental',...
        ['$F_m = ' sprintf('%.02f',fmax*1e-3) '$ kN'],...
        ['$F_m^{\mathrm{exp}} = ' sprintf('%.02f',fmax_exp) '$ kN'],...
        ['$F_c = ' sprintf('%.02f',fc*1e-3) '$ kN'],...
        ['$F_c^{\mathrm{exp}} = ' sprintf('%.02f',fc_exp) '$ kN'],...
        'Location','SouthEast','Interpreter',interpreter)
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    %% Display maximum damage-displacement curve
    figure('Name','Maximum damage vs displacement')
    clf
    plot(t*1e3,dmaxt,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Maximum damage','Interpreter',interpreter)
    mysaveas(pathname,'max_damage_displacement',formats);
    mymatlab2tikz(pathname,'max_damage_displacement.tex');
    
    %% Display energy-displacement curves
    figure('Name','Energies vs displacement')
    clf
    plot(t*1e3,Eut,'-b','LineWidth',linewidth)
    hold on
    plot(t*1e3,Edt,'-r','LineWidth',linewidth)
    plot(t*1e3,Eut+Edt,'-k','LineWidth',linewidth)
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Energy [J]','Interpreter',interpreter)
    legend('$\Psi_u$','$\Psi_c$','$\Psi_{\mathrm{tot}}$',...
        'Location','NorthWest','Interpreter',interpreter)
    mysaveas(pathname,'energies_displacement',formats);
    mymatlab2tikz(pathname,'energies_displacement.tex');
    
    %% Display outputs of iterative resolution
    figure('Name','Number of iterations vs displacement')
    clf
    plot(t*1e3,output.iteration,'-b','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Number of iterations','Interpreter',interpreter)
    mysaveas(pathname,'nb_iterations_displacement',formats);
    mymatlab2tikz(pathname,'nb_iterations_displacement.tex');
    
    figure('Name','Computing time vs displacement')
    clf
    plot(t*1e3,output.time,'-r','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Computing time [s]','Interpreter',interpreter)
    mysaveas(pathname,'cpu_time_displacement',formats);
    mymatlab2tikz(pathname,'cpu_time_displacement.tex');
    
    figure('Name','Error vs displacement')
    clf
    plot(t*1e3,output.error,'-k','LineWidth',linewidth)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Error','Interpreter',interpreter)
    mysaveas(pathname,'error_displacement',formats);
    mymatlab2tikz(pathname,'error_displacement.tex');
    
    %% Display solutions at different instants
    ampl = 0;
    tSnapshots = [0.27 0.28 0.29 0.30]*1e-3;
    rep = arrayfun(@(x) find(t>x-eps,1),tSnapshots);
    rep = [rep,length(T)];
    % tSnapshots = [tSnapshots,gett1(T)];
    % rep = arrayfun(@(x) find(t>x-eps,1),tSnapshots);
    
    for j=1:length(rep)
        dj = getmatrixatstep(dt,rep(j));
        uj = getmatrixatstep(ut,rep(j));
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            Hj = getmatrixatstep(Ht,rep(j));
        end
        
        plotSolution(S_phase,dj);
        mysaveas(pathname,['damage_t' num2str(rep(j))],formats,renderer);
        
        for i=1:Dim
            plotSolution(S,uj,'displ',i,'ampl',ampl);
            mysaveas(pathname,['displacement_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        end
        
        % for i=1:(Dim*(Dim+1)/2)
        %     plotSolution(S,uj,'epsilon',i,'ampl',ampl);
        %     mysaveas(pathname,['epsilon_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        %
        %     plotSolution(S,uj,'sigma',i,'ampl',ampl);
        %     mysaveas(pathname,['sigma_' num2str(i) '_t' num2str(rep(j))],formats,renderer);
        % end
        %
        % plotSolution(S,uj,'epsilon','mises','ampl',ampl);
        % mysaveas(pathname,['epsilon_von_mises_t' num2str(rep(j))],formats,renderer);
        %
        % plotSolution(S,uj,'sigma','mises','ampl',ampl);
        % mysaveas(pathname,['sigma_von_mises_t' num2str(rep(j))],formats,renderer);
        %
        % plotSolution(S,uj,'energyint','local','ampl',ampl);
        % mysaveas(pathname,['internal_energy_density_t' num2str(rep(j))],formats,renderer);
        %
        % if strcmpi(PFsolver,'historyfieldelem')
        %     figure('Name','Solution H')
        %     clf
        %     plot(Hj,S_phase);
        %     colorbar
        %     set(gca,'FontSize',fontsize)
        %     mysaveas(pathname,['internal_energy_density_history_t' num2str(rep(j))],formats,renderer);
        % elseif strcmpi(PFsolver,'historyfieldnode')
        %     plotSolution(S_phase,Hj,'ampl',ampl);
        %     mysaveas(pathname,['internal_energy_density_history_t' num2str(rep(j))],formats,renderer);
        % end
    end
end

%% Display evolution of solutions
if makeMovie
    ampl = 0;
    % ampl = getsize(S)/max(max(abs(getvalue(ut))))/20;
    
    options = {'plotiter',true,'plottime',false};
    framerate = 80;
    
    evolSolution(S_phase,dt,'FrameRate',framerate,'filename','damage','pathname',pathname,options{:});
    % for i=1:Dim
    %     evolSolution(S,ut,'displ',i,'ampl',ampl,'FrameRate',framerate,'filename',['displacement_' num2str(i)],'pathname',pathname,options{:});
    % end
    %
    % for i=1:(Dim*(Dim+1)/2)
    %     evolSolution(S,ut,'epsilon',i,'ampl',ampl,'FrameRate',framerate,'filename',['epsilon_' num2str(i)],'pathname',pathname,options{:});
    %     evolSolution(S,ut,'sigma',i,'ampl',ampl,'FrameRate',framerate,'filename',['sigma_' num2str(i)],'pathname',pathname,options{:});
    % end
    %
    % evolSolution(S,ut,'epsilon','mises','ampl',ampl,'FrameRate',framerate,'filename','epsilon_von_mises','pathname',pathname,options{:});
    % evolSolution(S,ut,'sigma','mises','ampl',ampl,'FrameRate',framerate,'filename','sigma_von_mises','pathname',pathname,options{:});
    % evolSolution(S,ut,'energyint','local','ampl',ampl,'FrameRate',framerate,'filename','internal_energy_density','pathname',pathname,options{:});
    % if strcmpi(PFsolver,'historyfieldelem')
    %     figure('Name','Solution H')
    %     clf
    %     T = setevolparam(T,'colorbar',true,'FontSize',fontsize,options{:});
    %     frame = evol(T,Ht,S_phase,'rescale',true);
    %     saveMovie(frame,'FrameRate',framerate,'filename','internal_energy_density_history','pathname',pathname);
    % elseif strcmpi(PFsolver,'historyfieldnode')
    %     evolSolution(S_phase,Ht,'ampl',ampl,'FrameRate',framerate,'filename','internal_energy_density_history','pathname',pathname,options{:});
    % end
end

%% Save solutions
if saveParaview
    [t,rep] = gettevol(T);
    for i=1:length(T)
        di = getmatrixatstep(dt,rep(i));
        ui = getmatrixatstep(ut,rep(i));
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            Hi = getmatrixatstep(Ht,rep(i));
        end
        
        switch lower(PFsolver)
            case 'historyfieldelem'
                write_vtk_mesh(S,{di,ui},{Hi},...
                    {'damage','displacement'},{'internal energy density history'},...
                    pathname,'solution',1,i-1);
            case 'historyfieldnode'
                write_vtk_mesh(S,{di,ui,Hi},[],...
                    {'damage','displacement','internal energy density history'},[],...
                    pathname,'solution',1,i-1);
            otherwise
                write_vtk_mesh(S,{di,ui},[],...
                    {'damage','displacement'},[],...
                    pathname,'solution',1,i-1);
        end
    end
    make_pvd_file(pathname,'solution',1,length(T));
end

% myparallel('stop');

% end
% end
% end
% end
% end
