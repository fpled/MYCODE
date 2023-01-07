%% Phase field fracture model - stochastic linear elasticity problem with single edge crack %%
%%------------------------------------------------------------------------------------------%%
% [Bourdin, Francfort, Marigo, 2000, JMPS] (isotropic phase field model with no split of Bourdin et al.)
% [Miehe, Welschinger, Hofacker, 2010 IJNME] (anisotropic phase field model of Miehe et al.)
% [Miehe, Hofacker, Welschinger, 2010, CMAME] (anisotropic phase field model of Miehe et al.)
% [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME] (anisotropic phase field model of Miehe et al.)
% [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM] (anisotropic phase field model of Miehe et al.)
% [Ambati, Gerasimov, De Lorenzis, 2015, CM] (hybrid isotropic-anisotropic phase field model of Ambati et al. compared with the isotropic one of Bourdin et al. and the anisotropic ones of Amor et al. and Miehe et al.)
% [Liu, Li, Msekh, Zuo, 2016, CMS] (anisotropic phase field model of Miehe et al.)
% [Wu, Nguyen, 2018, JMPS] (hybrid isotropic-anisotropic phase field model of Wu et al.)
% [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM] (anisotropic phase field model of Wu et al.)
% [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US] (anisotropic phase field model of Amor et al.)
% [Wu, Nguyen, Zhou, Huang, 2020, CMAME] (anisotropic phase field model of Wu et al.)
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME] (anisotropic phase field model of Nguyen et al.)
% [Hu, Guilleminot, Dolbow, 2020, CMAME] (anisotropic phase field model of Hu et al.)

% clc
clearvars
close all
% rng('default');

%% Input data
setProblem = true;
solveProblem = true;
displayModel = false;
displaySolution = false;

% test = true; % coarse mesh
test = false; % fine mesh

numWorkers = maxNumCompThreads;
% numWorkers = 1; maxNumCompThreads(1); % mono-thread computation

% Deterministic model parameters
Dim = 2; % space dimension Dim = 2, 3
symmetry = 'Isot'; % 'Isot', 'MeanIsot', 'Anisot'. Material symmetry
ang = 45; % clockwise material orientation angle around z-axis for anisotopic material [deg]
loading = 'Shear'; % 'Tension' or 'Shear'
PFmodel = 'Miehe'; % 'Bourdin', 'Amor', 'Miehe', 'He', 'Zhang', 'Spectral'
PFsplit = 'Strain'; % 'Strain' or 'Stress'
PFregularization = 'AT2'; % 'AT1' or 'AT2'
PFsolver = 'BoundConstrainedOptim'; % 'HistoryFieldElem', 'HistoryFieldNode' or 'BoundConstrainedOptim'
maxIter = 1; % maximum number of iterations at each loading increment
tolConv = 1e-2; % prescribed tolerance for convergence at each loading increment
initialCrack = 'GeometricCrack'; % 'GeometricCrack', 'GeometricNotch', 'InitialPhaseField'
FEmesh = 'Optim'; % 'Unif' or 'Optim'

% Random model parameters
randMat = struct('delta',0.2,'lcorr',1e-4); % random material parameters model
gc = 2.7e3;
aGc = 0.6*gc;
bGc = 1.4*gc;
% aGc = 0.9*gc;
% bGc = 1.1*gc;
% aGc = [0.7,1.2]*gc;
% bGc = [0.8,1.3]*gc;
randPF = struct('aGc',aGc,'bGc',bGc,'lcorr',Inf); % random phase field parameters model

numsamples = 1e4; % total number of samples
% numsamples = 500;
% sampleindices = 1:500;
switch lower(loading)
    case 'tension'
        sampleindices = 1:1120;
        % sampleindices = (1120+1):2240;
        % sampleindices = (2240+1):3360;
        % sampleindices = (3360+1):4320;
        % sampleindices = (4320+1):5280;
        % sampleindices = (5280+1):6240;
        % sampleindices = (6240+1):7200;
        % sampleindices = (7200+1):8160;
        % sampleindices = (8160+1):9120;
        % sampleindices = (9120+1):9600;
        % sampleindices = (9600+1):1e4;
    case 'shear'
        sampleindices = 1:480;
        % sampleindices = (480+1):960;
        % sampleindices = (960+1):1440;
        % sampleindices = (1440+1):1920;
        % sampleindices = (1920+1):2400;
        % sampleindices = (2400+1):2880;
        % sampleindices = (2880+1):3360;
        % sampleindices = (3360+1):3840;
        % sampleindices = (3840+1):4320;
        % sampleindices = (4320+1):4800;
        % sampleindices = (4800+1):5280;
        % sampleindices = (5280+1):5760;
        % sampleindices = (5760+1):6240;
        % sampleindices = (6240+1):6720;
        % sampleindices = (6720+1):7200;
        % sampleindices = (7200+1):7680;
        % sampleindices = (7680+1):8160;
        % sampleindices = (8160+1):8640;
        % sampleindices = (8640+1):9120;
        % sampleindices = (9120+1):9600;
        % sampleindices = (9600+1):1e4;
%         sampleindices = [9,13,55,61,101,109,127,198,226,282,304,339,345,353,365,378,408,434,463,479,492,...
%             512,548,588,636,652,653,660,682,695,710,729,741,765,767,792,832,845,864,879,891,966,968,975,992,...
%             1002,1008,1019,1102,1107,1115,1135,1152,1178,1187,1214,1217,1232,1363,1366,1391,1405,1409,1438,1473,1485,...
%             1532,1616,1620,1626,1639,1666,1702,1791,1795,1822,1827,1829,1858,1864,1870,1871,1917,1919,1938,...
%             2041,2049,2142,2151,2157,2200,2206,2237,2276,2332,2362,2365,2372,2421,2444,2456,2457,2465,2467,2476,2483,2493,...
%             2513,2519,2585,2595,2596,2627,2645,2657,2669,2678,2688,2697,2717,2737,2747,2752,2763,2824,2880,2892,2949,...
%             3017,3075,3085,3094,3130,3149,3158,3223,3268,3305,3313,3337,3354,3378,3383,3402,3444,...
%             3505,3508,3536,3538,3540,3570,3594,3598,3623,3655,3658,3674,3678,3721,3802,3854,3858,3859,3911,3937,3964,...
%             4014,4049,4138,4149,4165,4194,4211,4221,4281,4299,4311,4329,4333,4387,4391,4399,4406,4424,4453,4482,...
%             4500,4514,4518,4532,4627,4655,4678,4724,4741,4752,4783,4795,4853,4855,4884,4895,4995];
%         sampleindices = [5036,5052,5089,5104,5185,5197,5211,5242,5248,5298,5300,5304,5331,5357,5407,5410,5454,5466,5469,5496,5499,...
%             5507,5543,5551,5564,5609,5621,5640,5648,5657,5669,5688,5691,5718,5739,5757,5815,5817,5864,5874,5875,5900,5912,5931,5959,5992,...
%             6000,6018,6035,6072,6091,6098,6102,6140,6147,6148,6214,6215,6243,6260,6275,6315,6355,6386,6402,6413,6420,6434,6443,6473,6496,...
%             6511,6513,6523,6545,6551,6619,6620,6621,6623,6640,6654,6678,6682,6692,6714,6750,6771,6775,6792,6797,6819,6824,6851,6852,6861,6871,6943,6946,6947,6954,6992,...
%             7000,7042,7045,7068,7084,7090,7097,7111,7115,7142,7158,7170,7177,7292,7349,7372,7391,7392,7397,7399,7415,7419,7429,7435,7453,7499,...
%             7525,7578,7614,7626,7628,7637,7640,7655,7656,7663,7664,7683,7743,7747,7755,7761,7778,7902,7918,7943,7964,7993,...
%             8021,8059,8096,8103,8156,8159,8177,8182,8210,8258,8266,8283,8366,8397,8412,8416,...
%             8581,8589,8595,8623,8698,8747,8788,8789,8808,8824,8827,8861,8865,8866,8881,8895,8903,8926,8939,8942,8974,8995,...
%             9013,9030,9032,9090,9130,9144,9160,9200,9216,9220,9228,9277,9326,9346,9372,9384,9397,9405,9414,9446,9468,9483,9494,...
%             9504,9536,9541,9550,9560,9569,9608,9613,9665,9681,9711,9723,9724,9737,9777,9789,9806,9809,9874,9888,9895,9901,9929,9945,9950,9956];
%         sampleindices = [741,1822,3354,3854,3858,4724,5691,6035,6260,6619,8103,9013,9130];
        
%         sampleindices = 1:500;
%         sampleindices = [8,9,29,61,64,96,101,109,150,187,252,255,273,277,285,304,327,338,339,345,415,434,492,494];
%         sampleindices = [61,109,252,273,277,285,338,339,434];

%         sampleindices = [8,33,34,35,43,58,59,77,111,127,157,161,175,196,214,282,298,317,323,343,404,430,432,434];
%         sampleindices = [8,33,58,59,77,127,161,298,323,343,404];
%         sampleindices = [59,77,127,161,298,323,343,404];
    otherwise
        error('Wrong loading case');
end
N = length(sampleindices);

foldername = ['singleEdgeCrack' loading '_' num2str(Dim) 'D'];
filename = ['linElas' symmetry];
if strcmpi(symmetry,'anisot') % anisotropic material
    filename = [filename num2str(ang) 'deg'];
end
filename = [filename PFmodel PFsplit PFregularization PFsolver initialCrack...% 'MaxIter' num2str(maxIter) 'Tol' num2str(tolConv)
    'Mesh' FEmesh '_' num2str(numsamples) 'samples_from_' num2str(sampleindices(1)) '_to_' num2str(sampleindices(end))];
if any(randPF.aGc) && any(randPF.bGc)
    gcbounds = [randPF.aGc(:),randPF.bGc(:)]';
    filename = [filename '_RandPF_Gc' num2str(gcbounds(:)','_%g') '_Lcorr' num2str(randPF.lcorr,'_%g')];
end

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

%% Problem
if setProblem
    %% Domains and meshes
    L = 1e-3;
    a = L/2;
    if Dim==2
        e = 1;
        D = DOMAIN(2,[0.0,0.0],[L,L]);
        C = LIGNE([0.0,L/2],[a,L/2]);
    elseif Dim==3
        e = 0.1e-3;
        D = DOMAIN(3,[0.0,0.0,0.0],[L,L,e]);
        C = QUADRANGLE([0.0,L/2,0.0],[a,L/2,0.0],[a,L/2,e],[0.0,L/2,e]);
    end
    
    if Dim==2
        % clD = 6.25e-5; % [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME]
        % clD = 5e-5; % [Hu, Guilleminot, Dolbow, 2020, CMAME]
        % clD = 3e-5; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        % clD = 2e-5; % [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
        % clD = 3.9e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
        % clD = 2e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        % clD = 1e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        
        % clC = 3.906e-6; % [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME]
        % clC = 5e-6; % [Hu, Guilleminot, Dolbow, 2020, CMAME]
        % clC = 2.5e-6; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        % clC = 2e-6; % (shear test) [Miehe, Hofacker, Welschinger, 2010, CMAME]
        % clC = 1e-6; % (tension test) [Miehe, Welschinger, Hofacker, 2010 IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME]
        % clC = 6e-7; % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
        % clC = 3.9e-6; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
        % clC = 2e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        % clC = 1e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        switch lower(FEmesh)
            case 'unif'
                switch lower(symmetry)
                    case 'isot' % isotropic material
                        cl = 5e-6;
                    case 'anisot' % anisotropic material
                        cl = 4.25e-6;
                    otherwise
                        error('Wrong material symmetry class');
                end
                if test
                    cl = 1e-5;
                end
                clD = cl;
                clC = cl;
            case 'optim'
                clD = 2.5e-5;
                clC = 2.5e-6;
                if test
                    clD = 4e-5;
                    clC = 1e-5;
                end
            otherwise
                error('Wrong FE mesh')
        end
    elseif Dim==3
        switch lower(FEmesh)
            case 'unif'
                switch lower(symmetry)
                    case 'isot' % isotropic material
                        cl = 5e-6;
                    case 'anisot' % anisotropic material
                        cl = 4.25e-6;
                    otherwise
                        error('Wrong material symmetry class');
                end
                if test
                    cl = 2e-5;
                end
                clD = cl;
                clC = cl;
            case 'optim'
                clD = 4e-5;
                clC = 5e-6;
                if test
                    clD = 4e-5;
                    clC = 1e-5;
                end
            otherwise
                error('Wrong FE mesh')
        end
    end
    switch lower(initialCrack)
        case 'geometriccrack'
            S_phase = gmshdomainwithedgecrack(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,'duplicate',lower(loading),lower(symmetry),lower(PFmodel));
        case 'geometricnotch'
            c = 1e-5; % crack width
            S_phase = gmshdomainwithedgesmearedcrack(D,C,c,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,lower(loading),lower(symmetry),lower(PFmodel));
        case 'initialphasefield'
            S_phase = gmshdomainwithedgecrack(D,C,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'),Dim,lower(initialCrack),lower(loading),lower(symmetry),lower(PFmodel));
        otherwise
            error('Wrong model for initial crack');
    end
    S = S_phase;
    
    %% Phase field problem
    %% Material
    switch lower(symmetry)
        case {'isot','meanisot'} % almost surely or mean isotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 2.7e3; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
            % Regularization parameter (width of the smeared crack)
            % l = 3.75e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME]
            % l = 3e-5; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
            % l = 1.5e-5; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM], [Hu, Guilleminot, Dolbow, 2020, CMAME]
            l = 1e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
            % l = 7.5e-6; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            % l = 5e-6; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
            % l = 4e-6; % [Ambati, Gerasimov, De Lorenzis, 2015, CM]
            % eta = 0.052; w0 = 75.94; l = eta/sqrt(w0)*1e-3; % l = 6e-7; % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
        case 'anisot' % anisotropic material
            % Critical energy release rate (or fracture toughness)
            gc = 1e3; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
            % Regularization parameter (width of the smeared crack)
            l = 8.5e-6; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        otherwise
            error('Wrong material symmetry class');
    end
    % Small artificial residual stiffness
    % k = 1e-12;
    k = 0;
    
    % Material
    switch lower(PFregularization)
        case 'at1'
            % c0 = 8/3;
            K = 3/4*gc*l; % K = 2*(gc*l)/c0;
            R = 0;
            Qn = -3/8*gc/l; % Qn = -(gc/l)/c0;
        case 'at2'
            % c0 = 2;
            K = gc*l; % K = 2*(gc*l)/c0;
            R = gc/l; % R = 2*(gc/l)/c0;
            Qn = 0;
        otherwise
            error('Wrong regularization model');
    end
    mat_phase = FOUR_ISOT('k',K,'r',R,'qn',Qn,'PFregularization',PFregularization,'aGc',randPF.aGc,'bGc',randPF.bGc,'lcorr',randPF.lcorr);
    mat_phase = setnumber(mat_phase,1);
    S_phase = setmaterial(S_phase,mat_phase);
    
    %% Dirichlet boundary conditions
    switch lower(initialCrack)
        case 'geometriccrack'
            S_phase = final(S_phase,'duplicate');
        case 'geometricnotch'
            S_phase = final(S_phase);
        case 'initialphasefield'
            S_phase = final(S_phase);
            S_phase = addcl(S_phase,C,'T',1);
        otherwise
            error('Wrong model for initial crack');
    end
    
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
    option = 'DEFO'; % plane strain [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
    % option = 'CONT'; % plane stress [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
    switch lower(symmetry)
        case {'isot','meanisot'} % almost surely or mean isotropic material
            % Lame coefficients
            % lambda = 121.1538e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
            % mu = 80.7692e9; % [Miehe, Welschinger, Hofacker, 2010 IJNME]
            lambda = 121.15e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
            mu = 80.77e9; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
            if strcmpi(symmetry,'isot')
                % Young modulus and Poisson ratio
                if Dim==2
                    switch lower(option)
                        case 'defo'
                            E = mu*(3*lambda+2*mu)/(lambda+mu); % E = 210e9;
                            NU = lambda/(lambda+mu)/2; % NU = 0.3;
                        case 'cont'
                            E = 4*mu*(lambda+mu)/(lambda+2*mu);
                            NU = lambda/(lambda+2*mu);
                    end
                    % E = 210e9; NU = 0.2; % [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM]
                    % E = 210e9; NU = 0.3; % [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
                    % kappa = 121030e6; NU=0.227; lambda=3*kappa*NU/(1+NU); mu = 3*kappa*(1-2*NU)/(2*(1+NU)); E = 3*kappa*(1-2*NU); % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
                elseif Dim==3
                    E = mu*(3*lambda+2*mu)/(lambda+mu);
                    NU = lambda/(lambda+mu)/2;
                end
            elseif strcmpi(symmetry,'meanisot')
                % Elasticity matrix
                if Dim==2
                    Cmat = e*...
                        [lambda+2*mu,lambda,0;...
                        lambda,lambda+2*mu,0;...
                        0,0,mu];
                elseif Dim==3
                    Cmat = [lambda+2*mu,lambda,lambda,0,0,0;...
                        lambda,lambda+2*mu,lambda,0,0,0;...
                        lambda,lambda,lambda+2*mu,0,0,0;...
                        0,0,0,mu,0,0;...
                        0,0,0,0,mu,0;...
                        0,0,0,0,0,mu];
                end
            end
            
        case 'anisot' % anisotropic material
            if Dim==2
                switch lower(option)
                    case 'defo'
                        % [Nguyen, Yvonnet, Waldmann, He, 2020,IJNME]
                        % Elasticity matrix in reference material coordinate system [Pa]
                        Cmat = e*...
                            [65 20 0;
                            20 260 0;
                            0 0 30]*1e9;
                        theta = deg2rad(ang); % clockwise material orientation angle around z-axis [rad]
                        c = cos(theta);
                        s = sin(theta);
                        % Transition matrix for elasticity matrix from material coordinate system to global coordinate system
                        P = [c^2 s^2 -c*s;
                            s^2 c^2 c*s;
                            2*c*s -2*c*s c^2-s^2];
                        % Elasticity matrix in global coordinate system [Pa]
                        Cmat = P'*Cmat*P;
                    case 'cont'
                        error('Not implemented yet');
                end
            elseif Dim==3
                error('Not implemented yet');
            end
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
        case 'isot' % almost surely isotropic material
            mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit,'delta',randMat.delta,'lcorr',randMat.lcorr);
        case {'meanisot','anisot'} % mean isotropic or anisotropic material
            mat = ELAS_ANISOT('C',Cmat,'RHO',RHO,'DIM3',e,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel,'PFS',PFsplit,'delta',randMat.delta,'lcorr',randMat.lcorr);
        otherwise
            error('Wrong material symmetry class');
    end
    mat = setnumber(mat,1);
    S = setoption(S,option);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    if Dim==2
        BU = LIGNE([0.0,L],[L,L]);
        BL = LIGNE([0.0,0.0],[L,0.0]);
        BRight = LIGNE([L,0.0],[L,L]);
        BLeft = LIGNE([0.0,0.0],[0.0,L]);
        BFront = [];
        BBack = [];
    elseif Dim==3
        BU = PLAN([0.0,L,0.0],[L,L,0.0],[0.0,L,e]);
        BL = PLAN([0.0,0.0,0.0],[L,0.0,0.0],[0.0,0.0,e]);
        BRight = PLAN([L,0.0,0.0],[L,L,0.0],[L,0.0,e]);
        BLeft = PLAN([0.0,0.0,0.0],[0.0,L,0.0],[0.0,0.0,e]);
        BFront = PLAN([0.0,0.0,e],[L,0.0,e],[0.0,L,e]);
        BBack = PLAN([0.0,0.0,0.0],[L,0.0,0.0],[0.0,L,0.0]);
    end
    
    switch lower(initialCrack)
        case 'geometriccrack'
            S = final(S,'duplicate');
        otherwise
            S = final(S);
    end
    
    ud = 0;
    switch lower(loading)
        case 'tension'
            if Dim==2
                S = addcl(S,BU,{'UX','UY'},[0;ud]);
            elseif Dim==3
                S = addcl(S,BU,{'UX','UY','UZ'},[0;ud;0]);
            end
            S = addcl(S,BL,'UY');
        case 'shear'
            if Dim==2
                S = addcl(S,BU,{'UX','UY'},[ud;0]);
                S = addcl(S,BLeft,'UY');
                S = addcl(S,BRight,'UY');
            elseif Dim==3
                S = addcl(S,BU,{'UX','UY','UZ'},[ud;0;0]);
                S = addcl(S,BLeft,{'UY','UZ'});
                S = addcl(S,BRight,{'UY','UZ'});
                S = addcl(S,BFront,{'UY','UZ'});
                S = addcl(S,BBack,{'UY','UZ'});
            end
            S = addcl(S,BL);
        otherwise
            error('Wrong loading case');
    end
    
    %% Stiffness matrices and sollicitation vectors
    % [A,b] = calc_rigi(S);
    % b = -b;
    
    %% Time scheme
    switch lower(symmetry)
        case {'isot','meanisot'} % almost surely or mean isotropic material
            if Dim==2
                switch lower(loading)
                    case 'tension'
                        % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM]
                        % du = 1e-5 mm during the first 500 time steps (up to u = 5e-3 mm)
                        % du = 1e-6 mm during the last 1300 time steps (up to u = 6.3e-3 mm)
                        % dt0 = 1e-8;
                        % nt0 = 500;
                        % dt1 = 1e-9;
                        % nt1 = 1300;
                        % t0 = linspace(dt0,nt0*dt0,nt0);
                        % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        % t = [t0,t1];
                        
                        % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
                        % du = 1e-5 mm during 630 time steps (up to u = 6.3e-3 mm)
                        % dt = 1e-8;
                        % nt = 630;
                        % t = linspace(dt,nt*dt,nt);
                        
                        % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
                        % du = 1e-4 mm during 63 time steps (up to u = 6.3e-3 mm)
                        % dt = 1e-7;
                        % nt = 63;
                        % t = linspace(dt,nt*dt,nt);
                        
                        % [Liu, Li, Msekh, Zuo, 2016, CMS]
                        % du = 1e-4 mm during the first 50 time steps (up to u = 5e-3 mm)
                        % du = 1e-6 mm during the last 1300 time steps (up to u = 6.3e-3 mm)
                        % dt0 = 1e-7;
                        % nt0 = 50;
                        % dt1 = 1e-9;
                        % nt1 = 1300;
                        % t0 = linspace(dt0,nt0*dt0,nt0);
                        % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        % t = [t0,t1];
                        
                        % du = 1e-5 mm during the first 400 time steps (up to u = 4e-3 mm)
                        % du = 1e-6 mm during the last 4000 time steps (up to u = 8e-3 mm)
                        dt0 = 1e-8;
                        nt0 = 400;
                        dt1 = 1e-9;
                        nt1 = 4000;
                        if test
                            dt0 = 1e-7;
                            nt0 = 40;
                            dt1 = 1e-8;
                            nt1 = 400;
                        end
                        t0 = linspace(dt0,nt0*dt0,nt0);
                        t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        t = [t0,t1];
                    case 'shear'
                        % [Miehe, Welschinger, Hofacker, 2010 IJNME]
                        % du = 1e-4 mm during the first 100 time steps (up to u = 10e-3 mm)
                        % du = 1e-6 mm during the last 10 000 time steps (up to u = 20e-3 mm)
                        % dt0 = 1e-7;
                        % nt0 = 100;
                        % dt1 = 1e-9;
                        % nt1 = 10000;
                        % t0 = linspace(dt0,nt0*dt0,nt0);
                        % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        % t = [t0,t1];
                        
                        % [Liu, Li, Msekh, Zuo, 2016, CMS]
                        % du = 1e-4 mm during the first 50 time steps (up to u = 5e-3 mm)
                        % du = 1e-5 mm during the last 1500 time steps (up to u = 20e-3 mm)
                        % dt0 = 1e-7;
                        % nt0 = 50;
                        % dt1 = 1e-8;
                        % nt1 = 1500;
                        % t0 = linspace(dt0,nt0*dt0,nt0);
                        % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        % t = [t0,t1];
                        
                        % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
                        % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2019, AAM],
                        % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
                        % du = 1e-5 mm during 1500 time steps (up to u = 15e-3 mm)
                        % dt = 1e-8;
                        % nt = 1500;
                        
                        % du = 1e-5 mm during 2000 time steps (up to u = 20e-3 mm)
                        dt = 1e-8;
                        nt = 2000;
                        if test
                            dt = 5e-8;
                            nt = 400;
                        end
                        t = linspace(dt,nt*dt,nt);
                end
            elseif Dim==3
                % du = 1e-5 mm during 2500 time steps (up to u = 25e-3 mm)
                dt = 1e-8;
                nt = 2500;
                if test
                    dt = 1e-7;
                    nt = 250;
                end
                t = linspace(dt,nt*dt,nt);
            end
            T = TIMEMODEL(t);
            
        case 'anisot' % anisotropic material
            if Dim==2
                switch lower(loading)
                    case 'tension'
                        % du = 2e-5 mm (up to u = 10e-3 mm)
                        dt = 2e-8;
                        nt = 500;
                        if test
                            dt1 = 4e-8;
                            nt = 250;
                        end
                    case 'shear'
                        % du = 2e-5 mm (up to u = 20e-3 mm)
                        dt = 2e-8;
                        nt = 1000;
                        if test
                            dt = 4e-8;
                            nt = 500;
                        end
                end
                t = linspace(dt,nt*dt,nt);
            elseif Dim==3
                % du = 2e-5 mm (up to u = 20e-3 mm)
                dt = 2e-8;
                nt = 1000;
                if test
                    dt = 2e-7;
                    nt = 100;
                end
                t = linspace(dt,nt*dt,nt);
            end
            T = TIMEMODEL(t);
            
        otherwise
            error('Wrong material symmetry class');
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','BU','BL','BRight','BLeft','BFront','BBack','loading','symmetry','ang');
else
    load(fullfile(pathname,'problem.mat'),'T','S_phase','S','D','C','BU','BL','BRight','BLeft','BFront','BBack','loading','symmetry','ang');
end

%% Solution
if solveProblem
    myparallel('start',numWorkers);
    
    %% Solution
    tTotal = tic;
    
    fun = @(S_phase,S) solvePFDetLinElasSingleEdgeCrackForce(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,'maxiter',maxIter,'tol',tolConv);
    [ft,dmaxt,gc_sample] = solvePFStoLinElasForceGc(S_phase,S,T,fun,N,'numsamples',numsamples,'sampleindices',sampleindices);
    t = gettevol(T);
    idc = arrayfun(@(i) find(dmaxt(i,:)>=min(0.75,max(dmaxt(i,:))),1),1:N)';
    fc = arrayfun(@(i) ft(i,idc(i)),1:N)';
    udc = t(idc);
    [fmax,idmax] = max(ft,[],2);
    udmax = t(idmax);
    
    time = toc(tTotal);
    
    myparallel('stop');
    
    %% Statistical outputs of solution
    probs = [0.025 0.975];
    
    ft_mean = mean(ft);
    ft_std = std(ft);
    ft_ci = quantile(ft,probs);
    
    fmax_mean = mean(fmax);
    fmax_std = std(fmax);
    fmax_ci = quantile(fmax,probs);
    
    udmax_mean = mean(udmax);
    udmax_std = std(udmax);
    udmax_ci = quantile(udmax,probs);
    
    fc_mean = mean(fc);
    fc_std = std(fc);
    fc_ci = quantile(fc,probs);
    
    udc_mean = mean(udc);
    udc_std = std(udc);
    udc_ci = quantile(udc,probs);
    
    gc_mean = mean(gc_sample);
    gc_std = std(gc_sample);
    gc_ci = quantile(gc_sample,probs);
    
    npts = 100;
    [fmax_f,fmax_xi,fmax_bw] = ksdensity(fmax,'npoints',npts);
    [udmax_f,udmax_xi,udmax_bw] = ksdensity(udmax,'npoints',npts);
    [fc_f,fc_xi,fc_bw] = ksdensity(fc,'npoints',npts);
    [udc_f,udc_xi,udc_bw] = ksdensity(udc,'npoints',npts);
    [gc_f,gc_xi,gc_bw] = ksdensity(gc_sample,'npoints',npts);
    
    save(fullfile(pathname,'solution.mat'),'N','ft','dmaxt',...
        'ft_mean','ft_std','ft_ci','probs','time',...
        'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
        'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw',...
        'fc','fc_mean','fc_std','fc_ci','fc_f','fc_xi','fc_bw',...
        'udc','udc_mean','udc_std','udc_ci','udc_f','udc_xi','udc_bw',...
        'gc_sample','gc_mean','gc_std','gc_ci','gc_f','gc_xi','gc_bw');
else
    load(fullfile(pathname,'solution.mat'),'N','ft','dmaxt',...
        'ft_mean','ft_std','ft_ci','probs','time',...
        'fmax','fmax_mean','fmax_std','fmax_ci','fmax_f','fmax_xi','fmax_bw',...
        'udmax','udmax_mean','udmax_std','udmax_ci','udmax_f','udmax_xi','udmax_bw',...
        'fc','fc_mean','fc_std','fc_ci','fc_f','fc_xi','fc_bw',...
        'udc','udc_mean','udc_std','udc_ci','udc_f','udc_xi','udc_bw',...
        'gc_sample','gc_mean','gc_std','gc_ci','gc_f','gc_xi','gc_bw');
end

%% Outputs
fid = fopen(fullfile(pathname,'results.txt'),'w');
fprintf(fid,'Single edge crack\n');
fprintf(fid,'\n');
fprintf(fid,'dim      = %d\n',Dim);
fprintf(fid,'loading  = %s\n',loading);
fprintf(fid,'mat sym  = %s\n',symmetry);
if strcmpi(symmetry,'anisot')
    fprintf(fid,'angle    = %g deg\n',ang);
end
fprintf(fid,'PF model = %s\n',PFmodel);
fprintf(fid,'PF split = %s\n',PFsplit);
fprintf(fid,'PF regularization = %s\n',PFregularization);
fprintf(fid,'PF solver = %s\n',PFsolver);
fprintf(fid,'nb elements = %g\n',getnbelem(S));
fprintf(fid,'nb nodes    = %g\n',getnbnode(S));
fprintf(fid,'nb dofs     = %g\n',getnbddl(S));
fprintf(fid,'nb time dofs = %g\n',getnbtimedof(T));
fprintf(fid,'nb samples = %g\n',N);
fprintf(fid,'elapsed time = %f s\n',time);
fprintf(fid,'\n');

if Dim==2
    fprintf(fid,'mean(fmax)   = %g kN/mm\n',fmax_mean*1e-6);
    fprintf(fid,'std(fmax)    = %g kN/mm\n',fmax_std*1e-6);
elseif Dim==3
    fprintf(fid,'mean(fmax)   = %g kN\n',fmax_mean*1e-3);
    fprintf(fid,'std(fmax)    = %g kN\n',fmax_std*1e-3);
end
fprintf(fid,'disp(fmax)   = %g\n',fmax_std/fmax_mean);
if Dim==2
    fprintf(fid,'%d%% ci(fmax) = [%g,%g] kN/mm\n',(probs(2)-probs(1))*100,fmax_ci(1)*1e-6,fmax_ci(2)*1e-6);
elseif Dim==3
    fprintf(fid,'%d%% ci(fmax) = [%g,%g] kN\n',(probs(2)-probs(1))*100,fmax_ci(1)*1e-3,fmax_ci(2)*1e-3);
end
fprintf(fid,'\n');

if Dim==2
    fprintf(fid,'mean(fc)   = %g kN/mm\n',fc_mean*1e-6);
    fprintf(fid,'std(fc)    = %g kN/mm\n',fc_std*1e-6);
elseif Dim==3
    fprintf(fid,'mean(fc)   = %g kN\n',fc_mean*1e-3);
    fprintf(fid,'std(fc)    = %g kN\n',fc_std*1e-3);
end
fprintf(fid,'disp(fc)   = %g\n',fc_std/fc_mean);
if Dim==2
    fprintf(fid,'%d%% ci(fc) = [%g,%g] kN/mm\n',(probs(2)-probs(1))*100,fc_ci(1)*1e-6,fc_ci(2)*1e-6);
elseif Dim==3
    fprintf(fid,'%d%% ci(fc) = [%g,%g] kN\n',(probs(2)-probs(1))*100,fc_ci(1)*1e-3,fc_ci(2)*1e-3);
end
fprintf(fid,'\n');

fprintf(fid,'mean(udmax)   = %g mm\n',udmax_mean*1e3);
fprintf(fid,'std(udmax)    = %g mm\n',udmax_std*1e3);
fprintf(fid,'disp(udmax)   = %g\n',udmax_std/udmax_mean);
fprintf(fid,'%d%% ci(udmax) = [%g,%g] mm\n',(probs(2)-probs(1))*100,udmax_ci(1)*1e3,udmax_ci(2)*1e3);
fprintf(fid,'\n');

fprintf(fid,'mean(udc)   = %g mm\n',udc_mean*1e3);
fprintf(fid,'std(udc)    = %g mm\n',udc_std*1e3);
fprintf(fid,'disp(udc)   = %g\n',udc_std/udc_mean);
fprintf(fid,'%d%% ci(udc) = [%g,%g] mm\n',(probs(2)-probs(1))*100,udc_ci(1)*1e3,udc_ci(2)*1e3);
fprintf(fid,'\n');

gc_xiunif = linspace(min(gc_xi),max(gc_xi),1e3);
gc_funif = myunifpdf(gc_xiunif,aGc,bGc);
gc_ciunif = myunifinv(probs,aGc,bGc);

[gc_meanunif,gc_varunif] = myunifstat(aGc,bGc);
gc_stdunif = sqrt(gc_varunif);
fprintf(fid,'mean(gc)   = %g N/mm (estimate)\n',gc_mean*1e-3);
fprintf(fid,'           = %g N/mm (exact)\n',gc_meanunif*1e-3);
fprintf(fid,'std(gc)    = %g N/mm (estimate)\n',gc_std*1e-3);
fprintf(fid,'           = %g N/mm (exact)\n',gc_stdunif*1e-3);
fprintf(fid,'disp(gc)   = %g (estimate)\n',gc_std/gc_mean);
fprintf(fid,'           = %g (exact)\n',gc_stdunif/gc_meanunif);
fprintf(fid,'%d%% ci(gc) = [%g,%g] N/mm (estimate)\n',(probs(2)-probs(1))*100,gc_ci(1)*1e-3,gc_ci(2)*1e-3);
fprintf(fid,'           = [%g,%g] N/mm (exact)\n',gc_ciunif(1)*1e-3,gc_ciunif(2)*1e-3);
fclose(fid);

%% Display
if displayModel
    [t,rep] = gettevol(T);
    
    %% Display domains, boundary conditions and meshes
%     plotDomain({D,C},'legend',false);
%     mysaveas(pathname,'domain',formats,renderer);
%     mymatlab2tikz(pathname,'domain.tex');
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    ampl = 0.5;
    v = calc_init_dirichlet(S);
    [hN,legN] = vectorplot(S,'U',v,ampl,'r','LineWidth',1);
    % legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions_displacement',formats,renderer);
    
    [hD_phase,legD_phase] = plotBoundaryConditions(S_phase,'legend',false);
    % legend([hD_phase,hN_phase],[legD_phase,legN_phase],'Location','NorthEastOutside')
    mysaveas(pathname,'boundary_conditions_damage',formats,renderer);
    
    % plotModel(S,'legend',false);
    % mysaveas(pathname,'mesh',formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,'mesh',formats,renderer);
end

%% Display statistics of solutions
if displaySolution
    [t,~] = gettevol(T);
    
    %% Display force-displacement curve
    figure('Name','Force vs displacement')
    clf
    plot(t*1e3,ft_mean*((Dim==2)*1e-6+(Dim==3)*1e-3),'-b','Linewidth',linewidth)
    hold on
    ciplot(ft_ci(1,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),ft_ci(2,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),t*1e3,'b');
    alpha(0.2)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    l = legend('mean function',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'Location','NorthWest');
    set(l,'Interpreter','latex')
    mysaveas(pathname,'force_displacement',formats);
    mymatlab2tikz(pathname,'force_displacement.tex');
    
    figure('Name','Forces vs displacement')
    clf
    color = distinguishable_colors(N);
    for i=1:N
        plot(t*1e3,ft(i,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),'LineStyle','-','Color',color(i,:),'Linewidth',linewidth)
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Displacement [mm]','Interpreter',interpreter)
    ylabel('Force [kN]','Interpreter',interpreter)
    mysaveas(pathname,'forces_displacement',formats);
    mymatlab2tikz(pathname,'forces_displacement.tex');

%     color = distinguishable_colors(N);
%     for i=1:N
%         if mod(i,100)==0
%             close all
%         end
%         figure('Name',['Force-displacement #' num2str(sampleindices(i))])
%         clf
%         plot(t*1e3,ft(i,:)*((Dim==2)*1e-6+(Dim==3)*1e-3),'LineStyle','-','Color',color(i,:),'Linewidth',linewidth)
%         hold on
%         scatter(udc(i)*1e3,fc(i)*((Dim==2)*1e-6+(Dim==3)*1e-3),'Marker','+','MarkerEdgeColor','r','Linewidth',linewidth)
%         scatter(udmax(i)*1e3,fmax(i)*((Dim==2)*1e-6+(Dim==3)*1e-3),'Marker','+','MarkerEdgeColor',color(i,:),'Linewidth',linewidth)
%         hold off
%         grid on
%         box on
%         set(gca,'FontSize',fontsize)
%         xlabel('Displacement [mm]','Interpreter',interpreter)
%         ylabel('Force [kN]','Interpreter',interpreter)
%         % mysaveas(pathname,['force_displacement_' num2str(sampleindices(i))],{'epsc','png'});
%         % mymatlab2tikz(pathname,['force_displacement_' num2str(sampleindices(i)) '.tex']);
%     end
    
    %% Display pdf of maximum force
    figure('Name','Probability Density Estimate: Maximum force')
    clf
    plot(fmax_xi*((Dim==2)*1e-6+(Dim==3)*1e-3),fmax_f,'-b','LineWidth',linewidth)
    hold on
    ind_fmax = find(fmax_xi>=fmax_ci(1) & fmax_xi<fmax_ci(2));
    area(fmax_xi(ind_fmax)*((Dim==2)*1e-6+(Dim==3)*1e-3),fmax_f(ind_fmax),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(fmax_mean*((Dim==2)*1e-6+(Dim==3)*1e-3),0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$f$ [kN]','Interpreter',interpreter)
    ylabel('$p_{F_{\mathrm{max}}}(f)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_fmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_fmax.tex');
    
    %% Display pdf of critical force
    figure('Name','Probability Density Estimate: Critical force')
    clf
    plot(fc_xi*((Dim==2)*1e-6+(Dim==3)*1e-3),fc_f,'-b','LineWidth',linewidth)
    hold on
    ind_fc = find(fc_xi>=fc_ci(1) & fc_xi<fc_ci(2));
    area(fc_xi(ind_fc)*((Dim==2)*1e-6+(Dim==3)*1e-3),fc_f(ind_fc),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(fc_mean*((Dim==2)*1e-6+(Dim==3)*1e-3),0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$f$ [kN]','Interpreter',interpreter)
    ylabel('$p_{F_c}(f)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_fc',formats,renderer);
    mymatlab2tikz(pathname,'pdf_fc.tex');
    
    %% Display pdf of maximum displacement
    figure('Name','Probability Density Estimate: Maximum displacement')
    clf
    plot(udmax_xi*1e3,udmax_f,'-b','LineWidth',linewidth)
    hold on
    ind_udmax = find(udmax_xi>=udmax_ci(1) & udmax_xi<udmax_ci(2));
    area(udmax_xi(ind_udmax)*1e3,udmax_f(ind_udmax),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(udmax_mean*1e3,0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$u$ [mm]','Interpreter',interpreter)
    ylabel('$p_{U_{D,\mathrm{max}}}(u)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_udmax',formats,renderer);
    mymatlab2tikz(pathname,'pdf_udmax.tex');
    
    %% Display pdf of critical displacement
    figure('Name','Probability Density Estimate: Critical displacement')
    clf
    plot(udc_xi*1e3,udc_f,'-b','LineWidth',linewidth)
    hold on
    ind_udc = find(udc_xi>=udc_ci(1) & udc_xi<udc_ci(2));
    area(udc_xi(ind_udc)*1e3,udc_f(ind_udc),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(udc_mean*1e3,0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$u$ [mm]','Interpreter',interpreter)
    ylabel('$p_{U_{D,c}}(u)$','Interpreter',interpreter)
    l = legend('pdf',...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_udc',formats,renderer);
    mymatlab2tikz(pathname,'pdf_udc.tex');

    %% Display pdf of fracture toughness
    figure('Name','Probability Density Estimate: Fracture toughness')
    clf
    plot(gc_xi*1e-3,gc_f,'-b','LineWidth',linewidth)
    hold on
    ind_gc = find(gc_xi>=gc_ci(1) & gc_xi<gc_ci(2));
    area(gc_xi(ind_gc)*1e-3,gc_f(ind_gc),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2)
    scatter(gc_mean*1e-3,0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b')
    plot(gc_xiunif*1e-3,gc_funif,'-r','LineWidth',linewidth)
    ind_gcunif = find(gc_xiunif>=gc_ciunif(1) & gc_xiunif<gc_ciunif(2));
    area(gc_xiunif(ind_gcunif)*1e-3,gc_funif(ind_gcunif),'FaceColor','r','EdgeColor','none','FaceAlpha',0.2)
    scatter(gc_meanunif*1e-3,0,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','r')
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$g$ [N/mm]','Interpreter',interpreter)
    ylabel('$p_{G_c}(g)$','Interpreter',interpreter)
    l = legend('estimate pdf',...
        ['estimate $' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'estimate mean value',...
        'exact uniform pdf',...
        ['exact $' num2str((probs(2)-probs(1))*100) '\%$ confidence interval'],...
        'exact mean value');
    set(l,'Interpreter',interpreter)
    mysaveas(pathname,'pdf_gc',formats,renderer);
    mymatlab2tikz(pathname,'pdf_gc.tex');
end
