% clc
clearvars
close all
rng('default');

Dim = 2; % space dimension Dim = 2, 3
PFmodel = 'HeFreddi'; % 'Bourdin', 'Amor', 'Miehe', 'HeAmor', 'HeFreddi', 'Zhang'
PFsplit = 'Strain'; % 'Strain' or 'Stress'

%% Domains and meshes
L = 1; % domain size
e = 1; % thickness
RHO = 1; % density
switch Dim
    case 1
        D = DOMAIN(1,0.0,L);
    case 2
         D = DOMAIN(2,[0.0,0.0],[L,L]);
         % elemtype = 'TRI3';
         elemtype = 'QUA4';
         % elemtype = 'TRI6';
    case 3
        D = DOMAIN(3,[0.0,0.0,0.0],[L,L,L]);
        % elemtype = 'TET4';
        elemtype = 'CUB8';
        % elemtype = 'TET10';
end
% option = 'DEFO'; % plane strain
option = 'CONT'; % plane stress
nbelem = repmat(1,1,Dim);
S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'option',option);

S = final(S);
d = calc_init_dirichlet(S);

%% Materials
% Young modulus
E = 210; % [Pa]
% Poisson ratio
NU = 0.3;
% Material
mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',e,'d',d,'u',0,'PFM',PFmodel,'PFS',PFsplit);
mat = setnumber(mat,1);
S = setmaterial(S,mat);

S = final(S);

u = calc_init_dirichlet(S);
mat = setparam(mat,'u',u);
mats = MATERIALS(S);
for m=1:length(mats)
    mats{m} = setparam(mats{m},'u',u);
end
S = actualisematerials(S,mats);

% A = calc_rigi(S);

elem = getgroupelem(S,1);
node = getnode(S);
xnode = node(elem);
% gauss = calc_gauss(elem,'rigi');
gauss = calc_gauss(elem,1);
xgauss = gauss.coord;

B = calc_B(elem,xnode,xgauss);

u = evalparam(mat,'u',elem,xnode,xgauss); % displacement field
ue = localize(elem,u);
se = B*ue; % strain field in Voigt notation

se(:,:) = 2*rand(Dim*(Dim+1)/2,1)-1; % strain field in Voigt notation
se(2) = -se(2); % for negative trace

P = calc_proj_notation(elem);
se_KM = P\se; % strain field in Kelvin-Mandel notation

switch lower(PFmodel)
    case 'bourdin'
        Dp = calc_opmat(mat,elem,xnode,xgauss);
        Dm = 0;
    case 'amor'
        [Dp,Dm] = calc_opmat_Amor(mat,elem,xnode,xgauss,se,'check');
    case 'miehe'
        [Dp,Dm] = calc_opmat_Miehe(mat,elem,xnode,xgauss,se,'check');
    case {'heamor','hefreddi'}
        [Dp,Dm] = calc_opmat_He(mat,elem,xnode,xgauss,se,'check');
    case 'zhang'
        [Dp,Dm] = calc_opmat_Zhang(mat,elem,xnode,xgauss,se,'check');
    case 'spectral'
        [Dp,Dm] = calc_opmat_spectral(mat,elem,xnode,xgauss,se);
    case 'doublespectral'
        [Dp,Dm] = calc_opmat_doublespectral(mat,elem,xnode,xgauss);
    otherwise
        error(['Wrong phase field model ' model])
end
