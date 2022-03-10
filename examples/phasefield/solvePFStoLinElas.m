function [ft_sample,dt_mean,ut_mean,dt_var,ut_var,dt_sample,ut_sample] = solvePFStoLinElas(S_phase,S,T,fun,N,varargin)
% function [ft_sample,dt_mean,ut_mean,dt_var,ut_var,dt_sample,ut_sample] = solvePFStoLinElas(S_phase,S,T,fun,N,varargin)
% Solve stochastic Phase Field problem.

fun = fcnchk(fun);
nbSamples = getcharin('nbsamples',varargin,1);

sz_d = getnbddl(S_phase);
sz_u = getnbddl(S);

% Initialize samples
ft_sample = zeros(N,length(T));
dt_sample = zeros(nbSamples,sz_d,length(T));
ut_sample = zeros(nbSamples,sz_u,length(T));
% fmax_sample = zeros(N,1);

% Initialize statistical means and second-order moments
dt_mean = zeros(sz_d,length(T));
dt_moment2 = zeros(sz_d,length(T));
ut_mean = zeros(sz_u,length(T));
ut_moment2 = zeros(sz_u,length(T));

if ~verLessThan('matlab','9.2') % introduced in R2017a
    q = parallel.pool.DataQueue;
    afterEach(q, @nUpdateProgressBar);
else
    q = [];
end
textprogressbar('Solving problem: ');
j = 0;
parfor i=1:N
    if ~verLessThan('matlab','9.2') % introduced in R2017a
        send(q,i);
    end
    si = RandStream.create('mrg32k3a','NumStreams',N,'StreamIndices',i);

    % Generate random phase field parameters
    S_phasei = S_phase;
    mats_phase = MATERIALS(S_phase);
    node_phase = getnode(S_phase);
    for m=1:getnbgroupelem(S_phase)
        elem = getgroupelem(S_phase,m);
        mat = getmaterial(elem);
        if isparam(mat,'delta') && any(getparam(mat,'delta')>0) % random phase field parameters
            nbelem = getnbelem(elem);
            xnode = node_phase(elem);
            gauss = calc_gauss(elem,'mass');
            xgauss = gauss.coord;
            k = evalparam(mat,'k',elem,xnode,xgauss);
            r = evalparam(mat,'r',elem,xnode,xgauss);
            gc = sqrt(k.*r); % mean fracture toughness
            l = sqrt(k./r); % mean regularization parameter
            delta = getparam(mat,'delta'); % coefficients of variation for fracture toughness and regularization parameter
            if length(delta)==1
                deltaGc = delta; % coefficient of variation for fracture toughness
                deltaL = delta; % coefficient of variation for regularization parameter
            else
                deltaGc = delta(1); % coefficient of variation for fracture toughness
                deltaL = delta(2); % coefficient of variation for regularization parameter
            end
            aGc = 1/deltaGc^2;
            bGc = gc/aGc;
            aL = 1/deltaL^2;
            bL = l/aL;
            nU = 2;
            if deltaGc==0 || deltaL==0
                nU = 1;
            end
            if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
                lcorr = getparam(mat,'lcorr'); % spatial correlation length
                x = calc_x(elem,xnode,xgauss);
                x = getcoord(NODE(POINT(x(:,:,:))));
                Xi = shinozukaSample(si,x,lcorr,nU); % sample for bivariate Gaussian random field with statistically independent normalized Gaussian components
            else % random matrix model
                Xi = randn(si,1,nU); % sample for bivariate Gaussian random variable with statistically independent normalized Gaussian components
            end
            if deltaGc && deltaL
                rho = 0;
                if isparam(mat,'rcorr')
                    rho = getparam(mat,'rcorr'); % correlation coefficient between fracture toughness and regularization parameter
                end
                gc = gaminv(normcdf(Xi(:,1)),aGc,bGc); % sample for fracture toughness [N/m^2]
                l = gaminv(normcdf(rho*Xi(:,1) + sqrt(1-rho^2)*Xi(:,2)),aL,bL); % sample for regularization parameter [m]
            elseif deltaGc
                gc = gaminv(normcdf(Xi(:,1)),aGc,bGc); % sample for fracture toughness [N/m^2]
            else
                l = gaminv(normcdf(Xi(:,1)),aL,bL); % sample for regularization parameter [m]
            end
            if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
                if deltaGc
                    gc = reshape(gc,1,1,nbelem,gauss.nbgauss);
                    gc = MYDOUBLEND(gc);
                    gc = FEELEMFIELD({gc},'storage','gauss','type','scalar','ddl',DDL('gc'));
                end
                if deltaL
                    l = reshape(l,1,1,nbelem,gauss.nbgauss);
                    l = MYDOUBLEND(l);
                    l = FEELEMFIELD({l},'storage','gauss','type','scalar','ddl',DDL('l'));
                end
            end
            k = gc.*l;
            r = gc./l;
            mats_phase{m} = setparam(mats_phase{m},'k',k);
            mats_phase{m} = setparam(mats_phase{m},'r',r);
        end
    end
    S_phasei = actualisematerials(S_phasei,mats_phase);
    
    % Generate random material parameters
    Si = S;
    mats = MATERIALS(S);
    node = getnode(S);
    for m=1:getnbgroupelem(S)
        elem = getgroupelem(S,m);
        mat = getmaterial(elem);
        if isparam(mat,'delta') && getparam(mat,'delta')>0 % random material parameters
            nbelem = getnbelem(elem);
            xnode = node(elem);
            gauss = calc_gauss(elem,'rigi');
            xgauss = gauss.coord;
            if isa(mat,'ELAS_ISOT') % almost surely isotropic material
                E = evalparam(mat,'E',elem,xnode,xgauss);
                NU = evalparam(mat,'NU',elem,xnode,xgauss);
%                 % la = -24; % la < 1/5. Parameter controlling the level of statistical fluctuations
%                 % deltaC1 = 1/sqrt(1-la); % coefficient of variation for bulk modulus
%                 % deltaC2 = 1/sqrt(1-5*la); % coefficient of variation for shear modulus
%                 deltaC1 = getparam(mat,'delta'); % coefficient of variation for bulk modulus
%                 la = 1 - 1/deltaC1^2; % la < 1/5. Parameter controlling the level of statistical fluctuations
%                 deltaC2 = 1/sqrt(5/deltaC1^2 - 4); % coefficient of variation for shear modulus
%                 mC1 = E/3/(1-2*NU); % mean bulk modulus
%                 mC2 = E/(1+NU)/2; % mean shear modulus
%                 laC1 = (1-la)/mC1; % la1 > 0
%                 laC2 = (1-5*la)/mC2; % la2 > 0
%                 aC1 = 1-la; % a1 > 0
%                 bC1 = 1/laC1; % b1 > 0
%                 aC2 = 1-5*la; % a2 > 0
%                 bC2 = 1/laC2; % b2 > 0
%                 if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
%                     lcorr = getparam(mat,'lcorr');
%                     x = calc_x(elem,xnode,xgauss);
%                     x = getcoord(NODE(POINT(x(:,:,:))));
%                     Xi = shinozukaSample(si,x,lcorr,2); % sample for bivariate Gaussian random field with statistically independent normalized Gaussian components
%                 else % random matrix model
%                     Xi = randn(si,1,2); % sample for bivariate Gaussian random variable with statistically independent normalized Gaussian components
%                 end
%                 rho = 0;
%                 if isparam(mat,'rcorr')
%                     rho = getparam(mat,'rcorr'); % correlation coefficient between bulk and shear moduli
%                 end
%                 C1 = gaminv(normcdf(Xi(:,1)),aC1,bC1); % sample for bulk modulus [Pa]
%                 C2 = gaminv(normcdf(rho*Xi(:,1) + sqrt(1-rho^2)*Xi(:,2)),aC2,bC2); % sample for shear modulus [Pa]
%                 if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
%                     C1 = reshape(C1,1,1,nbelem,gauss.nbgauss);
%                     C2 = reshape(C2,1,1,nbelem,gauss.nbgauss);
%                     C1 = MYDOUBLEND(C1);
%                     C2 = MYDOUBLEND(C2);
%                     C1 = FEELEMFIELD({C1},'storage','gauss','type','scalar','ddl',DDL('C1'));
%                     C2 = FEELEMFIELD({C2},'storage','gauss','type','scalar','ddl',DDL('C2'));
%                 end
%                 % lambda = C1 - 2/3*C2; % [Pa]
%                 E = (9*C1.*C2)./(3*C1+C2); % [Pa]
%                 NU = (3*C1-2*C2)./(6*C1+2*C2);
                delta = getparam(mat,'delta'); % coefficients of variation for Young modulus and Poisson ratio
                if length(delta)==1
                    deltaE = delta; % 0 <= deltaE < 1/sqrt(2). coefficient of variation for Young modulus
                    deltaNU = delta; % coefficient of variation for Poisson ratio
                else
                    deltaE = delta(1); % 0 <= deltaE < 1/sqrt(2). coefficient of variation for Young modulus
                    deltaNU = delta(2); % coefficient of variation for Poisson ratio
                end
                if deltaE>=1/sqrt(2)
                    error(['Coefficient of variation for Young modulus must be < 1/sqrt(2) = ' num2str(1/sqrt(2))]);
                end
                aE = 1/deltaE^2; % aE > 2
                bE = E/aE; % 0 < bE = E/aE < E/2 since E > 0 and aE > 2
                m2NU = 2*NU; % 0 < m2NU < 1
                a2NU = (1-m2NU)/deltaNU^2-m2NU; % a2NU > 0
                b2NU = a2NU/m2NU-a2NU; % b2NU > 0
                nU = 2;
                if deltaE==0 || deltaNU==0
                    nU = 1;
                end
                if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
                    lcorr = getparam(mat,'lcorr');
                    x = calc_x(elem,xnode,xgauss);
                    x = getcoord(NODE(POINT(x(:,:,:))));
                    Xi = shinozukaSample(si,x,lcorr,nU); % sample for bivariate Gaussian random field with statistically independent normalized Gaussian components
                else % random matrix model
                    Xi = randn(si,1,nU); % sample for bivariate Gaussian random variable with statistically independent normalized Gaussian components
                end
                if deltaE && deltaNU
                    E = gaminv(normcdf(Xi(:,1)),aE,bE); % sample for Young modulus [Pa]
                    NU = betainv(normcdf(Xi(:,2)),a2NU,b2NU)/2; % sample for Poisson ratio
                elseif deltaE
                    E = gaminv(normcdf(Xi(:,1)),aE,bE); % sample for Young modulus [Pa]
                else
                    NU = betainv(normcdf(Xi(:,1)),a2NU,b2NU)/2; % sample for Poisson ratio
                end
                if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
                    if deltaE
                        E = reshape(E,1,1,nbelem,gauss.nbgauss);
                        E = MYDOUBLEND(E);
                        E = FEELEMFIELD({E},'storage','gauss','type','scalar','ddl',DDL('E'));
                    end
                    if deltaNU
                        NU = reshape(NU,1,1,nbelem,gauss.nbgauss);
                        NU = MYDOUBLEND(NU);
                        NU = FEELEMFIELD({NU},'storage','gauss','type','scalar','ddl',DDL('NU'));
                    end
                end
                mats{m} = setparam(mats{m},'E',E);
                mats{m} = setparam(mats{m},'NU',NU);
            elseif isa(mat,'ELAS_ANISOT') % anisotropic material
                Cmat = evalparam(mat,'C',elem,xnode,xgauss); % mean elasticity matrix
                delta = getparam(mat,'delta'); % coefficient of variation for elasticity matrix/field
                mL = chol(Cmat); % upper triangular matrix of the Cholesky factor of mean elasticity matrix
                n = size(Cmat,1);
                if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
                    lcorr = getparam(mat,'lcorr'); % spatial correlation length
                    x = calc_x(elem,xnode,xgauss);
                    x = getcoord(NODE(POINT(x(:,:,:))));
                    Xi = shinozukaSample(si,x,lcorr,n*(n+1)/2); % sample for multivariate Gaussian random field with statistically independent normalized Gaussian components
                    Cmat = randAnisotElasField(delta,mL,shiftdim(Xi,1)); % sample for non-Gaussian random elasticity field
                    Cmat = Cmat(:,:,:); % n-by-n-by-nx array
                    Cmat = reshape(Cmat,n,n,nbelem,gauss.nbgauss);
                    Cmat = MYDOUBLEND(Cmat);
                    syscoordgauss = getsyscoordlocal(elem);
                    fieldddl = DDL(DDLTENS4('C',syscoordgauss));
                    Cmat = FEELEMFIELD({Cmat},'storage','gauss','type','scalar','ddl',fieldddl);
                else % random matrix model
                    Xi = randn(si,n*(n+1)/2,1); % sample for multivariate Gaussian random variable with statistically independent normalized Gaussian components
                    Cmat = randAnisotElasMatrix(delta,mL,Xi); % sample for non-Gaussian random elasticity matrix
                end
                mats{m} = setparam(mats{m},'C',Cmat);
            else
                error('Wrong material symmetry class');
            end
        end
    end
    Si = actualisematerials(Si,mats);
    
    % Solve deterministic problem
    [dt,ut,ft] = fun(S_phasei,Si);
    
    % Compute second-order statistics
    dt_val = getvalue(dt);
    dt_mean = dt_mean + dt_val/N;
    dt_moment2 = dt_moment2 + dt_val.^2/N;
    ut_val = getvalue(ut);
    ut_mean = ut_mean + ut_val/N;
    ut_moment2 = ut_moment2 + ut_val.^2/N;
    
    ft_sample(i,:) = ft;
    if i<=nbSamples
        dt_sample(i,:,:) = getvalue(dt);
        ut_sample(i,:,:) = getvalue(ut);
    end
    % fmax_sample(i) = max(ft);
end
textprogressbar(' done');

% Compute unbiased variances
if N>1
    dt_var = (N/(N-1))*(dt_moment2 - dt_mean.^2);
    ut_var = (N/(N-1))*(ut_moment2 - ut_mean.^2);
else
    dt_var = zeros(sz_d,length(T));
    ut_var = zeros(sz_u,length(T));
end

function nUpdateProgressBar(~)
j = j+1;
textprogressbar(j/N*100,sprintf('(%d/%d)',j,N));
end

end
