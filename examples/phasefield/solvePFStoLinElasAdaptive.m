function [ft_sample,dt_sample,ut_sample,St_phase_sample,St_sample] = solvePFStoLinElasAdaptive(S_phase,S,T,fun,N,varargin)
% function [ft_sample,dt_sample,ut_sample,St_phase_sample,St_sample] = solvePFStoLinElasAdaptive(S_phase,S,T,fun,N,varargin)
% Solve stochastic Phase Field problem with mesh adaptation.

fun = fcnchk(fun);
filename = getcharin('filename',varargin,'filename');
pathname = getcharin('pathname',varargin,'.');
nbSamples = getcharin('nbsamples',varargin,1);

Dim = getdim(S);

% Initialize samples
ft_sample = zeros(N,length(T));
dt_sample = cell(nbSamples,length(T));
ut_sample = cell(nbSamples,length(T));
St_phase_sample = cell(nbSamples,length(T));
St_sample = cell(nbSamples,length(T));
% fmax_sample = zeros(N,1);

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
%         if isparam(mat,'delta') && any(getparam(mat,'delta')>0) % random phase field parameters
%             xnode = node_phase(elem);
%             gauss = calc_gauss(elem,'mass');
%             xgauss = gauss.coord;
%             delta = getparam(mat,'delta'); % coefficients of variation for fracture toughness and regularization parameter
%             if length(delta)==1
%                 deltaGc = delta; % 0 <= deltaGc < 1/sqrt(2). coefficient of variation for fracture toughness
%                 deltaL = delta; % 0 <= deltaL < 1/sqrt(2). coefficient of variation for regularization parameter
%             else
%                 deltaGc = delta(1); % 0 <= deltaGc < 1/sqrt(2). coefficient of variation for fracture toughness
%                 deltaL = delta(2); % 0 <= deltaL < 1/sqrt(2). coefficient of variation for regularization parameter
%             end
%             if deltaGc<0 || deltaGc>=1/sqrt(2)
%                 error('Coefficient of variation delta = %g for fracture toughness should be between 0 and %g',deltaGc,1/sqrt(2))
%             end
%             if deltaL<0 || deltaL>=1/sqrt(2)
%                 error('Coefficient of variation delta = %g for regularization parameter should be between 0 and %g',deltaL,1/sqrt(2))
%             end
%             nU = 2;
%             if deltaGc==0 || deltaL==0
%                 nU = 1;
%             end
%             if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
%                 lcorr = getparam(mat,'lcorr'); % spatial correlation length
%                 x = calc_x(elem,xnode,xgauss);
%                 x = getcoord(NODE(POINT(x(:,:,:))));
%                 [~,~,Z,Phi,k,c] = shinozukaSample(si,x,lcorr,nU); % sample for bivariate Gaussian random field with statistically independent normalized Gaussian components
%                 shinozukaPF = @(x) shinzukaFun(x,Z,Phi,k,c);
%                 mats_phase{m} = setparam(mats_phase{m},'shinozuka',shinozukaPF);
%             else % random vector model
%                 k = evalparam(mat,'k',elem,xnode,xgauss);
%                 r = evalparam(mat,'r',elem,xnode,xgauss);
%                 gc = sqrt(k.*r); % mean fracture toughness
%                 l = sqrt(k./r); % mean regularization parameter
%                 aGc = 1/deltaGc^2; % aGc > 2
%                 bGc = gc/aGc; % 0 < bGc = gc/aGc < gc/2 since gc > 0 and aGc > 2
%                 aL = 1/deltaL^2; % aL > 2
%                 bL = l/aL; % 0 < bL = l/aL < l/2 since l > 0 and aL > 2
%                 Xi = randn(si,1,nU); % sample for bivariate Gaussian random variable with statistically independent normalized Gaussian components
%                 if deltaGc && deltaL
%                     rho = 0;
%                     if isparam(mat,'rcorr')
%                         rho = getparam(mat,'rcorr'); % correlation coefficient between fracture toughness and regularization parameter
%                     end
%                     gc = gaminv(normcdf(Xi(:,1)),aGc,bGc); % sample for fracture toughness [N/m^2]
%                     l = gaminv(normcdf(rho*Xi(:,1) + sqrt(1-rho^2)*Xi(:,2)),aL,bL); % sample for regularization parameter [m]
%                 elseif deltaGc
%                     gc = gaminv(normcdf(Xi(:,1)),aGc,bGc); % sample for fracture toughness [N/m^2]
%                 else
%                     l = gaminv(normcdf(Xi(:,1)),aL,bL); % sample for regularization parameter [m]
%                 end
%                 k = gc.*l;
%                 r = gc./l;
%                 mats_phase{m} = setparam(mats_phase{m},'k',k);
%                 mats_phase{m} = setparam(mats_phase{m},'r',r);
%             end
%         end
        if isparam(mat,'aGc') && isparam(mat,'bGc') % random phase field parameters
            xnode = node_phase(elem);
            gauss = calc_gauss(elem,'mass');
            xgauss = gauss.coord;
            nU = 1;
            if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
                lcorr = getparam(mat,'lcorr'); % spatial correlation length
                x = calc_x(elem,xnode,xgauss);
                x = getcoord(NODE(POINT(x(:,:,:))));
                [~,~,Z,Phi,k,c] = shinozukaSample(si,x,lcorr,nU); % sample for univariate Gaussian random field
                shinozukaPF = @(x) shinzukaFun(x,Z,Phi,k,c);
                mats_phase{m} = setparam(mats_phase{m},'shinozuka',shinozukaPF);
            else % random variable model
                k = evalparam(mat,'k',elem,xnode,xgauss);
                r = evalparam(mat,'r',elem,xnode,xgauss);
                l = sqrt(k./r); % regularization parameter
                aGc = getparam(mat,'aGc'); % lower bound(s) for fracture toughness aGc > 0
                bGc = getparam(mat,'bGc'); % upper bound(s) for fracture toughness bGc > aGc > 0
                if ~isscalar(aGc) || ~isscalar(bGc)
                    imax = max(length(aGc),length(bGc));
                    iGc = randi(si,imax); % support index for univariate uniform distribution with multiple lower and upper endpoints aGc and bGc
                    if ~isscalar(aGc)
                        aGc = aGc(iGc); % lower bound for fracture toughness aGc > 0
                    end
                    if ~isscalar(bGc)
                        bGc = bGc(iGc); % upper bound for fracture toughness bGc > aGc > 0
                    end
                    mats_phase{m} = setparam(mats_phase{m},'aGc',aGc);
                    mats_phase{m} = setparam(mats_phase{m},'bGc',bGc);
                end
                if aGc<0 || bGc<0
                    error('Lower bound a = %g and upper bound b = %g for fracture toughness must be positive (superior to 0)',aGc,bGc)
                end
                if aGc>bGc
                    error('Lower bound a = %g must be inferior to upper bound b = %g for fracture toughness',aGc,bGc)
                end
                Xi = randn(si,1,nU); % sample for univariate Gaussian random variable
                gc = unifinv(normcdf(Xi),aGc,bGc); % sample for fracture toughness [N/m^2]
                k = gc.*l;
                r = gc./l;
                mats_phase{m} = setparam(mats_phase{m},'k',k);
                mats_phase{m} = setparam(mats_phase{m},'r',r);
            end
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
            xnode = node(elem);
            gauss = calc_gauss(elem,'rigi');
            xgauss = gauss.coord;
            if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
                lcorr = getparam(mat,'lcorr'); % spatial correlation length
                x = calc_x(elem,xnode,xgauss);
                x = getcoord(NODE(POINT(x(:,:,:))));
                n = Dim*(Dim+1)/2;
                nU = n*(n+1)/2; % multivariate Gaussian random field for anisotropic material
                if isa(mat,'ELAS_ISOT') % almost surely isotropic material
                    nU = 2; % bivariate Gaussian random field
                end
                [~,~,Z,Phi,k,c] = shinozukaSample(si,x,lcorr,nU); % sample for multivariate Gaussian random field with statistically independent normalized Gaussian components
                shinozukaMat = @(x) shinozukaFun(x,Z,Phi,k,c);
                mats{m} = setparam(mats{m},'shinozuka',shinozukaMat);
            else % random matrix model
                if isa(mat,'ELAS_ISOT') % almost surely isotropic material
                    E = evalparam(mat,'E',elem,xnode,xgauss); % mean Young modulus
                    NU = evalparam(mat,'NU',elem,xnode,xgauss); % mean Poisson ratio
%                     % la = -24; % la < 1/5. Parameter controlling the level of statistical fluctuations
%                     % deltaC1 = 1/sqrt(1-la); % coefficient of variation for bulk modulus
%                     % deltaC2 = 1/sqrt(1-5*la); % coefficient of variation for shear modulus
%                     deltaC1 = getparam(mat,'delta'); % coefficient of variation for bulk modulus
%                     if deltaC1<0
%                         error('Coefficient of variation delta = %g for bulk modulus must be positive (superior to 0)',deltaC1)
%                     end
%                     la = 1 - 1/deltaC1^2; % la < 1/5. Parameter controlling the level of statistical fluctuations
%                     deltaC2 = 1/sqrt(5/deltaC1^2 - 4); % coefficient of variation for shear modulus
%                     mC1 = E/3/(1-2*NU); % mean bulk modulus
%                     mC2 = E/(1+NU)/2; % mean shear modulus
%                     laC1 = (1-la)/mC1; % la1 > 0
%                     laC2 = (1-5*la)/mC2; % la2 > 0
%                     aC1 = 1-la; % a1 > 0
%                     bC1 = 1/laC1; % b1 > 0
%                     aC2 = 1-5*la; % a2 > 0
%                     bC2 = 1/laC2; % b2 > 0
%                     nU = 2;
%                     Xi = randn(si,1,nU); % sample for bivariate Gaussian random variable with statistically independent normalized Gaussian components
%                     rho = 0;
%                     if isparam(mat,'rcorr')
%                         rho = getparam(mat,'rcorr'); % correlation coefficient between bulk and shear moduli
%                     end
%                     C1 = gaminv(normcdf(Xi(:,1)),aC1,bC1); % sample for bulk modulus [Pa]
%                     C2 = gaminv(normcdf(rho*Xi(:,1) + sqrt(1-rho^2)*Xi(:,2)),aC2,bC2); % sample for shear modulus [Pa]
%                     % lambda = C1 - 2/3*C2; % [Pa]
%                     E = (9*C1.*C2)./(3*C1+C2); % [Pa]
%                     NU = (3*C1-2*C2)./(6*C1+2*C2);
                    delta = getparam(mat,'delta'); % coefficients of variation for Young modulus and Poisson ratio
                    if length(delta)==1
                        deltaE = delta; % 0 <= deltaE < 1/sqrt(2). coefficient of variation for Young modulus
                        deltaNU = delta; % coefficient of variation for Poisson ratio
                   else
                        deltaE = delta(1); % 0 <= deltaE < 1/sqrt(2). coefficient of variation for Young modulus
                        deltaNU = delta(2); % coefficient of variation for Poisson ratio
                    end
                    if deltaE<0 || deltaE>=1/sqrt(2)
                        error(['Coefficient of variation delta = %d for Young modulus should be between 0 and %d',deltaE,1/sqrt(2)]);
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
                    Xi = randn(si,1,nU); % sample for bivariate Gaussian random variable with statistically independent normalized Gaussian components
                    if deltaE && deltaNU
                        E = gaminv(normcdf(Xi(:,1)),aE,bE); % sample for Young modulus [Pa]
                        NU = betainv(normcdf(Xi(:,2)),a2NU,b2NU)/2; % sample for Poisson ratio
                    elseif deltaE
                        E = gaminv(normcdf(Xi(:,1)),aE,bE); % sample for Young modulus [Pa]
                    else
                        NU = betainv(normcdf(Xi(:,1)),a2NU,b2NU)/2; % sample for Poisson ratio
                    end
                    mats{m} = setparam(mats{m},'E',E);
                    mats{m} = setparam(mats{m},'NU',NU);
                elseif isa(mat,'ELAS_ANISOT') % anisotropic material
                    Cmat = evalparam(mat,'C',elem,xnode,xgauss); % mean elasticity matrix
                    delta = getparam(mat,'delta'); % coefficient of variation for elasticity matrix
                    mL = chol(Cmat); % upper triangular matrix of the Cholesky factor of mean elasticity matrix
                    n = size(Cmat,1);
                    nU = n*(n+1)/2;
                    Xi = randn(si,nU,1); % sample for multivariate Gaussian random variable with statistically independent normalized Gaussian components
                    Cmat = randAnisotElasMatrix(delta,mL,Xi); % sample for non-Gaussian random elasticity matrix
                    mats{m} = setparam(mats{m},'C',Cmat);
                else
                    error('Wrong material symmetry class');
                end
            end
        end
    end
    Si = actualisematerials(Si,mats);
    
    % Copy msh, mesh and sol files
    filenamei = [filename '_' num2str(i)];
    G = GMSHFILE(fullfile(pathname,filename));
    Gi = GMSHFILE(fullfile(pathname,filenamei));
    command = ['cp ' getfilemsh(G) ' ' getfilemsh(Gi) ';'...
        'cp ' getfilemesh(G) ' ' getfilemesh(Gi) ';'...
        'cp ' getfilesol(G) ' ' getfilesol(Gi)];
    dos(command);

    % Solve deterministic problem
    [dt,ut,ft,St_phase,St] = fun(S_phasei,Si,filenamei);
    ft_sample(i,:) = ft;
    if i<=nbSamples
        dt_sample(i,:) = dt;
        ut_sample(i,:) = ut;
        St_phase_sample(i,:) = St_phase;
        St_sample(i,:) = St;
    end
    % fmax_sample(i) = max(ft);
end
textprogressbar(' done');

function nUpdateProgressBar(~)
j = j+1;
textprogressbar(j/N*100,sprintf('(%d/%d)',j,N));
end

end
