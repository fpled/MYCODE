function S = setmaterialproperties(S,materials)
% function S = setmaterialproperties(S,materials)

node = getnode(S);
for m=1:length(materials)
    mat = materials{m};
    if isparam(mat,'delta') && getparam(mat,'delta')>0 && isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
        elem = getgroupelem(S,m);
        nbelem = getnbelem(elem);
        xnode = node(elem);
        gauss = calc_gauss(elem,'rigi');
        xgauss = gauss.coord;
        x = calc_x(elem,xnode,xgauss);
        x = getcoord(NODE(POINT(x(:,:,:))));
        shinozukaMat = getparam(mat,'shinozuka');
        Xi = shinozukaMat(x); % sample for multivariate Gaussian random field with statistically independent normalized Gaussian components
        if isa(mat,'ELAS_ISOT') % almost surely isotropic material
            E = evalparam(mat,'E',elem,xnode,xgauss);
            NU = evalparam(mat,'NU',elem,xnode,xgauss);
%             % la = -24; % la < 1/5. Parameter controlling the level of statistical fluctuations
%             % deltaC1 = 1/sqrt(1-la); % coefficient of variation for bulk modulus
%             % deltaC2 = 1/sqrt(1-5*la); % coefficient of variation for shear modulus
%             deltaC1 = getparam(mat,'delta'); % coefficient of variation for bulk modulus
%             la = 1 - 1/deltaC1^2; % la < 1/5. Parameter controlling the level of statistical fluctuations
%             deltaC2 = 1/sqrt(5/deltaC1^2 - 4); % coefficient of variation for shear modulus
%             mC1 = E/3/(1-2*NU); % mean bulk modulus
%             mC2 = E/(1+NU)/2; % mean shear modulus
%             laC1 = (1-la)/mC1; % la1 > 0
%             laC2 = (1-5*la)/mC2; % la2 > 0
%             aC1 = 1-la; % a1 > 0
%             bC1 = 1/laC1; % b1 > 0
%             aC2 = 1-5*la; % a2 > 0
%             bC2 = 1/laC2; % b2 > 0
%             rho = 0;
%             if isparam(mat,'rcorr')
%                 rho = getparam(mat,'rcorr'); % correlation coefficient between bulk and shear moduli
%             end
%             C1 = gaminv(normcdf(Xi(:,1)),aC1,bC1); % sample for bulk modulus [Pa]
%             C2 = gaminv(normcdf(rho*Xi(:,1) + sqrt(1-rho^2)*Xi(:,2)),aC2,bC2); % sample for shear modulus [Pa]
%             C1 = reshape(C1,1,1,nbelem,gauss.nbgauss);
%             C2 = reshape(C2,1,1,nbelem,gauss.nbgauss);
%             C1 = MYDOUBLEND(C1);
%             C2 = MYDOUBLEND(C2);
%             C1 = FEELEMFIELD({C1},'storage','gauss','type','scalar','ddl',DDL('C1'));
%             C2 = FEELEMFIELD({C2},'storage','gauss','type','scalar','ddl',DDL('C2'));
%             % lambda = C1 - 2/3*C2; % [Pa]
%             E = (9*C1.*C2)./(3*C1+C2); % [Pa]
%             NU = (3*C1-2*C2)./(6*C1+2*C2);
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
            if deltaE && deltaNU
                E = gaminv(normcdf(Xi(:,1)),aE,bE); % sample for Young modulus [Pa]
                NU = betainv(normcdf(Xi(:,2)),a2NU,b2NU)/2; % sample for Poisson ratio
            elseif deltaE
                E = gaminv(normcdf(Xi(:,1)),aE,bE); % sample for Young modulus [Pa]
            else
                NU = betainv(normcdf(Xi(:,1)),a2NU,b2NU)/2; % sample for Poisson ratio
            end
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
            mat = setparam(mat,'E',E);
            mat = setparam(mat,'NU',NU);
        elseif isa(mat,'ELAS_ANISOT') % anisotropic material
            Cmat = evalparam(mat,'C',elem,xnode,xgauss); % mean elasticity matrix
            delta = getparam(mat,'delta'); % coefficient of variation for elasticity matrix/field
            mL = chol(Cmat); % upper triangular matrix of the Cholesky factor of mean elasticity matrix
            Cmat = randAnisotElasField(delta,mL,shiftdim(Xi,1)); % sample for non-Gaussian random elasticity field
            Cmat = Cmat(:,:,:); % n-by-n-by-nx array
            Cmat = reshape(Cmat,n,n,nbelem,gauss.nbgauss);
            Cmat = MYDOUBLEND(Cmat);
            syscoordgauss = getsyscoordlocal(elem);
            fieldddl = DDL(DDLTENS4('C',syscoordgauss));
            Cmat = FEELEMFIELD({Cmat},'storage','gauss','type','scalar','ddl',fieldddl);
            mat = setparam(mat,'C',Cmat);
        else
            error('Wrong material symmetry class');
        end
    end
    S = setmaterial(S,mat,m);
end
