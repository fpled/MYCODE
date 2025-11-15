i=1;
si = RandStream.create('mrg32k3a','NumStreams',N,'StreamIndices',i);
Si=S;
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
            % % la = -24; % la < 1/5. Parameter controlling the level of statistical fluctuations
            % % deltaC1 = 1/sqrt(1-la); % coefficient of variation for bulk modulus
            % % deltaC2 = 1/sqrt(1-5*la); % coefficient of variation for shear modulus
            % deltaC1 = getparam(mat,'delta'); % coefficient of variation for bulk modulus
            % if deltaC1<0
            %     error('Coefficient of variation delta = %g for bulk modulus must be positive (superior to 0)',deltaC1)
            % end
            % la = 1 - 1/deltaC1^2; % la < 1/5. Parameter controlling the level of statistical fluctuations
            % deltaC2 = 1/sqrt(5/deltaC1^2 - 4); % coefficient of variation for shear modulus
            % mC1 = E/3/(1-2*NU); % mean bulk modulus
            % mC2 = E/(1+NU)/2; % mean shear modulus
            % laC1 = (1-la)/mC1; % la1 > 0
            % laC2 = (1-5*la)/mC2; % la2 > 0
            % aC1 = 1-la; % a1 > 0
            % bC1 = 1/laC1; % b1 > 0
            % aC2 = 1-5*la; % a2 > 0
            % bC2 = 1/laC2; % b2 > 0
            % nU = 2;
            % if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
            %     lcorr = getparam(mat,'lcorr'); % spatial correlation length
            %     x = calc_x(elem,xnode,xgauss);
            %     x = getcoord(NODE(POINT(x(:,:,:))));
            %     Xi = shinozukaSample(si,x,lcorr,nU); % sample for bivariate Gaussian random field with statistically independent normalized Gaussian components
            % else % random matrix model
            %     Xi = randn(si,1,nU); % sample for bivariate Gaussian random variable with statistically independent normalized Gaussian components
            % end
            % rho = 0;
            % if isparam(mat,'rcorr')
            %     rho = getparam(mat,'rcorr'); % correlation coefficient between bulk and shear moduli
            % end
            % C1 = gaminv(normcdf(Xi(:,1)),aC1,bC1); % sample for bulk modulus [Pa]
            % C2 = gaminv(normcdf(rho*Xi(:,1) + sqrt(1-rho^2)*Xi(:,2)),aC2,bC2); % sample for shear modulus [Pa]
            % if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
            %     C1 = reshape(C1,1,1,nbelem,gauss.nbgauss);
            %     C2 = reshape(C2,1,1,nbelem,gauss.nbgauss);
            %     C1 = MYDOUBLEND(C1);
            %     C2 = MYDOUBLEND(C2);
            %     C1 = FEELEMFIELD({C1},'storage','gauss','type','scalar','ddl',DDL('C1'));
            %     C2 = FEELEMFIELD({C2},'storage','gauss','type','scalar','ddl',DDL('C2'));
            % end
            % % lambda = C1 - 2/3*C2; % [Pa]
            % E = (9*C1.*C2)./(3*C1+C2); % [Pa]
            % NU = (3*C1-2*C2)./(6*C1+2*C2);
            delta = getparam(mat,'delta'); % coefficients of variation for Young modulus and Poisson ratio
            deltaNUsup = min(sqrt((1-2*NU)/(1+2*NU)),(1-2*NU)/2/sqrt(NU*(1-NU)));
            if isscalar(delta)
                deltaE = delta; % 0 <= deltaE < 1/sqrt(2). coefficient of variation for Young modulus
                deltaNU = delta; % 0 <= deltaNU < deltaNUsup. coefficient of variation for Poisson ratio
            else
                deltaE = delta(1); % 0 <= deltaE < 1/sqrt(2). coefficient of variation for Young modulus
                deltaNU = delta(2); % 0 <= deltaNU < deltaNUsup. coefficient of variation for Poisson ratio
            end
            if deltaE<0 || deltaE>=1/sqrt(2)
                error(['Coefficient of variation delta = %d for Young modulus should be between 0 and %d',deltaE,1/sqrt(2)]);
            end
            if deltaNU<0 || deltaNU>=deltaNUsup
                error(['Coefficient of variation delta = %d for Poisson ratio should be between 0 and %d',deltaNU,deltaNUsup]);
            end
            aE = 1/deltaE^2; % aE > 2
            bE = E/aE; % 0 < bE = E/aE < E/2 since E > 0 and aE > 2
            m2NU = 2*NU; % 0 < m2NU < 1
            a2NU = (1-m2NU)/deltaNU^2-m2NU; % a2NU > 1 since deltaNU < deltaNUsup < sqrt((1-m2NU)/(1+m2NU))
            b2NU = a2NU/m2NU-a2NU; % b2NU > 1 since deltaNU < deltaNUsup < (1-2*NU)/2/sqrt(NU*(1-NU)) = (1-m2NU)/sqrt(m2NU*(2-m2NU))
            nU = 2;
            if deltaE==0 || deltaNU==0
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
            nU = n*(n+1)/2;
            if isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
                lcorr = getparam(mat,'lcorr'); % spatial correlation length
                x = calc_x(elem,xnode,xgauss);
                x = getcoord(NODE(POINT(x(:,:,:))));
                Xi = shinozukaSample(si,x,lcorr,nU); % sample for multivariate Gaussian random field with statistically independent normalized Gaussian components
                Cmat = randAnisotElasField(delta,mL,shiftdim(Xi,1)); % sample for non-Gaussian random elasticity field
                Cmat = Cmat(:,:,:); % n-by-n-by-nx array
                Cmat = reshape(Cmat,n,n,nbelem,gauss.nbgauss);
                Cmat = MYDOUBLEND(Cmat);
                syscoordgauss = getsyscoordlocal(elem);
                fieldddl = DDL(DDLTENS4('C',syscoordgauss));
                Cmat = FEELEMFIELD({Cmat},'storage','gauss','type','scalar','ddl',fieldddl);
            else % random matrix model
                Xi = randn(si,nU,1); % sample for multivariate Gaussian random variable with statistically independent normalized Gaussian components
                Cmat = randAnisotElasMatrix(delta,mL,Xi); % sample for non-Gaussian random elasticity matrix
            end
            mats{m} = setparam(mats{m},'C',Cmat);
        else
            error('Wrong material symmetry class');
        end
    end
end
mats
Si = actualisematerials(Si,mats);
Eelem = calc_parammat(Si,'E');
Enode = calc_parammat(Si,'E','smooth');
NUelem = calc_parammat(Si,'NU');
NUnode = calc_parammat(Si,'NU','smooth');

figure('Name','Young''s modulus element field')
clf
plot(Eelem*1e-9,Si)
colorbar
set(gca,'FontSize',fontsize)
mysaveas(pathname,['Young_modulus_sample_' num2str(i) '_elem'],formats,renderer);


figure('Name','Young''s modulus nodal field')
clf
plot(Enode*1e-9,Si)
colorbar
set(gca,'FontSize',fontsize)
mysaveas(pathname,['Young_modulus_sample_' num2str(i) '_node'],formats,renderer);

figure('Name','Poisson''s ratio element field')
clf
plot(NUelem,Si)
colorbar
set(gca,'FontSize',fontsize)
mysaveas(pathname,['Poisson_ratio_sample_' num2str(i) '_elem'],formats,renderer);

figure('Name','Poisson''s ratio nodal field')
clf
plot(NUnode,Si)
colorbar
set(gca,'FontSize',fontsize)
mysaveas(pathname,['Poisson_ratio_sample_' num2str(i) '_node'],formats,renderer);