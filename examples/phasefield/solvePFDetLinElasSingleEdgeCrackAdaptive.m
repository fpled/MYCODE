function [dt,ut,ft,St_phase,St,Ht] = solvePFDetLinElasSingleEdgeCrackAdaptive(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,sizemap,varargin)
% function [dt,ut,ft,St_phase,St,Ht] = solvePFDetLinElasSingleEdgeCrackAdaptive(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,sizemap,varargin)
% Solve deterministic Phase Field problem with mesh adaptation.

display_ = ischarin('display',varargin);
filename = getcharin('filename',varargin,'gmsh_domain_single_edge_crack');
pathname = getcharin('pathname',varargin,'.');
gmshoptions = getcharin('gmshoptions',varargin,'-v 0');
mmgoptions = getcharin('mmgoptions',varargin,'-nomove -v -1');

Dim = getdim(S);

t = gett(T);

dt = cell(1,length(T));
ut = cell(1,length(T));
ft = zeros(1,length(T));
if nargout>=4
    St_phase = cell(1,length(T));
end
if nargout>=5
    St = cell(1,length(T));
end
if nargout>=6
    Ht = cell(1,length(T));
end
% dinct = cell(1,length(T)); % increment of phase field
% tol = 1e-12;

d = calc_init_dirichlet(S_phase);
u = calc_init_dirichlet(S);
if strcmpi(PFsolver,'historyfieldnode')
    H = FENODEFIELD(calc_energyint(S,u,'node','positive'));
else
    H = calc_energyint(S,u,'positive','intorder','mass');
end

if ~strcmpi(PFsolver,'historyfieldelem') && ~strcmpi(PFsolver,'historyfieldnode')
    optimFun = 'lsqlin'; % 'lsqlin' or 'lsqnonlin' or 'fmincon'
    % optimFun = 'lsqnonlin';
    % optimFun = 'fmincon';

    displayoptim = 'off';
    % displayoptim = 'iter';
    % displayoptim = 'iter-detailed';
    % displayoptim = 'final';
    % displayoptim = 'final-detailed';

    % tolX = 1e-6; % tolerance on the parameter value
    % tolFun = 1e-6; % tolerance on the function value
    % maxFunEvals = Inf; % maximum number of function evaluations

    % optimAlgo = 'interior-point';
    optimAlgo = 'trust-region-reflective';
    % optimAlgo = 'sqp';
    % optimAlgo = 'active-set';
    % optimAlgo = 'levenberg-marquardt';

    % options  = optimoptions(optimFun,'Display',displayoptim,'StepTolerance',tolX,'FunctionTolerance',tolFun,...
    %     'OptimalityTolerance',tolFun...%,'MaxFunctionEvaluations',maxFunEvals...%,'Algorithm',optimAlgo...
    %     ,'SpecifyObjectiveGradient',true...
    %     );
    % options  = optimoptions(optimFun,'Display',displayoptim,...
    %     'SpecifyObjectiveGradient',true);
    options = optimoptions(optimFun,'Display',displayoptim,'Algorithm',optimAlgo);
end

if display_
    fprintf('\n+-----------+-----------+-----------+----------+----------+------------+------------+\n');
    fprintf('|   Iter    |  u [mm]   |  f [kN]   | Nb nodes | Nb elems |  norm(d)   |  norm(u)   |\n');
    fprintf('+-----------+-----------+-----------+----------+----------+------------+------------+\n');
    fprintf('| %4d/%4d | %6.3e | %6.3e | %8d | %8d | %9.4e | %9.4e |\n',0,length(T),0,0,getnbnode(S),getnbelem(S),0,0);
end

materials_phase = MATERIALS(S_phase);
materials = MATERIALS(S);

for i=1:length(T)
    
    % Internal energy field
    switch lower(PFsolver)
        case 'historyfieldelem'
            h_old = getvalue(H);
            H = calc_energyint(S,u,'positive','intorder','mass');
            h = getvalue(H);
            for p=1:getnbgroupelem(S)
                he = double(h{p});
                he_old = double(h_old{p});
                rep = find(he <= he_old);
                he(rep) = he_old(rep);
                h{p} = MYDOUBLEND(he);
            end
            H = FEELEMFIELD(h,'storage',getstorage(H),'type',gettype(H),'ddl',getddl(H));
        case 'historyfieldnode'
            h_old = double(H);
            H = FENODEFIELD(calc_energyint(S,u,'node','positive'));
            h = double(H);
            rep = find(h <= h_old);
            h(rep) = h_old(rep);
            H = setvalue(H,h);
        otherwise
            H = calc_energyint(S,u,'positive','intorder','mass');
    end
    
    % Phase field
    mats_phase = MATERIALS(S_phase);
    for m=1:length(mats_phase)
        mat = mats_phase{m};
        if isparam(mat,'delta') && any(getparam(mat,'delta')>0) && isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
            r = getparam(mat,'r');
        else
            r = getparam(materials_phase{m},'r');
        end
        if strcmpi(PFsolver,'historyfieldnode')
            mats_phase{m} = setparam(mats_phase{m},'r',r+2*H);
        else
            mats_phase{m} = setparam(mats_phase{m},'r',r+2*H{m});
        end
    end
    S_phase = actualisematerials(S_phase,mats_phase);
    
    [A_phase,b_phase] = calc_rigi(S_phase);
    b_phase = -b_phase + bodyload(S_phase,[],'QN',2*H);
    
    % d_old = d;
    switch lower(PFsolver)
        case {'historyfieldelem','historyfieldnode'}
            d = A_phase\b_phase;
        otherwise
            d0 = freevector(S_phase,d);
            lb = d0;
            ub = ones(size(d0));
            switch optimFun
                case 'lsqlin'
                    [d,err,~,exitflag,output] = lsqlin(A_phase,b_phase,[],[],[],[],lb,ub,d0,options);
                case 'lsqnonlin'
                    fun = @(d) funlsqnonlinPF(d,A_phase,b_phase);
                    [d,err,~,exitflag,output] = lsqnonlin(fun,d0,lb,ub,options);
                case 'fmincon'
                    fun = @(d) funoptimPF(d,A_phase,b_phase);
                    [d,err,exitflag,output] = fmincon(fun,d0+eps,[],[],[],[],lb,ub,[],options);
            end
    end
    d = unfreevector(S_phase,d);
    
    % Mesh adaptation
    mats = MATERIALS(S);
    S_phase_old = S_phase;
    % S_old = S;
    cl = sizemap(d);
    S_phase = adaptmesh(S_phase_old,cl,fullfile(pathname,filename),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    S = S_phase;
    
    node_phase = getnode(S_phase);
    for m=1:length(mats_phase)
        mat = materials_phase{m};
        if isparam(mat,'delta') && any(getparam(mat,'delta')>0) && isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
            elem = getgroupelem(S_phase,m);
            nbelem = getnbelem(elem);
            xnode = node_phase(elem);
            gauss = calc_gauss(elem,'mass');
            xgauss = gauss.coord;
            x = calc_x(elem,xnode,xgauss);
            x = getcoord(NODE(POINT(x(:,:,:))));
            shinozukaPF = getparam(mat,'shinozuka');
            Xi = shinozukaPF(x); % sample for bivariate Gaussian random field with statistically independent normalized Gaussian components
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
            k = gc.*l;
            r = gc./l;
            mats_phase{m} = setparam(mats_phase{m},'k',k);
            mats_phase{m} = setparam(mats_phase{m},'r',r);
        end
        S_phase = setmaterial(S_phase,mats_phase{m},m);
    end
    S_phase = final(S_phase);
    
    % Update fields
    P_phase = calcProjection(S_phase,S_phase_old,[],'free',false,'full',true);
    d = P_phase'*d;
    % d_old = P_phase'*d_old;
    % dinc = d - d_old;
    % dincmin = min(dinc); if dincmin<-tol, dincmin, end
    
    numnodes = find(d>=1);
    if ~isempty(numnodes)
        S_phase = addcl(S_phase,numnodes,'T',1);
    end
    
    if strcmpi(PFsolver,'historyfieldnode')
        h = P_phase'*h;
        H = setvalue(H,h);
    end
    
    % P = calcProjection(S,S_old,[],'free',false,'full',true);
    P = kron(P_phase,eye(Dim));
    u = P'*u;
    
    % Displacement field
    node = getnode(S);
    for m=1:length(mats)
        mats{m} = setparam(mats{m},'d',d);
        mats{m} = setparam(mats{m},'u',u);
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
%                 rho = 0;
%                 if isparam(mat,'rcorr')
%                     rho = getparam(mat,'rcorr'); % correlation coefficient between bulk and shear moduli
%                 end
%                 C1 = gaminv(normcdf(Xi(:,1)),aC1,bC1); % sample for bulk modulus [Pa]
%                 C2 = gaminv(normcdf(rho*Xi(:,1) + sqrt(1-rho^2)*Xi(:,2)),aC2,bC2); % sample for shear modulus [Pa]
%                 C1 = reshape(C1,1,1,nbelem,gauss.nbgauss);
%                 C2 = reshape(C2,1,1,nbelem,gauss.nbgauss);
%                 C1 = MYDOUBLEND(C1);
%                 C2 = MYDOUBLEND(C2);
%                 C1 = FEELEMFIELD({C1},'storage','gauss','type','scalar','ddl',DDL('C1'));
%                 C2 = FEELEMFIELD({C2},'storage','gauss','type','scalar','ddl',DDL('C2'));
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
                mats{m} = setparam(mats{m},'E',E);
                mats{m} = setparam(mats{m},'NU',NU);
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
                mats{m} = setparam(mats{m},'C',Cmat);
            else
                error('Wrong material symmetry class');
            end
        end
        S = setmaterial(S,mats{m},m);
    end
    S = final(S);
    ud = t(i);
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
    
    if strcmpi(PFsolver,'historyfieldelem')
        H = calc_energyint(S,u,'positive','intorder','mass');
    end
    
    [A,b] = calc_rigi(S,'nofree');
    b = -b;
    
    u = freematrix(S,A)\b;
    u = unfreevector(S,u);
    
    switch lower(loading)
        case 'tension'
            numddl = findddl(S,'UY',BU);
        case 'shear'
            numddl = findddl(S,'UX',BU);
        otherwise
            error('Wrong loading case');
    end
    f = A(numddl,:)*u;
    f = sum(f);
    
    % Update fields
    dt{i} = d;
    ut{i} = u;
    ft(i) = f;
    if nargout>=4
        St_phase{i} = S_phase;
    end
    if nargout>=5
        St{i} = S;
    end
    if nargout>=6
        if strcmpi(PFsolver,'historyfieldnode')
            Ht{i} = double(H);
        else
            Ht{i} = reshape(double(mean(H,4)),[getnbelem(S),1]);
        end
    end
    % dinct{i} = dinc;
    
    if display_
        fprintf('| %4d/%4d | %6.3e | %6.3e | %8d | %8d | %9.4e | %9.4e |\n',i,length(T),t(i)*1e3,ft(i)*((Dim==2)*1e-6+(Dim==3)*1e-3),getnbnode(S),getnbelem(S),norm(dt{i}),norm(ut{i}));
    end
end

if display_
    fprintf('+-----------+-----------+-----------+----------+----------+------------+------------+\n');
end

end
