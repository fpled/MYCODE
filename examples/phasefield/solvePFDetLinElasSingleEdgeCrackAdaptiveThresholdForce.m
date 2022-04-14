function [ft,T] = solvePFDetLinElasSingleEdgeCrackAdaptiveThresholdForce(S_phase,S,T,PFsolver,C,BU,BL,BRight,BLeft,BFront,BBack,loading,sizemap,varargin)
% function [ft,T] = solvePFDetLinElasSingleEdgeCrackAdaptiveThresholdForce(S_phase,S,T,PFsolver,C,BU,BL,BRight,BLeft,BFront,BBack,loading,sizemap,varargin)
% Solve deterministic Phase Field problem with mesh adaptation.

display_ = ischarin('display',varargin);
filename = getcharin('filename',varargin,'gmsh_domain_single_edge_crack');
pathname = getcharin('pathname',varargin,'.');
gmshoptions = getcharin('gmshoptions',varargin,'-v 0');
mmgoptions = getcharin('mmgoptions',varargin,'-nomove -v -1');

Dim = getdim(S);

dt0 = T.dt0;
dt1 = T.dt1;
tf = T.tf;
dthreshold = T.dthreshold;

materials_phase = MATERIALS(S_phase);
materials = MATERIALS(S);
S_phase = setphasefieldproperties(S_phase,materials_phase);
S = setmaterialproperties(S,materials);

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
    fprintf('\n+------+-----------+-----------+----------+----------+------------+------------+\n');
    fprintf('| Iter |  u [mm]   |  f [kN]   | Nb nodes | Nb elems |  norm(d)   |  norm(u)   |\n');
    fprintf('+------+-----------+-----------+----------+----------+------------+------------+\n');
    fprintf('| %4d | %6.3e | %6.3e | %8d | %8d | %9.4e | %9.4e |\n',0,0,0,getnbnode(S),getnbelem(S),0,0);
end

dbthreshold = 0.99;
numddlb = findddl(S_phase,'T',BRight);
db = d(numddlb,:);

i = 0;
ti = 0;
dti = dt0;
while ti < tf
    i = i+1;
    
    if any(db > dbthreshold)
        ti = ti + dti;
        f = 0;
    else
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
                lb(lb==1) = 1-eps;
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
        if any(d > dthreshold)
            dti = dt1;
        end
        ti = ti + dti;
        
        d = unfreevector(S_phase,d);
        numddlb = findddl(S_phase,'T',BRight);
        db = d(numddlb,:);
        
        % Displacement field
        mats = MATERIALS(S);
        for m=1:length(mats)
            mats{m} = setparam(mats{m},'d',d);
            mats{m} = setparam(mats{m},'u',u);
        end
        S = actualisematerials(S,mats);
        S = removebc(S);
        ud = ti;
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
    end
    
    % Update fields
    t(i) = ti;
    ft(i) = f;
    
    if display_
        fprintf('| %4d | %6.3e | %6.3e | %8d | %8d | %9.4e | %9.4e |\n',i,t(i)*1e3,f*((Dim==2)*1e-6+(Dim==3)*1e-3),getnbnode(S),getnbelem(S),norm(d),norm(u));
    end
    
    if ti < tf && ~any(db > dbthreshold)
        % Mesh adaptation
        S_phase_old = S_phase;
        S_phase_ref = addcl(S_phase_old,C,'T',1);
        d_ref = freevector(S_phase_ref,d);
        d_ref = unfreevector(S_phase_ref,d_ref);
        % S_old = S;
        cl = sizemap(d_ref);
        S_phase = adaptmesh(S_phase_ref,cl,fullfile(pathname,filename),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
        S = S_phase;
        
        % Update phase field properties
        S_phase = setphasefieldproperties(S_phase,materials_phase);
        S_phase = final(S_phase,'duplicate');
        
        % Update material properties
        S = setmaterialproperties(S,materials);
        S = final(S,'duplicate');
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
        
        % Update fields
        P_phase = calcProjection(S_phase,S_phase_old,[],'free',false,'full',true);
        d = P_phase'*d;
        
        if strcmpi(PFsolver,'historyfieldnode')
            h = P_phase'*h;
            H = setvalue(H,h);
        elseif strcmpi(PFsolver,'historyfieldelem')
            H = calc_energyint(S,u,'positive','intorder','mass');
        end
        
        % P = calcProjection(S,S_old,[],'free',false,'full',true);
        P = kron(P_phase,eye(Dim));
        u = P'*u;
    end
end

if display_
    fprintf('+------+-----------+-----------+----------+----------+------------+------------+\n');
end

T = TIMEMODEL(t);

end
