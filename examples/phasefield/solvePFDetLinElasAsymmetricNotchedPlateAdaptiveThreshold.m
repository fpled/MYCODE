function [dt,ut,ft,T,St_phase,St,Ht] = solvePFDetLinElasAsymmetricNotchedPlateAdaptiveThreshold(S_phase,S,T,PFsolver,C,BU,BL,BR,H1,H2,H3,PU,PL,PR,sizemap,varargin)
% function [dt,ut,ft,T,St_phase,St,Ht] = solvePFDetLinElasAsymmetricNotchedPlateAdaptiveThreshold(S_phase,S,T,PFsolver,C,BU,BL,BR,H1,H2,H3,PU,PL,PR,sizemap,varargin)
% Solve deterministic Phase Field problem with mesh adaptation.

display_ = getcharin('display',varargin,true);
displayIter = getcharin('displayiter',varargin,false);
tolConv = getcharin('tol',varargin,1e-2);
maxIter = getcharin('maxiter',varargin,100);
filename = getcharin('filename',varargin,'gmsh_domain_asymmetric_notched_plate');
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
    fprintf('\n+------+---------+-----------+-----------+----------+----------+------------+------------+\n');
    fprintf('| Iter | Nb iter |  u [mm]   |  f [kN]   | Nb nodes | Nb elems |  norm(d)   |  norm(u)   |\n');
    fprintf('+------+---------+-----------+-----------+----------+----------+------------+------------+\n');
    fprintf('| %4d | %7d | %6.3e | %6.3e | %8d | %8d | %9.4e | %9.4e |\n',0,0,0,0,getnbnode(S),getnbelem(S),0,0);
end

i = 0;
ti = 0;
dti = dt0;
while ti < tf
    i = i+1;
    
    switch lower(PFsolver)
        case 'historyfieldelem'
            h_old = getvalue(H);
        case 'historyfieldnode'
            h_old = double(H);
    end
    d_old = d;
    errConv = Inf;
    nbIter = 0;
    while (errConv > tolConv) && (nbIter < maxIter)
        nbIter = nbIter+1;
        
        % Internal energy field
        switch lower(PFsolver)
            case 'historyfieldelem'
                H = calc_energyint(S,u,'positive','intorder','mass');
                h = getvalue(H);
                for p=1:getnbgroupelem(S)
                    he = double(h{p});
                    he_old = double(h_old{p});
                    rep = find(he <= he_old);
                    he(rep) = he_old(rep);
                    h{p} = MYDOUBLEND(he);
                end
                H = setvalue(H,h);
            case 'historyfieldnode'
                H = FENODEFIELD(calc_energyint(S,u,'node','positive'));
                h = double(H);
                rep = find(h <= h_old);
                h(rep) = h_old(rep);
                H = setvalue(H,h);
            otherwise
                H = calc_energyint(S,u,'positive','intorder','mass');
        end
        
        % Phase field
        if strcmpi(PFsolver,'historyfieldnode')
            r = FENODEFIELD(calc_parammat(S_phase,'r','node'));
            qn = FENODEFIELD(calc_parammat(S_phase,'qn','node'));
        else
            r = calc_parammat(S_phase,'r');
            qn = calc_parammat(S_phase,'qn');
        end
        mats_phase = MATERIALS(S_phase);
        R = r+2*H;
        for m=1:length(mats_phase)
            if strcmpi(PFsolver,'historyfieldnode')
                mats_phase{m} = setparam(mats_phase{m},'r',R);
            else
                mats_phase{m} = setparam(mats_phase{m},'r',R{m});
            end
        end
        S_phase = actualisematerials(S_phase,mats_phase);
        
        [A_phase,b_phase] = calc_rigi(S_phase);
        Q = 2*H+qn;
        if strcmpi(PFsolver,'historyfieldnode')
            q = double(Q);
            q = max(q,0);
            Q = setvalue(Q,q);
        else
            q = getvalue(Q);
            for p=1:getnbgroupelem(S)
                qe = double(q{p});
                qe = max(qe,0);
                q{p} = MYDOUBLEND(qe);
            end
            Q = setvalue(Q,q);
        end
        b_phase = -b_phase + bodyload(S_phase,[],'QN',Q);
        
        d_prev = d;
        switch lower(PFsolver)
            case {'historyfieldelem','historyfieldnode'}
                d = A_phase\b_phase;
            otherwise
                d0 = freevector(S_phase,d_old);
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
        d = unfreevector(S_phase,d);
        
        % Displacement field
        mats = MATERIALS(S);
        for m=1:length(mats)
            mats{m} = setparam(mats{m},'d',d);
            mats{m} = setparam(mats{m},'u',u);
        end
        S = actualisematerials(S,mats);
        if nbIter==1
            S = removebc(S);
            ti = ti + dti;
            ud = -ti;
            S = addcl(S,PU,'UY',ud);
            S = addcl(S,PL,{'UX','UY'});
            S = addcl(S,PR,'UY');
        end
        
        [A,b] = calc_rigi(S,'nofree');
        b = -b;
        
        u = freematrix(S,A)\b;
        u = unfreevector(S,u);
        
        numddl = findddl(S,'UY',PU);
        f = -A(numddl,:)*u;
        f = sum(f);
        
        errConv = norm(d-d_prev,'Inf');
        if displayIter
            fprintf('sub-iter #%2.d : error = %.3e\n',nbIter,errConv);
        end
    end

    % Update fields
    dt{i} = d;
    ut{i} = u;
    ft(i) = f;
    t(i) = ti;
    if nargout>=5
        St_phase{i} = S_phase;
    end
    if nargout>=6
        St{i} = S;
    end
    if nargout>=7
        if strcmpi(PFsolver,'historyfieldnode')
            Ht{i} = double(H);
        else
            Ht{i} = reshape(double(mean(H,4)),[getnbelem(S),1]);
        end
    end
    
    if display_
        fprintf('| %4d | %7d | %6.3e | %6.3e | %8d | %8d | %9.4e | %9.4e |\n',i,nbIter,t(i)*1e3,f*((Dim==2)*1e-6+(Dim==3)*1e-3),getnbnode(S),getnbelem(S),norm(d),norm(u));
    end
    
    if ti < tf
        % Mesh adaptation
        S_phase_old = S_phase;
        S_phase_ref = addcl(S_phase_old,C,'T',1);
        S_phase_ref = addcl(S_phase_ref,H1,'T',1);
        S_phase_ref = addcl(S_phase_ref,H2,'T',1);
        S_phase_ref = addcl(S_phase_ref,H3,'T',1);
        d_ref = freevector(S_phase_ref,d);
        d_ref = unfreevector(S_phase_ref,d_ref);
        % S_old = S;
        cl = sizemap(d_ref);
        S_phase = adaptmesh(S_phase_ref,cl,fullfile(pathname,filename),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
        S = S_phase;
        
        % Update phase field properties
        S_phase = setphasefieldproperties(S_phase,materials_phase);
        S_phase = final(S_phase,'duplicate');
        S_phase = addcl(S_phase,BU,'T');
        S_phase = addcl(S_phase,BL,'T');
        S_phase = addcl(S_phase,BR,'T');
        
        % Update material properties
        S = setmaterialproperties(S,materials);
        S = final(S,'duplicate');
        S = addcl(S,PU,'UY',ud);
        S = addcl(S,PL,{'UX','UY'});
        S = addcl(S,PR,'UY');
        
        % Update fields
        P_phase = calcProjection(S_phase,S_phase_old,[],'free',false,'full',true);
        d = P_phase'*d;
        
        % P = calcProjection(S,S_old,[],'free',false,'full',true);
        P = kron(P_phase,eye(Dim));
        u = P'*u;
        
        if strcmpi(PFsolver,'historyfieldnode')
            h = P_phase'*h;
            H = setvalue(H,h);
        elseif strcmpi(PFsolver,'historyfieldelem')
            H = calc_energyint(S,u,'positive','intorder','mass');
        end
    end
end

if display_
    fprintf('+------+---------+-----------+-----------+----------+----------+------------+------------+\n');
end

T = TIMEMODEL(t);

end
