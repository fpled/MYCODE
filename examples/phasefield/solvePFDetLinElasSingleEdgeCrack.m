function [dt,ut,ft,Ht] = solvePFDetLinElasSingleEdgeCrack(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,varargin)
% function [dt,ut,ft,Ht] = solvePFDetLinElasSingleEdgeCrack(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,varargin)
% Solve deterministic Phase Field problem.

display_ = getcharin('display',varargin,true);
displayIter = getcharin('displayiter',varargin,false);
tolConv = getcharin('tol',varargin,1e-2);
maxIter = getcharin('maxiter',varargin,100);
dbthreshold = getcharin('damageboundarythreshold',varargin,0.99);

Dim = getdim(S);

t = gett(T);

dt = cell(1,length(T));
ut = cell(1,length(T));
ft = zeros(1,length(T));
if nargout>=4
    Ht = cell(1,length(T));
end

d = calc_init_dirichlet(S_phase);
u = calc_init_dirichlet(S);
if strcmpi(PFsolver,'historyfieldnode')
    H = FENODEFIELD(calc_energyint(S,u,'node','positive'));
    r = FENODEFIELD(calc_parammat(S_phase,'r','node'));
    qn = FENODEFIELD(calc_parammat(S_phase,'qn','node'));
else
    H = calc_energyint(S,u,'positive','intorder','mass');
    r = calc_parammat(S_phase,'r');
    qn = calc_parammat(S_phase,'qn');
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
    fprintf('\n+-----------+---------+-----------+-----------+------------+------------+\n');
    fprintf('|   Iter    | Nb iter |  u [mm]   |  f [kN]   |  norm(d)   |  norm(u)   |\n');
    fprintf('+-----------+---------+-----------+-----------+------------+------------+\n');
end

numddlb = findddl(S_phase,'T',BRight);
db = d(numddlb,:);

for i=1:length(T)
    
    nbIter = 0;
    if any(db > dbthreshold)
        f = 0;
    else
        switch lower(PFsolver)
            case 'historyfieldelem'
                h_old = getvalue(H);
            case 'historyfieldnode'
                h_old = double(H);
        end
        d_old = d;
        errConv = Inf;
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
                    lb = freevector(S_phase,d_old);
                    lb(lb==1) = 1-eps;
                    ub = ones(size(lb));
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
            db = d(numddlb,:);
            
            % Displacement field
            mats = MATERIALS(S);
            for m=1:length(mats)
                mats{m} = setparam(mats{m},'d',d);
                mats{m} = setparam(mats{m},'u',u);
            end
            S = actualisematerials(S,mats);
            if nbIter==1
                S = removebc(S);
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
            
            errConv = norm(d-d_prev,'Inf');
            if displayIter
                fprintf('sub-iter #%2.d : error = %.3e\n',nbIter,errConv);
            end
        end
    end
    
    % Update fields
    dt{i} = d;
    ut{i} = u;
    ft(i) = f;
    if nargout>=4
        if strcmpi(PFsolver,'historyfieldnode')
            Ht{i} = double(H);
        else
            Ht{i} = reshape(double(mean(H,4)),[getnbelem(S),1]);
        end
    end
    
    if display_
        fprintf('| %4d/%4d | %7d | %6.3e | %6.3e | %9.4e | %9.4e |\n',i,length(T),nbIter,t(i)*1e3,f*((Dim==2)*1e-6+(Dim==3)*1e-3),norm(d),norm(u));
    end
end

if display_
    fprintf('+-----------+---------+-----------+-----------+------------+------------+\n');
end

dt = TIMEMATRIX(dt,T,size(d));
ut = TIMEMATRIX(ut,T,size(u));
if nargout>=4
    if strcmpi(PFsolver,'historyfieldnode')
        Ht = TIMEMATRIX(Ht,T,size(d));
    else
        Ht = TIMEMATRIX(Ht,T,[getnbelem(S),1]);
    end
end

end
