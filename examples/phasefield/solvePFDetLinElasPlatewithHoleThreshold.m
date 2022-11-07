function [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElasPlatewithHoleThreshold(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P0,varargin)
% function [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElasPlatewithHoleThreshold(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P0,varargin)
% Solve deterministic Phase Field problem.

display_ = getcharin('display',varargin,true);
displayIter = getcharin('displayiter',varargin,false);
tolConv = getcharin('tol',varargin,1e-2);
maxIter = getcharin('maxiter',varargin,100);
dbthreshold = getcharin('damageboundarythreshold',varargin,0.99);

Dim = getdim(S);

dt0 = T.dt0;
dt1 = T.dt1;
tf = T.tf;
dthreshold = T.dthreshold;

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
    fprintf('\n+------+---------+-----------+-----------+-----------+-----------+-----------+');
    fprintf('\n| Iter | Nb iter |  u [mm]   |  f [kN]   |  max(d)   |  Ed [J]   |  Eu [J]   |');
    fprintf('\n+------+---------+-----------+-----------+-----------+-----------+-----------+\n');
end

numddlbr = findddl(S_phase,'T',BRight);
numddlbl = findddl(S_phase,'T',BLeft);
dbr = d(numddlbr,:);
dbl = d(numddlbl,:);

i = 0;
ti = 0;
dti = dt0;
while ti < tf
    i = i+1;
    tIter = tic;
    nbIter = 0;
    if any(dbr > dbthreshold) && any(dbl > dbthreshold)
        ti = ti + dti;
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
                    d0 = freevector(S_phase,d_old);
                    lb = d0;
                    lb(lb==1) = 1-eps;
                    ub = ones(size(d0));
                    switch optimFun
                        case 'lsqlin'
                            d = lsqlin(A_phase,b_phase,[],[],[],[],lb,ub,d0,options);
                        case 'lsqnonlin'
                            fun = @(d) funlsqnonlinPF(d,A_phase,b_phase);
                            d = lsqnonlin(fun,d0,lb,ub,options);
                        case 'fmincon'
                            fun = @(d) funoptimPF(d,A_phase,b_phase);
                            d = fmincon(fun,d0+eps,[],[],[],[],lb,ub,[],options);
                    end
            end
            if any(d > dthreshold)
                dti = dt1;
            end
            d = unfreevector(S_phase,d);
            dbr = d(numddlbr,:);
            dbl = d(numddlbl,:);
            
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
                S = addcl(S,BU,'UY',ud);
                S = addcl(S,BL,'UY');
                S = addcl(S,P0,'UX');
            end
            
            [A,b] = calc_rigi(S,'nofree');
            b = -b;
            
            u = freematrix(S,A)\b;
            u = unfreevector(S,u);
            
            errConv = norm(d-d_prev,'Inf');
            if displayIter
                fprintf('sub-iter #%2.d : error = %.3e\n',nbIter,errConv);
            end
        end
        
        % Force
        numddl = findddl(S,'UY',BU);
        f = -A(numddl,:)*u;
        f = sum(f);
        
        % Energies
        for m=1:length(mats_phase)
            if strcmpi(PFsolver,'historyfieldnode')
                mats_phase{m} = setparam(mats_phase{m},'r',r);
            else
                mats_phase{m} = setparam(mats_phase{m},'r',r{m});
            end
        end
        S_phase = actualisematerials(S_phase,mats_phase);
        
        A_phase = calc_rigi(S_phase,'nofree');
        b_phase = bodyload(S_phase,[],'QN',qn,'nofree');
        
        Ed = 1/2*d'*A_phase*d - d'*b_phase;
        Eu = 1/2*u'*A*u;
    end
    
    % Update fields
    dt{i} = d;
    ut{i} = u;
    ft(i) = f;
    t(i) = ti;
    if nargout>=4
        if strcmpi(PFsolver,'historyfieldnode')
            Ht{i} = double(H);
        else
            Ht{i} = reshape(double(mean(H,4)),[getnbelem(S),1]);
        end
    end
    if nargout>=5
        Edt(i) = Ed;
    end
    if nargout>=6
        Eut(i) = Eu;
    end
    if nargout>=7
        iteration(i) = nbIter;
        time(i) = toc(tIter);
        err(i) = errConv;
    end
    
    if display_
        fprintf('| %4d | %7d | %9.3e | %9.3e | %9.3e | %9.3e | %9.3e |\n',i,nbIter,t(i)*1e3,f*((Dim==2)*1e-6+(Dim==3)*1e-3),max(d),Ed,Eu);
    end
end

if display_
    fprintf('+------+---------+-----------+-----------+-----------+-----------+-----------+\n');
end

T = TIMEMODEL(t);
dt = TIMEMATRIX(dt,T,size(d));
ut = TIMEMATRIX(ut,T,size(u));
if nargout>=4
    if strcmpi(PFsolver,'historyfieldnode')
        Ht = TIMEMATRIX(Ht,T,size(d));
    else
        Ht = TIMEMATRIX(Ht,T,[getnbelem(S),1]);
    end
end
if nargout>=7
    output.iteration = iteration;
    output.time = time;
    output.error = err;
end

end
