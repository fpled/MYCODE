function [dt,ut,ft,St_phase,St,Ht,Edt,Eut,output] = solvePFDetLinElasPlatewithHoleAdaptive(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P0,C,sizemap,varargin)
% function [dt,ut,ft,St_phase,St,Ht,Edt,Eut,output] = solvePFDetLinElasPlatewithHoleAdaptive(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,P0,C,sizemap,varargin)
% Solve deterministic Phase Field problem with mesh adaptation.

display_ = getcharin('display',varargin,true);
displayIter = getcharin('displayiter',varargin,false);
tolConv = getcharin('tol',varargin,1e-2);
maxIter = getcharin('maxiter',varargin,100);
dbthreshold = getcharin('damageboundarythreshold',varargin,0.99);
filename = getcharin('filename',varargin,'gmsh_domain_with_hole');
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
if nargout>=7
    Edt = zeros(1,length(T));
end
if nargout>=8
    Eut = zeros(1,length(T));
end
if nargout>=9
    iteration = zeros(1,length(T));
    time = zeros(1,length(T));
    err = zeros(1,length(T));
end

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
    fprintf('\n+-----------+---------+-----------+-----------+-----------+-----------+-----------+----------+----------+');
    fprintf('\n|   Iter    | Nb iter |  u [mm]   |  f [kN]   |  max(d)   |  Ed [J]   |  Eu [J]   | Nb nodes | Nb elems |');
    fprintf('\n+-----------+---------+-----------+-----------+-----------+-----------+-----------+----------+----------+');
    fprintf('\n| %4d/%4d | %4d | %7d | %9.3e | %9.3e | %9.3e | %9.3e | %9.3e | %8d | %8d |\n',0,length(T),0,0,0,0,0,0,getnbnode(S),getnbelem(S));
end

numddlbr = findddl(S_phase,'T',BRight);
numddlbl = findddl(S_phase,'T',BLeft);
dbr = d(numddlbr,:);
dbl = d(numddlbl,:);

for i=1:length(T)
    tIter = tic;
    nbIter = 0;
    if any(dbr > dbthreshold) && any(dbl > dbthreshold)
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
                            d = lsqlin(A_phase,b_phase,[],[],[],[],lb,ub,d0,options);
                        case 'lsqnonlin'
                            fun = @(d) funlsqnonlinPF(d,A_phase,b_phase);
                            d = lsqnonlin(fun,d0,lb,ub,options);
                        case 'fmincon'
                            fun = @(d) funoptimPF(d,A_phase,b_phase);
                            d = fmincon(fun,d0+eps,[],[],[],[],lb,ub,[],options);
                    end
            end
            d = unfreevector(S_phase,d);
            numddlbr = findddl(S_phase,'T',BRight);
            numddlbl = findddl(S_phase,'T',BLeft);
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
                ud = -t(i);
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
            if any(dbr > dbthreshold) && any(dbl > dbthreshold)
                break
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
    if nargout>=7
        Edt(i) = Ed;
    end
    if nargout>=8
        Eut(i) = Eu;
    end
    if nargout>=9
        iteration(i) = nbIter;
        time(i) = toc(tIter);
        err(i) = errConv;
    end
    
    if display_
        fprintf('| %4d/%4d | %4d | %7d | %9.3e | %9.3e | %9.3e | %9.3e | %9.3e | %8d | %8d |\n',i,length(T),nbIter,t(i)*1e3,f*((Dim==2)*1e-6+(Dim==3)*1e-3),max(d),Ed,Eu,getnbnode(S),getnbelem(S));
    end
    
    if i<length(T) && ~(any(dbr > dbthreshold) && any(dbl > dbthreshold))
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
        S_phase = final(S_phase);
        
        % Update material properties
        S = setmaterialproperties(S,materials);
        S = final(S);
        S = addcl(S,BU,'UY',ud);
        S = addcl(S,BL,'UY');
        S = addcl(S,P0,'UX');
        
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
    fprintf('+-----------+---------+-----------+-----------+----------+----------+------------+------------+\n');
end

if nargout>=9
    output.iteration = iteration;
    output.time = time;
    output.error = err;
end

end
