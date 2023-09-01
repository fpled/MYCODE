function [dt,ut,ft,T,St_phase,St,Ht,Edt,Eut,output] = solvePFDetLinElasAsymmetricNotchedPlateAdaptiveThreshold(S_phase,S,T,PFsolver,C,BU,BL,BR,H1,H2,H3,PU,PL,PR,initialCrack,sizemap,varargin)
% function [dt,ut,ft,T,St_phase,St,Ht,Edt,Eut,output] = solvePFDetLinElasAsymmetricNotchedPlateAdaptiveThreshold(S_phase,S,T,PFsolver,C,BU,BL,BR,H1,H2,H3,PU,PL,PR,initialCrack,sizemap,varargin)
% Solve deterministic Phase Field problem with mesh adaptation.

display_ = getcharin('display',varargin,true);
displayIter = getcharin('displayiter',varargin,false);
maxIter = getcharin('maxiter',varargin,100);
tolConv = getcharin('tol',varargin,1e-3);
critConv = getcharin('crit',varargin,'Energy');
filename = getcharin('filename',varargin,'gmsh_asymmetric_notched_plate');
pathname = getcharin('pathname',varargin,'.');
gmshoptions = getcharin('gmshoptions',varargin,'-v 0');
mmgoptions = getcharin('mmgoptions',varargin,'-nomove -v -1');

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    checkConvSol = ~isempty(strfind(lower(critConv),'solution'));
    checkConvRes = ~isempty(strfind(lower(critConv),'residual'));
    checkConvEnergy = ~isempty(strfind(lower(critConv),'energy'));
else
    checkConvSol = contains(critConv,'solution','IgnoreCase',true);
    checkConvRes = contains(critConv,'residual','IgnoreCase',true);
    checkConvEnergy = contains(critConv,'energy','IgnoreCase',true);
end

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
    qn = FENODEFIELD(calc_parammat(S_phase,'qn','node'));
else
    H = calc_energyint(S,u,'positive','intorder','mass');
    qn = calc_parammat(S_phase,'qn');
end
if checkConvEnergy
    Ae_phase = calc_rigi(S_phase,'nofree');
    be_phase = bodyload(S_phase,[],'QN',qn,'nofree');
    Ed = 1/2*d'*Ae_phase*d - d'*be_phase;
    A = calc_rigi(S,'nofree');
    Eu = 1/2*u'*A*u;
    E = Ed+Eu;
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
    fprintf('\n+------+---------+-----------+-----------+-----------+-----------+-----------+----------+----------+');
    fprintf('\n| Iter | Nb iter |  u [mm]   |  f [kN]   |  max(d)   |  Ed [J]   |  Eu [J]   | Nb nodes | Nb elems |');
    fprintf('\n+------+---------+-----------+-----------+-----------+-----------+-----------+----------+----------+');
    fprintf('\n| %4d | %7d | %9.3e | %9.3e | %9.3e | %9.3e | %9.3e | %8d | %8d |\n',0,0,0,0,0,0,0,getnbnode(S),getnbelem(S));
end

i = 0;
ti = 0;
dti = dt0;
while ti < tf-eps
    i = i+1;
    tIter = tic;
    nbIter = 0;
    if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
        H_old = H;
    end
    d_old = d;
    if strcmpi(PFsolver,'historyfieldnode')
        r = FENODEFIELD(calc_parammat(S_phase,'r','node'));
        qn = FENODEFIELD(calc_parammat(S_phase,'qn','node'));
    else
        r = calc_parammat(S_phase,'r');
        qn = calc_parammat(S_phase,'qn');
    end
    Ae_phase = calc_rigi(S_phase,'nofree');
    be_phase = bodyload(S_phase,[],'QN',qn,'nofree');
    if checkConvRes
        [S_phase,A_phase,b_phase] = calcphasefieldoperator(S_phase,r,qn,H);
    end
    errConv = Inf;
    errConvs = 0; errConvr = 0; errConve = 0;
    while (errConv > tolConv) && (nbIter < maxIter)
        nbIter = nbIter+1;
        if checkConvSol
            d_prev = d;
            u_prev = u;
        end
        if checkConvEnergy
            E_prev = E;
        end
        
        % Phase field
        if ~checkConvRes
            [S_phase,A_phase,b_phase] = calcphasefieldoperator(S_phase,r,qn,H);
        end
        
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
        dmax = max(d);
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
        
        % Internal energy field
        switch lower(PFsolver)
            case {'historyfieldelem','historyfieldnode'}
                H = calc_historyfield(S,u,H_old);
            otherwise
                H = calc_energyint(S,u,'positive','intorder','mass');
        end
        
        % Convergence
        if checkConvSol
            errConvd = norm(d-d_prev)/norm(d);
            errConvu = norm(u-u_prev)/norm(u);
            errConvs = max(errConvd,errConvu);
        end
        if checkConvRes
            % Phase field residual
            [S_phase,A_phase,b_phase] = calcphasefieldoperator(S_phase,r,qn,H);
            r_phase = A_phase*d - b_phase;
            errConvr = norm(r_phase)/norm(b_phase);
        end
        if checkConvEnergy
            % Energy
            Ed = 1/2*d'*Ae_phase*d - d'*be_phase;
            Eu = 1/2*u'*A*u;
            E = Ed+Eu;
            errConve = abs(E-E_prev)/abs(E);
        end
        if checkConvSol || checkConvRes || checkConvEnergy
            errConv = max(max(errConvs,errConvr),errConve);
            if displayIter
                fprintf('sub-iter #%2.d : error = %.3e',nbIter,errConv);
                if checkConvSol
                    fprintf(', err_s = %.3e',errConvs);
                end
                if checkConvRes
                    fprintf(', err_r = %.3e',errConvr);
                end
                if checkConvEnergy
                    fprintf(', err_e = %.3e',errConve);
                end
                fprintf('\n');
            end
        end
    end
    
    % Force
    numddl = findddl(S,'UY',PU);
    f = -A(numddl,:)*u;
    f = sum(f);
    
    % Energy
    if ~checkConvEnergy
        Ed = 1/2*d'*Ae_phase*d - d'*be_phase;
        Eu = 1/2*u'*A*u;
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
    if nargout>=8
        Edt(i) = Ed;
    end
    if nargout>=9
        Eut(i) = Eu;
    end
    if nargout>=10
        iteration(i) = nbIter;
        time(i) = toc(tIter);
        err(i) = errConv;
    end
    
    if display_
        fprintf('| %4d | %7d | %9.3e | %9.3e | %9.3e | %9.3e | %9.3e | %8d | %8d |\n',i,nbIter,t(i)*1e3,f*((Dim==2)*1e-6+(Dim==3)*1e-3),dmax,Ed,Eu,getnbnode(S),getnbelem(S));
    end
    
    if ti < tf-eps
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
        S_phase = final(S_phase);
        if strcmpi(initialCrack,'initialphasefield')
            S_phase = addcl(S_phase,C,'T',1);
        end
        S_phase = addcl(S_phase,BU,'T');
        S_phase = addcl(S_phase,BL,'T');
        S_phase = addcl(S_phase,BR,'T');
        
        % Update material properties
        S = setmaterialproperties(S,materials);
        S = final(S);
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
        else
            H = calc_energyint(S,u,'positive','intorder','mass');
        end
    end
end

if display_
    fprintf('+------+---------+-----------+-----------+-----------+-----------+-----------+----------+----------+\n');
end

T = TIMEMODEL(t);
if nargout>=10
    output.iteration = iteration;
    output.time = time;
    output.error = err;
end

end
