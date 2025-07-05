function [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElasLshapedPanel(S_phase,S,T,PFsolver,BL,BRight,BLeft,varargin)
% function [dt,ut,ft,Ht,Edt,Eut,output] = solvePFDetLinElasLshapedPanel(S_phase,S,T,PFsolver,BL,BRight,BLeft,varargin)
% Solve deterministic phase-field problem.

display_ = getcharin('display',varargin,true);
displayIter = getcharin('displayiter',varargin,false);
maxIter = getcharin('maxiter',varargin,100);
tolConv = getcharin('tol',varargin,1e-2);
critConv = getcharin('crit',varargin,'Energy');
dbthreshold = getcharin('damageboundarythreshold',varargin,0.999);

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    contain = @(str,pat) ~isempty(strfind(lower(str),pat));
else
    contain = @(str,pat) contains(str,pat,'IgnoreCase',true);
end
checkConvSol = contain(critConv,'solution');
checkConvRes = contain(critConv,'residual');
checkConvEnergy = contain(critConv,'energy');

Dim = getdim(S);

t = gett(T);

dt = cell(1,length(T));
ut = cell(1,length(T));
ft = zeros(1,length(T));
if nargout>=4
    Ht = cell(1,length(T));
end
if nargout>=5
    Edt = zeros(1,length(T));
end
if nargout>=6
    Eut = zeros(1,length(T));
end
if nargout>=7
    iteration = zeros(1,length(T));
    time = zeros(1,length(T));
    err = zeros(1,length(T));
end

d = calc_init_dirichlet(S_phase);
u = calc_init_dirichlet(S);
if strcmpi(PFsolver,'historyfieldnode')
    H = FENODEFIELD(calc_energyint(S,u,'node','positive','local'));
    r = FENODEFIELD(calc_parammat(S_phase,'r','node'));
    qn = FENODEFIELD(calc_parammat(S_phase,'qn','node'));
else
    H = calc_energyint(S,u,'intorder','mass','positive','local');
    r = calc_parammat(S_phase,'r');
    qn = calc_parammat(S_phase,'qn');
end
Ae_phase = calc_rigi(S_phase,'nofree');
be_phase = bodyload(S_phase,[],'QN',qn,'nofree');
if checkConvEnergy
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
    
    tolX = 100*eps; % tolerance on the parameter value
    tolFun = 100*eps; % tolerance on the function value
    maxFunEvals = Inf; % maximum number of function evaluations
    
    % optimAlgo = 'interior-point';
    optimAlgo = 'trust-region-reflective';
    % optimAlgo = 'sqp';
    % optimAlgo = 'active-set';
    % optimAlgo = 'levenberg-marquardt';
    
    switch optimFun
        case 'lsqlin'
            options = optimoptions(optimFun,'Display',displayoptim,'Algorithm',optimAlgo);
        case 'lsqnonlin'
            options  = optimoptions(optimFun,'Display',displayoptim,'Algorithm',optimAlgo,...
                'StepTolerance',tolX,'FunctionTolerance',tolFun,'OptimalityTolerance',tolFun,...
                'SpecifyObjectiveGradient',true);
        case 'fmincon'
            options  = optimoptions(optimFun,'Display',displayoptim,'Algorithm',optimAlgo,...
                'StepTolerance',tolX,'FunctionTolerance',tolFun,'OptimalityTolerance',tolFun,...
                'MaxFunctionEvaluations',maxFunEvals,...
                'SpecifyObjectiveGradient',true,'HessianFcn','objective');
    end
end

if display_
    fprintf('\n+-----------+---------+-----------+-----------+-----------+-----------+-----------+');
    fprintf('\n|   Iter    | Nb iter |  u [mm]   |  f [kN]   |  max(d)   |  Ed [J]   |  Eu [J]   |');
    fprintf('\n+-----------+---------+-----------+-----------+-----------+-----------+-----------+\n');
end

ismonotonic = ~any(diff(sign(t(t~=0))));
numddlb = findddl(S_phase,'T',BLeft);
db = d(numddlb,:);

for i=1:length(T)
    tIter = tic;
    nbIter = 0;
    if any(db > dbthreshold)
        f = 0;
    else
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            H_old = H;
        end
        d_old = d;
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
            dmax = max(d);
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
                % [Gerasimov, De Lorenzis, 2019, CMAME]
                % if Dim==2
                %     S = addcl(S,BRight,{'UX','UY'},[0;ud]);
                % elseif Dim==3
                %     S = addcl(S,BRight,{'UX','UY','UZ'},[0;ud;0]);
                % end
                S = addcl(S,BRight,'UY',ud);
                S = addcl(S,BL);
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
                    H = calc_energyint(S,u,'intorder','mass','positive','local');
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
            if any(db > dbthreshold)
                break
            end
        end
        
        % Force
        numddl = findddl(S,'UY',BRight);
        f = A(numddl,:)*u;
        f = sum(f);
        if ismonotonic
            f = abs(f);
        end
        
        % Energy
        if ~checkConvEnergy
            Ed = 1/2*d'*Ae_phase*d - d'*be_phase;
            Eu = 1/2*u'*A*u;
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
        fprintf('| %4d/%4d | %7d | %9.3e | %9.3e | %9.3e | %9.3e | %9.3e |\n',i,length(T),nbIter,t(i)*1e3,f*1e-3,dmax,Ed,Eu);
    end
end

if display_
    fprintf('+-----------+---------+-----------+-----------+-----------+-----------+-----------+\n');
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
if nargout>=7
    output.iteration = iteration;
    output.time = time;
    output.error = err;
end

end
