function [ft,dmaxt] = solvePFDetLinElasSingleEdgeCrackForce(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,varargin)
% function [ft,dmaxt] = solvePFDetLinElasSingleEdgeCrackForce(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,varargin)
% Solve deterministic phase-field problem.

display_ = getcharin('display',varargin,true);
displayIter = getcharin('displayiter',varargin,false);
displaySol = getcharin('displaysol',varargin,false);
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

ft = zeros(1,length(T));
dmaxt = zeros(1,length(T));

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
    fprintf('\n+-----------+---------+-----------+-----------+-----------+');
    fprintf('\n|   Iter    | Nb iter |  u [mm]   |  f [kN]   |  max(d)   |');
    fprintf('\n+-----------+---------+-----------+-----------+-----------+\n');
end

ismonotonic = ~any(diff(sign(t(t~=0))));
numddlb = findddl(S_phase,'T',BRight);
db = d(numddlb,:);

for i=1:length(T)
    
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
                S = addbcSingleEdgeCrack(S,ud,BU,BL,BLeft,BRight,BFront,BBack,loading);
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
            
            % Display solution fields
            if displaySol
                % Display phase field
                plotSolution(S_phase,d);
                
                % Display displacement field
                for j=1:Dim
                    plotSolution(S,u,'displ',j);
                end
                
                % Display internal energy field
                if strcmpi(PFsolver,'historyfieldnode')
                    plotSolution(S_phase,H);
                else
                    figure('Name','Solution H')
                    clf
                    plot(H,S_phase);
                    colorbar
                    set(gca,'FontSize',16)
                end
                
                % figure('Name','Solution H')
                % clf
                % subplot(1,2,1)
                % plot_sol(S,u,'energyint','local');
                % colorbar
                % title('Internal energy')
                % subplot(1,2,2)
                % plot_sol(S,u,'energyint',{'positive','local'});
                % colorbar
                % title('Positive internal energy')
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
        if ismonotonic
            f = abs(f);
        end
    end
    
    % Update fields
    ft(i) = f;
    dmaxt(i) = dmax;
    
    if display_
        fprintf('| %4d/%4d | %7d | %9.3e | %9.3e | %9.3e |\n',i,length(T),nbIter,t(i)*1e3,f*((Dim==2)*1e-6+(Dim==3)*1e-3),dmax);
    end
end

if display_
    fprintf('+-----------+---------+-----------+-----------+-----------+\n');
end

end
