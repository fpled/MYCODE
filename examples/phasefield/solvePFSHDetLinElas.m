function [dt,ht,ut,ft,Ht,Edt,Eht,Eut,output] = solvePFSHDetLinElas(S_phase,S_healing,S,T,PFsolver,addbc,findddlforce,findddlboundary,varargin)
% function [dt,ht,ut,ft,Ht,Edt,Eht,Eut,output] = solvePFSHDetLinElas(S_phase,S_healing,S,T,PFsolver,addbc,findddlforce,findddlboundary,varargin)
% Solve deterministic phase-field problem.

addbc = fcnchk(addbc);
findddlforce = fcnchk(findddlforce);
findddlboundary = fcnchk(findddlboundary);
display_ = getcharin('display',varargin,true);
displayIter = getcharin('displayiter',varargin,false);
displaySol = getcharin('displaysol',varargin,false);
maxIter = getcharin('maxiter',varargin,100);
tolConv = getcharin('tol',varargin,1e-2);
critConv = getcharin('crit',varargin,'Energy');
heff = getcharin('heff',varargin,0.5);
deff = getcharin('deff',varargin,@(d,h) d-heff*h);
dact = getcharin('dact',varargin,0.5);
fundact = getcharin('fundact',varargin,@(d,h) deff(d,h)>dact);
PFregularization = getcharin('PFregularization',varargin,'AT2');
dbth = getcharin('dbth',varargin,0.999);

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
ht = cell(1,length(T));
ut = cell(1,length(T));
ft = zeros(1,length(T));
if nargout>=5
    Ht = cell(1,length(T));
end
if nargout>=6
    Edt = zeros(1,length(T));
end
if nargout>=7
    Eht = zeros(1,length(T));
end
if nargout>=8
    Eut = zeros(1,length(T));
end
if nargout>=9
    iteration = zeros(1,length(T));
    time = zeros(1,length(T));
    err = zeros(1,length(T));
end

d = calc_init_dirichlet(S_phase);
h = calc_init_dirichlet(S_healing);
u = calc_init_dirichlet(S);
if strcmpi(PFsolver,'historyfieldnode')
    H = FENODEFIELD(calc_energyint(S,u,'node','positive','local'));
    r = FENODEFIELD(calc_parammat(S_phase,'r','node'));
    qn = FENODEFIELD(calc_parammat(S_phase,'qn','node'));
    rh = FENODEFIELD(calc_parammat(S_healing,'r','node'));
    qnh = FENODEFIELD(calc_parammat(S_healing,'qn','node'));
else
    H = calc_energyint(S,u,'intorder','mass','positive','local');
    r = calc_parammat(S_phase,'r');
    qn = calc_parammat(S_phase,'qn');
    rh = calc_parammat(S_healing,'r');
    qnh = calc_parammat(S_healing,'qn');
end
Ae_phase = calc_rigi(S_phase,'nofree');
be_phase = bodyload(S_phase,[],'QN',qn,'nofree');
Ae_healing = calc_rigi(S_healing,'nofree');
be_healing = bodyload(S_healing,[],'QN',qnh,'nofree');
if checkConvEnergy
    Ed = 1/2*d'*Ae_phase*d - d'*be_phase;
    Eh = 1/2*h'*Ae_healing*h - h'*be_healing;
    A = calc_rigi(S,'nofree');
    Eu = 1/2*u'*A*u;
    E = Ed+Eh+Eu;
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
            % options = optimoptions(optimFun,'Display',displayoptim,'Algorithm',optimAlgo,'MaxIterations',Inf);
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
    fprintf('\n+-----------+---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+');
    fprintf('\n|   Iter    | Nb iter |  u [mm]   |  f [kN]   |  max(d)   |  max(h)   |  max(D)   |  Ed [J]   |  Eh [J]   |  Eu [J]   |');
    fprintf('\n+-----------+---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+\n');
end

if displaySol
    fontsize = 16;
    fpos = get(groot,'DefaultFigurePosition');
    spos = get(groot,'ScreenSize');
    fd = figure('Name','Damage/Phase field',...
        'Position',[spos(1) fpos(2:4)]);
    clf
    fh = figure('Name','Healing field',...
        'Position',fpos);
    clf
    fD = figure('Name','Effective damage/phase field',...
        'Position',[spos(3)-fpos(3) fpos(2:4)]);
    clf
end

ismonotonic = ~any(diff(sign(t(t~=0))));
numddlb = findddlboundary(S_phase);
db = d(numddlb,:);
hb = h(numddlb,:);
Db = deff(db,hb);

for i=1:length(T)
    tIter = tic;
    nbIter = 0;
    if any(Db > dbth)
        f = 0;
    else
        if strcmpi(PFsolver,'historyfieldelem') || strcmpi(PFsolver,'historyfieldnode')
            H_old = H;
        end
        d_old = d;
        h_old = h;
        if checkConvRes
            [S_phase,A_phase,b_phase] = calcphasefieldoperator(S_phase,r,qn,H,h,heff,d,fundact);
            [S_healing,A_healing,b_healing] = calchealingfieldoperator(S_healing,rh,qnh,H,d,h,heff,fundact,PFregularization);
        end
        errConv = Inf;
        errConvs = 0; errConvr = 0; errConve = 0;
        while (errConv > tolConv) && (nbIter < maxIter)
            nbIter = nbIter+1;
            if checkConvSol
                d_prev = d;
                h_prev = h;
                u_prev = u;
            end
            if checkConvEnergy
                E_prev = E;
            end
            
            % Phase field
            if ~checkConvRes
                [S_phase,A_phase,b_phase] = calcphasefieldoperator(S_phase,r,qn,H,h,heff,d_old,fundact);
            end
            
            switch lower(PFsolver)
                case {'historyfieldelem','historyfieldnode'}
                    d = A_phase\b_phase;
                otherwise
                    d0 = freevector(S_phase,d_old);
                    h0 = freevector(S_phase,h);
                    lb = max([d0,heff*h0],[],2);
                    ub = ones(size(d0))+heff*h0;
                    % lb = d0;
                    % ub = [];
                    lb(lb==ub) = lb(lb==ub)-eps;
                    switch optimFun
                        case 'lsqlin'
                            d = lsqlin(A_phase,b_phase,[],[],[],[],lb,ub,d0,options);
                            % [d,resnormd,~,exitflagd,outputd] = lsqlin(A_phase,b_phase,[],[],[],[],lb,ub,d0,options);
                            % resnormd
                            % exitflagd
                            % outputd
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
            % db = d(numddlb,:);
            
            % Healing field
            if ~checkConvRes
                [S_healing,A_healing,b_healing] = calchealingfieldoperator(S_healing,rh,qnh,H,d,h_old,heff,fundact,PFregularization);
            end
            
            switch lower(PFsolver)
                case {'historyfieldelem','historyfieldnode'}
                    h = A_healing\b_healing;
                otherwise
                    h0 = freevector(S_healing,h_old);
                    d0 = freevector(S_healing,d);
                    lb = max([h0,(d0-1)/heff],[],2);
                    ub = min([ones(size(h0)),d0/heff],[],2);
                    % lb = h0;
                    % ub = ones(size(h0));
                    lb(lb==ub) = lb(lb==ub)-eps;
                    switch optimFun
                        case 'lsqlin'
                            h = lsqlin(A_healing,b_healing,[],[],[],[],lb,ub,h0,options);
                            % [h,resnormh,~,exitflagh,outputh] = lsqlin(A_healing,b_healing,[],[],[],[],lb,ub,h0,options);
                            % resnormh
                            % exitflagh
                            % outputh
                        case 'lsqnonlin'
                            fun = @(h) funlsqnonlinPF(h,A_healing,b_healing);
                            h = lsqnonlin(fun,h0,lb,ub,options);
                        case 'fmincon'
                            fun = @(h) funoptimPF(h,A_healing,b_healing);
                            h = fmincon(fun,h0+eps,[],[],[],[],lb,ub,[],options);
                    end
            end
            hmax = max(h);
            h = unfreevector(S_healing,h);
            % hb = h(numddlb,:);
            
            % Effective damage field
            D = deff(d,h);
            Dmax = max(D);
            D = unfreevector(S_phase,D);
            Db = D(numddlb,:);
            
            % Displacement field
            mats = MATERIALS(S);
            for m=1:length(mats)
                mats{m} = setparam(mats{m},'d',d);
                mats{m} = setparam(mats{m},'h',h);
                mats{m} = setparam(mats{m},'u',u);
            end
            S = actualisematerials(S,mats);
            if nbIter==1
                S = removebc(S);
                S = addbc(S,t(i));
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
                figure(fd)
                clf
                plot_sol(S_phase,d);
                colorbar
                set(gca,'FontSize',fontsize)
                figure(fh)
                clf
                plot_sol(S_healing,h);
                colorbar
                set(gca,'FontSize',fontsize)
                figure(fD)
                clf
                plot_sol(S_phase,D);
                colorbar
                set(gca,'FontSize',fontsize)
                
                % Display phase field
                % plotSolution(S_phase,d);
                
                % Display healing field
                % plotSolution(S_healing,h);
                
                % Display effective phase field
                % plotSolution(S_phase,D);
                
                % Display displacement field
                % for j=1:Dim
                %     plotSolution(S,u,'displ',j);
                % end
                
                % Display internal energy field
                % if strcmpi(PFsolver,'historyfieldnode')
                %     plotSolution(S_phase,H);
                % else
                %     figure('Name','Solution H')
                %     clf
                %     plot(H,S_phase);
                %     colorbar
                %     set(gca,'FontSize',fontsize)
                % end
                
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
                errConvh = norm(h-h_prev)/norm(h);
                errConvu = norm(u-u_prev)/norm(u);
                errConvs = max(errConvd,errConvh,errConvu);
            end
            if checkConvRes
                % Phase field residual
                [S_phase,A_phase,b_phase] = calcphasefieldoperator(S_phase,r,qn,H,h,heff,d_old,fundact);
                r_phase = A_phase*d - b_phase;
                errConvrd = norm(r_phase)/norm(b_phase);
                % Healing residual
                [S_healing,A_healing,b_healing] = calchealingfieldoperator(S_healing,rh,qnh,H,d,h_old,heff,fundact,PFregularization);
                r_healing = A_healing*h - b_healing;
                errConvrh = norm(r_healing)/norm(b_healing);
                errConvr = max(errConvrd,errConvrh);
            end
            if checkConvEnergy
                % Energy
                Ed = 1/2*d'*Ae_phase*d - d'*be_phase;
                Eh = 1/2*h'*Ae_healing*h - h'*be_healing;
                Eu = 1/2*u'*A*u;
                E = Ed+Eh+Eu;
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
            if any(Db > dbth)
                break
            end
        end
        
        % Force
        numddl = findddlforce(S);
        f = A(numddl,:)*u;
        f = sum(f);
        if ismonotonic
            f = abs(f);
        end
        
        % Energy
        if ~checkConvEnergy
            Ed = 1/2*d'*Ae_phase*d - d'*be_phase;
            Eh = 1/2*h'*Ae_healing*h - h'*be_healing;
            Eu = 1/2*u'*A*u;
        end
    end
    
    % Update fields
    dt{i} = d;
    ht{i} = h;
    ut{i} = u;
    ft(i) = f;
    if nargout>=5
        if strcmpi(PFsolver,'historyfieldnode')
            Ht{i} = double(H);
        else
            Ht{i} = reshape(double(mean(H,4)),[getnbelem(S),1]);
        end
    end
    if nargout>=6
        Edt(i) = Ed;
    end
    if nargout>=7
        Eht(i) = Eh;
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
        fprintf('| %4d/%4d | %7d | %9.3e | %9.3e | %9.3e | %9.3e | %9.3e | %9.3e | %9.3e | %9.3e |\n',i,length(T),nbIter,t(i)*1e3,f*1e-3,dmax,hmax,Dmax,Ed,Eh,Eu);
    end
end

% if displaySol
%     close(fd)
%     close(fh)
%     close(fD)
% end

if display_
    fprintf('+-----------+---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+\n');
end

dt = TIMEMATRIX(dt,T,size(d));
ht = TIMEMATRIX(ht,T,size(h));
ut = TIMEMATRIX(ut,T,size(u));
if nargout>=5
    if strcmpi(PFsolver,'historyfieldnode')
        Ht = TIMEMATRIX(Ht,T,size(d));
    else
        Ht = TIMEMATRIX(Ht,T,[getnbelem(S),1]);
    end
end
if nargout>=9
    output.iteration = iteration;
    output.time = time;
    output.error = err;
end

end
