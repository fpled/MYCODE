function [dt,ut,ft,Ht] = solvePFDetLinElasPlatewithHole(S_phase,S,T,PFsolver,BU,BL,P0,varargin)
% function [dt,ut,ft,Ht] = solvePFDetLinElasPlatewithHole(S_phase,S,T,PFsolver,BU,BL,P0,varargin)
% Solve deterministic Phase Field problem.

display_ = ischarin('display',varargin);

Dim = getdim(S);

t = gett(T);

dt = cell(1,length(T));
ut = cell(1,length(T));
ft = zeros(1,length(T));
if nargout>=4
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
    optimFun = 'lsqnonlin'; % 'fmincon' or 'lsqnonlin'
    % optimFun = 'fmincon';

    displayoptim = 'off';
    % displayoptim = 'iter';
    % displayoptim = 'iter-detailed';
    % displayoptim = 'final';
    % displayoptim = 'final-detailed';

    tolX = 1e-9; % tolerance on the parameter value
    tolFun = 1e-9; % tolerance on the function value
    % maxFunEvals = Inf; % maximum number of function evaluations

    % optimAlgo = 'interior-point';
    % optimAlgo = 'trust-region-reflective';
    % optimAlgo = 'sqp';
    % optimAlgo = 'active-set';
    % optimAlgo = 'levenberg-marquardt';

    % options  = optimoptions(optimFun,'Display',displayoptim,'StepTolerance',tolX,'FunctionTolerance',tolFun,...
    %     'OptimalityTolerance',tolFun...%,'MaxFunctionEvaluations',maxFunEvals...%,'Algorithm',optimAlgo...
    %     ,'SpecifyObjectiveGradient',true...
    %     );
    options  = optimoptions(optimFun,'Display',displayoptim,'StepTolerance',tolX,'FunctionTolerance',tolFun,'OptimalityTolerance',tolFun,...
        'SpecifyObjectiveGradient',true);
end

if display_
    fprintf('\n+-----------+-----------+-----------+------------+------------+\n');
    fprintf('|   Iter    |  u [mm]   |  f [kN]   |  norm(d)   |  norm(u)   |\n');
    fprintf('+-----------+-----------+-----------+------------+------------+\n');
end

mats_phase = MATERIALS(S_phase);
r = cell(length(mats_phase),1);
for m=1:length(mats_phase)
    r{m} = getparam(mats_phase{m},'r');
end

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
        if strcmpi(PFsolver,'historyfieldnode')
            mats_phase{m} = setparam(mats_phase{m},'r',r{m}+2*H);
        else
            mats_phase{m} = setparam(mats_phase{m},'r',r{m}+2*H{m});
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
            if i<=2
                d = A_phase\b_phase;
            else
                d0 = freevector(S_phase,d);
                lb = d0;
                ub = ones(size(d0));
                switch optimFun
                    case 'lsqnonlin'
                        fun = @(d) funlsqnonlinPF(d,A_phase,b_phase);
                        [d,err,~,exitflag,output] = lsqnonlin(fun,d0,lb,ub,options);
                    case 'fmincon'
                        fun = @(d) funoptimPF(d,A_phase,b_phase);
                        [d,err,exitflag,output] = fmincon(fun,d0+eps,[],[],[],[],lb,ub,[],options);
                end
            end
    end
    d = unfreevector(S_phase,d);
    % dinc = d - d_old;
    % dincmin = min(dinc); if dincmin<-tol, dincmin, end
    
    % Displacement field
    mats = MATERIALS(S);
    for m=1:length(mats)
        mats{m} = setparam(mats{m},'d',d);
        mats{m} = setparam(mats{m},'u',u);
    end
    S = actualisematerials(S,mats);
    S = removebc(S);
    ud = -t(i);
    if Dim==2
        S = addcl(S,BU,'UY',ud);
    elseif Dim==3
        S = addcl(S,BU,'UY',ud);
    end
    S = addcl(S,BL,'UY');
    S = addcl(S,P0,'UX');
    
    [A,b] = calc_rigi(S,'nofree');
    b = -b;
    
    u = freematrix(S,A)\b;
    u = unfreevector(S,u);
    
    numddl = findddl(S,'UY',BU);
    f = -A(numddl,:)*u;
    f = sum(f);
    
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
    % dinct{i} = dinc;
    
    if display_
        fprintf('| %4d/%4d | %6.3e | %6.3e | %9.4e | %9.4e |\n',i,length(T),t(i)*1e3,ft(i)*((Dim==2)*1e-6+(Dim==3)*1e-3),norm(dt{i}),norm(ut{i}));
    end
end

if display_
    fprintf('+-----------+-----------+-----------+------------+------------+\n');
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
% dinct = TIMEMATRIX(dinct,T,size(dinc));

end