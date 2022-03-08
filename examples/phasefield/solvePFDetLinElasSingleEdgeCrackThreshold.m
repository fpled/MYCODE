function [dt,ut,ft,Ht] = solvePFDetLinElasSingleEdgeCrackThreshold(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,varargin)
% function [dt,ut,ft,Ht] = solvePFDetLinElasSingleEdgeCrackThreshold(S_phase,S,T,PFsolver,BU,BL,BRight,BLeft,BFront,BBack,loading,varargin)
% Solve deterministic Phase Field problem.

display_ = ischarin('display',varargin);

Dim = getdim(S);

dt0 = T.dt0;
dt1 = T.dt1;
tf = T.tf;
dthreshold = T.dthreshold;

% dt = cell(1,length(T));
% ut = cell(1,length(T));
% ft = zeros(1,length(T));
% if nargout>=4
%     Ht = cell(1,length(T));
% end
% % dinct = cell(1,length(T)); % increment of phase field
% % tol = 1e-12;

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
    fprintf('\n+------+-----------+-----------+------------+------------+\n');
    fprintf('| Iter |  u [mm]   |  f [kN]   |  norm(d)   |  norm(u)   |\n');
    fprintf('+------+-----------+-----------+------------+------------+\n');
end

mats_phase = MATERIALS(S_phase);
r = cell(length(mats_phase),1);
for m=1:length(mats_phase)
    r{m} = getparam(mats_phase{m},'r');
end

i = 0;
ti = 0;
dti = dt0;
while ti < tf
    i = i+1;
    
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
            d0 = freevector(S_phase,d);
            ub = ones(size(d0));
            if isequal(d0,ub)
                d = d0;
            else
                lb = d0;
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
    end
    if any(d > dthreshold)
        dti = dt1;
    end
    ti = ti + dti;

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
    
    % Update fields
    t(i) = ti;
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
        fprintf('| %4d | %6.3e | %6.3e | %9.4e | %9.4e |\n',i,t(i)*1e3,ft(i)*((Dim==2)*1e-6+(Dim==3)*1e-3),norm(dt{i}),norm(ut{i}));
    end
end

if display_
    fprintf('+------+-----------+-----------+------------+------------+\n');
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
% dinct = TIMEMATRIX(dinct,T,size(dinc));

end
