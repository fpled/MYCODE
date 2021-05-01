function [dt,ut,ft] = solvePFDetLinElasAsymmetricNotchedPlate(S_phase,S,T,PU,PL,PR,varargin)
% function [dt,ut,ft] = solvePFDetLinElasAsymmetricNotchedPlate(S_phase,S,T,PU,PL,PR,varargin)
% Solve deterministic Phase Field problem.

display_ = ischarin('display',varargin);

t = gett(T);

dt = cell(1,length(T));
ut = cell(1,length(T));
ft = zeros(1,length(T));

sz_d = getnbddl(S_phase);
sz_u = getnbddl(S);
u = zeros(sz_u,1);
H = calc_energyint(S,u,'positive');

if display_
    fprintf('\n+----------+-----------+-----------+------------+------------+\n');
    fprintf('|   Iter   |  u [mm]   |  f [kN]   |  norm(d)   |  norm(u)   |\n');
    fprintf('+----------+-----------+-----------+------------+------------+\n');
end

mats_phase = MATERIALS(S_phase);
r = zeros(length(mats_phase),1);
for m=1:length(mats_phase)
    r(m) = getparam(mats_phase{m},'r');
end

for i=1:length(T)
    
    % Internal energy field
    h_old = getvalue(H);
    H = calc_energyint(S,u,'positive');
    h = getvalue(H);
    for p=1:getnbgroupelem(S)
        he = double(h{p});
        he_old = double(h_old{p});
        rep = find(he <= he_old);
        he(rep) = he_old(rep);
        h{p} = he;
    end
    H = FEELEMFIELD(h,'storage',getstorage(H),'type',gettype(H),'ddl',getddl(H));
    
    % Phase field
    mats_phase = MATERIALS(S_phase);
    for m=1:length(mats_phase)
        mats_phase{m} = setparam(mats_phase{m},'r',r(m)+2*H);
    end
    S_phase = actualisematerials(S_phase,mats_phase);
    
    [A_phase,b_phase] = calc_rigi(S_phase);
    b_phase = -b_phase + bodyload(S_phase,[],'QN',2*H);
    
    d = A_phase\b_phase;
    d = unfreevector(S_phase,d);
    
    % Displacement field
    mats = MATERIALS(S);
    for m=1:length(mats)
        mats{m} = setparam(mats{m},'d',d);
        mats{m} = setparam(mats{m},'u',u);
    end
    S = actualisematerials(S,mats);
    S = removebc(S);
    ud = -t(i);
    S = addcl(S,PU,'UY',ud);
    S = addcl(S,PL,{'UX','UY'});
    S = addcl(S,PR,'UY');
    
    [A,b] = calc_rigi(S,'nofree');
    b = -b;
    
    u = freematrix(S,A)\b;
    u = unfreevector(S,u);
    
    numddl = findddl(S,'UY',PU);
    f = -A(numddl,:)*u;
    % f = sum(f);
    
    % Update fields
    dt{i} = d;
    ut{i} = u;
    ft(i) = f;
    
    if display_
        fprintf('| %8d | %6.3e | %6.3e | %9.4e | %9.4e |\n',i,t(i)*1e3,ft(i)*1e-6,norm(dt{i}),norm(ut{i}));
    end
end

if display_
    fprintf('+----------+-----------+-----------+------------+------------+\n');
end

dt = TIMEMATRIX(dt,T,[sz_d,1]);
ut = TIMEMATRIX(ut,T,[sz_u,1]);

end
