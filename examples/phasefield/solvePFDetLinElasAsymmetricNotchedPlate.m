function [Ht,dt,ut,ft] = solvePFDetLinElasAsymmetricNotchedPlate(S,S_phase,T,PU,PL,PR,varargin)
% function [Ht,dt,ut,ft] = solvePFDetLinElasAsymmetricNotchedPlate(S,S_phase,T,PU,PL,PR,varargin)
% Solve deterministic Phase Field problem.

display_ = ischarin('display',varargin);

t = gett(T);

Ht = cell(1,length(T));
dt = cell(1,length(T));
ut = cell(1,length(T));
ft = zeros(1,length(T));

sz_phase = getnbddl(S_phase);
sz = getnbddl(S);
H = zeros(sz_phase,1);
u = zeros(sz,1);

if display_
    fprintf('\n+----------+-----------+-----------+------------+------------+------------+\n');
    fprintf('|   Iter   |  u [mm]   | f [kN/mm] |  norm(H)   |  norm(d)   |  norm(u)   |\n');
    fprintf('+----------+-----------+-----------+------------+------------+------------+\n');
end

mats_phase = MATERIALS(S_phase);
r = zeros(length(mats_phase),1);
for m=1:length(mats_phase)
    r(m) = getparam(mats_phase{m},'r');
end

for i=1:length(T)
    
    % Internal energy field
    h_old = double(H);
    H = FENODEFIELD(calc_energyint(S,u,'node','positive'));
    h = double(H);
    rep = find(h <= h_old);
    h(rep) = h_old(rep);
    H = setvalue(H,h);
    
    % Phase field
    mats_phase = MATERIALS(S_phase);
    for m=1:length(mats_phase)
        mats_phase{m} = setparam(mats_phase{m},'r',FENODEFIELD(r(m)+2*H));
    end
    S_phase = actualisematerials(S_phase,mats_phase);
    
    [A_phase,b_phase] = calc_rigi(S_phase);
    b_phase = -b_phase + bodyload(S_phase,[],'QN',FENODEFIELD(2*H));
    
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
    Ht{i} = double(H);
    dt{i} = d;
    ut{i} = u;
    ft(i) = f;
    
    if display_
        fprintf('| %8d | %6.3e | %6.3e | %9.4e | %9.4e | %9.4e |\n',i,t(i)*1e3,ft(i)*1e-6,norm(Ht{i}),norm(dt{i}),norm(ut{i}));
    end
end

if display_
    fprintf('+----------+-----------+-----------+------------+------------+------------+\n');
end

Ht = TIMEMATRIX(Ht,T,[sz_phase,1]);
dt = TIMEMATRIX(dt,T,[sz_phase,1]);
ut = TIMEMATRIX(ut,T,[sz,1]);

end