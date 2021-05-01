function [dt,ut,ft] = solvePFDetLinElasSingleEdgeCrackHnode(S_phase,S,T,BU,BL,BRight,BLeft,BFront,BBack,loading,varargin)
% function [dt,ut,ft] = solvePFDetLinElasSingleEdgeCrackHnode(S_phase,S,T,BU,BL,BRight,BLeft,BFront,BBack,loading,varargin)
% Solve deterministic Phase Field problem.

display_ = ischarin('display',varargin);
Dim = getdim(S);

t = gett(T);

dt = cell(1,length(T));
ut = cell(1,length(T));
ft = zeros(1,length(T));

sz_d = getnbddl(S_phase);
sz_u = getnbddl(S);
u = zeros(sz_u,1);
H = zeros(sz_d,1);

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
    ud = t(i);
    switch lower(loading)
        case 'tension'
            S = addcl(S,BU,'UY',ud);
            S = addcl(S,BL,'UY');
            if Dim==2
                S = addcl(S,POINT([0.0,0.0]),'UX');
            elseif Dim==3
                S = addcl(S,POINT([0.0,0.0,0.0]),{'UX','UZ'});
            end
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
            error('Wrong loading case')
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
            error('Wrong loading case')
    end
    f = A(numddl,:)*u;
    f = sum(f);
    
    % Update fields
    dt{i} = d;
    ut{i} = u;
    ft(i) = f;
    
    if display_
        fprintf('| %8d | %6.3e | %6.3e | %9.4e | %9.4e |\n',i,t(i)*1e3,ft(i)*((Dim==2)*1e-6+(Dim==3)*1e-3),norm(dt{i}),norm(ut{i}));
    end
end

if display_
    fprintf('+----------+-----------+-----------+------------+------------+\n');
end

dt = TIMEMATRIX(dt,T,[sz_d,1]);
ut = TIMEMATRIX(ut,T,[sz_u,1]);

end
