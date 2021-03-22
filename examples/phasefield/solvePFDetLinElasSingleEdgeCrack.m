function [Ht,dt,ut,ft] = solvePFDetLinElasSingleEdgeCrack(S,S_phase,T,gc,l,loading,L,e,varargin)
% function [Ht,dt,ut,ft] = solvePFDetLinElasSingleEdgeCrack(S,S_phase,T,gc,l,loading,L,e,varargin)
% Solving Phase Field problem.

display_ = getcharin('display',varargin,1);
Dim = getdim(S);
% Dirichlet boundary conditions
if Dim==2
    BU = LIGNE([0.0,L],[L,L]);
    BL = LIGNE([0.0,0.0],[L,0.0]);
    BRight = LIGNE([L,0.0],[L,L]);
    BLeft = LIGNE([0.0,0.0],[0.0,L]);
    BFront = [];
    BBack = [];
elseif Dim==3
    BU = PLAN([0.0,L,0.0],[L,L,0.0],[0.0,L,e]);
    BL = PLAN([0.0,0.0,0.0],[L,0.0,0.0],[0.0,0.0,e]);
    BRight = PLAN([L,0.0,0.0],[L,L,0.0],[L,0.0,e]);
    BLeft = PLAN([0.0,0.0,0.0],[0.0,L,0.0],[0.0,0.0,e]);
    BFront = PLAN([0.0,0.0,e],[L,0.0,e],[0.0,L,e]);
    BBack = PLAN([0.0,0.0,0.0],[L,0.0,0.0],[0.0,L,0.0]);
end

% Solution
t = gett(T);

Ht = cell(1,length(T));
dt = cell(1,length(T));
ut = cell(1,length(T));
ft = zeros(1,length(T));

sz_phase = getnbddl(S_phase);
sz = getnbddl(S);
H = zeros(sz_phase,1);
u = zeros(sz,1);
%
switch display_
    case 1
        fprintf('\n+----------+-----------+-----------+------------+------------+------------+\n');
        fprintf('|   Iter   |  u [mm]   |  f [kN]   |  norm(H)   |  norm(d)   |  norm(u)   |\n');
        fprintf('+----------+-----------+-----------+------------+------------+------------+\n');
    case 0
        fprintf('\nProgress\n');
end
%
lengthT = length(T);
for i=1:lengthT
    
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
        mats_phase{m} = setparam(mats_phase{m},'r',FENODEFIELD(gc/l+2*H));
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
    Ht{i} = double(H);
    dt{i} = d;
    ut{i} = u;
    ft(i) = f;
    
    switch display_
        case 1
            fprintf('| %8d | %6.3e | %6.3e | %9.4e | %9.4e | %9.4e |\n',i,t(i)*1e3,ft(i)*((Dim==2)*1e-6+(Dim==3)*1e-3),norm(Ht{i}),norm(dt{i}),norm(ut{i}));
        case 0
            fprintf([num2str(i),'/',num2str(lengthT),'\n']);
    end
end
%
switch display_
    case 1
        fprintf('+----------+-----------+-----------+------------+------------+------------+\n');
end
%
Ht = TIMEMATRIX(Ht,T,[sz_phase,1]);
dt = TIMEMATRIX(dt,T,[sz_phase,1]);
ut = TIMEMATRIX(ut,T,[sz,1]);