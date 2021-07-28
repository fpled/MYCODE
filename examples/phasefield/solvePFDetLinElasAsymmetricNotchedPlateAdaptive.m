function [dt,ut,ft,dinct,St_phase,St] = solvePFDetLinElasAsymmetricNotchedPlateAdaptive(S_phase,S,T,C,BU,BL,BR,H1,H2,H3,PU,PL,PR,sizemap,varargin)
% function [dt,ut,ft,dinct,St_phase,St] = solvePFDetLinElasAsymmetricNotchedPlateAdaptive(S_phase,S,T,C,BU,BL,BR,H1,H2,H3,PU,PL,PR,sizemap,varargin)
% Solve deterministic Phase Field problem with mesh adaptation.

display_ = ischarin('display',varargin);
filename = getcharin('filename',varargin,'gmsh_domain_asymmetric_notched_plate');
pathname = getcharin('pathname',varargin,'.');
gmshoptions = getcharin('gmshoptions',varargin,'-v 0');
mmgoptions = getcharin('mmgoptions',varargin,'-nomove -v -1');

t = gett(T);

dt = cell(1,length(T));
ut = cell(1,length(T));
ft = zeros(1,length(T));
dinct = cell(1,length(T)); % increment of phase field
tol = 1e-12;
St_phase = cell(1,length(T));
St = cell(1,length(T));

sz_d = getnbddl(S_phase);
sz_u = getnbddl(S);
sz_H = getnbelem(S);
u = zeros(sz_u,1);
d = zeros(sz_d,1);
H = FEELEMFIELD(zeros(sz_H,1),S);

if display_
    fprintf('\n+-----------+-----------+-----------+----------+----------+------------+------------+\n');
    fprintf('|   Iter    |  u [mm]   |  f [kN]   | Nb nodes | Nb elems |  norm(d)   |  norm(u)   |\n');
    fprintf('+-----------+-----------+-----------+----------+----------+------------+------------+\n');
    fprintf('| %4d/%4d | %6.3e | %6.3e | %8d | %8d | %9.4e | %9.4e |\n',0,length(T),0,0,getnbnode(S),getnbelem(S),0,0);
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
    
    dold = d;
    d = A_phase\b_phase;
    d = unfreevector(S_phase,d);
    dinc = d - dold;
    
    if ~isempty(find(dinc<-tol,1)), fprintf('Irreversibility error\n'), end
    
    % Mesh adaptation
    mats = MATERIALS(S);
    S_phase_old = S_phase;
    S_phase_old = addcl(S_phase_old,H1,'T',1);
    S_phase_old = addcl(S_phase_old,H2,'T',1);
    S_phase_old = addcl(S_phase_old,H3,'T',1);
    d_old = freevector(S_phase_old,d);
    d_old = unfreevector(S_phase_old,d_old);
    % S_old = S;
    cl = sizemap(d_old);
    S_phase = adaptmesh(S_phase_old,cl,fullfile(pathname,filename),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    S = S_phase;
    
    for m=1:length(mats_phase)
        S_phase = setmaterial(S_phase,mats_phase{m},m);
    end
    S_phase = final(S_phase,'duplicate');
    S_phase = addcl(S_phase,C,'T',1);
    S_phase = addcl(S_phase,BU,'T');
    S_phase = addcl(S_phase,BL,'T');
    S_phase = addcl(S_phase,BR,'T');
    
    P_phase = calcProjection(S_phase,S_phase_old,[],'free',false,'full',true);
    d = P_phase'*d;
    
    % Displacement field
    % P = calcProjection(S,S_old,[],'free',false,'full',true);
    P = kron(P_phase,eye(2));
    u = P'*u;
    
    for m=1:length(mats)
        mats{m} = setparam(mats{m},'d',d);
        mats{m} = setparam(mats{m},'u',u);
        S = setmaterial(S,mats{m},m);
    end
    S = final(S,'duplicate');
    ud = -t(i);
    S = addcl(S,PU,'UY',ud);
    S = addcl(S,PL,{'UX','UY'});
    S = addcl(S,PR,'UY');
    
    H = calc_energyint(S,u,'positive');
    
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
    dinct{i} = dinc;
    St_phase{i} = S_phase;
    St{i} = S;
    
    if display_
        fprintf('| %4d/%4d | %6.3e | %6.3e | %8d | %8d | %9.4e | %9.4e |\n',i,length(T),t(i)*1e3,ft(i)*1e-6,getnbnode(S),getnbelem(S),norm(dt{i}),norm(ut{i}));
    end
end

if display_
    fprintf('+-----------+-----------+-----------+----------+----------+------------+------------+\n');
end

% DO NOT WORK WITH MESH ADAPTATION
% dt = TIMEMATRIX(dt,T,[sz_d,1]);
% ut = TIMEMATRIX(ut,T,[sz_u,1]);

end
