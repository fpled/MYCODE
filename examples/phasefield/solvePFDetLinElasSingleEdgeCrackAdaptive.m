function [dt,ut,ft,dinct,St_phase,St] = solvePFDetLinElasSingleEdgeCrackAdaptive(S_phase,S,T,C,BU,BL,BRight,BLeft,BFront,BBack,loading,sizemap,varargin)
% function [dt,ut,ft,dinct,St_phase,St] = solvePFDetLinElasSingleEdgeCrackAdaptive(S_phase,S,T,C,BU,BL,BRight,BLeft,BFront,BBack,loading,sizemap,varargin)
% Solve deterministic Phase Field problem with mesh adaptation.

display_ = ischarin('display',varargin);
filename = getcharin('filename',varargin,'gmsh_domain_single_edge_crack');
pathname = getcharin('pathname',varargin,'.');
gmshoptions = getcharin('gmshoptions',varargin,'-v 0');
mmgoptions = getcharin('mmgoptions',varargin,'-nomove -v -1');

Dim = getdim(S);

t = gett(T);

dt = cell(1,length(T));
ut = cell(1,length(T));
ft = zeros(1,length(T));
dinct = cell(1,length(T)); % increment of phase field
tol = 1e-12;
St_phase = cell(1,length(T));
St = cell(1,length(T));

u = calc_init_dirichlet(S);
d = calc_init_dirichlet(S_phase);
H = calc_energyint(S,u,'positive');

if display_
    fprintf('\n+-----------+-----------+-----------+----------+----------+------------+------------+\n');
    fprintf('|   Iter    |  u [mm]   |  f [kN]   | Nb nodes | Nb elems |  norm(d)   |  norm(u)   |\n');
    fprintf('+-----------+-----------+-----------+----------+----------+------------+------------+\n');
    fprintf('| %4d/%4d | %6.3e | %6.3e | %8d | %8d | %9.4e | %9.4e |\n',0,0,0,0,getnbnode(S),getnbelem(S),0,0);
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
    
    % Mesh adaptation
    mats = MATERIALS(S);
    S_phase_old = S_phase;
    % S_old = S;
    cl = sizemap(d);
    S_phase = adaptmesh(S_phase_old,cl,fullfile(pathname,filename),'gmshoptions',gmshoptions,'mmgoptions',mmgoptions);
    S = S_phase;
    
    for m=1:length(mats_phase)
        S_phase = setmaterial(S_phase,mats_phase{m},m);
    end
    S_phase = final(S_phase,'duplicate');
    S_phase = addcl(S_phase,C,'T',1);
    
    P_phase = calcProjection(S_phase,S_phase_old,[],'free',false,'full',true);
    d = P_phase'*d;
    dold = P_phase'*dold;
    dinc = d - dold;
    % dincmin = min(dinc); if dincmin<-tol, dincmin, end
    
    % Displacement field
    % P = calcProjection(S,S_old,[],'free',false,'full',true);
    P = kron(P_phase,eye(Dim));
    u = P'*u;
    
    for m=1:length(mats)
        mats{m} = setparam(mats{m},'d',d);
        mats{m} = setparam(mats{m},'u',u);
        S = setmaterial(S,mats{m},m);
    end
    S = final(S,'duplicate');
    ud = t(i);
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
            error('Wrong loading case')
    end
    
    H = calc_energyint(S,u,'positive');
    
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
    dinct{i} = dinc;
    St_phase{i} = S_phase;
    St{i} = S;
    
    if display_
        fprintf('| %4d/%4d | %6.3e | %6.3e | %8d | %8d | %9.4e | %9.4e |\n',i,length(T),t(i)*1e3,ft(i)*((Dim==2)*1e-6+(Dim==3)*1e-3),getnbnode(S),getnbelem(S),norm(dt{i}),norm(ut{i}));
    end
end

if display_
    fprintf('+-----------+-----------+-----------+----------+----------+------------+------------+\n');
end

end
