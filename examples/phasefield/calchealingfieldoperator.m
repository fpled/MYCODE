function [S_healing,A_healing,b_healing] = calchealingfieldoperator(S_healing,r,qn,H,d,h,heff,fundact,PFregularization)
% function [S_healing,A_healing,b_healing] = calchealingfieldoperator(S_healing,r,qn,H,d,h,heff,fundact,PFregularization)

mats_healing = MATERIALS(S_healing);

R = r + 2*heff^2*H;
% if isa(R,'FENODEFIELD')
%     R = R + fundact(d,h).*r;
% else
%     node = getnode(S_healing);
%     for m=1:getnbgroupelem(S_healing)
%         elem = getgroupelem(S_healing,m);
%         xnode = node(elem);
%         gauss = calc_gauss(elem,'mass');
%         xgauss = gauss.coord;
%         % N = calc_N(elem,xnode,xgauss,'nbddlpernode',1);
%         N = calc_N(elem,xnode,xgauss);
%         % de = localize(elem,d,'scalar');
%         % he = localize(elem,h,'scalar');
%         de = localize(elem,d);
%         he = localize(elem,h);
%         de = N*de;
%         he = N*he;
%         R{m} = R{m} + fundact(de,he).*r;
%     end
%     % R = FEELEMFIELD(R,'storage','gauss','type','scalar','ddl',DDL('R'));
% end

% R = r;
% if isa(R,'FENODEFIELD')
%     R = R + 2*heff^2*fundact(d,h).^2.*H;
% else
%     node = getnode(S_healing);
%     for m=1:getnbgroupelem(S_healing)
%         elem = getgroupelem(S_healing,m);
%         xnode = node(elem);
%         gauss = calc_gauss(elem,'mass');
%         xgauss = gauss.coord;
%         % N = calc_N(elem,xnode,xgauss,'nbddlpernode',1);
%         N = calc_N(elem,xnode,xgauss);
%         % de = localize(elem,d,'scalar');
%         % he = localize(elem,h,'scalar');
%         de = localize(elem,d);
%         he = localize(elem,h);
%         de = N*de;
%         he = N*he;
%         R{m} = R{m} + 2*heff^2*fundact(de,he).^2.*H{m};
%     end
%     % R = FEELEMFIELD(R,'storage','gauss','type','scalar','ddl',DDL('R'));
% end
for m=1:length(mats_healing)
    if isa(R,'FENODEFIELD')
        mats_healing{m} = setparam(mats_healing{m},'r',R);
    else
        mats_healing{m} = setparam(mats_healing{m},'r',R{m});
    end
end
S_healing = actualisematerials(S_healing,mats_healing);

[A_healing,b_healing] = calc_rigi(S_healing);

Q = qn;
if isa(Q,'FENODEFIELD')
    switch lower(PFregularization)
        case 'at1'
            Q = Q - 2*fundact(d,h).*qn;
            % Q = Q - fundact(d,h).*qn;
            % Q = Q - 2*qn.*d;
        case 'at2'
            Q = Q + fundact(d,h).*r;
            % Q = Q + r.*d;
    end
    % Q = Q - 2*heff*fundact(d,h).*(1-d).*H;
    Q = Q - 2*heff*(1-d).*H;
else
    node = getnode(S_healing);
    for m=1:getnbgroupelem(S_healing)
        elem = getgroupelem(S_healing,m);
        xnode = node(elem);
        gauss = calc_gauss(elem,'mass');
        xgauss = gauss.coord;
        % N = calc_N(elem,xnode,xgauss,'nbddlpernode',1);
        N = calc_N(elem,xnode,xgauss);
        % de = localize(elem,d,'scalar');
        % he = localize(elem,h,'scalar');
        de = localize(elem,d);
        he = localize(elem,h);
        de = N*de;
        he = N*he;
        switch lower(PFregularization)
            case 'at1'
                Q{m} = Q{m} - 2*fundact(de,he).*qn{m};
                % Q{m} = Q{m} - fundact(de,he).*qn{m};
                % Q{m} = Q{m} - 2*qn{m}.*de;
            case 'at2'
                Q{m} = Q{m} + fundact(de,he).*r{m};
                % Q{m} = Q{m} + r{m}.*de;
        end
        % Q{m} = Q{m} - 2*heff*fundact(de,he).*(1-de).*H{m};
        Q{m} = Q{m} - 2*heff*(1-de).*H{m};
    end
    % Q = FEELEMFIELD(Q,'storage','gauss','type','scalar','ddl',DDL('Q'));
end
if isa(Q,'FENODEFIELD')
    q = double(Q);
    q = max(q,0);
else
    q = getvalue(Q);
    for p=1:length(q)
        qe = double(q{p});
        qe = max(qe,0);
        q{p} = MYDOUBLEND(qe);
    end
end
Q = setvalue(Q,q);
b_healing = -b_healing + bodyload(S_healing,[],'QN',Q);
