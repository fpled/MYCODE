function [S_phase,A_phase,b_phase] = calcphasefieldoperator(S_phase,r,qn,H,varargin)
% function [S_phase,A_phase,b_phase] = calcphasefieldoperator(S_phase,r,qn,H)
% function [S_phase,A_phase,b_phase] = calcphasefieldoperator(S_phase,r,qn,H,h,heff,d,fundact)

mats_phase = MATERIALS(S_phase);

R = r + 2*H;
% if nargin>=5
%     h = varargin{1};
%     heff = varargin{2};
%     d = varargin{3};
%     fundact = varargin{4};
%     if isa(R,'FENODEFIELD')
%         R = R + heff*(1-2*h);
%     else
%         node = getnode(S_phase);
%         for m=1:getnbgroupelem(S_phase)
%             elem = getgroupelem(S_phase,m);
%             xnode = node(elem);
%             gauss = calc_gauss(elem,'mass');
%             xgauss = gauss.coord;
%             % N = calc_N(elem,xnode,xgauss,'nbddlpernode',1);
%             N = calc_N(elem,xnode,xgauss);
%             % he = localize(elem,h,'scalar');
%             he = localize(elem,h);
%             he = N*he;
%             R{m} = R{m} + heff*(1-2*he);
%         end
%         % R = FEELEMFIELD(R,'storage','gauss','type','scalar','ddl',DDL('R'));
%     end
% end
for m=1:length(mats_phase)
    if isa(R,'FENODEFIELD')
        mats_phase{m} = setparam(mats_phase{m},'r',R);
    else
        mats_phase{m} = setparam(mats_phase{m},'r',R{m});
    end
end
S_phase = actualisematerials(S_phase,mats_phase);

[A_phase,b_phase] = calc_rigi(S_phase);

Q = qn + 2*H;
if nargin>=5
    h = varargin{1};
    heff = varargin{2};
    d = varargin{3};
    fundact = varargin{4};
    if isa(Q,'FENODEFIELD')
        % Q = Q + 2*heff*fundact(d,h).*h.*H;
        Q = Q + 2*heff*h.*H;
    else
        node = getnode(S_phase);
        for m=1:getnbgroupelem(S_phase)
            elem = getgroupelem(S_phase,m);
            xnode = node(elem);
            gauss = calc_gauss(elem,'mass');
            xgauss = gauss.coord;
            % N = calc_N(elem,xnode,xgauss,'nbddlpernode',1);
            N = calc_N(elem,xnode,xgauss);
            % he = localize(elem,h,'scalar');
            % de = localize(elem,d,'scalar');
            he = localize(elem,h);
            % de = localize(elem,d);
            he = N*he;
            % de = N*de;
            % Q{m} = Q{m} + 2*heff*fundact(de,he).*he.*H{m};
            Q{m} = Q{m} + 2*heff*he.*H{m};
        end
        % Q = FEELEMFIELD(Q,'storage','gauss','type','scalar','ddl',DDL('Q'));
    end
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
b_phase = -b_phase + bodyload(S_phase,[],'QN',Q);
