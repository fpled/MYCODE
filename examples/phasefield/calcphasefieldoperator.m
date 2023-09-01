function [S_phase,A_phase,b_phase] = calcphasefieldoperator(S_phase,r,qn,H)
% function [S_phase,A_phase,b_phase] = calcphasefieldoperator(S_phase,r,qn,H)

mats_phase = MATERIALS(S_phase);
R = r+2*H;
for m=1:length(mats_phase)
    if isa(R,'FENODEFIELD')
        mats_phase{m} = setparam(mats_phase{m},'r',R);
    else
        mats_phase{m} = setparam(mats_phase{m},'r',R{m});
    end
end
S_phase = actualisematerials(S_phase,mats_phase);

[A_phase,b_phase] = calc_rigi(S_phase);

Q = 2*H+qn;
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
