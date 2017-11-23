function [u_in,S] = solveTractionIsotTrans(param,S)
% function [u_in,S] = solveTractionIsotTrans(param,S)

EL = param(1);
NUL = param(2);
GL = param(3);

mats = MATERIALS(S);
for k=1:length(mats)
    mats{k} = setparam(mats{k},'EL',EL);
    mats{k} = setparam(mats{k},'NUL',NUL);
    mats{k} = setparam(mats{k},'GL',GL);
end
S = actualisematerials(S,mats);

[A,b] = calc_rigi(S);
b = -b;
u_in = A\b;

end
