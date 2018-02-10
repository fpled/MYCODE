function [u_in,S] = solveThreePointBendingNum(param,S)
% function [u_in,S] = solveThreePointBendingNum(param,S)

EL = param(1);
NUL = param(2);

mats = MATERIALS(S);
for k=1:length(mats)
    mats{k} = setparam(mats{k},'EL',EL);
    mats{k} = setparam(mats{k},'NUL',NUL);
end
S = actualisematerials(S,mats);

[A,b] = calc_rigi(S);
b = -b;
u_in = A\b;

end
