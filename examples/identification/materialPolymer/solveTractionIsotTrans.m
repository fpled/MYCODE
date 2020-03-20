function [u_in,S] = solveTractionIsotTrans(param,S)
% function [u_in,S] = solveTractionIsotTrans(param,S)

EL = param(1); % longitudinal Young modulus
NUL = param(2); % longitudinal Poisson ratio
GL = param(3); % longitudinal shear modulus

% Update materials
mats = MATERIALS(S);
for k=1:length(mats)
    mats{k} = setparam(mats{k},'EL',EL);
    mats{k} = setparam(mats{k},'NUL',NUL);
    mats{k} = setparam(mats{k},'GL',GL);
end
S = actualisematerials(S,mats);

% Solve
[A,b] = calc_rigi(S);
b = -b;
u_in = A\b;

end
