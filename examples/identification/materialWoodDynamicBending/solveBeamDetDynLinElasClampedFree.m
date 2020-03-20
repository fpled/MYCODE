function [ut,result,vt,at] = solveBeamDetDynLinElasClampedFree(param,S,N,M,b,u0,v0)
% function [ut,result,vt,at] = solveBeamDetDynLinElasClampedFree(param,S,N,M,b,u0,v0)

E = param(1)*1e9; % Young modulus [Pa]
alpha = param(2); % stiffness proportional Rayleigh (viscous) damping coefficient
beta = param(3); % mass proportional Rayleigh (viscous) damping coefficient

% Update materials
mats = MATERIALS(S);
for k=1:length(mats)
    mats{k} = setparam(mats{k},'E',E);
    % mats{k} = setparam(mats{k},'NU',NU);
end
S = actualisematerials(S,mats);

% Mass, stiffness and damping matrices and sollicitation vectors
K = calc_rigi(S);
C = alpha*K + beta*M;
    
[ut,result,vt,at] = ddsolve(N,b,M,K,C,u0,v0);

end
