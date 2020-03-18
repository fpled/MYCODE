function [ut,result,vt,at] = solveBeamDetDynLinElasClampedFree(param,S,N,M,b,funu0,v0)
% function [ut,result,vt,at] = solveBeamDetDynLinElasClampedFree(param,S,N,M,b,funu0,v0)

E = param(1)*1e9; % Young modulus [Pa]
NU = param(2); % Poisson ratio
delta = param(3)*1e-2; % initial static displacement [m]
alpha = param(4)*1e-4; % stiffness proportional Rayleigh (viscous) damping coefficient
beta = param(5)*1e-4; % mass proportional Rayleigh (viscous) damping coefficient

% Materials
mats = MATERIALS(S);
for k=1:length(mats)
    mats{k} = setparam(mats{k},'E',E);
    mats{k} = setparam(mats{k},'NU',NU);
end
S = actualisematerials(S,mats);
    
% Initial conditions
u0 = funu0(delta);
u0 = freevector(S,u0(:));

% Mass, stiffness and damping matrices and sollicitation vectors
K = calc_rigi(S);
C = alpha*K + beta*M;
    
[ut,result,vt,at] = ddsolve(N,b,M,K,C,u0,v0);

end
