function [ut,result,vt,at] = solveBeamDetDynLinElasClampedFree(param,S,N,M,b,funu0,v0,P1,varargin)
% function [ut,result,vt,at] = solveBeamDetDynLinElasClampedFree(param,S,N,M,b,funu0,v0,P1,varargin)

E = param(1)*1e9; % Young modulus [Pa]
alpha = param(2); % stiffness proportional Rayleigh (viscous) damping coefficient
beta = param(3); % mass proportional Rayleigh (viscous) damping coefficient
if length(param)==5
    c = param(4)*1e3; % junction rotational stiffness [N.m/rad]
    J = param(5); % moment of inertia [kg.m2/rad]=[N.m.s2/rad]
    u0 = funu0(E,c);
else
    u0 = funu0;
end
u0 = freevector(S,u0(:));

% Update materials
mats = MATERIALS(S);
for k=1:length(mats)
    mats{k} = setparam(mats{k},'E',E);
    % mats{k} = setparam(mats{k},'NU',NU);
end
S = actualisematerials(S,mats);

% Mass, stiffness and damping matrices and sollicitation vectors
% M = calc_mass(S);
K = calc_rigi(S);
if length(param)==5
    % [~,numnode,~] = intersect(S,P1,'strict',false);
    numnode = find(S.node==P1);
    numddl = findddl(S,'RZ',numnode,'free');
    K(numddl,numddl) = K(numddl,numddl) + c;
    M(numddl,numddl) = M(numddl,numddl) + J;
end
C = alpha*K + beta*M;

[ut,result,vt,at] = ddsolve(N,b,M,K,C,u0,v0);

end
