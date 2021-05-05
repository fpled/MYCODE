function [ut,result,vt,at] = solveBeamDetDynLinElasClampedFreeInitialTime(param,S,funN,M,funb,funu0,v0,P1,varargin)
% function [ut,result,vt,at] = solveBeamDetDynLinElasClampedFreeInitialTime(param,S,funN,M,funb,funu0,v0,P1,varargin)

tinit = param(1); % initial time [s]
E = param(2)*1e9; % Young modulus [Pa]
alpha = param(3); % stiffness proportional Rayleigh (viscous) damping coefficient
beta = param(4); % mass proportional Rayleigh (viscous) damping coefficient
if length(param)==6
    c = param(5)*1e3; % junction rotational stiffness [N.m/rad]
    J = param(6); % moment of inertia [kg.m2/rad]=[N.m.s2/rad]
    u0 = funu0(tinit,E,c);
else
    u0 = funu0(tinit);
end
u0 = freevector(S,u0(:));

% Time solver
N = funN(tinit);

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
if length(param)==6
    % [~,numnode,~] = intersect(S,P1,'strict',false);
    numnode = find(S.node==P1);
    numddl = findddl(S,'RZ',numnode,'free');
    K(numddl,numddl) = K(numddl,numddl) + c;
    M(numddl,numddl) = M(numddl,numddl) + J;
end
C = alpha*K + beta*M;
% b0 = zeros(getnbddlfree(S),1);
% loadFunction = @(N) zero(N);
% b = b0*loadFunction(N);
b = funb(tinit);

[ut,result,vt,at] = ddsolve(N,b,M,K,C,u0,v0);

end
