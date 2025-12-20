function [sigma,sigmaDebonding,epsilonDebonding,damage,...
    dissipatedEnergy,dissipatedEnergyDamage,dissipatedEnergyDebonding]...
    = solveModelRheo(param,epsilon,damageFun,tol,displayIter)
% function [sigma,sigmaDebonding,epsilonDebonding,damage,...
%     dissipatedEnergy,dissipatedEnergyDamage,dissipatedEnergyDebonding]...
%     = solveModelRheo(param,epsilon,damageFun,tol,displayIter)

if nargin<3 || isempty(damageFun)
    damageFun = 1;
end
if nargin<4 || isempty(tol)
    tol = eps;
end
if nargin<5 || isempty(displayIter)
    displayIter = false;
end

display = 'off';
% display = 'iter';
% display = 'iter-detailed';
% display = 'final'; % default
% display = 'final-detailed';

algo = 'trust-region-dogleg'; % default for fsolve
% algo = 'trust-region';
% algo = 'levenberg-marquardt';

tolX = eps; % tolerance on the parameter value
tolFun = eps; % tolerance on the function value
tolOpt = eps; % tolerance on the first-order optimality
maxIters = Inf; % maximum number of iterations
maxFunEvals = Inf; % maximum number of function evaluations
finDiffType = 'forward'; % finite diffferences, 'forward' or 'central' ('forward' by default)

options = optimoptions('fsolve','Display',display,'Algorithm',algo,...
            'StepTolerance',tolX,'FunctionTolerance',tolFun,'OptimalityTolerance',tolOpt,...
            'MaxIterations',maxIters,'MaxFunctionEvaluations',maxFunEvals,'FiniteDifferenceType',finDiffType);

ks = param(2); % elastic modulus of times steel proportion of steel in RC section [MPa]
kc = param(3); % elastic modulus of concrete times proportion of concrete in RC section [MPa]
sig0t = param(4); % first damage stress threshold for concrete under tension [MPa]
sig0c = param(5); % first damage stress threshold for concrete under compression [MPa]
Dsig0 = param(6); % first debonding stress threshold for concrete-steel bond [MPa]

switch damageFun
    case 1
        % damage parameters for concrete under tension
        alphaGt = param(7);
        alphadt = param(8);
        % damage parameters for concrete under compression
        alphaGc = param(9);
        alphadc = param(10);
        alphac = ((sig0c^2/(2*kc))^alphaGc)/((sig0t^2/(2*kc))^alphaGt); % dissymmetric behaviour between tension and compression
        
        zt = @(d) 1-d; % damage function under tension
        dzt = @(d) -1; % derivative of zt with respect to d
        zc = @(d) 1-d; % damage function under compression
        dzc = @(d) -1; % derivative of zc with respect to d
    case 2
        gammat = param(7); % damage parameter for concrete under tension
        gammac = param(8); % damage parameter for concrete under compression
        alphac = (1-gammac)*sig0c^2/((1-gammat)*sig0t^2); % dissymmetric behaviour between tension and compression
        
        zt = @(d) (1+gammat*d)./(1+d); % damage function under tension
        dzt = @(d) (gammat-1)./(1+d).^2; % derivative of zt with respect to d
        zc = @(d) (alphac+gammac*d)./(alphac+d); % damage function under compression
        dzc = @(d) alphac*(gammac-1)./(alphac+d).^2; % derivative of zc with respect to d
end

e0 = @(e,De) (ks+kc)*e+2*ks*De;
z = @(d,e,De) zt(d).*heaviside(e0(e,De)) + zc(d).*heaviside(-e0(e,De)); % damage function
dz = @(d,e,De) dzt(d).*heaviside(e0(e,De)) + dzc(d).*heaviside(-e0(e,De)); % derivative of z with respect to d

K = @(e,De,d) kc*z(d,e,De)/(1+z(d,e,De)); % equivalent stiffness for concrete alone
Ke = @(e,De,d) (ks+kc)*(ks+kc*z(d,e,De))./(2*ks+kc*(1+z(d,e,De))); % equivalent elastic stiffness
Kc = @(e,De,d) ks*kc*(z(d,e,De)-1)./(2*ks+kc*(1+z(d,e,De))); % equivalent coupling stiffness
KD = @(e,De,d) 2*ks*kc*(1+z(d,e,De))./(2*ks+kc*(1+z(d,e,De))); % equivalent debonding stiffness
dKe = @(e,De,d) kc*(ks+kc)^2*dz(d,e,De)./(2*ks+kc*(1+z(d,e,De))).^2; % derivative of Ke with respect to d
dKc = @(e,De,d) 2*ks*kc*(ks+kc)*dz(d,e,De)./(2*ks+kc*(1+z(d,e,De))).^2; % derivative of Kc with respect to d
dKD = @(e,De,d) 4*ks^2*kc*dz(d,e,De)./(2*ks+kc*(1+z(d,e,De))).^2; % derivative of KD with respect to d

W = @(e,De,d) 1/2*Ke(e,De,d).*e.^2 + Kc(e,De,d).*e.*De + 1/2*KD(e,De,d).*De.^2; % free energy
sig = @(e,De,d) Ke(e,De,d).*e + Kc(e,De,d).*De; % stress
Dsig = @(e,De,d) -(Kc(e,De,d).*e + KD(e,De,d).*De); % debonding stress
G = @(e,De,d) -(1/2*dKe(e,De,d).*e.^2 + dKc(e,De,d).*e.*De + 1/2*dKD(e,De,d).*De.^2); % energy restitution rate
% G = @(e,De,d) -kc/2*dz(d,e,De).*(((ks+kc).*e+2*ks*De)./(2*ks+kc*(1+z(d,e,De)))).^2; % energy restitution rate
Wdamage = @(e,De,d,d_old) G(e,De,d).*(d-d_old);
Wdebonding = @(e,De,d,De_old) Dsig(e,De,d).*(De-De_old);

ec = @(e,De,d) (sig(e,De,d)-ks*De)/(ks+kc); % strain for sound concrete
% ec = @(e,De,d) ((ks+kc*z(d,e,De))*e-2*ks*De)./(2*ks+kc*(1+z(d,e,De)));
es = @(e,De,d) ec(e,De,d)+De; % strain for steel
% es = @(e,De,d) (Ke(e,De,d)*e+(Kc(e,De,d)+kc)*De)/(ks+kc);

switch damageFun
    case 1
        % function of G for damage threshold function
        funG = @(e,De,d) (G(e,De,d).^alphaGt).*(1-d).^alphadt*heaviside(e0(e,De))...
                + 1/alphac*(G(e,De,d).^alphaGc).*(1-d).^alphadc*heaviside(-e0(e,De));
        % G0 = (sig0t^2/(2*kc))^alphaGt; % damage threshold value for G
        % G0 = ((sig0c^2/(2*kc))^alphaGc)/alphac;
    case 2
        % function of G for damage threshold function
        funG = @(e,De,d) G(e,De,d);
        % G0 = (1-gammat)*sig0t^2/(2*kc); % damage threshold value for G
        % G0 = (1-gammac)*sig0c^2/(2*kc*alphac);
end
G0 = funG(sig0t/K(0,0,0),0,0); % damage threshold value for G
funDamage = @(e,De,d) (funG(e,De,d)-G0)/G0; % damage threshold function
funDebonding = @(e,De,d) (Dsig(e,De,d).^2-Dsig0.^2)/Dsig0.^2; % debonding threshold function

if displayIter
    fprintf('\n+----------+----------+-----------+----------+----------+\n');
    fprintf('| e [1e-3] |     d    | De [1e-3] | Ds [MPa] |  s [MPa] |\n');
    fprintf('+----------+----------+-----------+----------+----------+\n');
end

damage = zeros(size(epsilon));
epsilonDebonding = zeros(size(epsilon));
sigma = zeros(size(epsilon));
sigmaDebonding = zeros(size(epsilon));
dissipatedEnergy = zeros(size(epsilon));
dissipatedEnergyDamage = zeros(size(epsilon));
dissipatedEnergyDebonding = zeros(size(epsilon));

sigma(1) = sig(epsilon(1),epsilonDebonding(1),damage(1));
for iter=2:length(epsilon)
    
    e = epsilon(iter);
    De = epsilonDebonding(iter-1);
    De_old = De;
    d = damage(iter-1);
    d_old = d;
    
    if funDamage(e,De,d)>tol || funDebonding(e,De,d)>tol % damage or debonding
        while true
            if funDamage(e,De,d)>tol % damage
                updateDamage = true;
                d = fsolve(@(d) funDamage(e,De,d),d,options);
            else
                updateDamage = false;
            end
            
            if funDebonding(e,De,d)>tol % debonding
                updateDebonding = true;
                De = fsolve(@(De) funDebonding(e,De,d),De,options);
            else
                updateDebonding = false;
            end
            
            if ~updateDebonding || ((updateDamage && abs(funDamage(e,De,d))<tol))
                break
            end
        end
        if d<d_old
            error('Damage cannot decrease');
        end
        if damageFun==1 && d>1
            break
        end
    end
    
    damage(iter) = d;
    epsilonDebonding(iter) = De;
    sigma(iter) = sig(e,De,d);
    sigmaDebonding(iter) = Dsig(e,De,d);
    dissipatedEnergyDamage(iter) = dissipatedEnergyDamage(iter-1) + Wdamage(e,De,d,d_old);
    dissipatedEnergyDebonding(iter) = dissipatedEnergyDebonding(iter-1) + Wdebonding(e,De,d,De_old);
    dissipatedEnergy(iter) = dissipatedEnergyDamage(iter) + dissipatedEnergyDebonding(iter);
    
    if displayIter
        fprintf('| %8g | %8g | %8g | %8g | %8g |\n',epsilon(iter),damage(iter),epsilonDebonding(iter)*1e3,sigmaDebonding(iter),sigma(iter));
    end
    
end

if displayIter
    fprintf('+----------+----------+----------+-----------+----------+\n');
end

end
