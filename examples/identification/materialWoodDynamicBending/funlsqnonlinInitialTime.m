function f = funlsqnonlinInitialTime(param,funuy_exp,S,funN,M,funb,funu0,v0,P1,varargin)
% function f = funlsqnonlinInitialTime(param,funuy_exp,S,funN,M,funb,funu0,v0,P1,varargin)

if ischarin('display',varargin)
    fprintf('tinit = %g s\n',param(1));
    fprintf('E     = %g GPa\n',param(2));
    fprintf('alpha = %g\n',param(3));
    fprintf('beta  = %g\n',param(4));
    if length(param)==6
        fprintf('c     = %g kN.m/rad\n',param(5));
        fprintf('J     = %g kg.m2/rad\n',param(6));
    end
    fprintf('\n');
end

% [ut,result,vt,at] = solveBeamDetDynLinElasClampedFreeInitialTime(param,S,funN,M,funb,funu0,v0,P1,varargin{:});
ut = solveBeamDetDynLinElasClampedFreeInitialTime(param,S,funN,M,funb,funu0,v0,P1,varargin{:});

ut = unfreevector(S,ut);
ut_val = getvalue(ut);
% Ut = ut_val(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
% Uxt = ut_val(findddl(S,'UX'),:);
Uyt = ut_val(findddl(S,'UY'),:);
% Rzt = ut_val(findddl(S,'RZ'),:);

uy = Uyt(end,:);

tinit = param(1); % [s]
uy_exp = funuy_exp(tinit); % [m]

f = uy(:) - uy_exp(:);

end