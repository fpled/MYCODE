function f = funlsqnonlin(param,uy_exp,S,N,M,b,funu0,v0,P1,varargin)
% function f = funlsqnonlin(param,uyexp,S,N,M,b,funu0,v0,P1,varargin)

if ischarin('display',varargin)
    fprintf('E     = %g GPa\n',param(1));
    fprintf('alpha = %g\n',param(2));
    fprintf('beta  = %g\n',param(3));
    if length(param)==5
        fprintf('c     = %g kN.m/rad\n',param(4));
        fprintf('J     = %g kg.m2/rad\n',param(5));
    end
    fprintf('\n');
end

% [ut,result,vt,at] = solveBeamDetDynLinElasClampedFree(param,S,N,M,b,funu0,v0,P1,varargin{:});
ut = solveBeamDetDynLinElasClampedFree(param,S,N,M,b,funu0,v0,P1,varargin{:});

ut = unfreevector(S,ut);
ut_val = getvalue(ut);
% Ut = ut_val(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
% Uxt = ut_val(findddl(S,'UX'),:);
Uyt = ut_val(findddl(S,'UY'),:);
% Rzt = ut_val(findddl(S,'RZ'),:);

uy = Uyt(end,:);

f = uy(:) - uy_exp(:);

end