function f = funlsqnonlin(param,uy_exp,S,N,M,b,u0,v0,varargin)
% function f = funlsqnonlin(param,uyexp,S,N,M,b,u0,v0,varargin)

fprintf('E  = %g GPa\n',param(1));
fprintf('alpha = %g\n',param(2));
fprintf('beta = %g\n',param(3));
fprintf('\n');

% [ut,result,vt,at] = solveBeamDetDynLinElasClampedFree(param,S,t,M,b,funu0,v0,varargin{:});
ut = solveBeamDetDynLinElasClampedFree(param,S,N,M,b,u0,v0,varargin{:});

ut = unfreevector(S,ut);
ut_val = getvalue(ut);
% Ut = ut_val(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
% Uxt = ut_val(findddl(S,'UX'),:);
Uyt = ut_val(findddl(S,'UY'),:);
% Rzt = ut_val(findddl(S,'RZ'),:);

uy = Uyt(end,:);

f = uy(:) - uy_exp(:);

end