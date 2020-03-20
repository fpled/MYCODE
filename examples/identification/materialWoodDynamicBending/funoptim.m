function f = funoptim(param,uy_exp,S,N,M,b,u0,v0,varargin)
% function f = funoptim(param,uy_exp,S,N,M,b,u0,v0,varargin)

% [ut,result,vt,at] = solveBeamDetDynLinElasClampedFree(param,S,N,M,b,u0,v0,varargin{:});
ut = solveBeamDetDynLinElasClampedFree(param,S,N,M,b,u0,v0,varargin{:});

ut = unfreevector(S,ut);
ut_val = getvalue(ut);
% Ut = ut_val(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
% Uxt = ut_val(findddl(S,'UX'),:);
Uyt = ut_val(findddl(S,'UY'),:);
% Rzt = ut_val(findddl(S,'RZ'),:);

uy = Uyt(end,:);

f = norm(uy - uy_exp(:))^2;

end