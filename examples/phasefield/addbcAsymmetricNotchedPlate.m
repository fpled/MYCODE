function S = addbcAsymmetricNotchedPlate(S,ud,PU,PL,PR)
% function S = addbcAsymmetricNotchedPlate(S,ud,PU,PL,PR)
% Add boundary conditions on displacement field for asymmetric notched plate problem.

S = addcl(S,PU,'UY',-ud);
S = addcl(S,PL,{'UX','UY'});
S = addcl(S,PR,'UY');

end