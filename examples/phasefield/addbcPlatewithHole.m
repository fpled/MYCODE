function S = addbcPlatewithHole(S,ud,BU,BL,P0)
% function S = addbcPlatewithHole(S,ud,BU,BL,P0)
% Add boundary conditions on displacement field for plate with hole problem.

S = addcl(S,BU,'UY',-ud);
S = addcl(S,BL,'UY');
S = addcl(S,P0,'UX');

end