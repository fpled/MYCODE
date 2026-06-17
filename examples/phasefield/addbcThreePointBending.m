function S = addbcThreePointBending(S,ud,BU,BL,BR,P0)
% function S = addbcThreePointBending(S,ud,BU,BL,BR,P0)
% Add boundary conditions on displacement field for three-point bending problem.

S = addcl(S,BU,'UY',-ud);
S = addcl(S,BL,'UY');
S = addcl(S,BR,'UY');
S = addcl(S,P0,'UX');

end