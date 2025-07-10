function S = addbcPlateWithHole(S,ud,BU,BL,P0)
% function S = addbcPlateWithHole(S,ud,BU,BL,P0)
% Add boundary conditions on displacement field for plate with hole problem.

% [Luo, Chen, Wang, Li, 2022, CM]
% S = addcl(S,BU,{'UX','UY'},[0;-ud]);
% S = addcl(S,BL);
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
S = addcl(S,BU,'UY',-ud);
S = addcl(S,BL,'UY');
S = addcl(S,P0,'UX');

end