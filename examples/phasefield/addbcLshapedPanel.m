function S = addbcLshapedPanel(S,ud,BL,BRight)
% function S = addbcLshapedPanel(S,ud,BL,BRright)
% Add boundary conditions on displacement field for L-shaped panel problem.

Dim = getdim(S);

% [Gerasimov, De Lorenzis, 2019, CMAME]
% if Dim==2
%     S = addcl(S,BRight,{'UX','UY'},[0;ud]);
% elseif Dim==3
%     S = addcl(S,BRight,{'UX','UY','UZ'},[0;ud;0]);
% end
S = addcl(S,BRight,'UY',ud);
S = addcl(S,BL);

end