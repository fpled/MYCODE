function S = addbcPlateWithHole(S,ud,BU,BL,P1,P2)
% function S = addbcPlateWithHole(S,ud,BU,BL,P1,P2)
% Add boundary conditions on displacement field for plate with hole problem.

Dim = getdim(S);

% [Luo, Chen, Wang, Li, 2022, CM]
% if Dim==2
%     S = addcl(S,BU,{'UX','UY'},[0;-ud]);
% elseif Dim==3
%     S = addcl(S,BU,{'UX','UY','UZ'},[0;-ud;0]);
% end
% S = addcl(S,BL);

if Dim==2
    % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Vu, Le Quang, He, 2024, AES]
    S = addcl(S,BU,'UY',-ud);
    S = addcl(S,BL,'UY');
    S = addcl(S,P1,'UX');
elseif Dim==3
    % [Nguyen, Yvonnet, Bornert, Chateau, Sab, Romani, Le Roy, 2016, IJF]
    S = addcl(S,BU,'UY',-ud);
    S = addcl(S,BL,'UY');
    S = addcl(S,P1,{'UX','UZ'});
    S = addcl(S,P2,'UZ');
end

end