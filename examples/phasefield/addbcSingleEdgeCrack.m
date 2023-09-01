function S = addbcSingleEdgeCrack(S,ud,BU,BL,BLeft,BRight,BFront,BBack,loading)
% S = addbcSingleEdgeCrack(S,ud,BU,BL,BLeft,BRight,BFront,BBack,loading)
% Add boundary conditions on displacement field for single edge crack problem.

Dim = getdim(S);

switch lower(loading)
    case 'tension'
        if Dim==2
            S = addcl(S,BU,{'UX','UY'},[0;ud]);
        elseif Dim==3
            S = addcl(S,BU,{'UX','UY','UZ'},[0;ud;0]);
        end
        S = addcl(S,BL,'UY');
    case 'shear'
        if Dim==2
            S = addcl(S,BU,{'UX','UY'},[ud;0]);
            S = addcl(S,BLeft,'UY');
            S = addcl(S,BRight,'UY');
        elseif Dim==3
            S = addcl(S,BU,{'UX','UY','UZ'},[ud;0;0]);
            S = addcl(S,BLeft,{'UY','UZ'});
            S = addcl(S,BRight,{'UY','UZ'});
            S = addcl(S,BFront,{'UY','UZ'});
            S = addcl(S,BBack,{'UY','UZ'});
        end
        S = addcl(S,BL);
    otherwise
        error('Wrong loading case');
end

end