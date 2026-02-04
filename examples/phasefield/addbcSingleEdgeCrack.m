function S = addbcSingleEdgeCrack(S,ud,BU,BL,BLeft,BRight,BFront,BBack,loading)
% S = addbcSingleEdgeCrack(S,ud,BU,BL,BLeft,BRight,BFront,BBack,loading)
% Add boundary conditions on displacement field for single edge crack problem.

Dim = getdim(S);

switch lower(loading)
    case 'tension'
        if Dim==2
            S = addcl(S,BU,{'UX','UY'},[0;ud]);
            % [Liu, Li, Msekh, Zuo, 2016, CMS], [Zhou, Rabczuk, Zhuang, 2018, AES]
            % S = addcl(S,BU,{'UX','UY'},[0;ud]);
            % S = addcl(S,BLeft,'UX');
            % S = addcl(S,BRight,'UX');
        elseif Dim==3
            S = addcl(S,BU,{'UX','UY','UZ'},[0;ud;0]);
            % [Liu, Li, Msekh, Zuo, 2016, CMS]
            % S = addcl(S,BU,{'UX','UY','UZ'},[0;ud;0]);
            % S = addcl(S,BLeft,{'UX','UZ'});
            % S = addcl(S,BRight,{'UX','UZ'});
            % S = addcl(S,BL,'UZ');
            
            % [Zhou, Rabczuk, Zhuang, 2018, AES]
            % S = addcl(S,BU,{'UX','UY','UZ'},[0;ud;0]);
            % S = addcl(S,BLeft,'UX');
            % S = addcl(S,BRight,'UX');
            % S = addcl(S,BFront,'UZ');
            % S = addcl(S,BBack,'UZ');
        end
        S = addcl(S,BL,'UY');
        % [Yu, Hou, Zheng, Xiao, Zhao, 2024, CM]
        % S = addcl(S,BL);
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
        % [Bourdin, Francfort, Marigo, 2000, JMPS], [Gmati, Mareau, Ammar, El_Arem, 2020, IJNME]
        % S = addcl(S,BL,{'UX','UY'},[-ud;0]);
    otherwise
        error('Wrong loading case');
end

end