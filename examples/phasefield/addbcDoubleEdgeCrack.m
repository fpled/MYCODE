function S = addbcDoubleEdgeCrack(S,ud,BU,BL,P0,setup)
% S = addbcDoubleEdgeCrack(S,ud,BU,BL,P0,setup)
% Add boundary conditions on displacement field for double edge crack problem (with two symmetric or asymmetric edge cracks).

Dim = getdim(S);

switch setup
    case 1
        if Dim==2
            S = addcl(S,BU,{'UX','UY'},[0;ud]);
        elseif Dim==3
            S = addcl(S,BU,{'UX','UY','UZ'},[0;ud;0]);
        end
        S = addcl(S,BL,'UY');
    case 2
        % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
        S = addcl(S,BU,'UY',ud);
        S = addcl(S,BL,'UY');
        if Dim==2
            S = addcl(S,P0,'UX');
        elseif Dim==3
            S = addcl(S,P0,{'UX','UZ'});
        end
    case 3
        % [Ribeiro Nogueira, Rastiello, Giry, Gatuingt, Callari, 2023, AJCE]
        % S = addcl(S,BU,'UY',ud);
        % S = addcl(S,BL,'UY');
        % if Dim==2
        %     S = addcl(S,P0,'UX');
        % elseif Dim==3
        %     S = addcl(S,P0,{'UX','UZ'});
        % end
        
        % [Tong, Shen, Shao, Chen, 2020, EFM], [Li, Lu, Huang, Yang, 2022, OE]
        % S = addcl(S,BU,'UY',ud);
        % S = addcl(S,BL);
        % if Dim==2
        %     S = addcl(S,P0,'UX');
        % elseif Dim==3
        %     S = addcl(S,P0,{'UX','UZ'});
        % end
        
        % [Shi, van Dam, van Mier, Sluys, 2000, MBS], [Alfaiate, Wells, Sluys, 2002, EFM], [Nguyen, 2005, PhD thesis], [Nguyen, Houlsby, 2007, IJNAMG], [Nguyen, Houlsby, 2008, IJNAMG], [Nguyen, 2008, IJSS], [Nguyen, Korsunsky, 2008, IJSS], [Nguyen, 2011, IJSS], [Galvez, Planas, Sancho, Reyes, Cendon, Casati, 2013, EFM], [Stefanou, Georgioudakis, Papadrakakis, 2014, MMUQMS], [Le, Nguyen, Bui, Sheikh, Kotousov, 2018, IJES], [Fang, Wu, Rabczuk, Wu, Sun, Li, 2020, CM], [Han, Li, Yu, Li, Zhang, 2022, JMPS], [Liu, Chen, Yuan, 2024, AAM]
        if Dim==2
            S = addcl(S,BU,{'UX','UY'},[0;ud]);
        elseif Dim==3
            S = addcl(S,BU,{'UX','UY','UZ'},[0;ud;0]);
        end
        S = addcl(S,BL);
    otherwise
        error('Wrong setup');
end

end