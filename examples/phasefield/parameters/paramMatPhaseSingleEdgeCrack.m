function [gc,l] = paramMatPhaseSingleEdgeCrack(symmetry,test)
% function [gc,l] = paramMatPhaseSingleEdgeCrack(symmetry,test)
%
% [gc,l] = paramMatPhaseSingleEdgeCrack(symmetry,test)
% Returns the fracture toughness 'gc' and regularization length 'l' given
% the material 'symmetry' and whether it is a 'test' case or not.

switch lower(symmetry)
    case 'isotropic' % isotropic material
        % Critical energy release rate (or fracture toughness)
        gc = 2.7e3; % [Miehe, Hofacker, Welschinger, 2010, CMAME]
        % Regularization parameter (width of the smeared crack)
        % l = 3.75e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME]
        % l = 3e-5; % [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
        % l = 1.5e-5; % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Hesch, Weinberg, 2014, IJNME], [Liu, Li, Msekh, Zuo, 2016, CMS], [Molnar, Gravouil, 2017, FEAD], [Zhou, Rabczuk, Zhuang, 2018, AES], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Hu, Guilleminot, Dolbow, 2020, CMAME]
        l = 1e-5; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Gerasimov, De Lorenzis, 2016, CMAME], [Wu, Nguyen, 2018, JMPS], [Gerasimov, De Lorenzis, 2019, CMAME], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        % l = 7.5e-6; % [Miehe, Welschinger, Hofacker, 2010, IJNME], [Miehe, Hofacker, Welschinger, 2010, CMAME], [Borden, Verhoosel, Scott, Hughes, Landis, 2012, CMAME], [Hesch, Weinberg, 2014, IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM], [Liu, Li, Msekh, Zuo, 2016, CMS], [Molnar, Gravouil, 2017, FEAD], [Zhou, Rabczuk, Zhuang, 2018, AES], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME], [Storvik, Both, Sargado, Nordbotten, Radu, 2021, CMAME]
        % l = 5e-6; % [Molnar, Gravouil, 2017, FEAD], [Wu, Nguyen, 2018, JMPS], [Wu, Nguyen, Zhou, Huang, 2020, CMAME]
        % l = 4e-6; % [Ambati, Gerasimov, De Lorenzis, 2015, CM]
        % eta = 0.052; w0 = 75.94; l = eta/sqrt(w0)*1e-3; % l = 6e-7; % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
    case 'anisotropic' % anisotropic material
        % Critical energy release rate (or fracture toughness)
        gc = 1e3; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
        % Regularization parameter (width of the smeared crack)
        l = 8.5e-6; % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
    otherwise
        error('Wrong material symmetry class');
end