function [E,NU] = paramMatIsot(Dim,lambda,mu,options)
% function [E,NU] = paramMatIsot(Dim,lambda,mu,options)
%
% [E,NU] = paramMatIsot(Dim,lambda,mu)
% Returns the Young modulus 'E' and Poisson coefficient 'NU' given the
% spatial dimension 'Dim' and the 3D Lame coefficients 'lambda' and 'mu'.
% All above properties can be double, MYDOUBLEND or FEELEMFIELD.
%
% [...] = hyperMatIsot(...,'PARAM1',val1,'PARAM2',val2,...) specifies
% parameter name/value pairs.
% Valid parameters are the following:
%
%  Parameter      Value
%  'option'       'DEFO' for plane strain or 'CONT' for plane stress in 2D.
%                 It is set to 'DEFO' by default

%% Inputs parsing
arguments
    % Mandatory positional arguments (no default value)
    Dim (1, 1) {mustBeInteger}
    lambda {mustBeReal, mustBePositive}
    mu {mustBeReal, mustBePositive}
    options.option {mustBeMember(options.option, {'DEFO' 'CONT'})} = 'DEFO'
end


if Dim==2
    switch lower(options.option)
        case 'defo'
            E = mu.*(3*lambda+2*mu)./(lambda+mu); %  E = 210e9;
            NU = lambda./(lambda+mu)/2; % NU = 0.3;
        case 'cont'
            E = 4*mu.*(lambda+mu)./(lambda+2*mu);
            NU = lambda./(lambda+2*mu);
    end
    % E = 210e9; NU = 0.2; % [Liu, Li, Msekh, Zuo, 2016, CMS], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM]
    % E = 210e9; NU = 0.3; % [Gerasimov, De Lorenzis, 2016, CMAME], [Molnar, Gravouil, 2017, FEAD], [Zhou, Rabczuk, Zhuang, 2018, AES], [Wu, Nguyen, 2018, JMPS], [Gerasimov, De Lorenzis, 2019, CMAME], [Wu, Nguyen, Zhou, Huang, 2020, CMAME], [Kristensen, Martinez-Paneda, 2020, TAFM]
    % kappa = 121030e6; NU=0.227; lambda=3*kappa*NU/(1+NU); mu = 3*kappa*(1-2*NU)/(2*(1+NU)); E = 3*kappa*(1-2*NU); % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
elseif Dim==3
    E = mu.*(3*lambda+2*mu)./(lambda+mu);
    NU = lambda./(lambda+mu)/2;
end