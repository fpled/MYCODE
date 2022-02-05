function [aC1,bC1,aC2,bC2] = hyperMatIsot(Dim,lambda,mu,deltaC1,options)
% function [aC1,bC1,aC2,bC2] = hyperMatIsot(Dim,lambda,mu,deltaC1,options)
%
% [aC1,bC1,aC2,bC2] = hyperMatIsot(Dim,lambda,mu,deltaC1)
% Returns the shape and scale parameters of the Gamma distribution
% (resp. 'a' and 'b') followed by random bulk and shear moduli
% (resp. 'C1' and 'C2').
%
% <strong>Inputs:</strong>
% 'Dim' is the spacial dimension.
% 'lambda' and 'mu' are the mean Lame coefficients.
% 'deltaC1' is the coefficient of variation of the random bulk modulus.
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
    deltaC1 {mustBeReal, mustBePositive}
    % Optional name-value arguments
    options.option {mustBeMember(options.option, {'DEFO' 'CONT'})} = 'DEFO'
end

%% Computation
switch Dim
    case 2
        switch lower(options.option)
            case 'defo'
                E = mu*(3*lambda+2*mu)/(lambda+mu); %  E = 210e9;
                NU = lambda/(lambda+mu)/2; % NU = 0.3;
            case 'cont'
                E = 4*mu*(lambda+mu)/(lambda+2*mu);
                NU = lambda/(lambda+2*mu);
        end
    case 3
        E = mu*(3*lambda+2*mu)/(lambda+mu);
        NU = lambda/(lambda+mu)/2;
end

la = 1 - 1/deltaC1^2; % la < 1/5. Parameter controlling the level of statistical fluctuations
% la = (1 - 1/deltaC2^2)/5;

mC1 = E/3/(1-2*NU); % mean bulk modulus
mC2 = mu; % mean shear modulus
laC1 = (1-la)/mC1;
laC2 = (1-5*la)/mC2;

aC1 = 1-la;
bC1 = 1/laC1;
aC2 = 1-5*la;
bC2 = 1/laC2;