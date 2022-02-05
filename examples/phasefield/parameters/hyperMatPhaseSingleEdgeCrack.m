function [aP1,bP1,aP2,bP2] = hyperMatPhaseSingleEdgeCrack(gc,l,deltaP1,deltaP2)
% function [aP1,bP1,aP2,bP2] = hyperMatPhaseSingleEdgeCrack(gc,l,deltaP1,deltaP2)
%
% [aP1,bP1,aP2,bP2] = hyperMatPhaseSingleEdgeCrack(gc,l,deltaP1,deltaP2)
% Returns the shape and scale parameters of the Gamma distribution
% (resp. 'a' and 'b') followed by random fracture toughness and
% regularization length (resp. 'P1' and 'P2').
%
% <strong>Inputs:</strong>
% 'gc' is the mean fracture toughness.
% 'l' is the mean regularization length.
% 'deltaP1' is the coefficient of variation of the random fracture toughness.
% 'deltaP2' is the coefficient of variation of the random regularization length.

%% Inputs parsing
arguments
    % Mandatory positional arguments (no default value)
    gc {mustBeReal, mustBePositive}
    l {mustBeReal, mustBePositive}
    deltaP1 {mustBeReal, mustBeNonnegative}
    deltaP2 {mustBeReal, mustBeNonnegative}
end

%% Computation

aP1 = 1/deltaP1^2;
bP1 = gc/aP1;
aP2 = 1/deltaP2^2;
bP2 = l/aP2;