function R = autocorrVec(V,ID)
% function R = autocorrVec(V,ID)
%
% R = autocorrVec(V,ID) returns the columns of the autocorrelation matrix
% of V specified by the indices in the vector ID. V is a matrix which
% columns represent realizations of a random vector.
% The full matrix is not computed, only the specified columns are.
%
% R = autocorrVec(V) returns the entire autocorrelation matrix.

arguments
    % Mandatory positional arguments (no default value)
    V (:, :) {mustBeNumeric, mustBeReal}
    % Optional positional arguments (given a default value)
    ID (1,:) {mustBeInteger, mustBePositive} = 1:size(V,1)
end

[nx,N] = size(V);
meanV = mean(V,2);

% Compute the columns of interest in the covariance matrix
covVec = V*V(ID,:)' - N*meanV*meanV(ID)';

% Compute the diagonal elements of the covariance matrix
covDiag = zeros(nx,1);
for kk = 1:nx
    covDiag(kk) = V(kk,:)*V(kk,:)';
end
covDiag = covDiag - N*meanV.^2;

% Compute the columns of interest in the autocorrelation matrix
R = covVec./sqrt(covDiag*covDiag(ID)');