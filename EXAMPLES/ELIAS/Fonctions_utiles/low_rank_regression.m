function [U,beta] = low_rank_regression(Psi,V,tolKL)
% function [U,beta] = low_rank_regression(Psi,V,tolKL)

%% Construction of the space
[U,s,~] = svd(V);
s = diag(s);
err = sqrt(1-cumsum(s.^2)/sum(s.^2));
m = find(err<tolKL);
if isempty(m)
    m = size(V,2);
else
    m = min(m);
end

U = U(:,1:m);

%% Regression
beta = (Psi'*Psi)\(Psi'*V'*U);

end
