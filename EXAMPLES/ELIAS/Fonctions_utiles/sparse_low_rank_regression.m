function [U,beta] = sparse_low_rank_regression(Psi,V,tolKL,k)
% function [U,beta] = sparse_low_rank_regression(Psi,V,tolKL,k)

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

%% LASSO
VU = V'*U;

beta = zeros(size(Psi,2),m);
for i=1:m
    fprintf('    LASSO : rank %d over %d\n',i,m);
    % beta(:,i) = l1qc_logbarrier(zeros(size(Psi,2),1), Psi, [], VU(:,i), 1e-10, 1e-3);
    
    lambdamax = find_lambdamax_l1_ls(Psi',VU(:,i));
    lambda = k*lambdamax;
    beta(:,i) = l1_ls(Psi,VU(:,i),lambda,1e-12,true);
end

end
