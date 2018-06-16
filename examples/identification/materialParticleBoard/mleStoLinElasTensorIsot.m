function lambda = mleStoLinElasTensorIsot(C1_data,C2_data)
% function lambda = mleStoLinElasTensorIsot(C1_data,C2_data)
% Maximum likelihood estimation for stochastic linear elastic tensor with
% isotropic symmetry
% C1_data: data for random coordinate C1
% C2_data: data for random coordinate C2
% lambda: parameters (lambda1,lambda2,lambda)
% lambda(1) = lambda_1 > 0
% lambda(2) = lambda_2 > 0
% lambda(3) = lambda < 1/5

n1_data = length(C1_data);
n2_data = length(C2_data);
data = [C1_data; C2_data];

nloglf = @(lambda,data,cens,freq) -n1_data*( (1-lambda(3))*log(lambda(1))...
    - gammaln(1-lambda(3)) ) - n2_data*( (1-5*lambda(3))*log(lambda(2))...
    - gammaln(1-5*lambda(3)) ) + lambda(3)*( sum(log(data(1:n1_data)))...
    + 5*sum(log(data(n1_data+1:end))) ) + lambda(1)*sum(data(1:n1_data))...
    + lambda(2)*sum(data(n1_data+1:end));
lambda = mle(data,'nloglf',nloglf,'start',[1 1 0],...
    'lowerbound',[0 0 -Inf],'upperbound',[Inf Inf 1/5]);

end

