function lambda = mleStoLinElasTensorIsot(C_data)
% function lambda = mleStoLinElasTensorIsot(C_data)
% Maximum likelihood estimation for stochastic linear elastic tensor with
% isotropic symmetry
% C_data(:,1): data for random coordinate C1
% C_data(:,2): data for random coordinate C2
% lambda: parameters (la1,la2,la)
% lambda(1) = la1 > 0
% lambda(2) = la2 > 0
% lambda(3) = la < 1/5

C1_data = C_data(:,1);
C2_data = C_data(:,2);
n1_data = length(C1_data);
n2_data = length(C2_data);
data = C_data(:); % data = [C1_data; C2_data];

mC = mean(C_data,1);

% initial guess
la  = 0; % la < 1/5
la1 = (1-la)/mC(1); % la1 > 0
la2 = (1-5*la)/mC(2); % la2 > 0

lambda0 = [la1 la2 la];
lb = [0 0 -Inf];
ub = [Inf Inf 1/5];

nloglf = @(lambda,data,cens,freq) -n1_data*( (1-lambda(3))*log(lambda(1)) - gammaln(1-lambda(3)) )...
    - n2_data*( (1-5*lambda(3))*log(lambda(2)) - gammaln(1-5*lambda(3)) )...
    + lambda(3)*( sum(log(data(1:n1_data))) + 5*sum(log(data(n1_data+1:end))) )...
    + lambda(1)*sum(data(1:n1_data))...
    + lambda(2)*sum(data(n1_data+1:end));
lambda = mle(data,'nloglf',nloglf,'Start',lambda0,...
    'LowerBound',lb,'UpperBound',ub);

end

