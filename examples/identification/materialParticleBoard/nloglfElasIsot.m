function nloglf = nloglfElasIsot(lambda,C_data,cens,freq)
% function nloglf = nloglfElasIsot(lambda,C_data,cens,freq)

if isvector(C_data)
    n_data = length(C_data)/2;
    C1_data = C_data(1:n_data);
    C2_data = C_data(n_data+1:end);
    % data = C_data;
else
    n_data = size(C_data,1);
    C1_data = C_data(:,1);
    C2_data = C_data(:,2);
    % data = C_data(:); % data = [C1_data; C2_data];
end

useRedParam = isscalar(lambda); % reduced parameterization
if useRedParam
    % reduced parameterization
    la  = lambda(1);     % la < 1/5
    a1  = 1-la;          % a1 > 0
    a2  = 1-5*la;        % a2 > 0
    mC1_data = mean(C1_data,1);
    mC2_data = mean(C2_data,1);
    la1 = a1/mC1_data; % la1 > 0
    la2 = a2/mC2_data; % la2 > 0
else
    % full parameterization
    la1 = lambda(1); % la1 > 0
    la2 = lambda(2); % la2 > 0
    la  = lambda(3); % la < 1/5
end

% Compute negative log-likelihood
nloglf = -n_data*( (1-la)*log(la1) - gammaln(1-la) )...
    - n_data*( (1-5*la)*log(la2) - gammaln(1-5*la) )...
    + la*( sum(log(C1_data)) + 5*sum(log(C2_data)) )...
    + la1*sum(C1_data)...
    + la2*sum(C2_data);

end
