function C_sample = mhsampleStoLinElasTensorIsotTrans_BUM(lambda,C_data,N)
% function C_sample = mhsampleStoLinElasTensorIsotTrans_BUM(lambda,C_data,N)
% Metropolis-Hastings Sampling using Blockwise Updating Method (BUM)
% for stochastic linear elastic tensor with transversely isotropic symmetry
% This method sometimes doesn't work, because the positive-definite of
% covariance matrix SIGMA isn't always guaranteed.
% If there is an error, rerun the code.
% See page 26 in 'Computational Statistics with Matlab'
% See also page 4 in 'Metropolis-Hastings sampling using multivariate gaussian tangents'
% 
% lambda = (la1,la2,la3,la4,la5,la)
% C_data: data set for random vector C=(C1,C2,C3)
% C_data(:,i): data for random coordinate Ci
% N: number of samples
% C_sample: sample set for random vector C=(C1,C2,C3)
% C_sample(:,i): data for random coordinate Ci

la1 = lambda(1);
la2 = lambda(2);
la3 = lambda(3);
% la4 = lambda(4);
% la5 = lambda(5);
la  = lambda(6);

mC = mean(C_data(:,1:3),1);

%% Sample generation
% Parameters of the trivariate normal distribution
Mu = mC; % Mu = [mC(1) mC(2) mC(3)];
Sigma = [-mC(1)^2/la -mC(3)^2/la -(mC(1)*mC(3))/la
         -mC(3)^2/la -mC(2)^2/la -(mC(2)*mC(3))/la
         -(mC(1)*mC(3))/la -(mC(2)*mC(3))/la -(mC(3)^2+mC(1)*mC(2))/(2*la)];
f = @(c) -la1*c(1)-la2*c(2)-la3*c(3)-la*log(c(1)*c(2)-c(3)^2);

C_sample = zeros(N,5);
n = 0;
while n<N
    n = n+1;
    c_old = mvnrnd(Mu,Sigma);
    f_old = f(c_old);
    g_old = [-la1-(c_old(2)*la)/(-c_old(3)^2+c_old(1)*c_old(2))
             -la2-(c_old(1)*la)/(-c_old(3)^2+c_old(1)*c_old(2))
             -la3+(2*c_old(3)*la)/(-c_old(3)^2+c_old(1)*c_old(2))];
    Sigma_old = [-c_old(1)^2/la -c_old(3)^2/la -(c_old(1)*c_old(3))/la
                 -c_old(3)^2/la -c_old(2)^2/la -(c_old(2)*c_old(3))/la
                 -(c_old(1)*c_old(3))/la -(c_old(2)*c_old(3))/la -(c_old(3)^2+c_old(1)*c_old(2))/(2*la)];
    Mu_old = (c_old' + Sigma_old*g_old)';
    q_old = mvnpdf(c_old,Mu_old,Sigma_old);
    %
    c_prop = mvnrnd(Mu_old,Sigma_old);
    log_q_prop = log(q_old);
    f_prop = f(c_prop);
    g_prop = [-la1-(c_prop(2)*la)/(-c_prop(3)^2+c_prop(1)*c_prop(2))
              -la2-(c_prop(1)*la)/(-c_prop(3)^2+c_prop(1)*c_prop(2))
              -la3+(2*c_prop(3)*la)/(-c_prop(3)^2+c_prop(1)*c_prop(2))];
    sigma_prop = [-c_prop(1)^2/la -c_prop(3)^2/la -(c_prop(1)*c_prop(3))/la
                  -c_prop(3)^2/la -c_prop(2)^2/la -(c_prop(2)*c_prop(3))/la
                  -(c_prop(1)*c_prop(3))/la -(c_prop(2)*c_prop(3))/la -(c_prop(3)^2+c_prop(1)*c_prop(2))/(2*la)];
    mu_prop = (c_prop' + sigma_prop*g_prop)';
    q_prop = mvnpdf(c_prop,mu_prop,sigma_prop);
    log_q_old = log(q_prop);
    r = exp( (f_prop-f_old)+(log_q_old-log_q_prop) );
    if r>1 && C_sample(n,1)>0 && C_sample(n,2)>0 && (C_sample(n,1)*C_sample(n,2)-C_sample(n,3)^2>0)
        C_sample(n,1:3) = c_prop;
    else
        s=rand;
        if s<r && C_sample(n,1)>0 && C_sample(n,2)>0 && (C_sample(n,1)*C_sample(n,2)-C_sample(n,3)^2>0)
            C_sample(n,1:3) = c_prop;
        else
            C_sample(n,1:3) = c_old;
        end
    end
end

end
