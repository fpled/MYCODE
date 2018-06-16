function [C_sample,lambda,lambda1,lambda2,lambda3,lambda4,lambda5] = mhsampleStoLinElasTensorIsotTrans_BUM(C1_data,C2_data,C3_data,C4_data,C5_data,N)
% function [C_sample,lambda,lambda1,lambda2,lambda3,lambda4,lambda5] = mhsampleStoLinElasTensorIsotTrans_BUM(C1_data,C2_data,C3_data,C4_data,C5_data,N)
% Metropolis-Hastings Sampling using Blockwise Updating Method (BUM)
% for stochastic linear elastic tensor with transversely isotropic symmetry
% This method sometimes doesn't work, because the positive-definite of
% covariance matrix SIGMA isn't always guaranteed.
% If there is an error, rerun the code.
% See page 26 in 'Computational Statistics with Matlab'
% See also page 4 in 'Metropolis-Hastings sampling using multivariate gaussian tangents'
% C1_data: data for random coordinate C1
% C2_data: data for random coordinate C2
% C3_data: data for random coordinate C3
% C4_data: data for random coordinate C4
% C5_data: data for random coordinate C5
% N: number of samples
% C_sample: sample set for random vector C=(C1,C2,C3,C4,C5)

mc1 = mean(C1_data);
mc2 = mean(C2_data);
mc3 = mean(C3_data);
mc4 = mean(C4_data);
mc5 = mean(C5_data);

%% Sample generation
% Method to calculate lambda which controls the level of fluctuations
lambda = -40; % negative number
lambda1 = -(mc2*lambda)/(-mc3^2+mc1*mc2);
lambda2 = -(mc1*lambda)/(-mc3^2+mc1*mc2);
lambda3 = (2*mc3*lambda)/(-mc3^2+mc1*mc2);
a = 1-2*lambda;
lambda4 = a/mc4;
lambda5 = a/mc5;
b4 = 1/lambda4;
b5 = 1/lambda5;

% Parameters of the trivariate normal distribution
Mu_mean = [mc1 mc2 mc3];
Sigma_mean = [-mc1^2/lambda -mc3^2/lambda -(mc1*mc3)/lambda
              -mc3^2/lambda -mc2^2/lambda -(mc2*mc3)/lambda
              -(mc1*mc3)/lambda -(mc2*mc3)/lambda -(mc3^2+mc1*mc2)/(2*lambda)];
f = @(c) -lambda1*c(1)-lambda2*c(2)-lambda3*c(3)-lambda*log(c(1)*c(2)-c(3)^2);

C_sample = zeros(5,N);
n = 0;
while n<N
    n = n+1;
    c_old = mvnrnd(Mu_mean,Sigma_mean);
    c1_old = c_old(1);
    c2_old = c_old(2);
    c3_old = c_old(3);
    f_old = f(c_old);
    g_old = [-lambda1-(c2_old*lambda)/(-c3_old^2+c1_old*c2_old)
             -lambda2-(c1_old*lambda)/(-c3_old^2+c1_old*c2_old)
             -lambda3+(2*c3_old*lambda)/(-c3_old^2+c1_old*c2_old)];
    sigma_old = [-c1_old^2/lambda -c3_old^2/lambda -(c1_old*c3_old)/lambda
                 -c3_old^2/lambda -c2_old^2/lambda -(c2_old*c3_old)/lambda
                 -(c1_old*c3_old)/lambda -(c2_old*c3_old)/lambda -(c3_old^2+c1_old*c2_old)/(2*lambda)];
    mu_old = (c_old' + sigma_old*g_old)';
    q_old = mvnpdf(c_old,mu_old,sigma_old);
    %
    c_prop = mvnrnd(mu_old,sigma_old);
    c1_prop = c_prop(1);
    c2_prop = c_prop(2);
    c3_prop = c_prop(3);
    log_q_prop = log(q_old);
    f_prop = f(c_prop);
    g_prop = [-lambda1-(c2_prop*lambda)/(-c3_prop^2+c1_prop*c2_prop)
              -lambda2-(c1_prop*lambda)/(-c3_prop^2+c1_prop*c2_prop)
              -lambda3+(2*c3_prop*lambda)/(-c3_prop^2+c1_prop*c2_prop)];
    sigma_prop = [-c1_prop^2/lambda -c3_prop^2/lambda -(c1_prop*c3_prop)/lambda
                  -c3_prop^2/lambda -c2_prop^2/lambda -(c2_prop*c3_prop)/lambda
                  -(c1_prop*c3_prop)/lambda -(c2_prop*c3_prop)/lambda -(c3_prop^2+c1_prop*c2_prop)/(2*lambda)];
    mu_prop = (c_prop' + sigma_prop*g_prop)';
    q_prop = mvnpdf(c_prop,mu_prop,sigma_prop);
    log_q_old = log(q_prop);
    r = exp( (f_prop-f_old)+(log_q_old-log_q_prop) );
    if r>1
        C_sample(1:3,n) = c_prop;
    else
        s=rand;
        if s<r
            C_sample(1:3,n) = c_prop;
        else
            C_sample(1:3,n) = c_old;
        end
    end
end
C_sample(4,:) = gamrnd(a,b4,N,1);
C_sample(5,:) = gamrnd(a,b5,N,1);
lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda];

end
