function C_sample = mhsampleStoLinElasTensorIsotTrans_BUM(param,C_data,N)
% function C_sample = mhsampleStoLinElasTensorIsotTrans_BUM(param,C_data,N)
% Metropolis-Hastings Sampling using Blockwise Updating Method (BUM)
% for stochastic linear elastic tensor with transversely isotropic symmetry
% This method sometimes doesn't work, because the positive-definite of
% covariance matrix SIGMA isn't always guaranteed.
% If there is an error, rerun the code.
% See page 26 in 'Computational Statistics with Matlab'
% See also page 4 in 'Metropolis-Hastings sampling using multivariate gaussian tangents'
% param = (lambda_1,lambda_2,lambda_3,lambda_4,lambda_5,lambda)
% C_data(:,1): data for random coordinate C1
% C_data(:,2): data for random coordinate C2
% C_data(:,3): data for random coordinate C3
% C_data(:,4): data for random coordinate C4
% C_data(:,5): data for random coordinate C5
% N: number of samples
% C_sample: sample set for random vector C=(C1,C2,C3,C4,C5)

lambda1 = param(1);
lambda2 = param(2);
lambda3 = param(3);
lambda4 = param(4);
lambda5 = param(5);
lambda  = param(6);

a = 1-2*lambda;
b4 = 1/lambda4;
b5 = 1/lambda5;

mc1 = mean(C_data(:,1));
mc2 = mean(C_data(:,2));
mc3 = mean(C_data(:,3));

%% Sample generation
% Parameters of the trivariate normal distribution
Mu = [mc1 mc2 mc3];
Sigma = [-mc1^2/lambda -mc3^2/lambda -(mc1*mc3)/lambda
         -mc3^2/lambda -mc2^2/lambda -(mc2*mc3)/lambda
         -(mc1*mc3)/lambda -(mc2*mc3)/lambda -(mc3^2+mc1*mc2)/(2*lambda)];
f = @(c) -lambda1*c(1)-lambda2*c(2)-lambda3*c(3)-lambda*log(c(1)*c(2)-c(3)^2);

C_sample = zeros(N,5);
n = 0;
while n<N
    n = n+1;
    c_old = mvnrnd(Mu,Sigma);
    c1_old = c_old(1);
    c2_old = c_old(2);
    c3_old = c_old(3);
    f_old = f(c_old);
    g_old = [-lambda1-(c2_old*lambda)/(-c3_old^2+c1_old*c2_old)
             -lambda2-(c1_old*lambda)/(-c3_old^2+c1_old*c2_old)
             -lambda3+(2*c3_old*lambda)/(-c3_old^2+c1_old*c2_old)];
    Sigma_old = [-c1_old^2/lambda -c3_old^2/lambda -(c1_old*c3_old)/lambda
                 -c3_old^2/lambda -c2_old^2/lambda -(c2_old*c3_old)/lambda
                 -(c1_old*c3_old)/lambda -(c2_old*c3_old)/lambda -(c3_old^2+c1_old*c2_old)/(2*lambda)];
    Mu_old = (c_old' + Sigma_old*g_old)';
    q_old = mvnpdf(c_old,Mu_old,Sigma_old);
    %
    c_prop = mvnrnd(Mu_old,Sigma_old);
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
C_sample(:,4) = gamrnd(a,b4,N,1);
C_sample(:,5) = gamrnd(a,b5,N,1);
lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda];

end
