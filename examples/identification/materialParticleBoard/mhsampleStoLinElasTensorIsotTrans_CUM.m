function [C_sample,lambda,lambda1,lambda2,lambda3,lambda4,lambda5] = mhsampleStoLinElasTensorIsotTrans_CUM(C1_data,C2_data,C3_data,C4_data,C5_data,N)
% function [C_sample,lambda,lambda1,lambda2,lambda3,lambda4,lambda5] = mhsampleStoLinElasTensorIsotTrans_CUM(C1_data,C2_data,C3_data,C4_data,C5_data,N)
% Metropolis-Hastings Sampling using Componentwise Updating Method (CUM)
% for stochastic linear elastic tensor with transversely isotropic symmetry
% The proposal distribution is symmetric with normal distribution for the 
% components c1, c2, c3.
% See page 31 in 'Computational Statistics with Matlab' 
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

sc1 = std(C1_data);
sc2 = std(C2_data);
sc3 = std(C3_data);

%% Sample generation
% Method to calculate lambda which controls the level of fluctuations
lambda = -200; % negative number
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

% Initial values for c1, c2, c3
c1 = unifrnd( min(C1_data), max(C1_data) );
c2 = unifrnd( min(C2_data), max(C2_data) );
c3 = unifrnd( min(C3_data), max(C3_data) );
n = 1;
C_sample = zeros(5,N);
C_sample(1:3,n) = [c1 c2 c3]';
while n<N
    n = n+1;
    
    new_c1 = normrnd(c1,sc1);
    pratio = mvnpdf([new_c1 c2 c3],Mu_mean,Sigma_mean)/...
             mvnpdf([    c1 c2 c3],Mu_mean,Sigma_mean);
    alpha = min([1 pratio]);
    u = rand;
    if u<alpha
        C_sample(1,n) = new_c1;
    else
        C_sample(1,n) = C_sample(1,n-1);
    end
    
    new_c2 = normrnd(c2,sc2);
    pratio = mvnpdf([c1 new_c2 c3],Mu_mean, Sigma_mean)/...
             mvnpdf([c1     c2 c3],Mu_mean, Sigma_mean );
    alpha = min([1 pratio]);
    u = rand;
    if u<alpha
        C_sample(2,n) = new_c2;
    else
        C_sample(2,n) = C_sample(2,n-1);
    end
    
    new_c3 = normrnd(c3,sc3);
    pratio = mvnpdf([c1 c2 new_c3],Mu_mean,Sigma_mean)/...
             mvnpdf([c1 c2     c3],Mu_mean,Sigma_mean);
    alpha = min([1 pratio]);
    u = rand;
    if u<alpha
        C_sample(3,n) = new_c3;
    else
        C_sample(3,n) = C_sample(3,n-1);
    end
    
end
C_sample(4,:) = gamrnd(a,b4,N,1);
C_sample(5,:) = gamrnd(a,b5,N,1);
lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda];

end
