function C_sample = mhsampleStoLinElasTensorIsotTrans_CUM(param,C_data,N)
% function C_sample = mhsampleStoLinElasTensorIsotTrans_CUM(param,C_data,N)
% Metropolis-Hastings Sampling using Componentwise Updating Method (CUM)
% for stochastic linear elastic tensor with transversely isotropic symmetry
% The proposal distribution is symmetric with normal distribution for the 
% components c1, c2, c3.
% See page 31 in 'Computational Statistics with Matlab' 
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

sc1 = std(C_data(:,1));
sc2 = std(C_data(:,2));
sc3 = std(C_data(:,3));

%% Sample generation
% Parameters of the trivariate normal distribution
Mu = [mc1 mc2 mc3];
Sigma = [-mc1^2/lambda -mc3^2/lambda -(mc1*mc3)/lambda
         -mc3^2/lambda -mc2^2/lambda -(mc2*mc3)/lambda
         -(mc1*mc3)/lambda -(mc2*mc3)/lambda -(mc3^2+mc1*mc2)/(2*lambda)];

% Initial values for c1, c2, c3
c1 = unifrnd( min(C1_data), max(C1_data) );
c2 = unifrnd( min(C2_data), max(C2_data) );
c3 = unifrnd( min(C3_data), max(C3_data) );
n = 1;
C_sample = zeros(N,5);
C_sample(n,1:3) = [c1 c2 c3];
while n<N
    n = n+1;
    
    new_c1 = normrnd(c1,sc1);
    pratio = mvnpdf([new_c1 c2 c3],Mu,Sigma)/...
             mvnpdf([    c1 c2 c3],Mu,Sigma);
    alpha = min([1 pratio]);
    u = rand;
    if u<alpha && new_c1>0
        C_sample(n,1) = new_c1;
    else
        C_sample(n,1) = C_sample(n-1,1);
    end
    
    new_c2 = normrnd(c2,sc2);
    pratio = mvnpdf([c1 new_c2 c3],Mu, Sigma)/...
             mvnpdf([c1     c2 c3],Mu, Sigma );
    alpha = min([1 pratio]);
    u = rand;
    if u<alpha && new_c2>0
        C_sample(n,2) = new_c2;
    else
        C_sample(n,2) = C_sample(n-1,2);
    end
    
    new_c3 = normrnd(c3,sc3);
    pratio = mvnpdf([c1 c2 new_c3],Mu,Sigma)/...
             mvnpdf([c1 c2     c3],Mu,Sigma);
    alpha = min([1 pratio]);
    u = rand;
    if u<alpha && (C_sample(n,1)*C_sample(n,2)-new_c3^2>0)
        C_sample(n,3) = new_c3;
    else
        C_sample(n,3) = C_sample(n-1,3);
    end
    
end
C_sample(:,4) = gamrnd(a,b4,N,1);
C_sample(:,5) = gamrnd(a,b5,N,1);
lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda];

end
