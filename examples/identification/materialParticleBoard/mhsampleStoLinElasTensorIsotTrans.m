function [C_sample,lambda,lambda1,lambda2,lambda3,lambda4,lambda5] = mhsampleStoLinElasTensorIsotTrans(C1_data,C2_data,C3_data,C4_data,C5_data,N)
% function [C_sample,lambda,lambda1,lambda2,lambda3,lambda4,lambda5] = mhsampleStoLinElasTensorIsotTrans(C1_data,C2_data,C3_data,C4_data,C5_data,N)
% Metropolis-Hastings Sampling for stochastic linear elastic tensor with
% transversely isotropic symmetry
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
proppdf = @(c) mvnpdf([c(1) c(2) c(3)],Mu_mean,Sigma_mean);
proprnd = @(c) mvnrnd(Mu_mean,Sigma_mean);
pdf = @(c) (c(1)*c(2)-c(3)^2)^(-lambda)*exp(-lambda1*c(1)-lambda2*c(2)-lambda3*c(3));
nsamples = N;
start = [0 0 0];
smpl = mhsample(start,nsamples,'pdf',pdf,'proppdf',proppdf, 'proprnd',proprnd);

C_sample(1:3,:) = smpl;
C_sample(4,:) = gamrnd(a,b4,N,1);
C_sample(5,:) = gamrnd(a,b5,N,1);
lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda];

end
