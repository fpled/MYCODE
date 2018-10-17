function C_sample = slicesampleStoLinElasTensorIsotTrans(param,C_data,N)
% function C_sample = slicesampleStoLinElasTensorIsotTrans(param,C_data,N)
% Slice Sampling for stochastic linear elastic tensor with
% transversely isotropic symmetry
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
pdf = @(c) (c(1)>0).*(c(2)>0).*(c(1)*c(2)-c(3)^2>0).*(c(1)*c(2)-c(3)^2)^(-lambda)*exp(-lambda1*c(1)-lambda2*c(2)-lambda3*c(3));
nsamples = N;
start = [mc1 mc2 mc3];
smpl = slicesample(start,nsamples,'pdf',pdf);

C_sample(:,1:3) = smpl;
C_sample(:,4) = gamrnd(a,b4,N,1);
C_sample(:,5) = gamrnd(a,b5,N,1);

end
