function C_sample = mhsampleStoLinElasTensorIsotTrans_CUM(lambda,C_data,N)
% function C_sample = mhsampleStoLinElasTensorIsotTrans_CUM(lambda,C_data,N)
% Metropolis-Hastings Sampling using Componentwise Updating Method (CUM)
% for stochastic linear elastic tensor with transversely isotropic symmetry
% The proposal distribution is symmetric with normal distribution for the 
% components c1, c2, c3.
% See page 31 in 'Computational Statistics with Matlab' 
% 
% lambda = (la1,la2,la3,la4,la5,la)
% C_data: data set for random vector C=(C1,C2,C3)
% C_data(:,1): data for random coordinate C1
% C_data(:,2): data for random coordinate C2
% C_data(:,3): data for random coordinate C3
% N: number of samples
% C_sample: sample set for random vector C=(C1,C2,C3)
% C_sample(:,1): data for random coordinate C1
% C_sample(:,2): data for random coordinate C2
% C_sample(:,3): data for random coordinate C3

la1 = lambda(1);
la2 = lambda(2);
la3 = lambda(3);
la4 = lambda(4);
la5 = lambda(5);
la  = lambda(6);

mC = mean(C_data(:,1:3),1);
sC = std(C_data(:,1:3),0,1);

%% Sample generation
% Parameters of the trivariate normal distribution
Mu = mC; % Mu = [mC(1) mC(2) mC(3)];
Sigma = [-mC(1)^2/la -mC(3)^2/la -(mC(1)*mC(3))/la
         -mC(3)^2/la -mC(2)^2/la -(mC(2)*mC(3))/la
         -(mC(1)*mC(3))/la -(mC(2)*mC(3))/la -(mC(3)^2+mC(1)*mC(2))/(2*la)];

% Initial values for c1, c2, c3
c1 = unifrnd( min(C_data(:,1)), max(C_data(:,1)) );
c2 = unifrnd( min(C_data(:,2)), max(C_data(:,2)) );
c3 = unifrnd( min(C_data(:,3)), max(C_data(:,3)) );
n = 1;
C_sample = zeros(N,5);
C_sample(n,1:3) = [c1 c2 c3];
while n<N
    n = n+1;
    
    new_c1 = normrnd(c1,sC(1));
    pratio = mvnpdf([new_c1 c2 c3],Mu,Sigma)/...
             mvnpdf([    c1 c2 c3],Mu,Sigma);
    alpha = min([1 pratio]);
    u = rand;
    if u<alpha && new_c1>0
        C_sample(n,1) = new_c1;
    else
        C_sample(n,1) = C_sample(n-1,1);
    end
    
    new_c2 = normrnd(c2,sC(2));
    pratio = mvnpdf([c1 new_c2 c3],Mu,Sigma)/...
             mvnpdf([c1     c2 c3],Mu,Sigma);
    alpha = min([1 pratio]);
    u = rand;
    if u<alpha && new_c2>0
        C_sample(n,2) = new_c2;
    else
        C_sample(n,2) = C_sample(n-1,2);
    end
    
    new_c3 = normrnd(c3,sC(3));
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

end
