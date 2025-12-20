function [C_sample,neval] = slicesampleStoLinElasTensorIsotTrans(lambda,C_data,N,varargin)
% function [C_sample,neval] = slicesampleStoLinElasTensorIsotTrans(lambda,C_data,N,varargin)
% Slice Sampling for stochastic linear elastic tensor with transversely isotropic symmetry
% lambda: parameter vector [la1, la2, la3, la4, la5, la] for full parameterization
%                          [la1, la2, la3, la]           for reduced parameterization
% C_data: data set for random vector C=(C1,C2,C3)
% C_data(:,i): data for random coordinate Ci
% N: number of samples
% C_sample: sample set for random vector C=(C1,C2,C3)
% C_sample(:,i): samples for random coordinate Ci
% neval: average number of function evaluations per sample that occurred in the slice sampling

useRedParam = (length(lambda)==4); % reduced parameterization

la1 = lambda(1); % la1 > 0
la2 = lambda(2); % la2 > 0
la3 = lambda(3); % la3 in R such that 2*sqrt(la1*la2)-la3 > 0
if useRedParam
    % reduced parameterization
    la  = lambda(4); % la < 1/2
else
    % full parameterization
    % la4 = lambda(4); % la4 > 0
    % la5 = lambda(5); % la5 > 0
    la  = lambda(6); % la < 1/2
end

mC_data = mean(C_data(:,1:3),1);

%% Sample generation
initial = getcharin('initial',varargin,mC_data); % initial point (start value) of the Markov Chain
nsamples = N; % number of samples to be generated
% pdf = @(c) (c(1)>0).*(c(2)>0).*(c(1)*c(2)-c(3)^2>0)...
%     .*(c(1)*c(2)-c(3)^2)^(-la).*exp(-la1*c(1)-la2*c(2)-la3*c(3)); % density function for target stationary distribution
% pdf = @(c) (c(1)>0).*(c(2)>0).*(c(1)*c(2)-c(3)^2>0)...
%     .*exp(-la*log(c(1)*c(2)-c(3)^2)-la1*c(1)-la2*c(2)-la3*c(3)); % density function for target stationary distribution
pdf = @(c) mypdf(c,la1,la2,la3,la); % density function for target stationary distribution
[C_sample,neval] = slicesample(initial,nsamples,'pdf',pdf);

end

function p = mypdf(c,la1,la2,la3,la)
supp = (c(1)>0).*(c(2)>0).*(c(1)*c(2)-c(3)^2>0);
% p = (c(1)*c(2)-c(3)^2)^(-la)*exp(-la1*c(1)-la2*c(2)-la3*c(3));
p = exp(-la*log(c(1)*c(2)-c(3)^2)-la1*c(1)-la2*c(2)-la3*c(3));
if p==0
    p = realmin;
end
p = supp.*p;
end
