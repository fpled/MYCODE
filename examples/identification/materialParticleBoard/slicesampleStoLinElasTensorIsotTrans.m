function C_sample = slicesampleStoLinElasTensorIsotTrans(lambda,C_data,N)
% function C_sample = slicesampleStoLinElasTensorIsotTrans(lambda,C_data,N)
% Slice Sampling for stochastic linear elastic tensor with
% transversely isotropic symmetry
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
% la4 = lambda(4);
% la5 = lambda(5);
la  = lambda(6);

mC = mean(C_data(:,1:3),1);

%% Sample generation
start = mC;
nsamples = N;
% pdf = @(c) (c(1)>0).*(c(2)>0).*(c(1)*c(2)-c(3)^2>0)...
%     .*(c(1)*c(2)-c(3)^2)^(-la).*exp(-la1*c(1)-la2*c(2)-la3*c(3));
% pdf = @(c) (c(1)>0).*(c(2)>0).*(c(1)*c(2)-c(3)^2>0)...
%     .*exp(-la*log(c(1)*c(2)-c(3)^2)-la1*c(1)-la2*c(2)-la3*c(3));
C_sample = slicesample(start,nsamples,'pdf',@(c) pdf(c,la1,la2,la3,la));

end

function p = pdf(c,la1,la2,la3,la)
supp = (c(1)>0).*(c(2)>0).*(c(1)*c(2)-c(3)^2>0);
p = exp(-la*log(c(1)*c(2)-c(3)^2)-la1*c(1)-la2*c(2)-la3*c(3));
if p==0
    p = exp(-745);
end
p = supp.*p;
end