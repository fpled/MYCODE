function f = funoptimlseIsotTrans(lambda,C_data,mC_data,dC_data,MCMCalg,varargin)
% function f = funoptimlseIsotTrans(lambda,C_data,mC_data,dC_data,MCMCalg,varargin)

la4 = lambda(4);
la5 = lambda(5);
la  = lambda(6);

a = 1-2*la;
b4 = 1/la4;
b5 = 1/la5;

N = 1e3; % number of samples

switch lower(MCMCalg)
    case 'mh'
        C_sample = mhsampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),N);
    case 'bum'
        C_sample = mhsampleStoLinElasTensorIsotTrans_BUM(lambda,C_data(:,1:3),N);
    case 'cum'
        C_sample = mhsampleStoLinElasTensorIsotTrans_CUM(lambda,C_data(:,1:3),N);
    case 'ss'
        C_sample = slicesampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),N);
    otherwise
        error(['MCMC algorithm ' MCMC ' not implemented'])
end

mC123 = mean(C_sample,1);
mC4 = a*b4;
mC5 = a*b5;
mC = [mC123 mC4 mC5];

vC123 = var(C_sample,0,1);
% vC123 = size(C_sample,1)/(size(C_sample,2)-1)*moment(C_sample,2,1);
vC4 = a*b4^2;
vC5 = a*b5^2;
vC = [vC123 vC4 vC5];

sC = sqrt(norm(vC));
dC = sC/norm(mC);

alpha = 1/2;

f = alpha * norm(mC - mC_data)^2/norm(mC_data)^2 + (1-alpha) * (dC - dC_data)^2/(dC_data)^2;

end
