function [f,C_sample] = funoptimlseIsotTrans(lambda,C_data,mC_data,nuC_data,MCMCalg,varargin)
% function [f,C_sample] = funoptimlseIsotTrans(lambda,C_data,mC_data,nuC_data,MCMCalg,varargin)

la4 = lambda(4);
la5 = lambda(5);
la  = lambda(6);

a = 1-2*la;
b4 = 1/la4;
b5 = 1/la5;

mC4 = a*b4;
mC5 = a*b5;
% phi_C = ln((C1*C2-C3^2) * C4^2 * C5^2) = ln(C1*C2-C3^2) + 2*ln(C4) + 2*ln(C5)
% nu_C = E{phi_C} = E{ln(C1*C2-C3^2)} + 2*E{ln(C4)} + 2*E{ln(C5)}
%                 = mphiC123          + mphiC4      + mphiC5
% with Gamma-distributed random variables C4~Gamma(a4,b4) and C5~Gamma(a5,b5)
% E(ln(X)) = psi(a)+log(b) for a Gamma-distributed random variable X~Gamma(a,b)
mphiC4 = 2*(psi(a)+log(b4));
mphiC5 = 2*(psi(a)+log(b5));

% convergence analysis
N = 100; % initial number of samples
addSamplesFactor = 0.1; % percentage of additional samples
tol = 1e-6; % prescribed stagnation tolerance
Nmax = 1e4; % maximal number of samples
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
mC = [mC123 mC4 mC5];

% vC123 = var(C_sample,0,1);
% % vC123 = size(C_sample,1)/(size(C_sample,2)-1)*moment(C_sample,2,1);
% vC4 = a*b4^2;
% vC5 = a*b5^2;
% vC = [vC123 vC4 vC5];
%
% sC = sqrt(norm(vC));
% dC = sC/norm(mC);

phiC123 = log(C_sample(:,1).*C_sample(:,2)-C_sample(:,3).^2);
mphiC123 = mean(phiC123,1);
nuC = mphiC123 + mphiC4 + mphiC5;

err = norm([mC nuC] - [mC_data nuC_data])^2/norm([mC_data nuC_data])^2;
err_stagn = 1;
while (err_stagn > tol) && (err > tol) && (N < Nmax)
    Nadd = ceil(addSamplesFactor*N);
    N = N + Nadd;
    switch lower(MCMCalg)
        case 'mh'
            C_sample_add = mhsampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),Nadd);
        case 'bum'
            C_sample_add = mhsampleStoLinElasTensorIsotTrans_BUM(lambda,C_data(:,1:3),Nadd);
        case 'cum'
            C_sample_add = mhsampleStoLinElasTensorIsotTrans_CUM(lambda,C_data(:,1:3),Nadd);
        case 'ss'
            C_sample_add = slicesampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),Nadd);
        otherwise
            error(['MCMC algorithm ' MCMC ' not implemented'])
    end
    C_sample = [C_sample;C_sample_add];
    mC123 = mean(C_sample,1);
    mC = [mC123 mC4 mC5];
    
    % vC123 = var(C_sample,0,1);
    % % vC123 = size(C_sample,1)/(size(C_sample,2)-1)*moment(C_sample,2,1);
    % vC4 = a*b4^2;
    % vC5 = a*b5^2;
    % vC = [vC123 vC4 vC5];
    %
    % sC = sqrt(norm(vC));
    % dC = sC/norm(mC);
    
    phiC123 = log(C_sample(:,1).*C_sample(:,2)-C_sample(:,3).^2);
    mphiC123 = mean(phiC123,1);
    nuC = mphiC123 + mphiC4 + mphiC5;
    
    err_old = err;
    err = norm([mC nuC] - [mC_data nuC_data])^2/norm([mC_data nuC_data])^2;
    err_stagn = abs(err-err_old);
end

f = norm([mC nuC] - [mC_data nuC_data])^2;
% f = norm([mC nuC] - [mC_data nuC_data])^2/norm([mC_data nuC_data])^2;

end
