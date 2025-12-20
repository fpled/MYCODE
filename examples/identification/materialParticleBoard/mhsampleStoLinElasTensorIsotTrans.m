function [C_sample,accept] = mhsampleStoLinElasTensorIsotTrans(lambda,C_data,N,MCMCalgo,varargin)
% function [C_sample,accept] = mhsampleStoLinElasTensorIsotTrans(lambda,C_data,N,MCMCalgo,varargin)
% Metropolis-Hastings Sampling for stochastic linear elastic tensor with transversely isotropic symmetry
% lambda: parameter vector [la1, la2, la3, la4, la5, la] for full parameterization
%                          [la1, la2, la3, la]           for reduced parameterization
% C_data: data set for random vector C=(C1,C2,C3)
% C_data(:,i): data for random coordinate Ci
% N: number of samples
% MCMCalgo: Metropolis-Hastings algorithm for Markov Chain Monte Carlo (MCMC) method
% MCMCalgo = 'IMH' or 'RWMH' ('RWMH' by default)
% C_sample: sample set for random vector C=(C1,C2,C3)
% C_sample(:,i): samples for random coordinate Ci
% accept: acceptance rate of the proposal distribution

if nargin<4 || isempty(MCMCalgo)
    MCMCalgo = 'RWMH'; % Random Walk Metropolis-Hastings algorithm
end

useRedParam = (length(lambda)==4); % reduced parameterization

la1 = lambda(1); % la1 > 0
la2 = lambda(2); % la2 > 0
la3 = lambda(3); % la3 in R such that 2*sqrt(la1*la2)-la3 > 0
if useRedParam
    % reduced parameterization
    la  = lambda(4); % la < 1/2 for integrability of target pdf, la < 0 for trivariate normal proposal pdf with positive-definite covariance matrix Sigma
else
    % full parameterization
    % la4 = lambda(4); % la4 > 0
    % la5 = lambda(5); % la5 > 0
    la  = lambda(6); % la < 1/2 for integrability of target pdf, la < 0 for trivariate normal proposal pdf with positive-definite covariance matrix Sigma
end

mC_data = mean(C_data(:,1:3),1);

%% Sample generation
% Covariance matrix for the trivariate normal distribution
Sigma = -1/la*[mC_data(1)^2,          mC_data(3)^2,           mC_data(1)*mC_data(3);
               mC_data(3)^2,          mC_data(2)^2,           mC_data(2)*mC_data(3);
               mC_data(1)*mC_data(3), mC_data(2)*mC_data(3), (mC_data(1)*mC_data(2)+mC_data(3)^2)/2];

start = getcharin('initial',varargin,mC_data); % initial point (start value) of the Markov Chain
nsamples = N; % number of samples to be generated
% pdf = @(c) (c(1)>0).*(c(2)>0).*(c(1)*c(2)-c(3)^2>0)...
%     .*(c(1)*c(2)-c(3)^2)^(-la)*exp(-la1*c(1)-la2*c(2)-la3*c(3)); % density function for target stationary distribution
% pdf = @(c) (c(1)>0).*(c(2)>0).*(c(1)*c(2)-c(3)^2>0)...
%     .*exp(-la*log(c(1)*c(2)-c(3)^2)-la1*c(1)-la2*c(2)-la3*c(3)); % density function for target stationary distribution
pdf = @(c) mypdf(c,la1,la2,la3,la); % density function for target stationary distribution
switch MCMCalgo
    case 'IMH'
        % Independent Metropolis-Hastings sampling
        mu = mC_data;                        % empirical mean value independent of current values
        proppdf = @(x,y) mvnpdf(x,mu,Sigma); % asymmetric density function for proposal distribution
        proprnd = @(x) mvnrnd(mu,Sigma);     % random number generator for proposal distribution
        [C_sample,accept] = mhsample(start,nsamples,'pdf',pdf,'proppdf',proppdf,'proprnd',proprnd); % asymmetric proppdf (works for both Independent and Random Walk Metropolis-Hastings sampling)
    case 'RWMH'
        % Random Walk Metropolis-Hastings sampling
        % proppdf = @(x,y) mvnpdf(x,y,Sigma); % symmetric density function for proposal distribution
        proprnd = @(x) mvnrnd(x,Sigma);       % random number generator for proposal distribution
        [C_sample,accept] = mhsample(start,nsamples,'pdf',pdf,'proprnd',proprnd,'symmetric',true); % symmetric proppdf (works only for Random Walk Metropolis-Hastings sampling)
end

% The invariant distribution π_J(x) of the jump chain is proportional to π(x)*(1-ρ(x)) 
% with ρ(x) = rejection probability at x = 1 - int alpha(x,y) q(y|x) dy,
% where q(y|x) = proposal pdf and alpha(x,y) = acceptance probability for a proposal (candidate) y from x.
% Only in the special and very rare case where 1−ρ(x) is (almost) constant
% in x would lead to π_J = π, otherwise π_J differs from π.
% The correct Monte Carlo approximation to the invariant distribution π(x) is given by the full chain, 
% with each state repeated according to how many iterations the chain stayed there
% if accept~=1
%     deltas = diff(C_sample); % (N-1)xd
%     isAccept = [false; any(deltas~=0,2)]; % new sample row is accepted if any coordinate changed
%     C_sample = C_sample(isAccept,:); % accepted samples
% end

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
