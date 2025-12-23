function [f,C_sample,varargout] = funoptimlseElasIsotTrans(lambda,C_data,mC_data,nuC_data,N,MCMCalgo,s,varargin)
% function [f,C_sample,varargout] = funoptimlseElasIsotTrans(lambda,C_data,mC_data,nuC_data,N,MCMCalgo,s,varargin)

if nargin<6 || isempty(MCMCalgo)
    MCMCalgo = 'RWMH'; % Random Walk Metropolis-Hastings algorithm
end
if nargin==7 && ~isempty(s)
    rng(s); % initialize the random number generator based in the settings contained in s
end

useRedParam = (length(lambda)==4); % reduced parameterization
if useRedParam
    % reduced parameterization
    la  = lambda(4); % la < 1/2, la < 0 for MH sampling using a trivariate normal proposal pdf with positive-definite covariance matrix Sigma
else
    % full parameterization
    la4 = lambda(4); % la4 > 0
    la5 = lambda(5); % la5 > 0
    la  = lambda(6); % la < 1/2, la < 0 for MH sampling using a trivariate normal proposal pdf with positive-definite covariance matrix Sigma
end
a = 1-2*la; % a > 0
if useRedParam
    % reduced parameterization
    la4 = a/mC_data(4); % la4 > 0
    la5 = a/mC_data(5); % la5 > 0
end
b4 = 1/la4; % b4 > 0
b5 = 1/la5; % b5 > 0

mC45 = gamstat(a,[b4,b5]);
% [mC45,vC45] = gamstat(a,[b4,b5]);

% mC4 = gamstat(a,b4); % mC4 = a*b4;
% mC5 = gamstat(a,b5); % mC5 = a*b5;
% [mC4,vC4] = gamstat(a,b4); % vC4 = a*b4^2;
% [mC5,vC5] = gamstat(a,b5); % vC5 = a*b5^2;
% mC45 = [mC4 mC5]; % mC45 = a*[b4,b5];
% vC45 = [vC4 vC5]; % vC45 = a*[b4,b5].^2;

% stdC45 = sqrt(v45); % stdC45 = sqrt(a)*[b4,b5];
% cvC45 = stdC45./mC45; % cvC45 = 1/sqrt(a)*[1,1];

% For a Gamma-distributed random variable X~Gamma(a,b), phiX = log(X),
% nuX = E{phiX} = E{log(X)} = psi(a)+log(b)
% For Gamma-distributed random variables C4~Gamma(a,b4) and C5~Gamma(a,b5),
% phiC4 = log(C4) and phiC5 = log(C5),
% nuC4 = E{phiC4} = E{log(C4)} = psi(a)+log(b4)
% nuC5 = E{phiC5} = E{log(C5)} = psi(a)+log(b5)
nuC4 = psi(a)+log(b4);
nuC5 = psi(a)+log(b5);

switch lower(MCMCalgo)
    case {'imh','rwmh'}
        [C_sample,accept] = mhsampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),N,MCMCalgo,varargin{:});
        if nargout==3
            varargout{1} = accept;
        end
    case 'ss'
        [C_sample,neval] = slicesampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),N,varargin{:});
        if nargout==3
            varargout{1} = neval;
        end
    otherwise
        error(['MCMC algorithm ' MCMC ' not implemented'])
end

mC123 = mean(C_sample,1);
% vC123 = var(C_sample,0,1); % vC123 = N/(N-1)*moment(C_sample,2,1);
% stdC123 = std(C_sample,0,1);
% cvC123 = stdC123./mC123;

mC = [mC123 mC45];
% vC = [vC123 vC45];
% stdC = [stdC123 stdC45];
% cvC = [cvC123 cvC45];
% stdC = sqrt(vC);
% cvC  = stdC./mC;
% sC = sqrt(norm(vC));
% mCnorm = norm(mC);
% dC = sC/mCnorm;

% phiC123 = log(det([C123]) = log(C1*C2-C3^2)
% phiC = log(det([C]) = log((C1*C2-C3^2) * C4^2 * C5^2)
%                     = log(C1*C2-C3^2) + 2*log(C4) + 2*log(C5)
%                     = phiC123         + 2*phiC4   + 2*phiC5
% nuC = E{phiC} = E{log(C1*C2-C3^2)} + 2*E{log(C4)} + 2*E{log(C5)}
%               = E{phiC123}         + 2*E{phiC4}   + 2*E{phiC5}
%               = nuC123             + 2*nuC4       + 2*nuC5
phiC123 = log(C_sample(:,1).*C_sample(:,2)-C_sample(:,3).^2);
nuC123 = mean(phiC123,1);
nuC = nuC123 + 2*(nuC4 + nuC5);

fC = [mC nuC];
fC_data = [mC_data nuC_data];

% f = norm(fC - fC_data)^2; % squared norm of residual
% f = norm(fC - fC_data)^2/norm(fC_data)^2; % relative squared norm of residual
% f = norm((fC - fC_data)./fC_data)^2; % squared norm of relative residual
f = norm(mC - mC_data)^2/norm(mC_data)^2 + abs(nuC - nuC_data)^2/abs(nuC_data)^2; % relative squared norm of residual

% if useRedParam
%     % fprintf('lambda = (%g, %g, %g, %g, %g, %g), fval = %g, resnorm = %g, error(norm(mean(C))) = %g, error = (%g, %g, %g, %g, %g, %g)\n',lambda(1:3),la4,la5,lambda(4),f,norm(fC - fC_data)^2,norm(mC - mC_data)/norm(mC_data),abs(fC - fC_data)./abs(fC_data));
%     fprintf('lambda = (%g, %g, %g, %g, %g, %g)\n',lambda(1:3),la4,la5,lambda(4));
% else
%     % fprintf('lambda = (%g, %g, %g, %g, %g, %g), fval = %g, resnorm = %g, error(norm(mean(C))) = %g, error = (%g, %g, %g, %g, %g, %g)\n',lambda,f,norm(fC - fC_data)^2,norm(mC - mC_data)/norm(mC_data),abs(fC - fC_data)./abs(fC_data));
%     fprintf('lambda = (%g, %g, %g, %g, %g, %g)\n',lambda);
% end
% fprintf('fval = %g, resnorm = %g, error(norm(mean(C))) = %g, error = (%g, %g, %g, %g, %g, %g)\n',f,norm(fC - fC_data)^2,norm(mC - mC_data)/norm(mC_data),abs(fC - fC_data)./abs(fC_data));

end
