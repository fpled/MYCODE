function f = funlsqnonlinElasIsot(lambda,mC_data,nuC_data)
% function f = funlsqnonlinElasIsot(lambda,mC_data,nuC_data)

useRedParam = isscalar(lambda); % reduced parameterization
if useRedParam
    % reduced parameterization
    la  = lambda(1);     % la < 1/5
    a1  = 1-la;          % a1 > 0
    a2  = 1-5*la;        % a2 > 0
    la1 = a1/mC_data(1); % la1 > 0
    la2 = a2/mC_data(2); % la2 > 0
else
    % full parameterization
    la1 = lambda(1); % la1 > 0
    la2 = lambda(2); % la2 > 0
    la  = lambda(3); % la < 1/5
    a1 = 1-la;       % a1 > 0
    a2 = 1-5*la;     % a2 > 0
end
b1 = 1/la1; % b1 > 0
b2 = 1/la2; % b2 > 0

mC = gamstat([a1,a2],[b1,b2]);
% [mC,vC] = gamstat([a1,a2],[b1,b2]);

% mC1 = gamstat(a1,b1); % mC1 = a1*b1;
% mC2 = gamstat(a2,b2); % mC2 = a2*b2;
% [mC1,vC1] = gamstat(a1,b1); % vC1 = a1*b1^2;
% [mC2,vC2] = gamstat(a2,b2); % vC2 = a2*b2^2;
% mC = [mC1 mC2]; % mC = [a1*b1,a2*b2];
% vC = [vC1 vC2]; % vC = [a1*b1^2,a2*b2^2];
% stdC = sqrt(vC); % stdC = [sqrt(a1)*b1,sqrt(a2)*b2];
% cvC = stdC./mC; % cvC = 1./[sqrt(a1),sqrt(a2)];
% sC = sqrt(norm(vC));
% mCnorm = norm(mC);
% dC = sC/mCnorm;

% phi_C = log(96*C1*C2^5) = log(96) + log(C1) + 5*log(C2)
% nu_C = E{phi_C}  = E{log(96*C1*C2^5)} = log(96) + E{log(C1)} + 5*E{log(C2)}
% with Gamma-distributed random variables C1~Gamma(a1,b1) and C2~Gamma(a2,b2)
% E(log(X)) = psi(a)+log(b) for a Gamma-distributed random variable X~Gamma(a,b)
nuC = log(96) + psi(a1)+log(b1) + 5*(psi(a2)+log(b2));

fC = [mC nuC];
fC_data = [mC_data nuC_data];

f = fC - fC_data; % residual
% f = (fC - fC_data)./fC_data; % relative residual

% if useRedParam
%     % fprintf('lambda = (%g, %g, %g), resnorm = %g, error(norm(mean(C))) = %g, error = (%g, %g, %g)\n',la1,la2,lambda,norm(f)^2,norm(fC - fC_data)^2,norm(mC - mC_data)/norm(mC_data),abs(fC - fC_data)./abs(fC_data));
%     fprintf('lambda = (%g, %g, %g)\n',la1,la2,lambda);
% else
%     % fprintf('lambda = (%g, %g, %g), resnorm = %g, error(norm(mean(C))) = %g, error = (%g, %g, %g)\n',lambda,norm(f)^2,norm(mC - mC_data)/norm(mC_data),abs(fC - fC_data)./abs(fC_data));
%     fprintf('lambda = (%g, %g, %g)\n',lambda);
% end
% fprintf('resnorm = %g, error(norm(mean(C))) = %g, error = (%g, %g, %g)\n',norm(f)^2,norm(mC - mC_data)/norm(mC_data),abs(fC - fC_data)./abs(fC_data));

end
