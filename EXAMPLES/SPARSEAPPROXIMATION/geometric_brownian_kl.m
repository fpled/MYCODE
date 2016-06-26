function X = geometric_brownian_kl(xi,c,s,T,n)
% function X = geometric_brownian_kl(xi,c,s,T,n)

t = linspace(0,T,n+1);
dt = T/n;

Ns = size(xi,1);
M = size(xi,2);
B = 0;
for i=1:M
    B = B+xi(:,i)*sqrt(2)/pi/(i-1/2)*sin(pi*(i-1/2)*t);
end
dB = B(:,2:end)-B(:,1:end-1);

% dB = sqrt(dt)*xi;
X = zeros(Ns,n+1);
% B = zeros(Ns,n+1);
x0 = 1;
X(:,1) = x0;
for i=1:n
    % B(:,i+1) = B(:,i)+dB(:,i);
    X(:,i+1) = X(:,i) + c*X(:,i)*dt + s*X(:,i).*dB(:,i);
    % X(i+1) = X(i) + c*X(i)*dt + s*dB(i);
end

end
