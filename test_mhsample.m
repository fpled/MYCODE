rng default  % For reproducibility
delta = 5;
pdf = @(x) normpdf(x);                      % target pdf
proppdf = @(x,y) unifpdf(y-x,-delta,delta); % proposal pdf
proprnd = @(x) x + rand*2*delta - delta;    % proposal random sampler
% nsamples = 15000;
nsamples = 1e6;
[x,accept] = mhsample(1,nsamples,'pdf',pdf,'proprnd',proprnd,'symmetric',1);
if accept~=1
    deltas = diff(x); % (N-1)xd
    isAccept = [false; any(deltas~=0,2)]; % new sample row is accepted if any coordinate changed
    x_accept = x(isAccept,:); % accepted samples
else
    x_accept = x;
end

figure
nbins = 50;
% h = histfit(x,nbins); % h = histogram(x,nbins);
% h = histfit(x); % h = histogram(x,ceil(sqrt(length(x))));
% h(1).FaceColor = [.8 .8 1];
h = histogram(x,'Normalization','pdf');
hold on
fplot(pdf,h.BinLimits,'-b','LineWidth',1)
[f,xi] = ksdensity(x);
plot(xi,f,'--r','LineWidth',1)
h_accept = histogram(x_accept,'Normalization','pdf');
[f_accept,xi_accept] = ksdensity(x_accept);
plot(xi_accept,f_accept,'--k','LineWidth',1)
xlabel('$x$','Interpreter','latex')
ylabel('$p_X(x)$','Interpreter','latex')
legend('hist','normpdf','estpdf','hist accepted','estpdf accepted')
hold off

rng default;  % For reproducibility
s = rng;
alpha = 2.43;
beta = 1;
pdf = @(x)gampdf(x,alpha,beta);                            % target pdf
proppdf = @(x,y)gampdf(x,floor(alpha),floor(alpha)/alpha); % proposal pdf
proprnd = @(x)sum(...                                      % proposal random sampler
              exprnd(floor(alpha)/alpha,floor(alpha),1));
% nsamples = 5000;
% smpl = mhsample(1,nsamples,'pdf',pdf,'proprnd',proprnd,...
%                 'proppdf',proppdf);
% xxhat = cumsum(smpl.^2)./(1:nsamples)';

rng(s)
smpl_5e3 = mhsample(1,5e3,'pdf',pdf,'proprnd',proprnd,...
                'proppdf',proppdf);
xxhat_5e3 = cumsum(smpl_5e3.^2)./(1:5e3)';

rng(s)
smpl_1e4 = mhsample(1,1e4,'pdf',pdf,'proprnd',proprnd,...
                'proppdf',proppdf);
xxhat_1e4 = cumsum(smpl_1e4.^2)./(1:1e4)';

rng(s)
smpl_1e5 = mhsample(1,1e5,'pdf',pdf,'proprnd',proprnd,...
                'proppdf',proppdf);
xxhat_1e5 = cumsum(smpl_1e5.^2)./(1:1e5)';

rng(s)
smpl_1e6 = mhsample(1,1e6,'pdf',pdf,'proprnd',proprnd,...
                'proppdf',proppdf);
xxhat_1e6 = cumsum(smpl_1e6.^2)./(1:1e6)';

[m,v] = gamstat(alpha,beta);
m2 = v + m^2;

figure
plot(1:5e3,xxhat_5e3,'-b','LineWidth',1)
hold on
plot(1:1e4,xxhat_1e4,'-m','LineWidth',1)
plot(1:1e5,xxhat_1e5,'-g','LineWidth',1)
plot(1:1e6,xxhat_1e6,'-c','LineWidth',1)
plot([1,1e6],[m2,m2],'-r','LineWidth',1)
hold off