rng default  % For reproducibility
delta = 5;
pdf = @(x) normpdf(x);
proppdf = @(x,y) unifpdf(y-x,-delta,delta);
proprnd = @(x) x + rand*2*delta - delta;   
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