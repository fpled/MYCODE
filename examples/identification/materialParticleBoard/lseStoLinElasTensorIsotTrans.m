function [lambda,err,exitflag,output] = lseStoLinElasTensorIsotTrans(C_data,MCMCalg)
% function [lambda,err,exitflag,output] = lseElasTensorIsot(C_data,MCMCalg)
% Least-squares estimation for stochastic linear elastic tensor with
% transversely isotropic symmetry
% C_data: data set for random vector C=(C1,C2,C3,C4,C5)
% C_data(:,i): data for random coordinate Ci
% MCMCalg: algorithm for Markov-Chain Monte Carlo (MCMC) method
% MCMCalg = 'MH', 'BUM', 'CUM' or 'SS' ('MH' by default)
% lambda: parameters (la1,la2,la3,la4,la5,la)
% lambda(1) = la1 > 0
% lambda(2) = la2 > 0
% lambda(3) = la3 in R such that 2*sqrt(la1*la2)-la3 > 0
% lambda(4) = la4 > 0
% lambda(5) = la5 > 0
% lambda(6) = la < 1/2

if nargin==1 || isempty(MCMCalg)
    MCMCalg = 'MH';
end

% empirical estimates
mC_data = mean(C_data,1);
% vC_data = var(C_data,0,1);
% % vC_data = length(C_data)/(length(C_data)-1)*moment(C_data,2,1);
% sC_data = sqrt(norm(vC_data));
% dC_data = sC_data/norm(mC_data);
phiC_data = log((C_data(:,1).*C_data(:,2)-C_data(:,3).^2).*(C_data(:,4).^2).*(C_data(:,5).^2));
nuC_data = mean(phiC_data,1);

% initial guess
la = -100; % la < 1/2
la1 = -(mC_data(2)*la)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la1 > 0
la2 = -(mC_data(1)*la)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la2 > 0
la3 = (2*mC_data(3)*la)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la3 in R such that 2*sqrt(la1*la2)-la3 > 0
a = 1-2*la; % a > 0
% la4 = a/mC_data(4); % la4 > 0
% la5 = a/mC_data(5); % la5 > 0
la4 = -2*la/mC_data(4); % la4 > 0
la5 = -2*la/mC_data(5); % la5 > 0
b4 = 1/la4; % b4 > 0
b5 = 1/la5; % b5 > 0

lambda0 = [la1 la2 la3 la4 la5 la];

% [f0,C_sample] = funoptimlseIsotTrans(lambda0,C_data,mC_data,nuC_data,MCMCalg);
% 
% N = size(C_sample,1);
% mC1s = arrayfun(@(x) mean(C_sample(1:x,1),1),1:N);
% mC2s = arrayfun(@(x) mean(C_sample(1:x,2),1),1:N);
% mC3s = arrayfun(@(x) mean(C_sample(1:x,3),1),1:N);
% mphiC4 = 2*(psi(a)+log(b4));
% mphiC5 = 2*(psi(a)+log(b5));
% nuCs = arrayfun(@(x) mean(log(C_sample(1:x,1).*C_sample(1:x,2)-C_sample(1:x,3).^2),1) + mphiC4 + mphiC5,1:N);
% 
% figure('Name','Convergence mean')
% clf
% plot(1:N,mC1s,'-b','LineWidth',1)
% hold on
% plot(1:N,mC2s,'-r','LineWidth',1)
% plot(1:N,mC3s,'-g','LineWidth',1)
% grid on
% box on
% set(gca,'FontSize',16)
% xlabel('Number of samples','Interpreter','latex')
% ylabel('Mean','Interpreter','latex')
% % xlabel('Nombre de r\''ealisations','Interpreter','latex')
% % ylabel('Moyenne','Interpreter','latex')
% l = legend('$C_1$','$C_2$','$C_3$');
% set(l,'Interpreter','latex')
% 
% figure('Name','Convergence logarithmic mean of det([C])')
% clf
% plot(1:N,nuCs,'-b','LineWidth',1)
% grid on
% box on
% set(gca,'FontSize',16)
% xlabel('Number of samples','Interpreter','latex')
% ylabel('Logarithmic mean of $\det([C])$','Interpreter','latex')
% % xlabel('Nombre de r\''ealisations','Interpreter','latex')
% % ylabel('Moyenne du logarithme de $\det([C])$','Interpreter','latex')

lb = [0 0 -Inf 0 0 -Inf];
ub = [Inf Inf Inf Inf Inf 1/2];
nonlcon = @funconIsotTrans;

% display = 'off';
% display = 'iter';
% display = 'iter-detailed';
display = 'final';
% display = 'final-detailed';

tolX = 1e-12; % tolerance on the parameter value
tolFun = 1e-12; % tolerance on the function value
tolCon = 1e-12; % tolerance on the constraint violation
% maxFunEvals = 3e3; % maximum number of function evaluations

% options  = optimoptions('fmincon','Display',display);
% options  = optimoptions('fmincon','Display',display,'TolX',tolX,'TolFun',tolFun,'TolCon',tolCon);
options  = optimoptions('fmincon','Display',display,'StepTolerance',tolX,'FunctionTolerance',tolFun,'ConstraintTolerance',tolCon);

fun = @(lambda) funoptimlseIsotTrans(lambda,C_data,mC_data,nuC_data,MCMCalg);
[lambda,err,exitflag,output] = fmincon(fun,lambda0,[],[],[],[],lb,ub,nonlcon,options);

end
