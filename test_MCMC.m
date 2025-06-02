% data
ET_data = mean_ET_data*1e-3; % [GPa]
GL_data = mean_GL_data*1e-3; % [GPa]
EL_data = mean_EL_data*1e-3; % [GPa]
NUL_data = mean_NUL_data;
NUT_data = 0.1+0.2*rand(length(mean_ET_data),1); % artificial data for NUT varying from 0.1 to 0.3
GT_data = ET_data./(2*(1+NUT_data)); % [GPa]
kT_data = (EL_data.*ET_data)./(2*(1-NUT_data).*EL_data-4*ET_data.*(NUL_data).^2); % [GPa]
C1_data = EL_data + 4*(NUL_data.^2).*kT_data; % [GPa]
C2_data = 2*kT_data; % [GPa]
C3_data = 2*sqrt(2)*kT_data.*NUL_data; % [GPa]
C4_data = 2*GT_data; % [GPa]
C5_data = 2*GL_data; % [GPa]
C_data = [C1_data(:) C2_data(:) C3_data(:) C4_data(:) C5_data(:)];

% empirical estimates
mC_data = mean(C_data,1);
% vC_data = var(C_data,0,1);
% % vC_data = length(C_data)/(length(C_data)-1)*moment(C_data,2,1);
% sC_data = sqrt(norm(vC_data));
% dC_data = sC/norm(mC_data);
phiC_data = log((C_data(:,1).*C_data(:,2)-C_data(:,3).^2).*(C_data(:,4).^2).*(C_data(:,5).^2));
nuC_data = mean(phiC_data,1);

% initial guess
la = -100; % la < 1/2
la1 = -(mC_data(2)*la)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la1 > 0
la2 = -(mC_data(1)*la)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la2 > 0
la3 = (2*mC_data(3)*la)/(mC_data(1)*mC_data(2)-mC_data(3)^2); % la3 in R such that 2*sqrt(la1*la2)-la3 > 0
a = 1-2*la; % a > 0
la4 = a/mC_data(4); % la4 > 0
la5 = a/mC_data(5); % la5 > 0

lambda = [la1 la2 la3 la4 la5 la];

b4 = 1/la4;
b5 = 1/la5;

mC4 = a*b4;
mC5 = a*b5;
mphiC4 = 2*(psi(a)+log(b4));
mphiC5 = 2*(psi(a)+log(b5));

% convergence analysis
N = 100; % initial number of samples
addSamplesFactor = 0.1; % percentage of additional samples
tol = 1e-6; % prescribed stagnation tolerance
Nmax = 1e4; % maximal number of samples
switch lower(MCMCalg)
    case 'mh'
        [C_sample,accept] = mhsampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),N);
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

err = norm([mC nuC] - [mC_data nuC_data])^2/norm([mC_data nuC_data])^2
err_stagn = 1;
while (err_stagn > tol) && (err > tol) && (N < Nmax)
    Nadd = ceil(addSamplesFactor*N);
    N = N + Nadd
    switch lower(MCMCalg)
        case 'mh'
            [C_sample_add,accept] = mhsampleStoLinElasTensorIsotTrans(lambda,C_data(:,1:3),Nadd);
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
    err = norm([mC nuC] - [mC_data nuC_data])^2/norm([mC_data nuC_data])^2
    err_stagn = abs(err-err_old)
end

mC1s = arrayfun(@(x) mean(C_sample(1:x,1),1),1:N);
mC2s = arrayfun(@(x) mean(C_sample(1:x,2),1),1:N);
mC3s = arrayfun(@(x) mean(C_sample(1:x,3),1),1:N);
nuCs = arrayfun(@(x) mean(phiC123(1:x,1),1) + mphiC4 + mphiC5,1:N);

figure('Name','Convergence mean of C1, C2, C3')
clf
plot(1:N,mC1s,'-b','LineWidth',linewidth)
hold on
plot(1:N,mC2s,'-r','LineWidth',linewidth)
plot(1:N,mC3s,'-g','LineWidth',linewidth)
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Number of samples','Interpreter',interpreter)
ylabel('Mean value','Interpreter',interpreter)
%xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
%ylabel('Moyenne','Interpreter',interpreter)
legend('$C_1$','$C_2$','$C_3$','Interpreter',interpreter)
% mysaveas(pathname,'convergence_mean_C1_C2_C3_lambda_init',formats);
% mymatlab2tikz(pathname,'convergence_mean_C1_C2_C3_lambda_init.tex');

figure('Name','Convergence logarithmic mean of det([C])')
clf
plot(1:N,nuCs,'-b','LineWidth',linewidth)
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('Number of samples','Interpreter',interpreter)
ylabel('Logarithmic mean of $\det([${\boldmath$C$}$])$','Interpreter',interpreter)
%xlabel('Nombre de r\''ealisations','Interpreter',interpreter)
%ylabel('Moyenne du logarithme de $\det([${\boldmath$C$}$])$','Interpreter',interpreter)
% mysaveas(pathname,'convergence_log_mean_detC_lambda_init',formats);
% mymatlab2tikz(pathname,'convergence_log_mean_detC_lambda_init.tex');

f = norm([mC nuC] - [mC_data nuC_data])^2
% f = norm([mC nuC] - [mC_data nuC_data])^2/norm([mC_data nuC_data])^2

% f0 = funoptimlseIsotTrans(lambda0,C_data,mC_data,nuC_data,MCMCalg)
