%%
% The purpose of this file is to obtain a reduced model
% using SVD for multiinclusion problem
% Observations
% 1. L1 regularization (alpha update) is  not used
% 2. Regularization in regression is optimal at 1e-02
% 3. No updates or alpha updates
%% Load Data
clear all
load ./svd_data/RV;
load ./svd_data/useppc;
load ./svd_data/S;
load ./svd_data/orefpc;
%% Sampling for regression
KL = 1; % check this value if you are using KL or normal SVD
gs = [];
nSamples = 1500;
[fn,rs] = random(useppc,nSamples);
gs(:,1:nSamples) = fn{:};
if KL == 1
    gsKL=[];
    meanMode = mean(gs,2);
    gsKL(:,1:nSamples) = gs(:,1:nSamples) - repmat(meanMode,1,nSamples);
    [U,Si,V] = svd(gsKL,0);
else
    [U,Si,V] = svd(gs,0);
end
%% Get the true error for modes
p=5;
modei=1;
modef=8;
SVDM=SEPARATEDSOLVER('order',p,'basistype',...
    1,'typebase',1,'regroup',0,'nmultelem',2,...
    'waveres',3,'modei',modei,'modef',modef);
[fsReg,rsReg,fsTest,rsTest] = cross_validation_svd(V,rs,modef);
[errorModes solution]=test_karloeve(SVDM,RV,fsReg,rsReg,fsTest,rsTest);
%% Plotting CV error for the first 8 modes
figure(110)
for i = 1:8
    ha(i) = subplot(2,4,i);
    nb=getm(solution{i});
    plot([1:nb], log10(errorModes{i}(1,[1:nb])),'--xb','LineWidth',2);
    %hold on
    %plot([1:nb], log10(errorModes_500{i}(1,[1:nb])),'-.or','LineWidth',2);
    %hold on
    %plot([1:nb], log10(errorModes_1500{i}(1,[1:nb])),'-*g','LineWidth',2);
    xlabel('Rank','fontsize',12);
    ylabel('log(Error)','fontsize',12);
    set(ha(i),'YLim',[-1.0 1.0]);
end
%% Plotting the cross validation error

semilogy([1:modes],min_cv_error,'-.sr','LineWidth',2);
xlabel('i','fontsize',14);
ylabel('Min CV error for W_{i}','fontsize',14);

%% Testing the reduced model
figure(3)
model_rs = random(RV);
model_val = randomeval(useppc,model_rs);
plot(model_val,S);
%% Making reduced model for integrated velocity on square domain
fx = LINFORM(0,1);
fx = setselgroup(fx,getnbgroupelem(S));
fx = fx{S}(:);
out = fx;

sol_fin = meanMode*one(PC);
for k=1:8
    sol_fin = sol_fin +U(:,k)*Si(k,k)*solution{k};
end
%% Monte Carlo Error Estimation

ntestsamples = 100000;
[qoi,sample_pts] = random(orefpc,ntestsamples);
orefni = randomeval(out'*sol_fin,sample_pts);
rel_error = abs(qoi - orefni)./qoi;
mean_rel_error = mean(rel_error);
var_i = rel_error - mean_rel_error*ones(ntestsamples,1)
var = (var_i'*var_i)/(ntestsamples-1);
std_dev = sqrt(var)/sqrt(ntestsamples)
%% Monte Carlo error values
error_val = [100,0.007801;500,0.003705;1000,0.001791];

%% Plotting error on integrated value
set(gca,'XTick',error_val(:,1),'ylim',[10^-3, 10^-2],'YTick',[10^-3,10^-2.75,10^-2.5,10^-2.25,10^-2]);
semilogy(error_val(:,1),(error_val(:,2)),'-.sr','LineWidth',2);
%axis([100 1000 1.0e-3 1.0e-02])
xlabel('Sample Size (Q)','fontsize',14);
ylabel('Error','fontsize',14);

%%
rmodel = zeros(size(model_val),1);
for i = 1:modes
    rmodel = rmodel+U(:,i)*Si(i,i)*randomeval(solution{i},model_rs);
    figure(4)
    plot(abs(model_val-rmodel)./abs(model_val),S);
    pause(1)
end
%% Plotting functions for different sample sizes in
% first 8 modes
figure(110)
for i = 1:8
    ha(i) = subplot(2,4,i);
    plot([1:nnnmax], log10(cv100([1:nnnmax],i)),'r-.+','LineWidth',2);
    hold on
    plot([1:nnnmax], log10(cv500([1:nnnmax],i)),'b-o','LineWidth',2);
    hold on
    plot([1:nnnmax], log10(cv1000([1:nnnmax],i)),'g--x','LineWidth',2);
    xlabel('Rank','fontsize',12);
    ylabel('log(\epsilon_{CV})','fontsize',12);
    set(ha(i),'YLim',[-1.0 1.0]);
end

