%%
% The purpose of this file is to obtain a reduced model
% using SVD for multiinclsion problem
% Observations
% 1. L1 regularization (alpha update) is  not used
% 2. Regularization in regression is optimal at 1e-02
% 3. No updates or alpha updates
%% Load Data
clear all
load RV;
load useppc;
load S;
load orefpc;
%% Sampling for regression
KL = 1; % check this value if you are using KL or normal SVD
gs = [];
nSamples = 100;
[fn,rs] = random(useppc,nSamples);
gs(:,1:nSamples) = fn{:};
if KL == 1
    gsKL=[];
    meanMode = mean(gs,2);
    %gsKL = gs - meanMode*ones(1,nSamples);
    gsKL(:,1:nSamples) = gs(:,1:nSamples) - repmat(meanMode,1,nSamples);
    [U,Si,V] = svd(gsKL,0);
else
    [U,Si,V] = svd(gs,0);
end
conver = [];
modes = 0;
m = size(nonzeros(gt(nonzeros(Si),1e-04)),1);
for N = 1:m
    Y = U(:,1:N)*Si(1:N,1:N)*(V(:,1:N))';
    conver(N) = norm(gs-(Y+repmat(meanMode,1,nSamples)),'fro');
    if conver(N)>1e-06;
        modes = modes+1;
    end
end

%% Plot the SVD convergence
figure(1)
plot([1:m],log10(conver),'LineWidth',2);
xlabel('i','fontsize',14);
ylabel('log(||U - V_{i} \sigma_{i} {W_{i}}^T||_{Fr})','fontsize',14);
sigma = nonzeros(Si);
figure(2);
semilogy([1:modes],sigma(1:modes),'LineWidth',2);
xlabel('i','fontsize',14);
ylabel('log(\sigma_{i})','fontsize',14);
%print(gcf,'-deps2c',svd_conv)

%% Plot different modes (no physical interpretation)
for i = 1:8
    subplot(2,4,i);
    plot(U(:,i),S);
end

%% Making th PC model and cross validation partition
solution = cell(modes,1);
min_cv_error = [];
p=2
[X,PC] = PCTPMODEL(RV,'order',p,'pcg','typebase',2);
nnnmax = 10;
CVO = cvpartition(size(gs,2),'kfold',10)
trIdx = CVO.training(1);
%%
cv_error = zeros(nnnmax,m);
for kk = 1:3
    fs = V(:,kk);
    rs_reg = [], fs_reg = [];
    rs_test = [], fs_test = [];
    
    k = 1; l = 1;
    for j = 1:length(trIdx);
        if trIdx(j) == 1
            rs_reg(k,:) = rs(j,:);
            fs_reg(k) = fs(j);
            k = k+1;
        else
            rs_test(l,:) = rs(j,:);
            fs_test(l) = fs(j);
            l = l+1;
        end
    end
    %for rr = 1:length(regParam)
    
    global u;
    
    %cv_error = zeros(nnnmax,m);
    errstd=[];
    errmean=[];
    w = 1./length(fs_reg);
    relnorm = sqrt(sum(w.*fs_reg.^2));
    u = PCTPMATRIXSUM(PC);
    stols = cell(1,getnbgroups(PC));
    Ns = length(fs);
    for nnn=1:nnnmax
        %for nnn=nnnmax:nnnmax
        error = zeros(1,size(fs_test));
        col = getcourbestyles(nnn,'nomarker');
        SVDM = SEPARATEDSOLVER('onebyone',true,...
            'nbfoncmax',nnn,'display',true,'pfixmax',10,'pfixtol',1e-4,...
            'pfixstagn',1e6,'tol',1e-12,'update',0,'alphaupdate',0);
        
        [u,result,fstemp,rstemp] = ...
            separated_decomposition_fun_regression_fast(SVDM,PC,fs_reg,rs_reg,'regularization',1e-02,'regType','l1','lambda',1e-04);
        usto{nnn}=u;
        
        
        for i = 1:size(rs_test)
            error(i) = fs_test(i) - randomeval(u,[rs_test(i,:)]);
        end
        
        cv_error(nnn,kk) = sqrt((error*error')/length(fs_test))/relnorm;
        
        if nnn==1
            solution{kk} = u;
            min_cv_error(kk) = cv_error(nnn,kk);
        else
            if cv_error(nnn,kk) < min_cv_error(kk)
                solution{kk} = u;
                min_cv_error(kk) = cv_error(nnn,kk);
            end
        end
    end
end
%end
%% Plotting CV error for the first 8 modes
figure(110)
for i = 1:8
    ha(i) = subplot(2,4,i);
    plot([1:nnnmax], log10(cv_error([1:nnnmax],i)),'-.sb','LineWidth',2);
    xlabel('Rank','fontsize',12);
    ylabel('log(\epsilon_{CV})','fontsize',12);
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

