function [f_sample,dt_mean,ut_mean,dt_var,ut_var] = solvePFStoLinElas(S_phase,S,T,fun,samples,varargin)
% function [f_sample,dt_mean,ut_mean,dt_var,ut_var] = solvePFStoLinElas(S_phase,S,T,fun,samples,varargin)
% Solve stochastic Phase Field problem.

fun = fcnchk(fun);

N = size(samples,1);
E_sample = samples(:,1);
NU_sample = samples(:,2);
gc_sample = samples(:,3);
l_sample = samples(:,4);

mats_phase = MATERIALS(S_phase);
mats = MATERIALS(S);

sz_d = getnbddl(S_phase);
sz_u = getnbddl(S);

f_sample = zeros(N,length(T));
% fmax_sample = zeros(N,1);

% Initialize statistical means and second-order moments
dt_mean = zeros(sz_d,length(T));
dt_moment2 = zeros(sz_d,length(T));
ut_mean = zeros(sz_u,length(T));
ut_moment2 = zeros(sz_u,length(T));

if ~verLessThan('matlab','9.2') % introduced in R2017a
    q = parallel.pool.DataQueue;
    afterEach(q, @nUpdateProgressBar);
else
    q = [];
end
textprogressbar('Solving problem: ');
j = 0;
parfor i=1:N
    if ~verLessThan('matlab','9.2') % introduced in R2017a
        send(q,i);
    end
    
    % Update phase field parameters
    gci = gc_sample(i);
    li = l_sample(i);
    S_phasei = S_phase;
    mats_phasei = mats_phase;
    for m=1:length(mats_phasei)
        mats_phasei{m} = setparam(mats_phasei{m},'k',gci*li);
        mats_phasei{m} = setparam(mats_phasei{m},'r',gci/li);
    end
    S_phasei = actualisematerials(S_phasei,mats_phasei);
    
    % Update material parameters
    Ei = E_sample(i);
    NUi = NU_sample(i);
    Si = S;
    matsi = mats;
    for m=1:length(matsi)
        matsi{m} = setparam(matsi{m},'E',Ei);
        matsi{m} = setparam(matsi{m},'NU',NUi);
    end
    Si = actualisematerials(Si,matsi);
    
    % Solve deterministic problem
    [dt,ut,ft] = fun(S_phasei,Si);
    
    % Compute second-order statistics
    dt_val = getvalue(dt);
    dt_mean = dt_mean + dt_val/N;
    dt_moment2 = dt_moment2 + dt_val.^2/N;
    ut_val = getvalue(ut);
    ut_mean = ut_mean + ut_val/N;
    ut_moment2 = ut_moment2 + ut_val.^2/N;
    
    f_sample(i,:) = ft;
    % fmax_sample(i) = max(ft);
end
textprogressbar(' done');

% Compute unbiased variances
dt_var = (N/(N-1))*(dt_moment2 - dt_mean.^2);
ut_var = (N/(N-1))*(ut_moment2 - ut_mean.^2);

function nUpdateProgressBar(~)
j = j+1;
textprogressbar(j/N*100,sprintf('(%d/%d)',j,N));
end

end
