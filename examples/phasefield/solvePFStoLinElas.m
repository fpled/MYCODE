function [f_sample,d_sample,u_sample] = solvePFStoLinElas(S_phase,S,T,fun,samples,varargin)
% function [f_sample,d_sample,u_sample] = solvePFStoLinElas(S_phase,S,T,fun,samples,varargin)
% Solve stochastic Phase Field problem.

fun = fcnchk(fun);
nbSamples = getcharin('nbsamples',varargin,3);

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
d_sample = zeros(nbSamples,sz_d,length(T));
u_sample = zeros(nbSamples,sz_u,length(T));
% fmax_sample = zeros(N,1);

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
    f_sample(i,:) = ft;
    if i<=nbSamples
        d_sample(i,:,:) = getvalue(dt);
        u_sample(i,:,:) = getvalue(ut);
    end
    % fmax_sample(i) = max(ft);
end
textprogressbar(' done');

function nUpdateProgressBar(~)
j = j+1;
textprogressbar(j/N*100,sprintf('(%d/%d)',j,N));
end

end
