function [f_sample,d_sample,u_sample,H_sample,S_phase_sample,S_sample] = solvePFStoLinElasAdaptive(S_phase,S,T,fun,samples,varargin)
% function [f_sample,d_sample,u_sample,H_sample,S_phase_sample,S_sample] = solvePFStoLinElasAdaptive(S_phase,S,T,fun,samples,varargin)
% Solve stochastic Phase Field problem with mesh adaptation.

fun = fcnchk(fun);
filename = getcharin('filename',varargin,'filename');
pathname = getcharin('pathname',varargin,'.');
nbSamples = getcharin('nbsamples',varargin,3);

N = size(samples,1);
E_sample = samples(:,1);
NU_sample = samples(:,2);
gc_sample = samples(:,3);
l_sample = samples(:,4);

mats_phase = MATERIALS(S_phase);
mats = MATERIALS(S);

f_sample = zeros(N,length(T));
d_sample = cell(nbSamples,length(T));
u_sample = cell(nbSamples,length(T));
H_sample = cell(nbSamples,length(T));
S_phase_sample = cell(N,length(T));
S_sample = cell(N,length(T));
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
    
    % Copy msh, mesh and sol files
    filenamei = [filename '_' num2str(i)];
    G = GMSHFILE(fullfile(pathname,filename));
    Gi = GMSHFILE(fullfile(pathname,filenamei));
    command = ['cp ' getfilemsh(G) ' ' getfilemsh(Gi) ';'...
        'cp ' getfilemesh(G) ' ' getfilemesh(Gi) ';'...
        'cp ' getfilesol(G) ' ' getfilesol(Gi)];
    dos(command);
    
    % Solve deterministic problem
    [dt,ut,ft,Ht,St_phase,St] = fun(S_phasei,Si,filenamei);
    f_sample(i,:) = ft;
    if i<=nbSamples
        d_sample(i,:) = dt;
        u_sample(i,:) = ut;
        H_sample(i,:) = Ht;
        S_phase_sample(i,:) = St_phase;
        S_sample(i,:) = St;
    end
    % fmax_sample(i) = max(ft);
end
textprogressbar(' done');

function nUpdateProgressBar(~)
j = j+1;
textprogressbar(j/N*100,sprintf('(%d/%d)',j,N));
end

end
