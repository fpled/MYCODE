function [uref,Aref,bref] = calcref_patch(Apatch,bpatch,P,tol)
% function [uref,Aref,bref] = calcref_patch(Apatch,bpatch,P,tol)
% Apatch,bpatch et P sont des cell

Aref = sparse(size(P{1},2),size(P{1},2));%PCRADIALMATRIX([getnbddlfree(S),getnbddlfree(S)],PC);
bref = sparse(size(P{1},2),1);%PCRADIALMATRIX([getnbddlfree(S),1],PC);
for k=1:length(Apatch)
    Aref = Aref + P{k}'*Apatch{k}*P{k} ;
    bref = bref + P{k}'*bpatch{k};
end

uref = solve(Aref,bref);

%uref = cgs(Aref,bref,tol);
