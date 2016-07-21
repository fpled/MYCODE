function [uref,Aref,bref] = calcref_patch_sep(Apatch,bpatch,P)
%function [uref,Aref,bref] = calcref_patch(Apatch,bpatch,P,tol)
% Apatch,bpatch et P sont des cell

Aref = mtimes(mtimes(P{1}',Apatch{1},1),P{1},1);
bref = mtimes(P{1}',bpatch{1},1);
for k=2:length(Apatch)
    Aref = Aref + mtimes(mtimes(P{k}',Apatch{k},1),P{k},1);
    bref = bref + mtimes(P{k}',bpatch{k},1);
end

solver = SEPSOLVER(getdim(Aref),'tol',1e-10,...
    'update',0,'updatedim',2:getdim(Aref),...
    'maxorder',100,'maxiter',5,'reference',[],...
    'errorindicator','none','itercrit',1e-3,...
    'righthandSD',false,'display',true);

uref = solve(Aref,bref,solver);

