function [udeter,Adeter,bdeter] = calcsoldeter_patch(S,Spart,P,kappa,r,varargin)
%function [udeter,Adeter,bdeter] = calcsoldeter_patch(S,Spart,P,kappa,varargin)

kr = cell(1,length(kappa));
kr{1} = 1;

for k = 2:length(kappa)
    kr{k} = kappa{k}(r);
end

Adeter = zeros(getnbddlfree(S),getnbddlfree(S));
adeter = MULTILINFORM([1,1,0]);
ldeter = LINFORM(0);
for k = 1:length(kappa)
    Atemp = adeter{Spart{k}}(:,:,kr{k});
    Adeter = Adeter + P{k}'*Atemp*P{k};
    bdeter = ldeter{S}(:);
end

udeter = solve(Adeter,bdeter);
udeter=unfreevector(S,udeter);
