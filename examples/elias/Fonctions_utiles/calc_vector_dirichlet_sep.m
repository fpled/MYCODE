function [bpatch,bpatchsep]=calc_vector_dirichlet_sep(psitp,PC,Spart,mas)
%function [bpatch{k},bpatchsep{k}]=calc_vector_dirichlet_sep(psitp{k},PC,Spart(fin){k},mas)

bpatch = setphi(PCTPMATRIX(PC,zeros(getnbddlfree(Spart),1)),one(PC));
Lb = PCTPMATRIX(PC,1);


% mas = getmasse(calc_masse(PC));
% mas = multisum(mas);
% mas = speye(size(mas)); 
psitp = getfuns(psitp);  
for i=1:length(psitp)

    Lbi = getphi(psitp{i});
%   for j=1 : getM(mas)
%       Lbi{j} = mas{j}*Lbi{j};
%   end
    Lbi = setphi(Lb,Lbi);  

    l = LINFORM(0,getphi0(psitp{i}),0);
    bpatch = bpatch + l{Spart}(:)*Lbi;
end
bpatchsep = SEPVECTOR(bpatch);


