function [A2patch,A2patchsep]=calc_matrix_dirichlet_dual_sep(psitp,PC,Spart)
%function [A2patch{k},A2patchsep{k}]=calc_matrix_dirichlet_dual_sep(psitp{k},PC,Spart(fin){k})
VS = getfuns(psitp);  
A2patch = setphi(PCTPMATRIX(PC,zeros(getnbddlfree(Spart),getnbddlfree(Spart))),calc_ximasse(one(PC)));
for j=1:getm(psitp)
    l1 = MULTILINFORM([1,1,0]);
    l2 = MULTILINFORM([0,1,1]);
    A2patch = A2patch+l1{Spart}(:,:,getphi0(VS{j}))*setphi(PCTPMATRIX(PC,1),getphi(VS{j}))+l2{Spart}(:,:,getphi0(VS{j}))*setphi(PCTPMATRIX(PC,1),getphi(VS{j}));
end
A2patch = calc_ximasse(A2patch);
A2patchsep = SEPMATRIX(A2patch);



