function [Apatch,Apatchsep] = calc_matrix_dirichlet_sep(psitp,PC,Spart)
%function [Apatch{k},Apatchsep{k}]=calc_matrix_dirichlet_sep(psitp{k},PC,Spart(fin){k})


Apatch = setphi(PCTPMATRIX(PC,zeros(getnbddlfree(Spart),getnbddlfree(Spart))),calc_ximasse(one(PC)));

VS = getfuns(psitp);
for i=1:getm(psitp)
    
    a1 = MULTILINFORM([1,1,0,0]);
    a2 = MULTILINFORM([0,1,0,1]);
    a3 = MULTILINFORM([1,0,1,0]);
    a4 = MULTILINFORM([0,0,1,1]);
    
    %Li = getphi(VS{i},i);
    %Li = PCMATRIX(Li,[size(Li,2),1],getpcgroup(PCphi,i));
    Li = getphi(VS{i});
    Li = setphi(PCTPMATRIX(PC,1),Li);
    for j=1:getm(psitp)
        
        L2=PCTPMATRIX(PC,1);
        Lj = getphi(VS{j});
        Lj = setphi(L2,Lj);
        lilj = calc_ximasse(Li)* calc_ximasse(Lj);
        
        
        
        Apatch = Apatch + a1{Spart}(:,:,getphi0(VS{i}),getphi0(VS{j}))*lilj+...
            a2{Spart}(:,:,getphi0(VS{i}),getphi0(VS{j}))*lilj+...
            a3{Spart}(:,:,getphi0(VS{i}),getphi0(VS{j}))*lilj+...
            a4{Spart}(:,:,getphi0(VS{i}),getphi0(VS{j}))*lilj;
    end
end
Apatchsep = SEPMATRIX(calc_ximasse(Apatch));


