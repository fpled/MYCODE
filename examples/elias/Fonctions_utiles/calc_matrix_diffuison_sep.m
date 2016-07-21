function [Apatch,Apatchsep]=calc_matrix_diffuison_sep(ktp,PC,Spart)
%function [Apatch{k},Apatchsep{k}]=calc_matrix_diffuison_sep(ktp{k},PC,Spart(fin){k})

L1 = setphi(PCTPMATRIX(PC,1),getphi(ktp{1}));
ak = BILINFORM(1,1,getphi0(ktp{1}),0);
Apatch = ak{Spart}(:,:)*L1;
for i=2:getm(ktp)
    Li = setphi(PCTPMATRIX(PC,1),getphi(ktp{i}));
    ak = BILINFORM(1,1,getphi0(ktp{i}),0);
    Apatch = Apatch+ak{Spart}(:,:)*Li;
end
Apatch = calc_ximasse(Apatch);
Apatchsep = SEPMATRIX(Apatch);
