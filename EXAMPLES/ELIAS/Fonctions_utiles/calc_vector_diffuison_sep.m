function [bpatch,bpatchsep]=calc_vector_diffuison_sep(PC,Spart)
%function [bpatch{k},bpatchsep{k}]=calc_vector_diffuison_sep(PC,Spart(fin){k})
lk = LINFORM(0);
bpatch = lk{Spart}(:)*one(PC);
bpatchsep = SEPVECTOR(bpatch);
