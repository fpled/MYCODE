function bpc = calc_pcvector_dirichlet(S,PC,psi)

b = PCRADIALMATRIX([getnbddlfree(S),1],PC);
psi=project(psi,PC);
for i=1:getm(psi)
    l = LINFORM(0,getV(psi,i),0);
    b = b+l{S}(:)*getL(psi,i);
end
mas = getmasse(calc_masse(PC));
mas = multisum(mas);
mas = speye(size(mas));

bpc = expand(b);
bpc = double(bpc);
bpc = bpc*mas;
bpc = PCMATRIX(bpc,[size(b,1),1],PC);
%bpc = spectral_decomposition(bpc,'tol',1e-8);
end
