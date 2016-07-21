function A = calc_pcmatrix_dirichlet(S,PC,psi)
%function A1 = calc_pcmatrix_dirichlet(S1,PCphi1,psi)
%S : Model, PC:= Polychaos , psi := fonctions caract�ristique.  

A = PCRADIALMATRIX([getnbddlfree(S),getnbddlfree(S)],PC);
L = calc_masse(getL(project(psi,PC)));
%L = calc_masse(getL(project(psi,PC)));      
for i=1:getm(psi)
%for i=1:getm(psi)

    a1 = MULTILINFORM([1,1,0,0]);
    a2 = MULTILINFORM([0,1,0,1]);
    a3 = MULTILINFORM([1,0,1,0]);
    a4 = MULTILINFORM([0,0,1,1]);


    for j=1:getm(psi)

        lilj = (L(i))*L(j);

        A = A + a1{S}(:,:,getV(psi,i),getV(psi,j))*lilj+...
            a2{S}(:,:,getV(psi,i),getV(psi,j))*lilj+...
            a3{S}(:,:,getV(psi,i),getV(psi,j))*lilj+...
            a4{S}(:,:,getV(psi,i),getV(psi,j))*lilj;
    end
end

A = calc_ximasse(A);
