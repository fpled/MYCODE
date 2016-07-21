function Msep = construire_sepmatrix(M,PC,pos)
%function Msep=construire_sepmatrix(Matrice,PC globale,position)
onek = cell(1,getM(PC{1})+1);
Msep = SEPMATRIX(getM(PC{1})+1);

for k=2:getM(PC{1})+1  
    onek{k} = get_ximasse(calc_ximasse(one(PC{k-1})));
    onek{k} = {onek{k}{1}};
end

if (getM(PC{1})==4)
    for ii=1:getm(M)
        pourcentage(ii,getm(M));
%
        Mi = getximasse(calc_ximasse(getL(M,ii),PC{pos-1}));
        Mi = {Mi{1}};
%  
        if pos==1 
            fprintf('Erreur : la position un est pour la fonction d�terministe, essayez 2')
            break;
        elseif pos==2
            Msep = Msep + SEPMATRIX([{getV(M,ii)},Mi,onek{3},onek{4},onek{5}]); 
        elseif pos==3
            Msep = Msep + SEPMATRIX([{getV(M,ii)},onek{2},Mi,onek{4},onek{5}]); 
        elseif pos==4
            Msep = Msep + SEPMATRIX([{getV(M,ii)},onek{2},onek{3},Mi,onek{5}]); 
        elseif pos==5
            Msep = Msep + SEPMATRIX([{getV(M,ii)},onek{2},onek{3},onek{4},Mi]); 
        end 
    end 
else fprintf('Pas programm�e')
end   
end
