function Bsep = construire_sepvector(b,PC,pos)
%function Msep=construire_sepmatrix(Vecteur,PC globale,position)
onek = cell(1,getM(PC{1})+1);
Bsep = SEPMATRIX(getM(PC{1})+1);

for k=2:getM(PC{1})+1
    onek{k} = double(one(PC{k-1}))';
    onek{k} = {onek{k}};
end


if getM(PC{1})==4
    for ii=1:getm(b)
        pourcentage(ii,getm(b));
        %
        bi = double(getL(b,ii))';
        bi={bi};
        %
        if pos==1
            fprintf('Erreur : la position 1 est pour la fonction d�terministe, essayez 2')
            break;
        elseif pos==2
            Bsep = Bsep + SEPMATRIX([{getV(b,ii)},bi,onek{3},onek{4},onek{5}]);
        elseif pos==3
            Bsep = Bsep + SEPMATRIX([{getV(b,ii)},onek{2},bi,onek{4},onek{5}]);
        elseif pos==4
            Bsep = Bsep + SEPMATRIX([{getV(b,ii)},onek{2},onek{3},bi,onek{5}]);
        elseif pos==5
            Bsep = Bsep + SEPMATRIX([{getV(b,ii)},onek{2},onek{3},onek{4},bi]);
        end
    end
else fprintf('Pas programm�e')
end
end
