function [ur,Si] = calc_sol_dirichlet_deter(S,lsi)
% function [ur] = calc_sol_dirichlet_deter(S,ls)
%S : Model
%lsi: levelset(r)

Si = LSMODEL(setmaterial(S,FOUR_ISOT()),LEVELSET(lsi));
Si = final(Si);
Si = addcl(Si,[],'T');

lsval =  getvalue(lseval(lsi,Si));
psir = -lsval;
%psir2 = psir.^2;

%a1 = BILINFORM(1,1,psi2,0);
a1 = MULTILINFORM([1,1,0,0]);
%a2 = BILINFORM(0,1,(psi2)./2,1);
a2 = MULTILINFORM([0,1,1,0]);
%a3 = BILINFORM(1,0,(psi2)./2,1);
a3 = MULTILINFORM([1,0,0,1]);
a4 = MULTILINFORM([0,0,1,1]);

%Ai = a1{S}(:,:) + a2{S}(:,:) + a3{S}(:,:) + a4{S}(:,:,psi,psi);
Ai = a1{Si}(:,:,psir,psir) + a2{Si}(:,:,psir,psir) + a3{Si}(:,:,psir,psir) + a4{Si}(:,:,psir,psir);
%
li = LINFORM(0,psir,0);
bi = li{Si}(:);
%
udetr = solve(Ai,bi);
udetr=unfreevector(Si,udetr);
ur = udetr.*psir;
return 
