function f = calc_coefdif_exp(c1,c2,xi,yi,lc,A)
% f = calc_coefdif_exp(c1,c2,xi,yi,lc,A)
% Calcul la partie exponnentielle d'un coefficient de difusion 
% (c1,c2) est le centre du domaine.
% (xi,yi) les coordonnees des noeuds.
% lc : longeur carac
% A : Amplitude

f = exp(-A*((xi-c1).^2 + (yi-c2).^2)/lc^2);
end
