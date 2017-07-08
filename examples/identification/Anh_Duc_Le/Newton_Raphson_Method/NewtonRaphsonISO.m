function [E_ident,nu_ident]=NewtonRaphsonISO(E_o,nu_o,U_exp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dans ce travail,on utilise la méthode Newton-Raphson pour               %
% identifier les coefficients élastique à partir des résultats            %
% expérimentaux.                                                          %
% Pour ce la, on utilise l'approximation de Taylor pour définir  et       %
% calculer numériquement les dérivés partielles de la fonction objectif   %
% f(E,nu)                                                                 %
%                 f(E,nu)=||U_ef(E,nu)-U_exp||^2                          % 
% Puis,on utilise la méthode Newton-Raphson pour résoudre le probleme     %
% non-linéaire df(E,nu)=0;                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_E=1e-3;
d_nu=1e-5;
Alpha_o=[E_o
    nu_o];
Alpha_t=ones(2,1);
k=0;
Phi_t=ones(2,1);

[Phi_o,Phi_p_o]=calculderivationISO(E_o,nu_o,U_exp,d_E,d_nu);  

while abs(Phi_t(1)) > 1e-9|| abs(Phi_t(2)) > 1e-10
    
Alpha_t= Alpha_o-Phi_p_o\Phi_o

[Phi_t,Phi_p_t]=calculderivationISO(Alpha_t(1),Alpha_t(2),U_exp,d_E,d_nu);

    if abs(Phi_t(1)) > 1e-9 || abs(Phi_t(2)) > 1e-10
        Phi_o=Phi_t;
        Phi_p_o=Phi_p_t;
        Alpha_o=Alpha_t;
    else
    end
    k=k+1  %calcul no de boucle
end
E_ident=Alpha_t(1);
nu_ident=Alpha_t(2);
end