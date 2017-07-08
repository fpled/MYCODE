function [E_l_ident,E_t_ident,nu_lt_ident,G_lt_ident]...
    =NewtonRaphsonISO_TR(E_l_o,E_t_o,nu_lt_o,G_lt_o,U_exp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dans ce travail,on utilise la méthode Newton-Raphson pour               %
% identifier les coefficients élastique à partir des résultats            %
% expérimentaux.                                                          %
% Pour ce la, on utilise l'approximation de Taylor pour définir  et       %
% calculer numériquement les dérivés partielles de la fonction objectif:  %
%     f(E_l,E_t,nu_lt,G_lt)=||U_ef(E_l,E_t,nu_lt,G_lt)-U_exp||^2          %
% Puis,on utilise la méthode Newton-Raphson pour résoudre le probleme     %
% non-linéaire df(E_l,E_t,nu_lt,G_ltu)=0;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dE_l=1e-3;
dE_t=1e-5;
dnu_lt=1e-5;
dG_lt=1e-5;

Alpha_o=[E_l_o
    E_t_o
    nu_lt_o
    G_lt_o];
Alpha_t=ones(4,1);
k=0;
Phi_t=ones(4,1);

[Phi_o,Phi_p_o]=calculderivationISO_TR(E_l_o,E_t_o,nu_lt_o,G_lt_o,U_exp,dE_l,dE_t,dnu_lt,dG_lt);  

while abs(Phi_t(1)) > 1e-10 || abs(Phi_t(2)) > 1e-11|| abs(Phi_t(3)) > 1e-11 || abs(Phi_t(4)) > 1e-11
    
Alpha_t= Alpha_o-inv(Phi_p_o)*Phi_o

[Phi_t,Phi_p_t]=calculderivationISO_TR(Alpha_t(1),Alpha_t(2),Alpha_t(3),Alpha_t(4),U_exp,dE_l,dE_t,dnu_lt,dG_lt);

    if abs(Phi_t(1)) > 1e-10 || abs(Phi_t(2)) > 1e-11|| abs(Phi_t(3)) > 1e-11 || abs(Phi_t(4)) > 1e-11
        Phi_o=Phi_t;
        Phi_p_o=Phi_p_t;
        Alpha_o=Alpha_t;
    else
    end
    k=k+1  %calcul no de boucle
end
E_l_ident=Alpha_t(1);
E_t_ident=Alpha_t(2);
nu_lt_ident=Alpha_t(3);
G_lt_ident=Alpha_t(4);
end