%% NEWTON_RAPHSON POUR LE PROBLEME ISOTROPE TRANVERSE

% clc
clearvars
close all

%% Generer un champ de deplacement experimental

E_l_exp=11.550;
E_t_exp=0.500;
nu_lt_exp=0.4;
G_lt_exp=0.550;
U_vec=flexion3pISO_TR(E_l_exp,E_t_exp,nu_lt_exp,G_lt_exp);
p=1;
Amplitude=(p/100)*U_vec;
Bruit=Amplitude.*(2*rand(length(U_vec),1)-1);
U_exp=U_vec+Bruit;

%% Identifier les coefficients elastique par la methode Newton-Raphson

E_l_o=10.000;
E_t_o=0.3;
nu_lt_o=0.3;
G_lt_o=0.3;

[E_l_ident,E_t_ident,nu_lt_ident,G_lt_ident]= NewtonRaphsonISO_TR(E_l_o,E_t_o,nu_lt_o,G_lt_o,U_exp)
