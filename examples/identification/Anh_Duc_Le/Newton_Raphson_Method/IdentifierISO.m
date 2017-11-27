%% NEWTON_RAPHSON POUR LE PROBLEME ISOTROPE

% clc
clearvars
close all 

%% Generer un champ de deplacement experimental

E_exp=11.550;
nu_exp=0.4;
U_vec=flexion3pISO(E_exp,nu_exp);
p=1;
Amplitude=(p/100)*U_vec;
Bruit=Amplitude.*(2*rand(length(U_vec),1)-1);
U_exp=U_vec+Bruit;


%% Identifier les coefficients elastique par la methode Newton-Raphson

E_o=10.000;
nu_o=0.2;

[E_ident,nu_ident]=NewtonRaphsonISO(E_o,nu_o,U_exp)
