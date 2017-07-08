function [Phi_o,Phi_p_o]=calculderivationISO(E_o,nu_o,U_exp,d_E,d_nu)

f_E_o_nu_o=sum(((flexion3pISO(E_o,nu_o)-U_exp).^2),1);
f_dE=sum(((flexion3pISO(E_o+d_E,nu_o)-U_exp).^2),1);
f_dnu=sum(((flexion3pISO(E_o,nu_o+d_nu)-U_exp).^2),1);

f_m_dE=sum(((flexion3pISO(E_o-d_E,nu_o)-U_exp).^2),1);
f_m_dnu=sum(((flexion3pISO(E_o,nu_o-d_nu)-U_exp).^2),1);

f_dE_dnu=sum(((flexion3pISO(E_o+d_E,nu_o+d_nu)-U_exp).^2),1);

Derivation_dE=(f_dE-f_E_o_nu_o)/d_E;
Derivation_dnu=(f_dnu-f_E_o_nu_o)/d_nu;
Derivation_2_dE_2=(f_m_dE-2*f_E_o_nu_o+f_dE)/(d_E^2);
Derivation_2_dnu_2=(f_m_dnu-2*f_E_o_nu_o+f_dnu)/(d_nu^2);
Derivation_2_dE_dnu=((f_dE_dnu-f_dnu)-(f_dE-f_E_o_nu_o))/(d_E*d_nu);

Phi_o=[Derivation_dE;Derivation_dnu];
Phi_p_o=[Derivation_2_dE_2 Derivation_2_dE_dnu
    Derivation_2_dE_dnu Derivation_2_dnu_2];
end