function [Phi_o,Phi_p_o]=calculderivationISO_TR(E_l_o,E_t_o,nu_lt_o,G_lt_o,U_exp,dE_l,dE_t,dnu_lt,dG_lt)

f_o=sum(((flexion3pISO_TR(E_l_o,E_t_o,nu_lt_o,G_lt_o)-U_exp).^2),1);
f_dE_l=sum(((flexion3pISO_TR(E_l_o+dE_l,E_t_o,nu_lt_o,G_lt_o)-U_exp).^2),1);
f_dE_t=sum(((flexion3pISO_TR(E_l_o,E_t_o+dE_t,nu_lt_o,G_lt_o)-U_exp).^2),1);
f_dnu_lt=sum(((flexion3pISO_TR(E_l_o,E_t_o,nu_lt_o+dnu_lt,G_lt_o)-U_exp).^2),1);
f_dG_lt=sum(((flexion3pISO_TR(E_l_o,E_t_o,nu_lt_o,G_lt_o+dG_lt)-U_exp).^2),1);

f_m_dE_l=sum(((flexion3pISO_TR(E_l_o-dE_l,E_t_o,nu_lt_o,G_lt_o)-U_exp).^2),1);
f_m_dE_t=sum(((flexion3pISO_TR(E_l_o,E_t_o-dE_t,nu_lt_o,G_lt_o)-U_exp).^2),1);
f_m_dnu_lt=sum(((flexion3pISO_TR(E_l_o,E_t_o,nu_lt_o-dnu_lt,G_lt_o)-U_exp).^2),1);
f_m_dG_lt=sum(((flexion3pISO_TR(E_l_o,E_t_o,nu_lt_o,G_lt_o-dG_lt)-U_exp).^2),1);


f_dE_l_dE_t=sum(((flexion3pISO_TR(E_l_o+dE_l,E_t_o+dE_t,nu_lt_o,G_lt_o)-U_exp).^2),1);
f_dE_l_dnu_lt=sum(((flexion3pISO_TR(E_l_o+dE_l,E_t_o,nu_lt_o+dnu_lt,G_lt_o)-U_exp).^2),1);
f_dE_l_dG_lt=sum(((flexion3pISO_TR(E_l_o+dE_l,E_t_o,nu_lt_o,G_lt_o+dG_lt)-U_exp).^2),1);
f_dE_t_dnu_lt=sum(((flexion3pISO_TR(E_l_o,E_t_o+dE_t,nu_lt_o+dnu_lt,G_lt_o)-U_exp).^2),1);
f_dE_t_dG_lt=sum(((flexion3pISO_TR(E_l_o,E_t_o+dE_t,nu_lt_o,G_lt_o+dG_lt)-U_exp).^2),1);
f_dnu_lt_dG_lt=sum(((flexion3pISO_TR(E_l_o,E_t_o,nu_lt_o+dnu_lt,G_lt_o+dG_lt)-U_exp).^2),1);

Derivation_dE_l=(f_dE_l-f_o)/dE_l;
Derivation_dE_t=(f_dE_t-f_o)/dE_t;
Derivation_dnu_lt=(f_dnu_lt-f_o)/dnu_lt;
Derivation_dG_lt=(f_dG_lt-f_o)/dG_lt;

Derivation_2_dE_l_2=(f_m_dE_l-2*f_o+f_dE_l)/(dE_l^2);
Derivation_2_dE_t_2=(f_m_dE_t-2*f_o+f_dE_t)/(dE_t^2);
Derivation_2_dnu_lt_2=(f_m_dnu_lt-2*f_o+f_dnu_lt)/(dnu_lt^2);
Derivation_2_dG_lt_2=(f_m_dG_lt-2*f_o+f_dG_lt)/(dG_lt^2);

Derivation_2_dE_l_dE_t=((f_dE_l_dE_t-f_dE_t)-(f_dE_l-f_o))/(dE_l*dE_t);
Derivation_2_dE_l_dnu_lt=((f_dE_l_dnu_lt-f_dnu_lt)-(f_dE_l-f_o))/(dE_l*dnu_lt);
Derivation_2_dE_l_dG_lt=((f_dE_l_dG_lt-f_dG_lt)-(f_dE_l-f_o))/(dE_l*dG_lt);
Derivation_2_dE_t_dnu_lt=((f_dE_t_dnu_lt-f_dnu_lt)-(f_dE_t-f_o))/(dE_t*dnu_lt);
Derivation_2_dE_t_dG_lt=((f_dE_t_dG_lt-f_dG_lt)-(f_dE_t-f_o))/(dE_t*dG_lt);
Derivation_2_dnu_lt_dG_lt=((f_dnu_lt_dG_lt-f_dG_lt)-(f_dnu_lt-f_o))/(dnu_lt*dG_lt);


Phi_o=[Derivation_dE_l;Derivation_dE_t;Derivation_dnu_lt;Derivation_dG_lt];

Phi_p_o=[Derivation_2_dE_l_2 Derivation_2_dE_l_dE_t Derivation_2_dE_l_dnu_lt Derivation_2_dE_l_dG_lt
    Derivation_2_dE_l_dE_t Derivation_2_dE_t_2 Derivation_2_dE_t_dnu_lt Derivation_2_dE_t_dG_lt
    Derivation_2_dE_l_dnu_lt Derivation_2_dE_t_dnu_lt Derivation_2_dnu_lt_2 Derivation_2_dnu_lt_dG_lt
    Derivation_2_dE_l_dG_lt Derivation_2_dE_t_dG_lt Derivation_2_dnu_lt_dG_lt Derivation_2_dG_lt_2];

end