function S_phase = addbcdamageAsymmetricNotchedPlateAdaptive(S_phase,C,H1,H2,H3)
% function S_phase = addbcdamageAsymmetricNotchedPlateAdaptive(S_phase,C,H1,H2,H3)
% Add boundary conditions on damage/phase field for mesh adaptation for asymmetric notched plate problem.

S_phase = addcl(S_phase,C,'T',1);
S_phase = addcl(S_phase,H1,'T',1);
S_phase = addcl(S_phase,H2,'T',1);
S_phase = addcl(S_phase,H3,'T',1);

end