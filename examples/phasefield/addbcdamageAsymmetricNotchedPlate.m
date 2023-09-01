function S_phase = addbcdamageAsymmetricNotchedPlate(S_phase,C,BU,BL,BR,initialCrack)
% function S_phase = addbcdamageAsymmetricNotchedPlate(S_phase,C,BU,BL,BR,initialCrack)
% Add boundary conditions on damage/phase field for asymmetric notched plate problem.

if strcmpi(initialCrack,'initialphasefield')
    S_phase = addcl(S_phase,C,'T',1);
end
S_phase = addcl(S_phase,BU,'T');
S_phase = addcl(S_phase,BL,'T');
S_phase = addcl(S_phase,BR,'T');

end