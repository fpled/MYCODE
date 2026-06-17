function S_phase = addbcdamageThreePointBendingAdaptive(S_phase,C,BU,BL,BR)
% function S_phase = addbcdamageThreePointBendingAdaptive(S_phase,C,BU,BL,BR)
% Add boundary conditions on damage/phase field for mesh adaptation for three-point bending  problem.

S_phase = addcl(S_phase,C,'T',1);
S_phase = addcl(S_phase,BU,'T',1);
S_phase = addcl(S_phase,BL,'T',1);
S_phase = addcl(S_phase,BR,'T',1);

end