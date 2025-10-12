function S_phase = addbcdamageDoubleEdgeCrackAdaptive(S_phase,Ca,Cb)
% function S_phase = addbcdamageDoubleEdgeCrackAdaptive(S_phase,Ca,Cb)
% Add boundary conditions on damage/phase field for mesh adaptation for double edge crack problem.

S_phase = addcl(S_phase,Ca,'T',1);
S_phase = addcl(S_phase,Cb,'T',1);

end