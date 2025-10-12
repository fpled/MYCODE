function S_phase = addbcdamageDoubleEdgeCrack(S_phase,Ca,Cb,initialCrack)
% function S_phase = addbcdamageDoubleEdgeCrack(S_phase,Ca,Cb,initialCrack)
% Add boundary conditions on damage/phase field for double edge crack problem.

if strcmpi(initialCrack,'initialphasefield')
    S_phase = addcl(S_phase,Ca,'T',1);
    S_phase = addcl(S_phase,Cb,'T',1);
end

end