function S_phase = addbcdamageSingleEdgeCrack(S_phase,C,initialCrack)
% function S_phase = addbcdamageSingleEdgeCrack(S_phase,C,initialCrack)
% Add boundary conditions on damage/phase field for single edge crack problem.

if strcmpi(initialCrack,'initialphasefield')
    S_phase = addcl(S_phase,C,'T',1);
end

end