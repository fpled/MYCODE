function S_phase = addbcdamageThreePointBending(S_phase,C,initialCrack)
% function S_phase = addbcdamageThreePointBending(S_phase,C,initialCrack)
% Add boundary conditions on damage/phase field for three-point bending problem.

if strcmpi(initialCrack,'initialphasefield')
    S_phase = addcl(S_phase,C,'T',1);
end

end