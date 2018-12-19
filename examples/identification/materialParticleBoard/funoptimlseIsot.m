function f = funoptimlseIsot(lambda,mC_data,nuC_data,varargin)
% function f = funoptimlseIsot(lambda,mC_data,nuC_data,varargin)

la1 = lambda(1);
la2 = lambda(2);
la  = lambda(3);

a1 = 1-la;
b1 = 1/la1;
a2 = 1-5*la;
b2 = 1/la2;

mC1 = a1*b1;
mC2 = a2*b2;
mC = [mC1 mC2];

% vC1 = a1*b1^2;
% vC2 = a2*b2^2;
% vC = [vC1 vC2];
% sC = sqrt(norm(vC));
% dC = sC/norm(mC);

% phi_C = ln(96*C1*C2^5) = ln(96) + ln(C1) + 5*ln(C2)
% nu_C = E{phi_C}  = E{ln(96*C1*C2^5)} = ln(96) + E{ln(C1)} + 5*E{ln(C2)}
% with Gamma-distributed random variables C1~Gamma(a1,b1) and C2~Gamma(a2,b2)
nuC = log(96) + psi(a1)+log(b1) + 5*(psi(a2)+log(b2));

f = [mC nuC] - [mC_data nuC_data];

end
