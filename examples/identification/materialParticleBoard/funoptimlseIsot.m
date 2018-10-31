function f = funoptimlseIsot(lambda,mC_data,dC_data,varargin)
% function f = funoptimlseIsot(lambda,mC_data,dC_data,varargin)

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

vC1 = a1*b1^2;
vC2 = a2*b2^2;
vC = [vC1 vC2];
sC = sqrt(norm(vC));
dC = sC/norm(mC);

alpha = 1/2;
f = alpha * norm(mC - mC_data)^2/norm(mC_data)^2 + (1-alpha) * (dC - dC_data)^2/(dC_data)^2;

end
