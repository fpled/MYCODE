function x = functionIdentificationIsot(U_exp,noise)
E0 = 10;
nu0 = 0.2;
x0 = [E0 nu0];
% lb = [0 -1];
lb = [0 0];
ub = [Inf 0.5];

tolX = 1e-14;
tolFun = 1e-14;
display = 'off';
optionslsqnonlin  = optimoptions('lsqnonlin','Display',display,'TolX',tolX,'TolFun',tolFun);

funlsqnonlin = @(x) funlsqnonlinIsot(x,U_exp,noise);

[x,~,~,~,~] = lsqnonlin(@(x) funlsqnonlin(x),x0,lb,ub,optionslsqnonlin);

end
