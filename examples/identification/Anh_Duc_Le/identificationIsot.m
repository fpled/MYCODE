%% Identification of isotropic elastic properties %%
%%------------------------------------------------%%

% clc
clearvars
close all
% rng('default');

%% Experimental values
E_exp = 11.55;
nu_exp = 0.4;
x_exp = [E_exp nu_exp];
U_exp = solveThreePointBendingIsot(x_exp);
fprintf('Exp values\n');
fprintf('----------\n');
fprintf('E_exp  = %g GPa\n',x_exp(1));
fprintf('nu_exp = %g\n',x_exp(2));

%% Identification
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
optionsfminsearch = optimset('Display',display,'TolX',tolX,'TolFun',tolFun);
optionsfminunc    = optimoptions('fminunc','Display',display,'TolX',tolX,'TolFun',tolFun);
optionsfmincon    = optimoptions('fmincon','Display',display,'TolX',tolX,'TolFun',tolFun);

ampl = 1e-2;
noise = ampl.*(2*rand(length(U_exp),1)-1).*U_exp;
% noise = 0;

funlsqnonlin = @(x) funlsqnonlinIsot(x,U_exp,noise);
funoptim = @(x) funoptimIsot(x,U_exp,noise);

t1 = tic;
[x1,err1,~,exitflag1,output1] = lsqnonlin(funlsqnonlin,x0,lb,ub,optionslsqnonlin);
fprintf('\nlsqnonlin');
fprintf('\n---------\n');
fprintf('E   = %g GPa\n',x1(1));
fprintf('nu  = %g\n',x1(2));
fprintf('err = %g\n',sqrt(err1)/norm(U_exp));
% fprintf('exitflag = %g\n',exitflag1);
% disp(output1);
toc(t1)

t2 = tic;
[x2,err2,exitflag2,output2] = fminsearch(funoptim,x0,optionsfminsearch);
fprintf('\nfminsearch');
fprintf('\n----------\n');
fprintf('E   = %g GPa\n',x2(1));
fprintf('nu  = %g\n',x2(2));
fprintf('err = %g\n',sqrt(err2)/norm(U_exp));
% fprintf('exitflag = %g\n',exitflag2);
% disp(output2);
toc(t2)

t3 = tic;
[x3,err3,exitflag3,output3] = fminunc(funoptim,x0,optionsfminunc);
fprintf('\nfminunc');
fprintf('\n-------\n');
fprintf('E   = %g GPa\n',x3(1));
fprintf('nu  = %g\n',x3(2));
fprintf('err = %g\n',sqrt(err3)/norm(U_exp));
% fprintf('exitflag = %g\n',exitflag3);
% disp(output3);
toc(t3)

t4 = tic;
[x4,err4,exitflag4,output4] = fmincon(funoptim,x0,[],[],[],[],lb,ub,[],optionsfmincon);
fprintf('\nfmincon');
fprintf('\n-------\n');
fprintf('E   = %g GPa\n',x4(1));
fprintf('nu  = %g\n',x4(2));
fprintf('err = %g\n',sqrt(err4)/norm(U_exp));
% fprintf('exitflag = %g\n',exitflag4);
% disp(output4);
toc(t4)
