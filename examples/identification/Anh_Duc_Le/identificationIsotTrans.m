%% Identification of transversely isotropic elastic properties %%
%%-------------------------------------------------------------%%

% clc
clear all
close all
% set(0,'DefaultFigureVisible','off');
% rng('default');

%% Experimental values
EL_exp = 11.55;
ET_exp = 0.5;
nuL_exp = 0.4;
GL_exp = 0.55;
x_exp = [EL_exp ET_exp nuL_exp GL_exp];
U_exp = ThreePointsBendingIsotTrans(x_exp);
fprintf('Exp values\n');
fprintf('----------\n');
fprintf('EL_exp  = %g GPa\n',x_exp(1));
fprintf('ET_exp  = %g\n',x_exp(2));
fprintf('nuL_exp = %g\n',x_exp(3));
fprintf('GL_exp  = %g GPa\n',x_exp(4));

%% Identification
EL0 = 10.0;
ET0 = 0.2;
nuL0 = 0.2;
GL0 = 0.2;
x0 = [EL0 ET0 nuL0 GL0];
% lb = [0 0 -1 0];
lb = [0 0 0 0];
ub = [Inf Inf 0.5 Inf];

tolX = 1e-14;
tolFun = 1e-14;
display = 'off';
optionslsqnonlin  = optimoptions('lsqnonlin','Display',display,'TolX',tolX,'TolFun',tolFun);
optionsfminsearch = optimset('Display',display,'TolX',tolX,'TolFun',tolFun);
optionsfminunc    = optimoptions('fminunc','Display',display,'TolX',tolX,'TolFun',tolFun,'Algorithm','quasi-newton');
optionsfmincon    = optimoptions('fmincon','Display',display,'TolX',tolX,'TolFun',tolFun);

ampl = 1e-2;
noise = ampl.*(2*rand(length(U_exp),1)-1).*U_exp;
% noise = 0;

funlsqnonlin = @(x) funlsqnonlinIsotTrans(x,U_exp,noise);
funoptim = @(x) funoptimIsotTrans(x,U_exp,noise);

t1 = tic;
[x1,err1,~,exitflag1,output1] = lsqnonlin(@(x) funlsqnonlin(x),x0,lb,ub,optionslsqnonlin);
fprintf('\nlsqnonlin');
fprintf('\n---------\n');
fprintf('EL  = %g GPa\n',x1(1));
fprintf('ET  = %g GPa\n',x1(2));
fprintf('nuL = %g\n',x1(3));
fprintf('GL  = %g GPa\n',x1(4));
fprintf('err = %g\n',sqrt(err1)/norm(U_exp));
% fprintf('exitflag = %g\n',exitflag1);
% disp(output1);
toc(t1)

t2 = tic;
[x2,err2,exitflag2,output2] = fminsearch(@(x) funoptim(x),x0,optionsfminsearch);
fprintf('\nfminsearch');
fprintf('\n----------\n');
fprintf('EL  = %g GPa\n',x2(1));
fprintf('ET  = %g GPa\n',x2(2));
fprintf('nuL = %g\n',x2(3));
fprintf('GL  = %g GPa\n',x2(4));
fprintf('err = %g\n',sqrt(err2)/norm(U_exp));
% fprintf('exitflag = %g\n',exitflag2);
% disp(output2);
toc(t2)

t3 = tic;
[x3,err3,exitflag3,output3] = fminunc(@(x) funoptim(x),x0,optionsfminunc);
fprintf('\nfminunc');
fprintf('\n-------\n');
fprintf('EL  = %g GPa\n',x3(1));
fprintf('ET  = %g GPa\n',x3(2));
fprintf('nuL = %g\n',x3(3));
fprintf('GL  = %g GPa\n',x3(4));
fprintf('err = %g\n',sqrt(err3)/norm(U_exp));
% fprintf('exitflag = %g\n',exitflag3);
% disp(output3);
toc(t3)

t4 = tic;
[x4,err4,exitflag4,output4] = fmincon(@(x) funoptim(x),x0,[],[],[],[],lb,ub,[],optionsfmincon);
fprintf('\nfmincon');
fprintf('\n-------\n');
fprintf('EL  = %g GPa\n',x4(1));
fprintf('ET  = %g GPa\n',x4(2));
fprintf('nuL = %g\n',x4(3));
fprintf('GL  = %g GPa\n',x4(4));
fprintf('err = %g\n',sqrt(err4)/norm(U_exp));
% fprintf('exitflag = %g\n',exitflag4);
% disp(output4);
toc(t4)
