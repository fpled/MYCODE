%clear all; close all;
Dom=DOMAIN(2,[0,0],[1,1]);
[Sini,Xini,Yini] = mesh(Dom,50,50);
Sini = convertelem(Sini,'TRI3');
%%
Snew = Sini;
myboxes = {DOMAIN(2,[0.5,0.5],[0.9,0.9])};

Snew = setparamgroupelem(Snew,'partition',0);
[Stemp,node1,numelem1] = intersect(Snew,myboxes{1});
[Snew,newgroup] = separateelemwithnum(Snew,numelem1);
Snew = setparamgroupelem(Snew,'partition',1,newgroup);
figure(3)
clf
plotparamelem(Snew,'partition')

cltype = 2;
%%
S = createddlnode(Sini,DDL('u'));
switch cltype
    case 1
        S = addcl(S,POINT([0,0]),'u',0);
    case 2
        for k=1:4
            S = addcl(S,getedge(Dom,k),'u',0);
        end
end

Spart=cell(1,length(myboxes)+1);
Spart{1} = getnumgroupelemwithparam(Snew,'partition',0);
for k=1:length(myboxes)
    Spart{k+1} = getnumgroupelemwithparam(Snew,'partition',k);
end

for k=1:(length(myboxes)+1)
    Spart{k} = keepgroupelem(Snew,Spart{k});
    Spart{k} = removenodewithoutelem(Spart{k});
end

% surfaces fines et grossieres
for k = 2:length(myboxes)+1
    surface{k} = createddlnode(sortnodenumber(create_boundary(Spart{k})));
    xsurface{k} = getcoord(getnode(surface{k}));
    Msurface{k} =  calc_matrix(BILINFORM(0,0),surface{k});
end
surface{1} = createddlnode(sortnodenumber(setdiff(create_boundary(Spart{1}),create_boundary(S))));
%
xsurface{1} = getcoord(getnode(surface{1}));
Msurface{1} =  calc_matrix(BILINFORM(0,0),surface{1});

%%
a = BILINFORM(1,1);
S0 = createddlnode(Spart{1},DDL('u'));
switch cltype
    case 1
        S0 = addcl(S0,POINT([0,0]),'u',0);
    case 2
        for k=1:4
            S0 = addcl(S0,getedge(Dom,k),'u',0);
        end
end

A0 = calc_matrix(a,S0);
P0 = calc_P(S,S0);
P0 = freevector(S,P0,2);
P0 = freevector(S0,P0,1);

S1 = createddlnode(Spart{2},DDL('u'));
S1 = final(S1);
x1 = getcoord(getnode(S1));
y1 = x1(:,2);
x1 = x1(:,1);

switch cltype
    case 1
        %        kappa = 1+...%(x1>=.4).*(x1<=.6).*(y1>=.4).*(y1<=.6).*...
        %    90*exp(-10*(x1-.7+y1-.7).^2/.1^2-(y1-x1-.7+.7).^2/.1^2);%.*cos(x1*5);
        kappa = 1+...%(x1>=.4).*(x1<=.6).*(y1>=.4).*(y1<=.6).*...
            10*exp(-5*(x1-.7+y1-.7).^2/.1^2-5*(y1-x1-.7+.7).^2/.1^2);%.*cos(x1*5);
    case 2
        kappa = 1+...%(x1>=.4).*(x1<=.6).*(y1>=.4).*(y1<=.6).*...
            100*exp(-3*(x1-.7).^2/.1^2-3*(y1-.7).^2/.1^2);%.*cos(x1*5);
%         kappa = 1.;
end

figure(66)
clf
plot(FENODEFIELD(kappa),S1)%,'surface')
a1 = BILINFORM(1,1,kappa,0);
A1 = calc_matrix(a1,S1);
a1tilde = BILINFORM(1,1,1,0);
A1tilde = calc_matrix(a1tilde,S1);
P1 = calc_P(S,S1);
P1 = freevector(S,P1,2);
P1 = freevector(S1,P1,1);

Sglobal = Sini;
switch cltype
    case 1
        f = LINFORM(0,0);
    case 2
        f = LINFORM(0,1);
end
b0 = f{S0}(:);
b1 = f{S1}(:);
AG0 = P0'*A0*P0;
AG1tilde = P1'*A1tilde*P1;
AG1 = P1'*A1*P1;
AGtilde = AG0+AG1tilde;
AG = AG0+AG1;
fG0 = P0'*b0;
fG1 = P1'*b1;
switch cltype
    case 1
        F1 = getedge(Dom,2);
        F1 = intersect(create_boundary(Sini),F1);
        F1 = createddlnode(F1,DDL('u'));
        PF1 = freevector(S,calc_P(S,F1),2);
        f1 = LINFORM(0,1);
        fG0 = fG0+ PF1'*f1{F1}(:);
        F2 = getedge(Dom,4);
        F2 = intersect(create_boundary(Sini),F2);
        F2 = createddlnode(F2,DDL('u'));
        PF2 = freevector(S,calc_P(S,F2),2);
        f2 = LINFORM(0,-1);
        fG0 = fG0+ PF2'*f2{F2}(:);
end
fG = fG0+fG1;

U = solvesingular(AG,fG);
%U = unfreevector(S,U);
figure
%plot_sol(S,U,'epsilon',1)
plot(FENODEFIELD(U),S)


%%
Gamma = create_boundary(S1);
Gamma = createddlnode(Gamma,DDL('u'));
P1Gamma = calc_P(S1,Gamma);
P0Gamma = freevector(S0,calc_P(S0,Gamma),2);
PSGamma = freevector(S,calc_P(S,Gamma),2);

overlap = 1; % si 1 : prolongement du pb global sous le patch
%

Un1 = zeros(getnbddlfree(S),1);
Un = zeros(getnbddlfree(S),1);
wn1 = zeros(getnbddlfree(S1),1);
ln1 = zeros(getnbddlfree(Gamma),1);
wn = zeros(getnbddlfree(S1),1);
ln = zeros(getnbddlfree(Gamma),1);

clear errUn
rho=1; % param de relaxation

delta = 1e-3;
deltaw = delta;
approx = 0;
errUnstagn=[];
for k=1:50
    if overlap
        Un1 = AGtilde\(fG0-PSGamma'*Msurface{1}*ln+AG1tilde*Un);
    else
        Un1 = A0\(P0*fG0-P0Gamma'*Msurface{1}*ln);
        Un1 = P0'*Un1;
    end
    Un1 = rho*Un1+(1-rho)*Un;
    errUnstagn(k) = norm(Un-Un1)/norm(Un1);
    if approx
        Uneps = reshape(unfreevector(S,Un1),size(Xini));
        Uneps = svdtruncate(Uneps,delta);
        Uneps = Uneps(:);
        Un1 = freevector(S,Uneps);
    end
    
    Unwn = P0'*P0*Un + P1'*(wn-P1Gamma'*P1Gamma*wn);
    
    errUn(k) = norm(Unwn-U)/norm(U);
    Un = Un1;
    T = [A1,-P1Gamma'*Msurface{1};-Msurface{1}*P1Gamma,sparse(length(ln),length(ln))]\...
        [b1;-Msurface{1}*PSGamma*Un];
    wn = full(T(1:length(wn)));
    ln = full(T(length(wn)+1:end));
    
    if approx
        wn = reshape(wn,sqrt(length(wn)),sqrt(length(wn)));
        wn = svdtruncate(wn,deltaw);
        wn = wn(:);
        ln = P1Gamma*(A1*wn-b1);
    end
%     if errUn(k)<1e-13
%         break
%     end
end
figure(45)
if overlap
    semilogy(errUn,'b')
else
    semilogy(errUn,'r')
end
hold on
%semilogy(errUnstagn,'k--')

%% test spectral radius
% methode avec overlap
PSI = [sparse(length(ln),length(wn)),speye(length(ln),length(ln))]*...
    ([A1,-P1Gamma';-P1Gamma,sparse(length(ln),length(ln))]\...
    [sparse(length(wn),length(Un));-PSGamma]);
opt.disp=0;
Aoverlap = AGtilde\(PSGamma'*PSI + AG0);
eigs(Aoverlap,10,'LM',opt)
% methode sans overlap
Anonoverlap = A0\(P0Gamma'*PSI*P0' + A0);
eigs(Anonoverlap,10,'LM',opt)

%%
figure(56)
clf
subplot(1,3,1)

plot(FENODEFIELD(Unwn),S)%,'surface')
colorbar
subplot(1,3,2)
plot(FENODEFIELD(Un),S)%,'surface')
colorbar
subplot(1,3,3)
plot(FENODEFIELD(wn),S1)%,'surface')
colorbar
%%
figure(57)
clf
subplot(1,3,1)
plot_sol(S,(Unwn),'epsilon',1)
colorbar
ca=caxis;
subplot(1,3,2)
plot_sol(S,(Un),'epsilon',1)
colorbar
caxis(ca)
subplot(1,3,3)
plot_sol(S1,(wn),'epsilon',1)
colorbar
caxis(ca)

%%
% svdfin=svds(reshape(unfreevector(S,Unwn),size(Xini)),30)
% svdgros=svds(reshape(unfreevector(S,Un),size(Xini)),30)
% svdfinfin=svds(reshape(unfreevector(S1,wn),sqrt(length(wn)),sqrt(length(wn))),30)
% figure(676)
% semilogy(svdgros,'b')
% hold on
% semilogy(svdfin,'r')
% semilogy(svdfinfin,'g')
