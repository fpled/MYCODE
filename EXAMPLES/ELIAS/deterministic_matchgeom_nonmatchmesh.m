clear all
close all
D = DOMAIN(2,[0,0],[5,5]);
P1=[2,2];
P2=[4,4];

patchbox = DOMAIN(2,P1,P2);

ng = 30;
rf = .05;
% domaine grossier Sg
Sg = mesh(D,ng,ng);
Sg = convertelem(Sg,'TRI3');
Sg = createddlnode(sortnodenumber(Sg));
Sg = addcl(Sg,[],'u');

% patch fin Sfin
Sfin = gmsh(patchbox,rf);
Sfin = createddlnode(sortnodenumber(Sfin));
Sdfin = addcl(Sfin,[],'u');

% partition des elements de Sg
Sg = setparamgroupelem(Sg,'partition',0);
[temp,node1,numelem1] = intersect(Sg,patchbox);
[Sg,newgroup] = separateelemwithnum(Sg,numelem1);
Sg = setparamgroupelem(Sg,'partition',1,newgroup);

% domaines grossiers Sgout et Sgin
numelem = getnumgroupelemwithparam(Sg,'partition',1);
Sgin = removenodewithoutelem(keepgroupelem(Sg,numelem));
Sgin = createddlnode(sortnodenumber(Sgin));

numelem = getnumgroupelemwithparam(Sg,'partition',0);
Sgout = removenodewithoutelem(keepgroupelem(Sg,numelem));
Sgout = createddlnode(sortnodenumber(Sgout));
for kk=1:4
    Sgout = addcl(Sgout,getedge(D,kk),'u');
end

% surface fine et surface grossi�re
Gfin = createddlnode(sortnodenumber(create_boundary(Sfin)));
Ggout = createddlnode(sortnodenumber(setdiff(create_boundary(Sgout),create_boundary(Sg))));
figure(100);plot(Ggout,'numnode')
Ggin = createddlnode(sortnodenumber(create_boundary(Sgin)));
figure(101);plot(Ggin,'numnode') % identique � Ggout
%Ggin=Ggout;
MGfin =  calc_matrix(BILINFORM(0,0),Gfin);
MGgout =  calc_matrix(BILINFORM(0,0),Ggout);

PSgout2Ggout = freevector(Sgout,calc_P_transfer(Sgout,Ggout),2);
PSgin2Ggin = calc_P_transfer(Sgin,Ggin);
PSfin2Gfin = calc_P_transfer(Sfin,Gfin);
%PGfin2Ggin = calc_P_transfer(Gfin,Ggin);

figure(1)
clf
plot(Sg)
plot(Sfin,'color','r')
plot(Sgin,'color','c')
plot(Sgout,'color','b')

%PSgout2Ggin = freevector(Sgout,calc_P_transfer(Sgout,Ggin),2);%=PSgout2Ggout
%PSgin2Ggout = calc_P_transfer(Sgin,Ggout);%=PSgin2Ggin



%% projecteur interface grossiere et interface fine
%%%% cas d'un maillage compatible
%PGgin2Ggout = calc_P_transfer(Ggin,Ggout);
%PGgout2Gfin = PGfin2Ggout';

%%%% geometrie compatible mais maillage non compatible
xGfin = getcoord(getnode(Gfin));xGgout = getcoord(getnode(Ggout));

% pour aller du grossier vers le fin
PGgout2Gfin = eval_sol(Ggout,speye(getnbddl(Ggout),getnbddl(Ggout)),POINT(xGfin),'u');
PGgout2Gfin = sparse(squeeze(PGgout2Gfin)');
% pour aller du fin vers le grossier
PGfin2Ggout = (PGgout2Gfin'*MGfin*PGgout2Gfin)\(PGgout2Gfin'*MGfin)
%ou
%PGfin2Ggout = eval_sol(Gfin,speye(getnbddl(Gfin),getnbddl(Gfin)),POINT(xGgout),'u');
%PGfin2Ggout = sparse(squeeze(PGfin2Ggoutbis)');

couplage = 2;

if couplage ==1 % couplage fin
    nG = size(MGfin,1);
else % couplage grossier
    nG = size(MGgout,1);
end


%%
A = calc_matrix(BILINFORM(1,1),Sg);
b = calc_vector(LINFORM(0),Sg);
u = A\b;
figure(5)
plot(u,Sg);
title('REFERENCE')
%% sans multiplicateur

Agout = calc_matrix(BILINFORM(1,1),Sgout);
%Agin  = calc_matrix(BILINFORM(1,1),Sgin);
Afin  = calc_matrix(BILINFORM(1,1),Sfin);
bgout = calc_vector(LINFORM(0),Sgout);
%bgin = calc_vector(LINFORM(0),Sgin);
bfin = calc_vector(LINFORM(0),Sfin);

PSgout2Sfin = PSfin2Gfin'*PGgout2Gfin*PSgout2Ggout;
%PSfin2Sgout = PSgout2Ggout'*PGfin2Ggout*PSfin2Gfin;

A11 = Agout+PSgout2Sfin'*Afin*PSgout2Sfin;
A12 = PSgout2Sfin'*freevector(Sdfin,Afin,1)';
A22 = freematrix(Sdfin,Afin);
A21 = freevector(Sdfin,Afin,1)*PSgout2Sfin;
b1 =  bgout + PSgout2Sfin'*bfin;
b2 =  freevector(Sdfin,bfin);

A = [A11,A12; A21,A22];
b = [b1;b2];

Uzref = solve(A,b);
Uref = Uzref(1:size(Agout,1));
zref = Uzref(size(Agout,1)+1:end);
wref = PSgout2Sfin*Uref + unfreevector(Sdfin,zref);

figure(15)
clf
subplot(1,2,1)
plot_sol(Sgout,Uref,'edgecolor','none');%caxis(cax);
ax=axis;
title('U')
colorbar
subplot(1,2,2)
plot_sol(Sfin,wref,'edgecolor','none');cax=caxis;
axis(ax)
title('w')
colorbar
figure(16)
clf
plot_sol(Sgout,Uref,'edgecolor','none');%caxis(cax);
plot_sol(Sfin,wref,'edgecolor','none');
title('U+w')
colorbar

%% avec multiplicateur

Agout = calc_matrix(BILINFORM(1,1),Sgout);
%Agin  = calc_matrix(BILINFORM(1,1),Sgin);
Afin  = calc_matrix(BILINFORM(1,1),Sfin);
bgout = calc_vector(LINFORM(0),Sgout);
%bgin = calc_vector(LINFORM(0),Sgin);
bfin = calc_vector(LINFORM(0),Sfin);

ngout = size(Agout,1);nfin = size(Afin,1);
%PSfin2Sgout = PSgout2Ggout'*PGfin2Ggout*PSfin2Gfin;
A11 = Agout;
A22 = Afin;
b1 =  bgout;
b2 =  bfin;
if couplage==1
    A13 = PSgout2Ggout'*PGgout2Gfin'*MGfin;
    A23 = -PSfin2Gfin'*MGfin;
    A31 = MGfin*PGgout2Gfin*PSgout2Ggout;
    A32 = -MGfin*PSfin2Gfin;
else
    A13 = PSgout2Ggout'*MGgout;
    A23 = -PSfin2Gfin'*PGfin2Ggout'*MGgout;
    A31 = A13';%A31 = MGgout*PSgout2Ggout;
    A32 = A23';%A32 = -MGgout*PGfin2Ggout*PSfin2Gfin;
end
A = [A11,sparse(ngout,nfin),A13; sparse(nfin,ngout),A22,A23;A31,A32,sparse(nG,nG)];
b = [b1;b2;zeros(nG,1)];

Uzref = solve(A,b);
Uref = Uzref(1:ngout);
wref = Uzref(ngout+1:ngout+nfin);

figure(17)
clf
subplot(1,2,1)
plot_sol(Sgout,Uref,'edgecolor','none');%caxis(cax);
ax=axis;
title('U')
subplot(1,2,2)
plot_sol(Sfin,wref,'edgecolor','none');cax=caxis;
axis(ax)
title('w')
figure(18)
clf
plot_sol(Sgout,Uref,'edgecolor','none');%caxis(cax);
plot_sol(Sfin,wref,'edgecolor','none');
title('U+w')

