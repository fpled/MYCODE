clear all
close all
D = DOMAIN(2,[0,0],[5,5]);
P1=[0.,3.1];
P2=[2.1,5];
P1=[1.,2.1];
P2=[3.1,4];

patchbox = DOMAIN(2,P1,P2);

ng = 30;
rg = 0.5;
rf = .05;
compatgeom = 1;
% domaine grossier Sg
if compatgeom
    Sg = mesh(D,ng,ng);
    Sg = convertelem(Sg,'TRI3');
else
    Sg = gmsh(D,rg);
end
Sg = createddlnode(sortnodenumber(Sg));

eD=1:3; % edges with homogenous Dirichlet BD
eN=4; % % edges with homogenous Neumann BD
for kk=eD
    Sg = addcl(Sg,getedge(D,kk),'u');
end
%Sg = addcl(Sg,[],'u');

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
%%
Sgin = Sg;
%Sgin = clearfaces(Sg);
Sgin = keepeleminnode(removenodewithoutelem((keepgroupelem(Sgin,numelem))));
Sgin = createddlnode(sortnodenumber(Sgin));

numelem = getnumgroupelemwithparam(Sg,'partition',0);
Sgout = Sg;
%Sgout = clearfaces(Sg);
Sgout = keepeleminnode(removenodewithoutelem((keepgroupelem(Sgout,numelem))));
Sgout = createddlnode(sortnodenumber(Sgout));
for kk=1:3
    Sgout = addcl(Sgout,getedge(D,kk),'u');
end

% surface fine
Gfin = createddlnode(sortnodenumber(create_boundary(Sfin)));
Gfin=create_boundary(Sfin);
for kk=1:4
    Gfin=(setdiff(Gfin,getedge(D,kk)));
end
Gfin = createddlnode(sortnodenumber(Gfin));

MGfin =  calc_matrix(BILINFORM(0,0),Gfin);

PSfin2Gfin = calc_P_transfer(Sfin,Gfin);

figure(1)
clf
plot(Sg)
plot(Sfin,'color','r')
plot(Sgin,'color','c')
plot(Sgout,'color','b')


%% projecteur domaine grossier et interface fine
versionM=1;
%%%% geometrie et maillage non compatibles
xGfin = getcoord(getnode(Gfin));

lightcoupling = 1;
if lightcoupling==1
    ls=LSRECTANGLE(P1(1),P1(2),P2(1),P2(2));
    Sgoutsplit = LSMODEL(Sgout,ls);
    Sgoutsplit = changeelemnumber(lssplitelem(Sgoutsplit));
    Sgoutsplit = calc_connec(Sgoutsplit);
    groupcut = find(cellfun(@(C) strcmp(getlstype(C),'cut') || strcmp(getlstype(C),'in'),Sgoutsplit.groupelem));
    elemcut = getnumelem(Sgoutsplit,groupcut);
    nodecoupling=find(sum(Sgoutsplit.connec.elem2node(:,elemcut),2));
    if versionM==1
        epb=getedges(patchbox);
        for i=1:size(epb,2)
            nodecoupling=union(nodecoupling,getnumnodeelem(intersect(create_boundary(Sgoutsplit),epb{i})));
        end
    end
% pour avoir une couche de noeud en plus, il faut faire la ligne suivante
    nodecoupling = find(sum(Sgoutsplit.connec.node2node(:,nodecoupling),2));

else
    nodecoupling = 1:getnbddl(Sgout);
end
% pour aller du grossier vers le fin
I = speye(getnbddl(Sgout),getnbddl(Sgout));
PSgout2Gfin = sparse(getnbddl(Gfin),getnbddl(Sgout));
temp = eval_sol(Sgout,I(:,nodecoupling),POINT(xGfin),'u');
temp = sparse(squeeze(temp)');
PSgout2Gfin(:,nodecoupling)= temp;
PSgout2Gfin= freevector(Sgout,PSgout2Gfin,2);

figure(32)
clf
plot(Sgout)
plot(Sgout.node(nodecoupling),'*')
plot(Sfin,'color','r')



couplage = 1; % avec compatgeom=0, couplage=1

if couplage ==1 %couplage fin
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
%% A(U+w,dU + dw) = L(dU+dw)

Agout = calc_matrix(BILINFORM(1,1),Sgout);
Afin  = calc_matrix(BILINFORM(1,1),Sfin);
bgout = calc_vector(LINFORM(0),Sgout);
bfin = calc_vector(LINFORM(0),Sfin);

PSgout2Sfin = PSfin2Gfin'*PSgout2Gfin;

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
Afin  = calc_matrix(BILINFORM(1,1),Sfin);
bgout = calc_vector(LINFORM(0),Sgout);
bfin = calc_vector(LINFORM(0),Sfin);

ngout = size(Agout,1);nfin = size(Afin,1);
A11 = Agout;
A22 = Afin;
b1 =  bgout;
b2 =  bfin;
if couplage==1
    A13 = PSgout2Gfin'*MGfin;
    A23 = -PSfin2Gfin'*MGfin;
    A31 = MGfin*PSgout2Gfin;
    A32 = -MGfin*PSfin2Gfin;
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

