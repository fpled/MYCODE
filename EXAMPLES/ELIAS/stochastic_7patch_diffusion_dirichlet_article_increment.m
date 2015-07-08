tolsvd = 1e-15;
tolmultisvd = 1e-3;
tolcgs = 1e-10;
tol = 1e-15;
toliter = 1e-15;
nbiter = 200;
cgsiter = 1000;
tolPGD = 1e-1;
PGDrank = 300;
%%
D = DOMAIN(2,[0,0],[10,10]);
% D = DOMAIN(2,[0,0],[5,5]);
myboxes = { DOMAIN(2,[1,1],[3,3])  DOMAIN(2,[4,1],[6,3]) DOMAIN(2,[7,1],[9,3])...
    DOMAIN(2,[1,5],[3,7])  DOMAIN(2,[4,7],[6,9]) DOMAIN(2,[7,5],[9,7])...
    DOMAIN(2,[0,8],[2,10])...
    }
% myboxes = { DOMAIN(2,[1,1],[3,3])   DOMAIN(2,[7,1],[9,3])...
%             DOMAIN(2,[4,5],[6,7])...
%             DOMAIN(2,[0,8],[2,10])...
%              }

% myboxes = {DOMAIN(2,[2,2],[4,4])}
ng = 20;
rf = .05;
% domaine grossier Sg
S = mesh(D,ng,ng);
S = convertelem(S,'TRI3');
S = createddlnode(sortnodenumber(S));
S = addcl(S,[],'u');
S = setparamgroupelem(S,'partition',0,1:getnbgroupelem(S));
% patch fin Sfin

for k = 2:length(myboxes)+1
    % Spartf{k} = gmsh(myboxes{k-1},rf);
    % Spartf{k} = gmsh2femobject(2,'C:\PROGS\gmsh-2.5.0-Windows\domfin.msh');
    Spartf{k} = mesh(myboxes{k-1},ng+9,ng+9);
    Spartf{k} = convertelem(Spartf{k},'TRI3');
    Spartf{k} = createddlnode(sortnodenumber(Spartf{k}));
    Spartfd{k} = addcl(Spartf{k},[],'u');
end
% partition des elements de Sg
for k = 1:length(myboxes)
    [temp,node1,numelem1] = intersect(S,myboxes{k});
    [S,newgroup] = separateelemwithnum(S,numelem1);
    S = setparamgroupelem(S,'partition',k,newgroup);
    
    % domaines grossiers Sgout et Sgin
    numelem = getnumgroupelemwithparam(S,'partition',k);
    Spart{k+1} = removenodewithoutelem(keepgroupelem(S,numelem));
    Spart{k+1} = createddlnode(sortnodenumber(Spart{k+1}));
end
numelem = getnumgroupelemwithparam(S,'partition',0);
Spart{1} = removenodewithoutelem(keepgroupelem(S,numelem));
Spart{1} = createddlnode(sortnodenumber(Spart{1}));

for kk=1:4
    Spart{1} = addcl(Spart{1},getedge(D,kk),'u');
end
%
% surfaces fines et grossieres
for k = 2:length(myboxes)+1
    surfacef{k} = createddlnode(sortnodenumber(create_boundary(Spartf{k})));
    xsurfacef{k} = getcoord(getnode(surfacef{k}));
    surface{k} = createddlnode(sortnodenumber(create_boundary(Spart{k})));
    xsurface{k} = getcoord(getnode(surface{k}));
    Msurfacef{k} =  calc_matrix(BILINFORM(0,0),surfacef{k});
end
surface{1} = createddlnode(sortnodenumber(setdiff(create_boundary(Spart{1}),create_boundary(S))));
%
xsurface{1} = getcoord(getnode(surface{1}));
Msurface{1} =  calc_matrix(BILINFORM(0,0),surface{1});

for k =2:length(myboxes)+1
    coord = getcoord(getnode(Spartf{k}));
    xi{k} = coord(:,1);
    yi{k} = coord(:,2);
    nbelem{k} = getnbelem(Spartf{k})
    helempatchf{k} = sqrt(2*getvolume(myboxes{k-1})/getnbelem(Spartf{k}))
end
nbelemS = getnbelem(S)
helemS = sqrt(2*2/nbelemS)
for k =2:length(myboxes)+1
    coordt = getcoord(getnode(Spart{k}));
    xit{k} = coordt(:,1);
    yit{k} = coordt(:,2);
    helempatch{k} = sqrt(2*getvolume(myboxes{k-1})/getnbelem(Spartf{k}))
end

numpatchdiff = [1 3 5];
numpatchdir = [2 4 6 7];
% numpatchdir=1;
%%
figure(1)
clf
plot(Spart{1})
hold on
for k=2:length(myboxes)+1
    plot(Spartf{k},'color','r');
    hold on
end
%% Projecteurs
P{1} = calc_P_transfer(S,Spart{1});
P{1} = freevector(S,freevector(Spart{1},P{1},1),2);
for k = 2:(length(myboxes)+1)
    P{k} = calc_P_transfer(S,Spart{k});
    P{k} = freevector(S,P{k},2);
    Pgamma{k} = calc_P_transfer(S,surface{k});
    Pgamma{k} = freevector(S,Pgamma{k},2);
    Ppatchgamma{k} = calc_P_transfer(Spart{k},surface{k});
    Phorspatchgamma{k} = calc_P_transfer(Spart{1},surface{k});
    Phorspatchgamma{k} = freevector(Spart{1},Phorspatchgamma{k},2);
end
%% Projecteurs fins
for k = 2:(length(myboxes)+1)
    Ppatchfgammaf{k} = calc_P_transfer(Spartf{k},surfacef{k});
end
%% projecteur interface grossiere et interface fine
%%%% geometrie compatible mais maillage non compatible

% pour aller du fin vers le grossier
for k=2:length(myboxes)+1
    Pgammaoutgammaf{k} = eval_sol(surface{k},speye(getnbddl(surface{k}),getnbddl(surface{k})),POINT(xsurfacef{k}),'u');
    Pgammaoutgammaf{k} = sparse(squeeze(Pgammaoutgammaf{k})');
    % pour aller du grossier vers le fin
    Pgammafgammaout{k} = (Pgammaoutgammaf{k}'*Msurfacef{k}*Pgammaoutgammaf{k})\(Pgammaoutgammaf{k}'*Msurfacef{k});
    % Pgammafgammaout{k} = eval_sol(surfacef{k},speye(getnbddl(surfacef{k}),getnbddl(surfacef{k})),POINT(xsurface{k}),'u');
    % Pgammafgammaout{k} = sparse(squeeze(Pgammafgammaout{k})');
end

%%  Donnees stochastiques

rrr = 8;
ppp = 2;
nbva = 15;
lg = POLYFE(linspace(0,1,rrr+1),ppp);
% PC = POLYCHAOSTP(RANDPOLYS(lg,nbva),ppp,'typebase',1,'groups',{1});
% [X,PC] = PCTPMODEL(RANDVARS(RVUNIFORM(0,1),nbva),'randpolys',RANDPOLYS(lg,nbva),'order',ppp,'groups',{1:5,6:7,8:12,13:14,15,16:17,18});
% [X,PC] = PCTPMODEL(RANDVARS(RVUNIFORM(0,1),nbva),'randpolys',RANDPOLYS(lg,nbva),'order',ppp,'groups',{1:5,6:7,8:12,13:14,15,16:17,18});
% [X,PC] = PCTPMODEL(RANDVARS(RVUNIFORM(0,1),nbva),'randpolys',RANDPOLYS(lg,nbva),'order',ppp,'groups',{1:5,6,7:11,12,13,14,15});
% fedim = 1:15;%[6,12,14,15]
% femesh = ones(1,nbva);femesh([6,12,14,15])=6;
[X,PC] = PCTPMODEL(RANDVARS(RVUNIFORM(0,1),nbva),'order',[3*ones(1,5),2,3*ones(1,5),2,3,2,2],...
    'typebase',1,'groups',{1:5,6,7:11,12,13,14,15},'fedim',[6,12,14,15],'femesh',{rrr,rrr,rrr,rrr},'pcg');
PC = calc_massegroups(PC);


%%
for rr=1:getnbgroups(PC)
    RANDVARS(getpcgroup(PC,rr))
end
%% calcul des masses
mas = PCTPMATRIX(PC,1);
for i=1 : getnbgroups(PC)
    masi = getmasse(calc_masse(getpcgroup(PC,i)));
    masi = multisum(masi);
    mas{i} = masi;
end
%% Definition des coef diff
beta=0;
for i=1:length(numpatchdiff)
    num = numpatchdiff(i)+1;
    xl1 = getvertex(myboxes{num-1},1);
    xl2 = getvertex(myboxes{num-1},3);
    
    if i<length(numpatchdiff)
        %
        c1 = (xl1(1)+xl2(1))/2;
        c2 = (xl1(2)+xl2(2))/2;
        %
        c3 = (xl1(1)+xl2(1))/2.1;
        c4 = c2;
        %
        c5 = c1;
        c6 = (xl1(2)+xl2(2))/2.4;
        %
        c7 = c3;
        c8 = c6;
        %
        c9 = (xl1(1)+xl2(1))/1.9;
        c10 = (xl1(2)+xl2(2))/2.6;
        
        lc1=0.25;lc2=0.2;lc3=0.2;A=2;
        f1{num} = calc_coefdif_exp(c1,c2,xi{num},yi{num},lc1,A);
        f2{num} = calc_coefdif_exp(c3,c4,xi{num},yi{num},lc1,A);
        f3{num} = calc_coefdif_exp(c5,c6,xi{num},yi{num},lc2,A);
        f4{num} = calc_coefdif_exp(c7,c8,xi{num},yi{num},lc2,A);
        f5{num} = calc_coefdif_exp(c9,c10,xi{num},yi{num},lc3,A);
        
        kappa{num} = @(x) 1 + 10*(f1{num}*(x(:,1)+1)/2+f2{num}*(x(:,2)+1)/2+f3{num}*(x(:,3)+1)/2+f4{num}*(x(:,4)+1)/2+f5{num}*(x(:,5)+1)/2);
        kappatilde{num} = 1 + beta*10*(f1{num}+f2{num}+f3{num}+f4{num}+f5{num});
    elseif i == length(numpatchdiff)
        cend1 = (xl1(1)+xl2(1))/2;
        cend2 = (xl1(2)+xl2(2))/2;
        kappa{num} = @(x) 1 + 3*(abs(xi{num}-cend1)<0.5).*(abs(yi{num}-cend2)<0.5)*(x(:,1)+1)/2;
        kappatilde{num} =  1 + 3*(abs(xi{num}-cend1)<0.5).*(abs(yi{num}-cend2)<0.5)*beta;
    end
end
%% decomp sur le chaos
for ii=1:length(numpatchdiff)
    k = numpatchdiff(ii)+1
    kp = kappa{k};
    kappapc{k} = decompmatrixiter(getpcgroup(PC,k-1),[],[],@(x) kp(x));
end
%%
for ii=1:length(numpatchdiff)
    k = numpatchdiff(ii)+1
    kr = random(kappapc{k});
    % kdiag = eye(size(kr,1),size(kr,1));
%     for jj=1:size(kr,1)
%         kdiag(jj,jj)=kr(jj);
%     end
    figure(487+k)
    clf
    % surf(xi{k},yi{k},kdiag)
    plot(kr,Spartf{k},'surface');
    % plot3(xi{k},yi{k},kr)
end
%% decomposition du coefficient de diffusion
%
km{1} = 1;
for ii=1:length(numpatchdiff)
    k = numpatchdiff(ii)+1
    [km{k},result{k}] =spectral_decomposition(kappapc{k},'display',true);
    %
    ktp{k} = setphi(PCTPMATRIX(PC,getV(km{k},1)),double(getL(km{k},1))',k-1);
    for i=2:getm(km{k})
        ktp{k} = ktp{k} + setphi(PCTPMATRIX(PC,getV(km{k},i)),double(getL(km{k},i))',k-1);
    end
    ktpsep{k} = SEPVECTOR(ktp{k});
end
%% definition des LEVELSET

for i=1:length(numpatchdir)
    num = numpatchdir(i)+1;
    xl1 = getvertex(myboxes{num-1},1);
    xl2 = getvertex(myboxes{num-1},3);
    
    if i<3
        %
        c1 = (xl1(1)+xl2(1))/2;
        c2 =  (xl1(2)+xl2(2))/2;
        %
        % ls{num} = @(x) -LSCIRCLE(c1+0.2*(x(:,1)),c2,.3+0.3*x(:,2));
        ls{num} = @(x) -LSCIRCLE(c1,c2,.3+0.1*x(:,1));
        % ls{num} = @(x) -LSCIRCLE(3,3,.5+0.1*x(:,1));
    elseif i==3
        c1 = (xl1(1)+xl2(1))/2;
        c2 =  (xl1(2)+xl2(2))/2;
        ls{num} = @(x)  -LSCIRCLE(c1+0.1*x(:,1),c2,.3);
        
        %
    elseif i==length(numpatchdir)
        ls1=@(x) LSCIRCLE(x(:,1),10-x(:,1),x(:,1));
        ls2 = @(x) -LSHYPERPLAN(2,0,10-x(:,1),1,-1);
        ls{num} = @(x) union(ls1(0.5*x+1),ls2(0.5*x+1));
    end
end

%% decomp sur le chaos
for ii=1:length(numpatchdir)
    k = numpatchdir(ii)+1
    fun{k} = @(x) getvalue(lseval(ls{k}(x),Spartf{k}));
    phipc{k} = decompmatrixiter(getpcgroup(PC,k-1),[],[],@(x) fun{k}(x));
    figure(447)
    clf
    % plotparamelem(Sls,'lsnumber')
    plot(Spart{1})
    [lsr{k},r] = random(phipc{k});
    lsplot(lssplitelem(LSMODEL(Spartf{k},LEVELSET(lsr{k}))));
    contourplot(LEVELSET(lsr{k}),Spartf{k});
end

%% decompositon spectrale de la fonction caracteristique
for ii=1:length(numpatchdir)
    k = numpatchdir(ii)+1
    [psim{k},resultpsi{k}] =spectral_decomposition(-phipc{k},'display',true,'tol',6e-3,'reference',-phipc{k});
    psitp{k} =  setphi(PCTPMATRIX(PC,getV(psim{k},1)),double(getL(psim{k},1))',k-1);
    for i=2:getm(psim{k})
        psitp{k} = psitp{k} + setphi(PCTPMATRIX(PC,getV(psim{k},i)),double(getL(psim{k},i))',k-1);
    end
    psitpsep{k} = SEPVECTOR(psitp{k});
    
end
%% coef de diff fictif avec mesure de domaine
kappatilde{1}=1;
for ii=1:length(numpatchdir)
    num = numpatchdir(ii)+1
    kappatilde{num}=1;
end

% if length(myboxes)==1
%     kappatilde{2} = (getvolume(myboxes{1})-(pi*expect(r^2)))/getvolume(myboxes{1})
% end

%% Formulation du prob globale
AS=zeros(getnbddlfree(S));
Atilde = cell(1,length(myboxes)+1);

for k=1:length(myboxes)+1
    Atilde{k} = calc_matrix(BILINFORM(1,1,kappatilde{k},0),Spart{k});
    Atilde{k} = P{k}'*Atilde{k}*P{k};
    AS = AS + Atilde{k};
end

Atilde{1} = calc_matrix(BILINFORM(1,1,kappatilde{1},0),Spart{1});
%%
for k=1:length(myboxes)+1
    bltilde = LINFORM(0);
    btilde{k} = bltilde{Spart{k}}(:);
end

%% Formulation du probleme sur le patch.

Apatch{1} = Atilde{1}*calc_ximasse(one(PC));
Apatchsep{1} = SEPMATRIX(Apatch{1});
%%
for ii=1:length(numpatchdiff)
    k = numpatchdiff(ii)+1
    [Apatch{k},Apatchsep{k}] = calc_matrix_diffuison_sep(ktp{k},PC,Spartf{k});
end
%%
for ii=1:length(numpatchdir)
    k = numpatchdir(ii)+1
    [Apatch{k},Apatchsep{k}] = calc_matrix_dirichlet_sep(psitp{k},PC,Spartf{k});
    [A2patch{k},A2patchsep{k}] = calc_matrix_dirichlet_dual_sep(psitp{k},PC,Spartf{k});
end
%%
bpatch{1} = btilde{1}*one(PC);
bpatchsep{1} = SEPMATRIX([{btilde{1}} getphi(one(PC))]);
%%
for ii=1:length(numpatchdiff)
    k = numpatchdiff(ii)+1
    [bpatch{k},bpatchsep{k}] = calc_vector_diffuison_sep(PC,Spartf{k})
end
%%
for ii=1:length(numpatchdir)
    k = numpatchdir(ii)+1
    [bpatch{k},bpatchsep{k}] = calc_vector_dirichlet_sep(psitp{k},PC,Spartf{k},mas)
end
%%
Dl = cell(1,length(myboxes)+1);
Dlsep = cell(1,length(myboxes)+1);
dint = BILINFORM(0,0);

for k=2:length(myboxes)+1
    Dl{k} = dint{surface{k}}(:,:)*one(PC);
    Dl{k} = calc_ximasse(Dl{k});
    % Dl{k} = speye(getnbddlfree(surface{k}));
    Dlsep{k} = SEPMATRIX(Dl{k});
end


Dlf = cell(1,length(myboxes)+1);
Dlfsep = cell(1,length(myboxes)+1);
for k=2:length(myboxes)+1
    Dlf{k} = calc_matrix(dint,surfacef{k});
    Dlf{k} = Dlf{k}*one(PC);
    Dlf{k} = calc_ximasse(Dlf{k});
    Dlfsep{k} = SEPMATRIX(Dlf{k});
    % Dl{k} = speye(getnbddlfree(surface{k}));
end

% Dlpsi = cell(1,length(myboxes)+1);
% 
% for k=2:length(myboxes)+1
%     Dlpsi{k} = PCRADIALMATRIX([getnbddl(surfacef{k}),getnbddl(surfacef{k})],PC);
%     for i=1:getm(psim{k})
%         dint2 = MULTILINFORM([0,0,0]);
%         Dlpsi{k} = Dlpsi{k} + dint2{surfacef{k}}(:,:,Ppatchfgammaf{k}*getV(psim{k},i))*getL(psim{k},i);
%     end
%     Dlpsi{k} = calc_ximasse(Dlpsi{k});
% end

% % Solution de reference
% couplage=1;
% multiplicateur=0;
% 
% tempref1 = sparse(getnbddlfree(Spart{1}),getnbddlfree(Spart{1}));
% tempref2 = sparse(getnbddlfree(Spart{1}),getnbddlfree(Spartfd{2}));
% tempref3 = sparse(getnbddlfree(Spartfd{2}),getnbddlfree(Spart{1}));
% for k=2:length(myboxes)+1
%     tempref1 = tempref1 + Phorspatchgamma{2}'*Pgammaoutgammaf{2}'*Ppatchfgammaf{k}*calc_matrix(BILINFORM(1,1),Spartf{k})*Ppatchfgammaf{k}'*Pgammaoutgammaf{2}*Phorspatchgamma{2};
%     tempref2 = tempref2 + Phorspatchgamma{2}'*Pgammaoutgammaf{2}'*Ppatchfgammaf{k}*freevector(Spartfd{k},A2patch{k},1)';
%     tempref3 = tempref3 + freevector(Spartfd{k},A2patch{k},1)*Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Phorspatchgamma{k};
% end
% 
% Aref = [Atilde{1}+tempref1,tempref2; tempref3 ,freematrix(Spartfd{k},Apatch{k})];
% Bref = [btilde{1}+ Phorspatchgamma{2}'*Pgammaoutgammaf{2}'*Ppatchfgammaf{2}*calc_vector(LINFORM(0),Spartf{k});freevector(Spartfd{2},bpatch{2})];
% 
% Uzref = solve(Aref,Bref);
% Uref = Uzref(1:size(Atilde{1},1));
% 
% zref = Uzref(size(Atilde{1},1)+1:end);
% %wref = Ppatchgammaf{k}'*Pfgsurface{k}*Phorspatchgamma{k}*Uref + unfreevector(Sdf{k},zref);
% wref = Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Phorspatchgamma{k}*Uref + psim{2}.*unfreevector(Spartfd{k},zref);
% 
% [Ur,r] = random(Uref);
% wr = randomeval(wref,r);
% 
% 
% figure(16)
% clf
% 
% plot_sol(Spart{1},Ur,'edgecolor','none');cax=caxis;
% hold on;
% plot_sol(Spartf{2},wr,'edgecolor','none');%caxis(cax);
% colorbar
%%
%%
% lmax = calc_lmax_patch_dirichlet(S,Spartf,Spartfd,AS,Atilde,Apatch,A2patch,psim,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf);
% lmin = calc_lmin_patch_dirichlet(lmax,S,Spartf,Spartfd,AS,Atilde,Apatch,A2patch,psim,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf);
%% algo patch

% rho=2/(lmax+lmin);


rho = 0.4%/lmax;

U_table={};
w_table={};
w1_table={};
erreur_table={};

for k=2:length(myboxes)+1
    lambdaprec{k} = SEPMATRIX([{zeros(getnbddlfree(surfacef{k}),1)} getphi(one(PC))],1);
    wprec{k} = SEPMATRIX([{zeros(getnbddlfree(Spartf{k}),1)} getphi(one(PC))],1);
    w1prec{k} = SEPMATRIX([{zeros(getnbddlfree(Spartfd{k}),1)} getphi(one(PC))],1);
end


%
Uprec = SEPMATRIX([{zeros(getnbddlfree(S),1)} getphi(one(PC))],1);
err=[];
errref=[];
fprintf('\n');

for i=1:nbiter
    temp = SEPMATRIX([{zeros(getnbddlfree(S),1)} getphi(one(PC))],1);
    
    
    for k=2:length(myboxes)+1
        temp = temp-mtimes(Pgamma{k}'*Pgammaoutgammaf{k}',Dlfsep{k},1)*lambdaprec{k};
        % temp = temp - Pgamma{k}'*(Pgammaoutgammaf{k}'*Dlf{k}*lambdaprec{k});
        
        temp = temp + (mtimes(Atilde{k},Uprec,1));
        % temp = temp +  Atilde{k}*Uprec;
    end
    bglob = mtimes(P{1}',bpatchsep{1},1)+temp;
    % Uhat = cgs(AS,  P{1}'*bpatch{1} + temp,tolcgs,cgsiter);
    Uhat = mldivide(AS,bglob,1);
    U = rho * Uhat + (1-rho)*Uprec;
    % U = spectral_decomposition(U);
    %
    if (i ~=1)
        U = multisvd(U,'tol',tolmultisvd,'display',false,'maxorder',100);
    end
    %
    for ii=1:length(numpatchdiff)
        k = numpatchdiff(ii)+1
        
        Bpatchrsep{k} = freevector(bpatchsep{k},1,Spartfd{k})- freevector(Apatchsep{k}*mtimes(Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k},U,1),1,Spartfd{k});
        % Bpatchr{k} = freevector(Spartfd{k},bpatch{k})- freevector(Spartfd{k},Apatch{k}*(Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k}*U));
        % Apatchr{k} = freematrix(Spartfd{k},Apatch{k});
        Apatchrsep{k} = freematrix(Apatchsep{k},1,Spartfd{k});
        
        
        solver = SEPSOLVER(getdim(Apatchrsep{k}),'tol',tolPGD,...
            'update',0,'updatedim',2:getdim(Apatchrsep{k}),...
            'maxorder',PGDrank,'maxiter',5,'reference',[],...
            'errorindicator','none','itercrit',1e-3,...
            'righthandSD',false,'display',true);
        
        [dw1{k},resultsep{k}] = solve(Apatchrsep{k},Bpatchrsep{k}-Apatchrsep{k}*w1prec{k},solver);
        w1{k} = w1prec{k} + dw1{k};
        % w1{k} = unfreevector(w1{k},1,Spartfd{k});
        %
        % w1{k} = dw1{k}+unfreevector(Spartfd{k},w1prec{k});
        % w{k} =  w1{k} + mtimes(Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k},U,1);
        w{k} = calc_ximasse(PCTPMATRIX(unfreevector(w1{k},1,Spartfd{k}),PC))+ ...
            Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k}*PCTPMATRIX(U,PC);
        %
        blambda= mtimes(Ppatchfgammaf{k},Apatchsep{k},1)*mtimes(Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k},U,1)+mtimes(Ppatchfgammaf{k},Apatchsep{k}',1)*unfreevector(w1{k},1,Spartfd{k})-mtimes(Ppatchfgammaf{k},SEPMATRIX([{calc_vector(LINFORM(0),Spartf{k})} getphi(one(PC))],1),1);
        
        %
        lambda{k} = mldivide(random(Dlf{k}),blambda,1);
        
        if (i ~=1)
            lambda{k} = multisvd(lambda{k},'tol',tolmultisvd,'display',false,'maxorder',100);
        end
    end
    
    for ii=1:length(numpatchdir)
        k= numpatchdir(ii)+1;
        Bpatchrsep{k} = freevector(bpatchsep{k},1,Spartfd{k})- freevector(A2patchsep{k}*mtimes(Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k},U,1),1,Spartfd{k});
        Apatchrsep{k} = freematrix(Apatchsep{k},1,Spartfd{k});
        
        solver = SEPSOLVER(getdim(Apatchrsep{k}),'tol',tolPGD,...
            'update',0,'updatedim',2:getdim(Apatchrsep{k}),...
            'maxorder',PGDrank,'maxiter',5,'reference',[],...
            'errorindicator','none','itercrit',1e-3,...
            'righthandSD',false,'display',true);
        
        [dw1{k},resultsep{k}] = solve(Apatchrsep{k},Bpatchrsep{k}-Apatchrsep{k}*w1prec{k},solver);
        w1{k} = w1prec{k} + dw1{k};
        w{k} = psitp{k}.*calc_ximasse((unfreevector(Spartfd{k},PCTPMATRIX(w1{k},PC))))+ ...
            Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k}*PCTPMATRIX(U,PC);
        
        %
        blambdadir = mtimes(Ppatchfgammaf{k}*calc_matrix(BILINFORM(1,1),Spartf{k})*Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k},U,1)+mtimes(Ppatchfgammaf{k},A2patchsep{k}',1)*unfreevector(w1{k},1,Spartfd{k})-SEPMATRIX([{Ppatchfgammaf{k}*calc_vector(LINFORM(0),Spartf{k})} getphi(one(PC))]);
        lambda{k} = mldivide(random(Dlf{k}),blambdadir,1);
        if (i ~=1)
            lambda{k} = multisvd(lambda{k},'tol',tolmultisvd,'display',false,'maxorder',100);
        end
        %
    end
    err(i)=0;
    Q=10;
    for rr=1:Q
        [Ur,r] = random(PCTPMATRIX(U,PC));
        err(i) = err(i) + (norm(Ur-randomeval(PCTPMATRIX(Uprec,PC),r))/norm(Ur));
    end
    err(i) = (1/Q)*err(i)
    for k=2:length(myboxes)+1
        wprec{k} = w{k};
        lambdaprec{k} = lambda{k};
    end
    Uprec = U;
    w1prec=w1;
    fprintf('iteration #%d : ',i)
    fprintf(' ,err stagn = %d\n',err(i))
    
    U_table{i} = {U};
    w_table{i} = {w};
    w1_table{i} = {w1};
    erreur_table{i} = {err(i)};
end
%%
