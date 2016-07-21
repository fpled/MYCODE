D = DOMAIN(2,[0,0],[5,5]);
myboxes = {DOMAIN(2,[3,3],[4,4]),DOMAIN(2,[1,.5],[2,1.5]),DOMAIN(2,[1,2],[2,3]),DOMAIN(2,[3,.5],[4,1.5])};
% myboxes = {DOMAIN(2,[3,3],[4,4])};
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
    Spartf{k} = mesh(myboxes{k-1},ng+20,ng+20);
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
mmm=1;
rrr=1;
ppp=3;

PCphi=cell(1,length(myboxes)+1);
for k=2:(length(myboxes)+1)
    PCphi{k} = POLYCHAOS(setnumber(RANDPOLYS(POLYFE(rrr,ppp),mmm),k-1),ppp,'typebase',1);
    PCphi{k} = calc_masse(PCphi{k});
    if(k==2)
        PCphi{1} = PCphi{2};
    else
        PCphi{1} = union(PCphi{1},PCphi{k});
    end
end
PC = PCphi{1};
% PCphi{1} = union(PCphi{2:end});

%%
erreur_alpha = {};
alpha=[10];
aa=1
% for aa=1:length(alpha)

Kappa = cell(1,length(myboxes)+1);
Kappapc = cell(1,length(myboxes)+1);
kappa{1}=1;
kappapc{1}=1;

for k =2:length(myboxes)+1
    coord = getcoord(getnode(Spartf{k}));
    xi{k} = coord(:,1);
    yi{k} = coord(:,2);
end
%
if length(myboxes)==1
    lc1 = .25;
    kappa{2} = @(x) 1+...
        alpha(aa)*exp(-2*(xi{2}-3.5).^2/lc1^2-2*(yi{2}-3.5).^2/lc1^2)*x(:,1);
    
elseif length(myboxes)==4
    lc1 = .25;
    kappa{2} = @(x) 1+...
        alpha(aa)*exp(-2*(xi{2}-3.5).^2/lc1^2-2*(yi{2}-3.5).^2/lc1^2)*x(:,1);
    
    
    lc2 = .25;
    kappa{3} = @(x) 1+...
        alpha(aa)*exp(-2*(xi{3}-1.5).^2/lc2^2-2*(yi{3}-1).^2/lc2^2)*x(:,1);
    
    lc3 = .25;
    kappa{4} = @(x) 1+...
        alpha(aa)*exp(-2*(xi{4}-1.5).^2/lc3^2-2*(yi{4}-2.5).^2/lc3^2)*x(:,1);
    
    lc4 = .25;
    kappa{5} = @(x) 1+...
        alpha(aa)*exp(-2*(xi{5}-3.5).^2/lc4^2-2*(yi{5}-1).^2/lc4^2)*x(:,1);
end
%  lc4 = .25;
%  kappa{2} = @(x) 1+...
%  10*exp(-2*(xi{2}-3.5).^2/lc4^2-2*(yi{2}-1).^2/lc4^2)*x(:,1);

for k = 2:(length(myboxes)+1)
    kp = kappa{k};
    kappapc{k} = decompmatrixiter(PCphi{k},[],[],@(x) kp(x));
end
%% decomposition du coefficient de diffusion
Km = cell(1,length(myboxes)+1);
result = cell(1,length(myboxes)+1);
km{1} = 1;
for k = 2:(length(myboxes)+1)
    [km{k},result{k}] = spectral_decomposition(kappapc{k},'display',true,'tol',1e-20);
end
%%
% beta = [0 0.3 0.5 ]%0.7 0.9 1 2 10];
beta=0;
erreur_beta = {};
bb=1
% for bb=1:length(beta)

Kappatilde = cell(1,length(myboxes)+1);
kappatilde{1}=1;
for k =2:length(myboxes)+1
    coordt = getcoord(getnode(Spartf{k}));
    xit{k} = coordt(:,1);
    yit{k} = coordt(:,2);
end

if length(myboxes)==1
    lc1 = .25;
    kappatilde{2} =  1+...
        alpha(aa)*exp(-2*(xit{2}-3.5).^2/lc1^2-2*(yit{2}-3.5).^2/lc1^2)*beta(bb);
elseif length(myboxes)==4
    lc1 = .25;
    kappatilde{2} =  1+...
        alpha(aa)*exp(-2*(xit{2}-3.5).^2/lc1^2-2*(yit{2}-3.5).^2/lc1^2)*beta(bb);
    
    
    lc2 = .25;
    kappatilde{3} =  1+...
        alpha(aa)*exp(-2*(xit{3}-1.5).^2/lc2^2-2*(yit{3}-1).^2/lc2^2)*beta(bb);
    %
    lc3 = .25;
    kappatilde{4} =  1+...
        alpha(aa)*exp(-2*(xit{4}-1.5).^2/lc3^2-2*(yit{4}-2.5).^2/lc3^2)*beta(bb);
    %
    lc4 = .25;
    kappatilde{5} =  1+...
        alpha(aa)*exp(-2*(xit{5}-3.5).^2/lc4^2-2*(yit{5}-1).^2/lc4^2)*beta(bb);
end

%% Formulation du prob sans patch sur tout le domaine (Omega)
AS=zeros(getnbddlfree(S));
Atilde = cell(1,length(myboxes)+1);

for k=1:length(myboxes)+1
    Atilde{k} = calc_matrix(BILINFORM(1,1,kappatilde{k},0),Spart{k});
    Atilde{k} = P{k}'*Atilde{k}*P{k};
    AS = AS + Atilde{k};
end

Atilde{1} = calc_matrix(BILINFORM(1,1,kappatilde{1},0),Spart{1});
btilde{1} = calc_vector(LINFORM(0),Spart{1});

%% Formulation du probleme sur le patch.
Apatch = cell(1,length(myboxes)+1);
Apatch{1} = Atilde{1};

for k=2:length(myboxes)+1
    Apatch{k} = PCRADIALMATRIX([getnbddl(Spartf{k}),getnbddl(Spartf{k})],PC);
    Lk = getL(km{k});
    % Lk = decompmatrixiter(PC,[],[],@(x) randomeval(Lk,x(:,k-1)));
    Lk = project(Lk,PC);
    for i=1:getm(km{k})
        ak = BILINFORM(1,1,getV(km{k},i),0);
        Apatch{k} = Apatch{k}+ak{Spartf{k}}(:,:)*Lk(i);
    end
end

bpatch = cell(1,length(myboxes)+1);
lk = LINFORM(0);
bpatch{1} = PCRADIALMATRIX([getnbddlfree(Spart{1}),1],PC);
bpatch{1} = calc_vector(lk,Spart{1});
bpatch{1} = bpatch{1};
% bpatch{1} = P{1}'*bpatch{1};
for k=2:length(myboxes)+1
    bpatch{k} = PCRADIALMATRIX([getnbddlfree(Spartf{k}),1],PC);
    bpatch{k} = bpatch{k}+lk{Spartf{k}}(:)*one(PC);
end

Dl = cell(1,length(myboxes)+1);
dint = BILINFORM(0,0);

for k=2:length(myboxes)+1
    Dl{k} = calc_matrix(dint,surface{k});
    Dl{k} = Dl{k}*one(PC);
    Dl{k} = calc_ximasse(Dl{k});
    % Dl{k} =speye(getnbddlfree(surface{k}));
end

for k=2:length(myboxes)+1
    Dlf{k} = calc_matrix(dint,surfacef{k});
    Dlf{k} = Dlf{k}*one(PC);
    Dlf{k} = calc_ximasse(Dlf{k});
    % Dl{k} =speye(getnbddlfree(surface{k}));
end
% M = cell(1,length(myboxes)+1);
% %
% for k=2:length(myboxes)+1
%     M{k} = [Apatch{k} , -Ppatchfgammaf{k}'*Dl{k};-Dl{k}*Ppatchfgammaf{k} , sparse(getnbddlfree(surfacef{k}),getnbddlfree(surfacef{k}))*one(PC)];
%     M{k} = calc_ximasse(M{k});
% end
%% Solution de reference
% Uref = calcref_patch(Apatch,bpatch,P,1e-15);
couplage=1;

if length(myboxes)==1
    
    if (couplage==1)
        A13 = Phorspatchgamma{2}'*Pgammaoutgammaf{2}'*Dlf{2};
        A23 = -Ppatchfgammaf{2}'*Dlf{2};
        
        A11 = Atilde{1};
        A22 = Apatch{2};
        
        
        Aref = [A11,sparse(getnbddlfree(Spart{1}),getnbddlfree(Spartf{2})),A13;sparse(getnbddlfree(Spartf{2}),getnbddlfree(Spart{1})),A22, A23;A13',A23',sparse(getnbddlfree(surfacef{2}),getnbddlfree(surfacef{2}))]
        bref = [bpatch{1};bpatch{2};zeros(getnbddlfree(surfacef{2}),1)];
        Uzref = solve(Aref,bref);
        
        Uref = Uzref(1:getnbddlfree(Spart{1}));
        wref = Uzref(getnbddlfree(Spart{1})+1:getnbddlfree(Spartf{2})+getnbddlfree(Spart{1}));
    elseif(couplage==2)%grossier
        A13 = Phorspatchgamma{2}'*Dl{2};
        A23 = -Ppatchfgammaf{2}'*Pgammafgammaout{2}'*Dl{2};
        
        A11 = Atilde{1};
        A22 = Apatch{2};
        
        
        Aref = [A11,sparse(getnbddlfree(Spart{1}),getnbddlfree(Spartf{2})),A13;sparse(getnbddlfree(Spartf{2}),getnbddlfree(Spart{1})),A22, A23;A13',A23',sparse(getnbddlfree(surface{2}),getnbddlfree(surface{2}))]
        bref = [bpatch{1};bpatch{2};zeros(getnbddlfree(surface{2}),1)];
        Uzref = solve(Aref,bref);
        
        Uref = Uzref(1:getnbddlfree(Spart{1}));
        wref = Uzref(getnbddlfree(Spart{1})+1:getnbddlfree(Spartf{2})+getnbddlfree(Spart{1}));
        
    elseif couplage==3
        
        
        tempref1 = sparse(getnbddlfree(Spart{1}),getnbddlfree(Spart{1}));
        tempref2 = sparse(getnbddlfree(Spart{1}),getnbddlfree(Spartfd{2}));
        tempref3 = sparse(getnbddlfree(Spartfd{2}),getnbddlfree(Spart{1}));
        for k=2:length(myboxes)+1
            tempref1 = tempref1 + Phorspatchgamma{2}'*Pgammafgammaout{2}*Ppatchfgammaf{k}*calc_matrix(BILINFORM(1,1),Spartf{k})*Ppatchfgammaf{k}'*Pgammafgammaout{2}'*Phorspatchgamma{2};
            tempref2 = tempref2 + Phorspatchgamma{2}'*Pgammafgammaout{2}*Ppatchfgammaf{k}*freevector(Spartfd{k},Apatch{k},1)';
            tempref3 = tempref3 + freevector(Spartfd{k},Apatch{k},1)*Ppatchfgammaf{k}'*Pgammafgammaout{k}'*Phorspatchgamma{k};
        end
        
        Aref = [Atilde{1}+tempref1,tempref2; tempref3,freematrix(Spartfd{k},Apatch{k})];
        Bref = [btilde{1}+ Phorspatchgamma{2}'*Pgammafgammaout{2}*Ppatchfgammaf{2}*calc_vector(LINFORM(0),Spartf{k});freevector(Spartfd{2},bpatch{2})];
        
        Uzref = solve(Aref,Bref);
        Uref = Uzref(1:size(Atilde{1},1));
        
        zref = Uzref(size(Atilde{1},1)+1:end);
        % wref = Ppatchgammaf{k}'*Pfgsurface{k}*Phorspatchgamma{k}*Uref + unfreevector(Sdf{k},zref);
        wref = Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Phorspatchgamma{k}*Uref + unfreevector(Spartfd{k},zref);
        
    end
    
    [Ur,r] = random(Uref);
    wr = randomeval(wref,r);
    
    
    figure(16)
    clf
    
    plot_sol(Spart{1},Ur,'edgecolor','none');cax=caxis;
    hold on;
    plot_sol(Spartf{2},wr,'edgecolor','none');%caxis(cax);
    %
elseif (length(myboxes)==4)
    %
    if couplage==1
        for k=2:length(myboxes)+1
            A13{k} = Phorspatchgamma{k}'*Pgammaoutgammaf{k}'*Dlf{k};
            A23{k} = -Ppatchfgammaf{k}'*Dlf{k};
            A11{1} = Atilde{1};
            A22{k} = Apatch{k};
            sp1{k} = sparse(getnbddlfree(Spart{1}),getnbddlfree(Spartf{k}));
            sp2{k} = sparse(getnbddlfree(Spartf{2}),getnbddlfree(Spartf{k}));
            spg2{k} = sparse(getnbddlfree(Spartf{2}),getnbddlfree(surfacef{k}));
            sp3{k} = sparse(getnbddlfree(Spartf{3}),getnbddlfree(Spartf{k}));
            spg3{k} = sparse(getnbddlfree(Spartf{3}),getnbddlfree(surfacef{k}));
            sp4{k} = sparse(getnbddlfree(Spartf{4}),getnbddlfree(Spartf{k}));
            spg4{k} = sparse(getnbddlfree(Spartf{4}),getnbddlfree(surfacef{k}));
            sp5{k} = sparse(getnbddlfree(Spartf{5}),getnbddlfree(Spartf{k}));
            spg5{k} = sparse(getnbddlfree(Spartf{5}),getnbddlfree(surfacef{k}));
            spf2{k} = sparse(getnbddlfree(surfacef{2}),getnbddlfree(surfacef{k}));
            spf3{k} = sparse(getnbddlfree(surfacef{3}),getnbddlfree(surfacef{k}));
            spf4{k} = sparse(getnbddlfree(surfacef{4}),getnbddlfree(surfacef{k}));
            spf5{k} = sparse(getnbddlfree(surfacef{5}),getnbddlfree(surfacef{k}));
            bzeros{k} = zeros(getnbddlfree(surfacef{k}),1);
            sgp2{k} = sparse(getnbddlfree(surfacef{2}),getnbddlfree(Spartf{k}));
            sgp3{k} = sparse(getnbddlfree(surfacef{3}),getnbddlfree(Spartf{k}));
            sgp4{k} = sparse(getnbddlfree(surfacef{4}),getnbddlfree(Spartf{k}));
            sgp5{k} = sparse(getnbddlfree(surfacef{5}),getnbddlfree(Spartf{k}));
            n{k} = getnbddlfree(Spartf{k});
        end
        n{1} = getnbddlfree(Spart{1});
        Aref = [A11{1},sp1{2},sp1{3},sp1{4},sp1{5},A13{2},A13{3},A13{4},A13{5};...
            sp1{2}',A22{2},sp2{3},sp2{4},sp2{5},A23{2},spg2{3},spg2{4},spg2{5};...
            sp1{3}',sp3{2},A22{3},sp3{4},sp3{5},spg3{2},A23{3},spg3{4},spg3{5};...
            sp1{4}',sp4{2},sp4{3},A22{4},sp4{5},spg4{2},spg4{3},A23{4},spg4{5};...
            sp1{5}',sp5{2},sp5{3},sp5{4},A22{5},spg5{2},spg5{3},spg5{4},A23{5};...
            A13{2}',A23{2}',sgp2{3},sgp2{4},sgp2{5},spf2{2},spf2{3},spf2{4},spf2{5};...
            A13{3}',sgp3{2},A23{3}',sgp3{4},sgp3{5},spf3{2},spf3{3},spf3{4},spf3{5};...
            A13{4}',sgp4{2},sgp4{3},A23{4}',sgp4{5},spf4{2},spf4{3},spf4{4},spf4{5};...
            A13{5}',sgp5{2},sgp5{3},sgp5{4},A23{5}',spf5{2},spf5{3},spf5{4},spf5{5};...
            ];
        bref = [bpatch{1};bpatch{2};bpatch{3};bpatch{4};bpatch{5};...
            bzeros{2};bzeros{3};bzeros{4};bzeros{5}];
        
        Uzref = solve(Aref,bref);
        
        Uref = Uzref(1:n{1});
        wref{2} = Uzref(n{1}+1:n{1}+n{2});
        wref{3} = Uzref(n{1}+n{2}+1:n{1}+n{2}+n{3});
        wref{4} = Uzref(n{1}+n{2}+n{3}+1:n{1}+n{2}+n{3}+n{4});
        wref{5} = Uzref(n{1}+n{2}+n{3}+n{4}+1:n{1}+n{2}+n{3}+n{4}+n{5});
        %
    elseif couplage==2
        
        for k=2:length(myboxes)+1
            A13{k} = Phorspatchgamma{k}'*Dl{k};
            A23{k} = -Ppatchfgammaf{k}'*Pgammafgammaout{k}'*Dl{k};
            A11{1} = Atilde{1};
            A22{k} = Apatch{k};
            sp1{k} = sparse(getnbddlfree(Spart{1}),getnbddlfree(Spartf{k}));
            sp2{k} = sparse(getnbddlfree(Spartf{2}),getnbddlfree(Spartf{k}));
            spg2{k} = sparse(getnbddlfree(Spartf{2}),getnbddlfree(surface{k}));
            sp3{k} = sparse(getnbddlfree(Spartf{3}),getnbddlfree(Spartf{k}));
            spg3{k} = sparse(getnbddlfree(Spartf{3}),getnbddlfree(surface{k}));
            sp4{k} = sparse(getnbddlfree(Spartf{4}),getnbddlfree(Spartf{k}));
            spg4{k} = sparse(getnbddlfree(Spartf{4}),getnbddlfree(surface{k}));
            sp5{k} = sparse(getnbddlfree(Spartf{5}),getnbddlfree(Spartf{k}));
            spg5{k} = sparse(getnbddlfree(Spartf{5}),getnbddlfree(surface{k}));
            spf2{k} = sparse(getnbddlfree(surface{2}),getnbddlfree(surface{k}));
            spf3{k} = sparse(getnbddlfree(surface{3}),getnbddlfree(surface{k}));
            spf4{k} = sparse(getnbddlfree(surface{4}),getnbddlfree(surface{k}));
            spf5{k} = sparse(getnbddlfree(surface{5}),getnbddlfree(surface{k}));
            bzeros{k} = zeros(getnbddlfree(surface{k}),1);
            sgp2{k} = sparse(getnbddlfree(surface{2}),getnbddlfree(Spartf{k}));
            sgp3{k} = sparse(getnbddlfree(surface{3}),getnbddlfree(Spartf{k}));
            sgp4{k} = sparse(getnbddlfree(surface{4}),getnbddlfree(Spartf{k}));
            sgp5{k} = sparse(getnbddlfree(surface{5}),getnbddlfree(Spartf{k}));
            n{k} = getnbddlfree(Spartf{k});
        end
        n{1} = getnbddlfree(Spart{1});
        Aref = [A11{1},sp1{2},sp1{3},sp1{4},sp1{5},A13{2},A13{3},A13{4},A13{5};...
            sp1{2}',A22{2},sp2{3},sp2{4},sp2{5},A23{2},spg2{3},spg2{4},spg2{5};...
            sp1{3}',sp3{2},A22{3},sp3{4},sp3{5},spg3{2},A23{3},spg3{4},spg3{5};...
            sp1{4}',sp4{2},sp4{3},A22{4},sp4{5},spg4{2},spg4{3},A23{4},spg4{5};...
            sp1{5}',sp5{2},sp5{3},sp5{4},A22{5},spg5{2},spg5{3},spg5{4},A23{5};...
            A13{2}',A23{2}',sgp2{3},sgp2{4},sgp2{5},spf2{2},spf2{3},spf2{4},spf2{5};...
            A13{3}',sgp3{2},A23{3}',sgp3{4},sgp3{5},spf3{2},spf3{3},spf3{4},spf3{5};...
            A13{4}',sgp4{2},sgp4{3},A23{4}',sgp4{5},spf4{2},spf4{3},spf4{4},spf4{5};...
            A13{5}',sgp5{2},sgp5{3},sgp5{4},A23{5}',spf5{2},spf5{3},spf5{4},spf5{5};...
            ];
        bref = [bpatch{1};bpatch{2};bpatch{3};bpatch{4};bpatch{5};...
            bzeros{2};bzeros{3};bzeros{4};bzeros{5}];
        
        Uzref = solve(Aref,bref);
        
        Uref = Uzref(1:n{1});
        wref{2} = Uzref(n{1}+1:n{1}+n{2});
        wref{3} = Uzref(n{1}+n{2}+1:n{1}+n{2}+n{3});
        wref{4} = Uzref(n{1}+n{2}+n{3}+1:n{1}+n{2}+n{3}+n{4});
        wref{5} = Uzref(n{1}+n{2}+n{3}+n{4}+1:n{1}+n{2}+n{3}+n{4}+n{5});
        
    end
    figure(14)
    clf
    [Ur,r]=random(Uref);
    plot_sol(Spart{1},Ur,'edgecolor','none');cax=caxis;
    hold on;
    for k=2:length(myboxes)+1
        plot_sol(Spartf{k},randomeval(wref{k},r),'edgecolor','none');%caxis(cax);
    end
    colorbar
end
%% calcul de rho optimal
% lmax = calc_lmax_patch(S,Spartfd,AS,Atilde,Apatch,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf);
% l{dd}=lmax;

%% algo patch
errcgs = {};
% tolcgs = [1e-3 1e-5 1e-9 1e-13];
tolcgs = [1e-1 1e-3 1e-5];
for cg=1:length(tolcgs)
    rho=1%/lmax;
    for k=2:length(myboxes)+1
        if couplage==1
            lambdaprec{k} = zeros(getnbddlfree(surfacef{k}),1);
        elseif couplage==2
            lambdaprec{k} = zeros(getnbddlfree(surface{k}),1);
        end
        wprec{k} =  zeros(getnbddlfree(surfacef{k}),1);
    end
    %
    Uprec = zeros(getnbddlfree(S),1);
    
    err=[];
    errref=[];
    fprintf('\n');
    
    for i=1:30
        temp =zeros(1,getnbddlfree(S))';
        
        
        for k=2:length(myboxes)+1
            
            if couplage==1
                temp = temp - Pgamma{k}'*(Pgammaoutgammaf{k}'*Dlf{k}*lambdaprec{k});
            elseif couplage==2 %grossier
                temp = temp - Pgamma{k}'*Dl{k}*lambdaprec{k};
            end
            temp = temp +  Atilde{k}*Uprec;
        end
        
        Uhat = solve(AS,  P{1}'*btilde{1} + temp);
        U = rho * Uhat + (1-rho)*Uprec;
        %
        for k=2:length(myboxes)+1
            Apatchr{k} = freematrix(Spartfd{k},Apatch{k});
            Bpatchr{k} = freevector(Spartfd{k},bpatch{k}) - freevector(Spartfd{k},Apatch{k}*(Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k}*U));
            %
            w1{k} = solve(Apatchr{k},Bpatchr{k});
            % w1{k} = cgs(Apatchr{k},Bpatchr{k},tolcgs(cg));
            w1{k} = unfreevector(Spartfd{k},w1{k});
            %
            w{k} =  w1{k} + Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k}*U;
            %
            % blambda = Ppatchfgammaf{k}*Apatch{k}*w{k}-Ppatchfgammaf{k}*calc_vector(LINFORM(0),Spartf{k});
            blambda= Ppatchfgammaf{k}*Apatch{k}*Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k}*U+(Ppatchfgammaf{k}*(Apatch{k}'*unfreevector(Spartfd{k},w1{k}))-Ppatchfgammaf{k}*calc_vector(LINFORM(0),Spartf{k}));
            %
            if couplage==1
                lambda{k} = solve(random(Dlf{k}),blambda) ;
            elseif couplage==2
                lambda{k} = solve((random(Dl{k})),Pgammafgammaout{k}*blambda) ;
            end
            %
        end
        err(i) = norm(U-Uprec)/norm(U);
        errref(i) = norm(P{1}*U-Uref)/norm(Uref);
        for k=2:length(myboxes)+1
            wprec{k} = w{k};
            lambdaprec{k} = lambda{k};
        end
        Uprec = U;
        fprintf('iteration #%d : ',i)
        fprintf('err ref = %d , err stagn = %d\n',errref(i),err(i))
    end
    errcgs{cg} = errref;
end
erreur_beta{bb} = errref;
% end
erreur_alpha{aa} = erreur_beta
% end
%%
% save('C:\PROGS\test\erreur_4patchalpha_beta_rho_geoconf_meshnoconf_alpha10_beta0_ktildenew.mat','errref','Uzref','Uref','wref','rho','PC','S','Spart','w','U','lambda','P','kappa')
%%
figure(451)
[Ur,r] = random(U);
plot_sol(Spart{1},P{1}*Ur,'edgecolor','none')
hold on
for k=2:length(myboxes)+1
    plot_sol(Spartf{k},randomeval(w{k},r),'edgecolor','none')
end
colorbar
%%
% load('erreur_4patchalpha_beta_rho_geoconf_meshnoconf_alpha10_beta0_rho1.mat')
chemin=0;
% chemin = 'C:\Documents and Settings\Anthony\Mes documents\My Dropbox\ELIAS-MATHILDE-ANTHONY\CURRENT_PGD_PATCH_STOCHASTIC\figures\';
%%
close all;
maxe=[];
maxv = max([max(mavar(U)) max(mavar(w{2})) max(mavar(w{3})) max(mavar(w{4})) max(mavar(w{5}))])
for k=1:length(myboxes)
    maxe(k) = max([max(var_expect_cond(P{1}*U,k)) max(var_expect_cond(w{2},k)) max(var_expect_cond(w{3},k)) max(var_expect_cond(w{4},k)) max(var_expect_cond(w{5},k))])
end
carte_stat_patch_nonconf(S,Spart,Spartf,surface,U,w,Uref,Uzref,wref,P,maxv,maxe,'chemin',chemin)

%%
close all;plot_patch_nonconf(S,Spartf,Spart,U,w,Uref,wref,P,'grad',1,'chemin',chemin)
%%
close all;ploterrorUwU(w,U,Uzref)

