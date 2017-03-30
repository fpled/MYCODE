clc
clear all
close all
format long
% myparallel('start');

% VOption = 'Chol';
VOption = 'LDL';

% parametrage
% d=2;                                       % dimension de l'espace
n=3;                                           % dimension de la matrice
lx=1;      ly=1;                             % dimension du domaine
Lx=lx/50;     Ly=ly/50;                 % spatial correlation length
nelx=20;   nely=20;                     % spatial discretization
delta=0.5;
if delta < 0
    error('valeur de delta < 0');
elseif delta > sqrt((n+1)/(n+5))
    error('valeur de delta > ((n+1)(n+5)^-1)^0.5')
end
NumSimpower = 1:0.2:3;
NumSim = round(10.^NumSimpower);
NumSimMax=max(NumSim);

% information disponible : la moyenne de [A]
R=randn(n,n);
an=R*R';
LAnMat=chol(an);
deltaAn = delta/sqrt(n+1) * sqrt(1+(trace(an))^2/(trace(an^2)));

dx = lx/nelx;   dy = ly/nely;
PIx1 = dx*(1-1/sqrt(3))/2;   PIx2 = dx*(1+1/sqrt(3))/2;
PIy1 = dy*(1-1/sqrt(3))/2;   PIy2 = dy*(1+1/sqrt(3))/2;
PIx1Vec = PIx1:dx:(nelx-1)*dx+PIx1;   PIx2Vec = PIx2:dx:(nelx-1)*dx+PIx2;
PIy1Vec = PIy1:dy:(nely-1)*dy+PIy1;   PIy2Vec = PIy2:dy:(nely-1)*dy+PIy2;
etaxVec = union(PIx1Vec,PIx2Vec);     etayVec = union(PIy1Vec,PIy2Vec);
[PIx,PIy]=meshgrid(etaxVec,etayVec);
PIx=PIx';     PIy=PIy';
PIxVec = PIx(:);   PIyVec = PIy(:); 
PIx=PIx';     PIy=PIy';
N = length(PIxVec);   Nx = size(PIx,1);   Ny = size(PIx,2);

CVorg=zeros(N,N);
for i=1:N
        etax=PIxVec-PIxVec(i);
        etay=PIyVec-PIyVec(i);
        
        RU = (sinc(etax/(2*Lx)).^2) .* (sinc(etay/(2*Ly)).^2);
        CVorg(:,i)=RU;
end
CVorg(1:N+1:end)=1;
CVorg=0.5*(CVorg+CVorg');     % symetrisation

switch VOption
    case 'Chol'
        tic
        [T,num] = cholcov(CVorg,0);
        if num==0
            LV=chol(CVorg);
        else 
            [v,d] = eig(CVorg);
            D=diag(d);
            DD = find(D<=1e-12);
            D(DD) = 1e-10;
            d = diag(D);
            CVrcd=v*d*v';
            LV = chol(CVrcd);
        end
        toc
        VV = randn(N,n*(n+1)/2*NumSimMax);
        tic
        V = LV * VV;
        toc
        VTrp = V';
        VPrm = permute(reshape(VTrp,n*(n+1)/2,NumSimMax,N),[1 3 2]);
        
    case 'LDL'
        tic
        [L,D]=ldl(CVorg);
        [v0,d0]=eig(D);
        dd = diag(d0);     dd(dd<10e-10) = [];
        DD = diag(dd);    VV = v0(:,end-size(DD,1)+1:end);
        Dd = VV*sqrt(DD)*VV';
        toc
        VV = randn(N,n*(n+1)/2*NumSimMax);
        tic
        V = L*Dd*VV;
        toc
        VTrp = V';
        VPrm = permute(reshape(VTrp,n*(n+1)/2,NumSimMax,N),[1 3 2]);
end

% tic
% moyV=zeros(length(Num_Sim),1);
% errCV=zeros(length(Num_Sim),1);
% for i=1:length(Num_Sim)
%     Vsim=V(:,1:Num_Sim(i));
%     moyVsim=mean(Vsim,2);
%     moyV(i) = norm(moyVsim,inf);
%     
%     moyVsimMat = repmat(moyVsim,1,size(Vsim,2));
%     RVVec = Vsim-moyVsimMat;
%     
%     CVcal = cov(RVVec');
% %     CVcal = corrcoef(RVVec');
%     
%     errCV(i) = norm(CVcal-CVorg,'fro')/norm(CVorg,'fro')*100;
% end
% toc

% figure
% semilogx(Num_Sim,moyV);
% figure
% semilogx(Num_Sim,errCV);
% figure
% loglog(Num_Sim,moyV,Num_Sim,1./sqrt(Num_Sim));
% figure
% loglog(Num_Sim,errCV,Num_Sim,1./sqrt(Num_Sim));

% CVorg1D = reshape(CVorg(size(CVorg,1)/2-Ny/2,:),Nx,Ny);
% CVcal1D = reshape(CVcal(size(CVcal,1)/2-Ny/2,:),Nx,Ny);
% 
% figure
% surf(PIx,PIy,CVorg1D,'EdgeColor','none');
% view(3);
% axis equal; axis tight;
% shading interp
% figure
% surf(PIx,PIy,CVcal1D,'EdgeColor','none');
% view(3);
% axis equal; axis tight;
% shading interp

tic
VMatGauss = triu(ones(n,n),1);   VMatGamma = diag(ones(n,1));
VMatGauss = repmat(VMatGauss,1,1,N*NumSimMax);   VMatGamma = repmat(VMatGamma,1,1,N*NumSimMax);
VMatGauss(VMatGauss==1) = VPrm(1:n*(n-1)/2,:,:);   VMatGamma(VMatGamma==1) = VPrm(n*(n-1)/2+1:end,:,:);

sigma_n = delta*(n+1)^(-0.5);   % standard_deviation
LnGauss = sigma_n * VMatGauss;
FU = normcdf(VMatGamma);
alphaVec = (n+1)/(2*delta^2) + (1-(1:n)')/2;
alphaMat = diag(alphaVec);
alphaTen =  repmat(alphaMat,1,1,N*NumSimMax);
BVec=ones(n,1);
BMat= diag(BVec);
BTen =  repmat(BMat,1,1,N*NumSimMax);
h = gaminv(FU,alphaTen,BTen);
LnGamma = sigma_n * sqrt(2*h);
Ln = LnGauss + LnGamma;

LnTrp = permute(Ln,[2 1 3]);
aa = 1:n*N*NumSimMax;
ai = reshape(aa,[],N*NumSimMax);
bi = repmat(ai,n,1);
bj=aa(ones(n,1),:);
ii=reshape(bi,1,numel(bi));
jj=reshape(bj,1,numel(bj));
LnTrpMat = sparse(ii,jj,LnTrp(:));
LnPrm = permute(Ln,[1,3,2]);
LnMat = reshape(LnPrm,[],size(Ln,2),1);

GnMat = LnTrpMat * LnMat;
Gn = permute(reshape(GnMat, n, N*NumSimMax, n), [1 3 2]);

LAn = repmat(LAnMat,1,1,N*NumSimMax);
LAnPrm = permute(LAn,[1,3,2]);
LAnPrmMat = reshape(LAnPrm,[],size(LAn,2),1);
LAnTrp = repmat(LAnMat',1,1,N*NumSimMax);
LAnTrpMat = sparse(ii,jj,LAnTrp(:));
BnMat = LAnTrpMat * GnMat;
Bn = permute(reshape(BnMat, n, N*NumSimMax, n), [1 3 2]);
BnPIMat = sparse(ii,jj,Bn(:));
AnMat = BnPIMat * LAnPrmMat;
An = permute(reshape(AnMat, n, N*NumSimMax, n), [1 3 2]);

An11 = An(1,1,:);
toc

randSim = round(rand(5,1)*NumSimMax);
for i=1:length(randSim)
    indSim = randSim(i);
    An11Mat = reshape(An11(:,:,(indSim-1)*N+1:indSim*N),Nx,Ny);
    An11Mat = An11Mat';
    
    % figure
    % surface(PI1,PI2,A_n11Mat);
    figure
    hold on
    surf(PIx,PIy,An11Mat,'EdgeColor','none');
    % view(3)
    axis equal; axis tight;
    % axis off;
    tt=title(['Random Field A_n_,_1_1(x)']);
    set(tt,'FontName','Bell MT','Fontsize',12);
    shading interp
end

% compute spatial correlation lengths relative to x for Gn and An
GAaax = 1:n*Nx;
GAaix = reshape(GAaax,[],Nx);
GAbix = repmat(GAaix,n,1);
GAbjx=GAaax(ones(n,1),:);
GAiix=reshape(GAbix,1,numel(GAbix));
GAjjx=reshape(GAbjx,1,numel(GAbjx));
rGnx=zeros(Ny,Nx);
rAnx=zeros(Ny,Nx);
LGnxVec=zeros(Ny,1);
LAnxVec=zeros(Ny,1);
for i=1:Ny
    trGGnnx=zeros(NumSimMax,Nx);
    trAAnnx=zeros(NumSimMax,Nx);
    
    parfor j=1:NumSimMax
        Gnn = Gn;      Ann = An;
        Gnx = Gnn(:,:,((j-1)*N+(i-1)*Nx+1:(j-1)*N+(i-1)*Nx+Nx));      Anx = Ann(:,:,((j-1)*N+(i-1)*Nx+1:(j-1)*N+(i-1)*Nx+Nx));
        GnxMat = sparse(GAiix,GAjjx,Gnx(:));      AnxMat = sparse(GAiix,GAjjx,Anx(:));
        Gnx1Prm = permute(Gnx(:,:,1),[1,3,2]);      Anx1Prm = permute(Anx(:,:,1),[1,3,2]);
        Gnx1Mat = reshape(Gnx1Prm,[],size(Gnx,2),1);       Anx1Mat = reshape(Anx1Prm,[],size(Anx,2),1);
        Gnx1Mat = repmat(Gnx1Mat,Nx,1);      Anx1Mat = repmat(Anx1Mat,Nx,1);
        GGnx = GnxMat * Gnx1Mat;      AAnx = AnxMat * Anx1Mat;
        for k=1:Nx
            GGnnx = GGnx((k-1)*n+1:(k-1)*n+n,:);      AAnnx = AAnx((k-1)*n+1:(k-1)*n+n,:);
            trGGnnx(j,k) = trace(GGnnx);      trAAnnx(j,k) = trace(AAnnx);
        end
    end
    
    moytrGGnnx = mean(trGGnnx);      moytrAAnnx = mean(trAAnnx);
    rGnx(i,:) = 1/delta^2 * (1/n*moytrGGnnx-1);      rAnx(i,:) = (moytrAAnnx-trace(an^2)) / (deltaAn^2*trace(an^2));
    LGnxVec(i) = trapz(PIx(i,:),rGnx(i,:));      LAnxVec(i) = trapz(PIx(i,:),rAnx(i,:));
end
moyrGnx=mean(rGnx);
moyrAnx=mean(rAnx);
set(figure,'defaulttextinterpreter','latex'); clf
plot(PIx(1,:),moyrGnx);
tx = xlabel('x'); ty = ylabel('rGn(x)');
set(tx,'Fontsize',10); set(ty,'Fontsize',10);
tt=title(['Graph of function x $\mapsto$ rGn(x)']) ;
set(tt,'FontName','Bell MT','Fontsize',12); set(gca,'FontName','Times','FontSize',10) ;

set(figure,'defaulttextinterpreter','latex'); clf;
plot(PIx(1,:),moyrAnx);
tx = xlabel('x'); ty = ylabel('rAn(x)');
set(tx,'Fontsize',10); set(ty,'Fontsize',10);
tt=title(['Graph of function x $\mapsto$ rAn(x)']) ;
set(tt,'FontName','Bell MT','Fontsize',12); set(gca,'FontName','Times','FontSize',10) ;
LGnx = mean(LGnxVec);
LAnx = mean(LAnxVec);

% compute spatial correlation lengths relative to y for Gn and An
GAaay = 1:n*Ny;
GAaiy = reshape(GAaay,[],Ny);
GAbiy = repmat(GAaiy,n,1);
GAbjy=GAaay(ones(n,1),:);
GAiiy=reshape(GAbiy,1,numel(GAbiy));
GAjjy=reshape(GAbjy,1,numel(GAbjy));
rGny=zeros(Nx,Ny);
rAny=zeros(Nx,Ny);
LGnyVec=zeros(Nx,1);
LAnyVec=zeros(Nx,1);
for i=1:Nx
    trGGnny=zeros(NumSimMax,Ny);
    trAAnny=zeros(NumSimMax,Ny);
    
    parfor j=1:NumSimMax
        Gnn = Gn;      Ann = An;
        Gny = Gnn(:,:,((j-1)*N+(i-1)*Ny+1:(j-1)*N+(i-1)*Ny+Ny));      Any = Ann(:,:,((j-1)*N+(i-1)*Ny+1:(j-1)*N+(i-1)*Ny+Ny));
        GnyMat = sparse(GAiiy,GAjjy,Gny(:));      AnyMat = sparse(GAiiy,GAjjy,Any(:));
        Gny1Prm = permute(Gny(:,:,1),[1,3,2]);      Any1Prm = permute(Any(:,:,1),[1,3,2]);
        Gny1Mat = reshape(Gny1Prm,[],size(Gny,2),1);       Any1Mat = reshape(Any1Prm,[],size(Any,2),1);
        Gny1Mat = repmat(Gny1Mat,Ny,1);      Any1Mat = repmat(Any1Mat,Ny,1);
        GGny = GnyMat * Gny1Mat;      AAny = AnyMat * Any1Mat;
        for k=1:Ny
            GGnny = GGny((k-1)*n+1:(k-1)*n+n,:);      AAnny = AAny((k-1)*n+1:(k-1)*n+n,:);
            trGGnny(j,k) = trace(GGnny);      trAAnny(j,k) = trace(AAnny);
        end
    end
    
    moytrGGnny = mean(trGGnny);      moytrAAnny = mean(trAAnny);
    rGny(i,:) = 1/delta^2 * (1/n*moytrGGnny-1);      rAny(i,:) = (moytrAAnny-trace(an^2)) / (deltaAn^2*trace(an^2));
    LGnyVec(i) = trapz(PIy(:,i),rGny(i,:));      LAnyVec(i) = trapz(PIy(:,i),rAny(i,:));
end
moyrGny=mean(rGny);
moyrAny=mean(rAny);
set(figure,'defaulttextinterpreter','latex'); clf;
plot(PIy(:,1),moyrGny);
tx = xlabel('y'); ty = ylabel('rGn(y)');
set(tx,'Fontsize',10); set(ty,'Fontsize',10);
tt=title(['Graph of function y $\mapsto$ rGn(y)']) ;
set(tt,'FontName','Bell MT','Fontsize',12); set(gca,'FontName','Times','FontSize',10) ;

set(figure,'defaulttextinterpreter','latex'); clf;
plot(PIy(:,1),moyrAny);
tx = xlabel('y'); ty = ylabel('rAn(y)');
set(tx,'Fontsize',10); set(ty,'Fontsize',10);
tt=title(['Graph of function y $\mapsto$ rAn(y)']) ;
set(tt,'FontName','Bell MT','Fontsize',12); set(gca,'FontName','Times','FontSize',10) ;
LGny = mean(LGnyVec);
LAny = mean(LAnyVec);

% myparallel('stop');