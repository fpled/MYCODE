clear all
clc
%% number of simulation
N=30; % number of simulations
%% Select  values
lamda = -19; % model parameter controlling the level of statistical fluctiations (choice)
% E_bar_exp = 11.55;
% nu_bar_exp = 0.4;
k_bar = 19.25;
m_bar = 4.125;
phat0(1)=1-lamda;
phat0(2)=k_bar/(1-lamda);
phat0(3)=m_bar/(1-5*lamda);

k_0 = gamrnd(1-lamda,k_bar/(1-lamda),N,1);
m_0 = gamrnd(1-5*lamda,m_bar/(1-5*lamda),N,1);
e_0 = (9*k_0.*m_0)./(m_0+3*k_0);
n_0 = (3*k_0-2*m_0)./(6*k_0+2*m_0);

%% p.d.f reference
pdf_C0 = @(x,y) 1/(phat0(2)^phat0(1)*gamma(phat0(1)))*(x^(phat0(1)-1))*exp(-x/phat0(2))*...
    1/(phat0(3)^(5*phat0(1)-4)*gamma(5*phat0(1)-4))*(y^(5*phat0(1)-5))*exp(-y/phat0(3));

%% Reverse identification
ResulIdent=[];
x_exp=[e_0 n_0];
for i=1:N
    U_exp(:,i)=solveThreePointBendingIsot(x_exp(i,:));
    ampl = 1e-2;
    noise = ampl*(2*rand(length(U_exp(:,i)),1)-1).*U_exp(:,i);    
    x=functionIdentificationIsot(U_exp(:,i),noise);
    ResulIdent = [ResulIdent;x];
end
e = ResulIdent(:,1);
n = ResulIdent(:,2);
k = e./(-6*n+3);
m = e./(2*n+2);

%% Plot samples
fontsize = 16;
linewidth = 1;
markersize = 36;
interpreter = 'latex';

figure('Name','Experimental data of Bulk modulus')
clf
bar(1:length(k),k)
set(gca,'FontSize',fontsize)
set(gca,'XLim',[0,length(k)+1])
xlabel('\''Echantillon','Interpreter',interpreter);
ylabel('Bulk modulus (GPa)','Interpreter',interpreter);

figure('Name','Experimental data of Shear modulus')
clf
bar(1:length(m),m)
set(gca,'FontSize',fontsize)
set(gca,'XLim',[0,length(m)+1])
xlabel('\''Echantillon','Interpreter',interpreter);
ylabel('Shear modulus (GPa)','Interpreter',interpreter);


%% Maximum likelihood estimation
opts = statset('TolFun',eps,'TolX',eps,'FunValCheck','off','MaxIter',1e3,'MaxFunEvals',1e3);
data=[k;m];
nloglf = @(phat,data,cens,freq) length(k)*log(gamma(phat(1)))+length(k)*phat(1)*log(phat(2))...
      +(1-phat(1))*sum(log(data(1:length(k))))+1/phat(2)*sum(data(1:length(k)))...
      +length(m)*log(gamma(5*phat(1)-4))+length(m)*(5*phat(1)-4)*log(phat(3))...
      +5*(1-phat(1))*sum(log(data(length(k)+1:end)))+1/phat(3)*sum(data(length(k)+1:end));
phat = mle(data,'nloglf',nloglf,'start',[10 2 1],'options',opts)

%% Plot pdf 

pdf_C = @(x,y) 1/(phat(2)^phat(1)*gamma(phat(1)))*(x^(phat(1)-1))*exp(-x/phat(2))*...
    1/(phat(3)^(5*phat(1)-4)*gamma(5*phat(1)-4))*(y^(5*phat(1)-5))*exp(-y/phat(3));

mk=mean(k);
sk=std(k);
mm=mean(m);
sm=std(m);
xmin = max(0,mk-5*sk);
xmax = mk+5*sk;
x = linspace(xmin,xmax,1e2);
ymin = max(0,mm-5*sm);
ymax = mm+5*sm;
y = linspace(ymin,ymax,1e2);

[X2,Y2]=meshgrid(x,y);
for i=1:length(x)
    for j=1:length(x)
        Z2(i,j)= pdf_C(X2(i,j),Y2(i,j));
    end
end
figure
mesh(X2,Y2,Z2)
set(gca,'FontSize',fontsize)
xlabel('$c_1$','Interpreter',interpreter);
ylabel('$c_2$','Interpreter',interpreter);
zlabel('$p_C(c)$','Interpreter',interpreter);



pdf_Enu = @(x,y) (x/(-6*y+3))^(phat(1)-1)*(x/(2*y+2))^(5*phat(1)-5)...
    *x/(2*(1+y)^2*(1-2*y)^2)*exp(-1/phat(2)*x/(-6*y+3))*exp(-1/phat(3)*x/(2*y+2))*...
    1/(gamma(phat(1))*phat(2)^phat(1))*1/1/(gamma((5*phat(1)-4))*phat(3)^(5*phat(1)-4));
me=mean(e);
se=std(e);
mn=mean(n);
sn=std(n);
xmin = max(0,me-5*se);
xmax = me+5*se;
x = linspace(xmin,xmax,1e2);
ymin = max(0,mn-5*sn);
ymax = 0.5-sn;
y = linspace(ymin,ymax,1e2);
% [X,Y]=meshgrid(e,n);
% for i=1:length(e)
%     for j=1:length(e)
%         Z(i,j)= pdf_Enu(X(i,j),Y(i,j));
%     end
% end
[X2,Y2]=meshgrid(x,y);
for i=1:length(x)
    for j=1:length(y)
        Z2(i,j)= pdf_Enu(X2(i,j),Y2(i,j));
    end
end
figure
mesh(X2,Y2,Z2)
set(gca,'FontSize',fontsize)
xlabel('$e$','Interpreter',interpreter);
ylabel('$\nu$','Interpreter',interpreter);
zlabel('$p_{E,\nu}(e,n)$','Interpreter',interpreter);




