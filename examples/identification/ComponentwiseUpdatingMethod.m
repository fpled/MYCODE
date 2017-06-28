%Componentwise updating method;
%Page 31 in 'Computational Statistics with Matlab';
%Componentwise Metropolis sampler (the proposal distribution is symmetric)
%with normal distrubution for the components c1 c2 c3;

clc
clear all
close all

load('resultats_ET_GL_EL_nuL_all_samples')
fontsize = 16;
linewidth = 1;
markersize = 36;
interpreter = 'latex';
formats = {'fig','epsc2'};

Num_sample = 20;
for i = 1:Num_sample
   num_ech = ['C' num2str(i)];
   eval( ['ET_data(i) = mean_ET_' num_ech '_rest;'] ); %Gpa
   eval( ['GL_data(i) = mean_GL_' num_ech '_rest/1000;'] ); %Gpa 
   eval( ['EL_data(i) = mean_EL_' num_ech '_rest/1000;'] ); %Gpa
   eval( ['nuL_data(i) = mean_nuL_' num_ech '_rest;'] );  
end
nuT_data = 0.1+0.2*rand(1,20); %20 artificial datas for nuT vary from 0.1 to 0.3
GT_data = ET_data./(2*(1+nuT_data));
kT_data = (EL_data.*ET_data)./(2*(1-nuT_data).*EL_data-4*ET_data.*(nuL_data).^2);
C1_data = EL_data + 4*(nuL_data.^2).*kT_data;
C2_data = 2*kT_data; 
C3_data = 2*sqrt(2)*kT_data.*nuL_data;
C4_data = 2*GT_data;
C5_data = 2*GL_data;  
mc1 = mean(C1_data);
mc2 = mean(C2_data);
mc3 = mean(C3_data);
mc4 = mean(C4_data);
mc5 = mean(C5_data);
sc1 = std(C1_data);
sc2 = std(C2_data);
sc3 = std(C3_data);

%%We should proposer a method to calculate lambda which
%control the level of fluctuation
lambda = -20; % negative number
lambda1 = -(mc2*lambda)/(-mc3^2+mc1*mc2);
lambda2 = -(mc1*lambda)/(-mc3^2+mc1*mc2);
lambda3 = (2*mc3*lambda)/(-mc3^2+mc1*mc2);  
lambda4 = (1-2*lambda)/mc4;
lambda5 = (1-2*lambda)/mc5; %gamma parameters

%%Parameters of the trivariate normal distribution
Mu_mean = [mc1 mc2 mc3];
Sigma_mean = [ -mc1^2/lambda -mc3^2/lambda -(mc1*mc3)/lambda
               -mc3^2/lambda -mc2^2/lambda -(mc2*mc3)/lambda
               -(mc1*mc3)/lambda -(mc2*mc3)/lambda -(mc3^2+mc1*mc2)/(2*lambda) ];

%%Start values for c1 c2 c3
c1 = unifrnd( min(C1_data), max(C1_data) );
c2 = unifrnd( min(C2_data), max(C2_data) );
c3 = unifrnd( min(C3_data), max(C3_data) );
T = 5000;
t = 1;
Samples_C = zeros(5,T);
Samples_C(1:3,t) = [c1 c2 c3]'; 
while t<T
 t = t+1;
 
 new_c1 = normrnd( c1, sc1 );
 pratio = mvnpdf( [new_c1 c2 c3],Mu_mean, Sigma_mean )/...
          mvnpdf( [    c1 c2 c3],Mu_mean, Sigma_mean );
 alpha = min([1 pratio]);
 u = rand;
 if u<alpha
 Samples_C(1,t) = new_c1;
 else 
 Samples_C(1,t) = Samples_C(1,t-1);   
 end
 
 new_c2 = normrnd( c2, sc2 );
 pratio = mvnpdf( [c1 new_c2 c3],Mu_mean, Sigma_mean )/...
          mvnpdf( [c1     c2 c3],Mu_mean, Sigma_mean );
 alpha = min([1 pratio]);
 u = rand;
 if u<alpha
 Samples_C(2,t) = new_c2;
 else 
 Samples_C(2,t) = Samples_C(2,t-1);
 end
 
 new_c3 = normrnd( c3, sc3 );
 pratio = mvnpdf( [c1 c2 new_c3],Mu_mean, Sigma_mean )/...
          mvnpdf( [c1 c2     c3],Mu_mean, Sigma_mean );
 alpha = min([1 pratio]);
 u = rand;
 if u<alpha
 Samples_C(3,t) = new_c3;
 else 
 Samples_C(3,t) = Samples_C(3,t-1); 
 end
 
end

Samples_C(4,:) = gamrnd(1-2*lambda,1/lambda4,T,1);
Samples_C(5,:) = gamrnd(1-2*lambda,1/lambda5,T,1);
kT_Model= Samples_C(2,:)/2;
nuL_Model = (Samples_C(3,:)./kT_Model)/(2*sqrt(2));
EL_Model = Samples_C(1,:) - 4*(nuL_Model.^2).*kT_Model;
GT_Model = Samples_C(4,:)/2;
GL_Model = Samples_C(5,:)/2;
cons1 = (kT_Model.^(-1))+4*(nuL_Model.^2)./EL_Model;
cons2 =  GT_Model.^(-1);
ET_Model = 4*(cons1+cons2).^(-1);
nuT_Model = (ET_Model./GT_Model)/2-1;


figure('name','random samples of (ET,GL,nuT)')
clf
scatter3(ET_Model,GL_Model,nuT_Model,'.')
hold on
scatter3(ET_data,GL_data,nuT_data,'r+')
hold off
xlabel('Young modulus $E^T$ (GPa)','Fontsize',14,'Interpreter','latex');
ylabel('Shear modulus $G^L$ (GPa)','Fontsize',14,'Interpreter','latex');
zlabel('Poisson''s ratio $\nu^T$','Fontsize',14,'Interpreter','latex');

%% Plot pdfs
pdf_C4 = @(c4) gampdf(c4,1-2*lambda,1/lambda4); % Gamma probability density function of C4
pdf_C5 = @(c5) gampdf(c5,1-2*lambda,1/lambda5); % Gamma probability density function of C5
cdf_C4 = @(c4) gamcdf(c4,1-2*lambda,1/lambda4); % Gam
ma cumulative density function of C1
cdf_C5 = @(c5) gamcdf(c5,1-2*lambda,1/lambda5); % Gamma cumulative density function of C2

N = 100;
x4_min = max(0,mean(C4_data)-8*std(C4_data));
x4_max = mean(C4_data)+8*std(C4_data);
x5_min = max(0,mean(C5_data)-8*std(C5_data));
x5_max = mean(C5_data)+8*std(C5_data);

% Plot pdf of C1
figure('Name','Probability density function of C4')
clf
x4 = linspace(x4_min,x4_max,N);
y4 = pdf_C4(x4);
plot(x4,y4,'-b');
hold on
plot(C4_data,pdf_C4(C4_data),'r+');
hold off
grid on
box on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[min(x4),max(x4)])
xlabel('$c_4$ (Gpa)','Interpreter',interpreter)
ylabel('$p_{C_4}(c_4)$','Interpreter',interpreter)

% Plot pdf of C2
figure('Name','Probability density function of C5')
clf
x5 = linspace(x5_min,x5_max,N);
y5 = pdf_C5(x5);
plot(x5,y5,'-b');
hold on
plot(C5_data,pdf_C5(C5_data),'r+');
hold off
grid on
box on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[min(x5),max(x5)])
xlabel('$c_5$ (GPa)','Interpreter',interpreter)
ylabel('$p_{C_5}(c_5)$','Interpreter',interpreter)

%% Plot cdfs
% Plot cdf of C4
figure('name','cumulative distribution function of C4')
clf
z4 = cdf_C4(x4);
plot(x4,z4,'-b');
hold on
plot(C4_data,cdf_C4(C4_data),'r+');
hold off
grid on
box on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[min(x4),max(x4)])
xlabel('$c_4$ (GPa)','Interpreter',interpreter)
ylabel('$F_{C_4}(c_4)$','Interpreter',interpreter)

% Plot cdf of C2
figure('name','cumulative distribution function of C5')
clf
z5 = cdf_C5(x5);
plot(x5,z5,'-b');
hold on
plot(C5_data,cdf_C5(C5_data),'r+');
hold off
grid on
box on
set(gca,'FontSize',fontsize)
set(gca,'XLim',[min(x5),max(x5)])
xlabel('$c_5$ (GPa)','Interpreter',interpreter)
ylabel('$F_{C_5}(c_5)$','Interpreter',interpreter)


