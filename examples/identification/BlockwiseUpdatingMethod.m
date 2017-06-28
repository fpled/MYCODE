%Blockwise updating Method; 
%Page 26 in 'Computational Statistics with Matlab';
%See also page 4 in 'Metropolis-Hastings sampling using multivariate 
%gaussian tangents';
%This method sometimes doesn't work, because the Positive-definite of
%covariance matrix SIGMA isn't always guaranteed;
%If there is a error, rerun the code

clc
clear all
close all

% syms c1 c2 c3 lambda lambda1 lambda2 lambda3
% PHI = -lambda1*c1 - lambda2*c2 - lambda3*c3 - lambda*log(c1*c2-c3^2);
% SIGMA = -hessian(PHI,[c1 c2 c3])^(-1)
% gx = gradient(PHI,c3)

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
lambda = -40; % negative number
lambda1 = -(mc2*lambda)/(-mc3^2+mc1*mc2);
lambda2 = -(mc1*lambda)/(-mc3^2+mc1*mc2);
lambda3 = (2*mc3*lambda)/(-mc3^2+mc1*mc2);  
lambda4 = (1-2*lambda)/mc4;
lambda5 = (1-2*lambda)/mc5; %gamma parameters
Mu_mean = [mc1 mc2 mc3];
Sigma_mean = [ -mc1^2/lambda -mc3^2/lambda -(mc1*mc3)/lambda
               -mc3^2/lambda -mc2^2/lambda -(mc2*mc3)/lambda
               -(mc1*mc3)/lambda -(mc2*mc3)/lambda -(mc3^2+mc1*mc2)/(2*lambda) ];
f = @(c) ( -lambda1*c(1)-lambda2*c(2)-lambda3*c(3)-lambda*log(c(1)*c(2)-c(3)^2) );        
%
T = 5000;
Samples_C = zeros(5,T);
t = 0;
while t<T
 t = t+1;
    c_old = mvnrnd(Mu_mean,Sigma_mean);
    c1_old = c_old(1);
    c2_old = c_old(2);
    c3_old = c_old(3);
    f_old = f(c_old);
    g_old = [-lambda1-(c2_old*lambda)/(-c3_old^2+c1_old*c2_old)
             -lambda2-(c1_old*lambda)/(-c3_old^2+c1_old*c2_old)
             -lambda3+(2*c3_old*lambda)/(-c3_old^2+c1_old*c2_old)];
    sigma_old = [-c1_old^2/lambda -c3_old^2/lambda -(c1_old*c3_old)/lambda
                 -c3_old^2/lambda -c2_old^2/lambda -(c2_old*c3_old)/lambda
                 -(c1_old*c3_old)/lambda -(c2_old*c3_old)/lambda -(c3_old^2+c1_old*c2_old)/(2*lambda)];
    mu_old = (c_old' + sigma_old*g_old)';
    q_old = mvnpdf(c_old, mu_old, sigma_old);
    %
    c_prop = mvnrnd(mu_old,sigma_old);
    c1_prop = c_prop(1);
    c2_prop = c_prop(2);
    c3_prop = c_prop(3);
    log_q_prop = log(q_old);
    f_prop = f(c_prop);
    g_prop = [-lambda1-(c2_prop*lambda)/(-c3_prop^2+c1_prop*c2_prop)
              -lambda2-(c1_prop*lambda)/(-c3_prop^2+c1_prop*c2_prop)
              -lambda3+(2*c3_prop*lambda)/(-c3_prop^2+c1_prop*c2_prop)];
    sigma_prop = [-c1_prop^2/lambda -c3_prop^2/lambda -(c1_prop*c3_prop)/lambda
                  -c3_prop^2/lambda -c2_prop^2/lambda -(c2_prop*c3_prop)/lambda
                  -(c1_prop*c3_prop)/lambda -(c2_prop*c3_prop)/lambda -(c3_prop^2+c1_prop*c2_prop)/(2*lambda)];
    mu_prop = (c_prop' + sigma_prop*g_prop)';           
    q_prop = mvnpdf(c_prop, mu_prop, sigma_prop);
    log_q_old = log(q_prop);
    r = exp( (f_prop-f_old)+(log_q_old-log_q_prop) );
    if r>1
    Samples_C(1:3,t) = c_prop;
    else 
        s=rand;
        if s<r
        Samples_C(1:3,t) = c_prop;    
        else
        Samples_C(1:3,t) = c_old;
        end
    end
end
Samples_C(4,:) = gamrnd(1-2*lambda,1/lambda4,T,1);
Samples_C(5,:) = gamrnd(1-2*lambda,1/lambda5,T,1);
kT_Model= Samples_C(2,:)/2;
nuL_Model = (Samples_C(3,:)./kT_Model)/(2*sqrt(2));
EL_Model = Samples_C(1
,:) - 4*(nuL_Model.^2).*kT_Model;
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
cdf_C4 = @(c4) gamcdf(c4,1-2*lambda,1/lambda4); % Gamma cumulative density function of C1
cdf_C5 = @(c5) gamcdf(c5,1-2*lambda,1/lambda5); % Gamma cumulative density function of C2

N = 100;
x4_min = max(0,mean(C4_data)-8*std(C4_data));
x4_max = mean(C4_data)+8*std(C4_data);
x5_min = max(0,mean(C5_data)-8*std(C5_data));
x5_max = mean(C5_data)+8*std(C5_data);

% Plot pdf of C4
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

% Plot pdf of C5
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

% Plot cdf of C5
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
