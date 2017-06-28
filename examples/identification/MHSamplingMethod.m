clc
close all
clear all

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
proppdf = @(c) mvnpdf([c(1) c(2) c(3)],Mu_mean,Sigma_mean);
proprnd = @(c) mvnrnd(Mu_mean,Sigma_mean);
pdf = @(c) (c(1)*c(2)-c(3)^2)^(-lambda)*exp(-lambda1*c(1)-lambda2*c(2)-lambda3*c(3));
nsamples = 5000;
start = [0 0 0];
smpl = mhsample(start,nsamples,'pdf',pdf,'proppdf',proppdf, 'proprnd',proprnd);
