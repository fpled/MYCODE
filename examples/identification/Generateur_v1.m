clc
clear all
close all

n=4;    % dimension de la matrice

% information disponible : la moyenne de [A]
R=randn(n,n);
A_moy=R*R';
U_A_moy=chol(A_moy);

% parametrage
gamma_para=2;      % ordre du moment
delta_para=0.31;
lambda_para=(1-delta_para^2*(n-1)+(trace(A_moy))^2/trace(A_moy^2))/(2*delta_para^2);    
lambda_mini=max(0,gamma_para/n+(1-n)/2+1-1/n);
lambda_para=max(lambda_para,lambda_mini);
lambda_para_integer=round(lambda_para);

mean_Gauss = 0;     
sigma_Gauss = (2/(n-1+2*lambda_para))^0.5;   % standard_deviation
mA = n-1+2*lambda_para_integer;    % degrees of freedom of Wishart distribution
sigma_Wishart = (1/mA)*A_moy;    % covariance matrix

Nom_Sim_power=1:0.2:5;
Nom_Sim=round(10.^Nom_Sim_power);
Conv_moy_A_Wishart=zeros(length(Nom_Sim),1);
Conv_moy_A_X=zeros(length(Nom_Sim),1);
Conv_moy_A_U=zeros(length(Nom_Sim),1);
Conv_moy_A=zeros(length(Nom_Sim),1);
delta_A_Wishart=zeros(length(Nom_Sim),1);
delta_A_X=zeros(length(Nom_Sim),1);
delta_A_U=zeros(length(Nom_Sim),1);
delta_A=zeros(length(Nom_Sim),1);
EM_Ainv2_Wishart=zeros(length(Nom_Sim),1);
EM_Ainv2_X=zeros(length(Nom_Sim),1);
EM_Ainv2_U=zeros(length(Nom_Sim),1);
EM_Ainv2=zeros(length(Nom_Sim),1);
    

for sim=1:length(Nom_Sim)
    
    NS=Nom_Sim(sim);
    
    L_t_ll = zeros(1,n*NS);
    L_t_llp = zeros(1,n*(n-1)/2*NS);
    
    A_rand_Wishart = zeros(n,n,NS);
    A_rand_X = zeros(n,n,NS);
    A_rand_U = zeros(n,n,NS);
    A_rand = zeros(n,n,NS);
    
    rng('default');    % r¨¦initialisation des NS premiers MC simulations
    
    for k=1:NS
        
        % M¨¦thode 1 : Wishart distribution
        A_rand_Wishart(:,:,k) = wishrnd(sigma_Wishart,mA);
        
        % M¨¦thode 2 : d¨¦composition en loi normale centr¨¦e
        X_rand = mvnrnd(zeros(n,1),A_moy/mA,mA);
        A_rand_X(:,:,k) = X_rand'*X_rand;
        
        % M¨¦thode 3 : d¨¦composition en loi normale centr¨¦e r¨¦duite
        U_rand = randn(n,mA);
        P_rand_sum=zeros(n,n);
        for m=1:mA
            P_rand = U_A_moy' * U_rand(:,m) * (U_A_moy' * U_rand(:,m))';
            P_rand_sum = P_rand_sum + P_rand;
        end
        A_rand_U(:,:,k) = P_rand_sum/mA;
       
        % M¨¦thode 4 : cas g¨¦n¨¦rale (lambda not an integer)
        L_t_llp(1,(k-1)*n*(n-1)/2+1:(k-1)*n*(n-1)/2+n*(n-1)/2) = sigma_Gauss*randn(n*(n-1)/2,1) + mean_Gauss;
        
        L_rand = triu(ones(n,n),1);
        L_rand(L_rand==1)=L_t_llp(1,(k-1)*n*(n-1)/2+1:(k-1)*n*(n-1)/2+n*(n-1)/2)*2^(-0.5);
        
        B=1;
        for i=1:n
            A=(n-i+2*lambda_para)/2;
            Y_l = gamrnd(A,B,1,1);
            L_t_ll(1,(k-1)*n+i) = sigma_Gauss*(Y_l)^0.5;
        end
        
        L_rand_diag = diag(L_t_ll(1,(k-1)*n+1:(k-1)*n+n));
        L_rand = L_rand + L_rand_diag;
        
        G_rand=L_rand'*L_rand;
        
        A_rand(:,:,k) = U_A_moy'*G_rand*U_A_moy;
        
    end
    
    if NS == max(Nom_Sim)
        [pdf_llp,xi_llp] = ksdensity(L_t_llp);
        pdf_llp_analy = (2*pi)^(-0.5)*sigma_Gauss^(-1)*exp(-xi_llp.^2/(2*sigma_Gauss^2));
        figure
        % plot(xi_llp,pdf_llp,xi_llp,pdf_llp_analy)
        semilogy(xi_llp,pdf_llp,xi_llp,pdf_llp_analy)
        tx = xlabel('x'); ty = ylabel('log(p_X)');
        set(tx,'Fontsize',10); set(ty,'Fontsize',10);
        tt=title(['P.D.F pour les termes non diagonaux']) ;
        set(tt,'FontName','Bell MT','Fontsize',12); set(gca,'FontName','Times','FontSize',10) ;
        
        
        for i=1:n
            ind=i:n:n*NS;
            [pdf_ll,xi_ll] = ksdensity(L_t_ll(ind));
            
            pdf_ll_analy = (n-1+2*lambda_para)*xi_ll.*gampdf((n-1+2*lambda_para)*xi_ll.^2/2, (n-i+2*lambda_para)/2,1);
            
            figure
            %     plot(xi_ll,pdf_ll,xi_ll,pdf_ll_analy)
            semilogy(xi_ll,pdf_ll,xi_ll,pdf_ll_analy)
            tx = xlabel('x'); ty = ylabel('log(p_X)');
            set(tx,'Fontsize',10); set(ty,'Fontsize',10);
            tt=title(['P.D.F pour le ',num2str(i),'th terme diagonal']) ;
            set(tt,'FontName','Bell MT','Fontsize',12); set(gca,'FontName','Times','FontSize',10) ;
        end    
    end

    % Etude de convergence de l'esp¨¦rance math¨¦matique vers la moyenne
    EM_A_Wishart = mean(A_rand_Wishart,3);
    EM_A_X = mean(A_rand_X,3);
    EM_A_U = mean(A_rand_U,3);
    EM_A = mean(A_rand,3);
    
    Conv_moy_A_Wishart(sim,1) = (norm(EM_A_Wishart - A_moy,'fro'))^2/(norm(A_moy,'fro'))^2;
    Conv_moy_A_X(sim,1) = (norm(EM_A_X - A_moy,'fro'))^2/(norm(A_moy,'fro'))^2;
    Conv_moy_A_U(sim,1) = (norm(EM_A_U - A_moy,'fro'))^2/(norm(A_moy,'fro'))^2;
    Conv_moy_A(sim,1) = (norm(EM_A - A_moy,'fro'))^2/(norm(A_moy,'fro'))^2;
    
    % Etude de convergence du coefficient de variation et du moment d'ordre
    % deux de A^-1
    AAmoy2_Wishart=zeros(NS,1);
    AAmoy2_X=zeros(NS,1);
    AAmoy2_U=zeros(NS,1);
    AAmoy2=zeros(NS,1);
    Ainv2_Wishart=zeros(NS,1);
    Ainv2_X=zeros(NS,1);
    Ainv2_U=zeros(NS,1);
    Ainv2=zeros(NS,1);
    for k=1:NS
        AAmoy2_Wishart(k,1) = (norm(A_rand_Wishart(:,:,k) - A_moy,'fro'))^2;
        AAmoy2_X(k,1) = (norm(A_rand_X(:,:,k) - A_moy,'fro'))^2;
        AAmoy2_U(k,1) = (norm(A_rand_U(:,:,k) - A_moy,'fro'))^2;
        AAmoy2(k,1) = (norm(A_rand(:,:,k) - A_moy,'fro'))^2;
        Ainv2_Wishart(k,1) = (norm(inv(A_rand_Wishart(:,:,k)),'fro'))^2;
        Ainv2_X(k,1) = (norm(inv(A_rand_X(:,:,k)),'fro'))^2;
        Ainv2_U(k,1) = (norm(inv(A_rand_U(:,:,k)),'fro'))^2;
        Ainv2(k,1) = (norm(inv(A_rand(:,:,k)),'fro'))^2;
    end
    
    EM_AAmoy2_Wishart = mean(AAmoy2_Wishart);
    EM_AAmoy2_X = mean(AAmoy2_X);
    EM_AAmoy2_U = mean(AAmoy2_U);
    EM_AAmoy2 = mean(AAmoy2);
    delta_A_Wishart(sim,1)=(EM_AAmoy2_Wishart/(norm(A_moy,'fro'))^2)^0.5;
    delta_A_X(sim,1)=(EM_AAmoy2_X/(norm(A_moy,'fro'))^2)^0.5;
    delta_A_U(sim,1)=(EM_AAmoy2_U/(norm(A_moy,'fro'))^2)^0.5;
    delta_A(sim,1)=(EM_AAmoy2/(norm(A_moy,'fro'))^2)^0.5;
    
    EM_Ainv2_Wishart (sim,1)=mean(Ainv2_Wishart);
    EM_Ainv2_X (sim,1)=mean(Ainv2_X);
    EM_Ainv2_U (sim,1)=mean(Ainv2_U);
    EM_Ainv2 (sim,1)=mean(Ainv2);
    
    
%     loglog(Nom_Sim,1./(Nom_Sim).^0.5,Nom_Sim,Conv_moy_A_Wishart,Nom_Sim,Conv_moy_A_X,Nom_Sim,Conv_moy_A_U,Nom_Sim,Conv_moy_A)
%     semilogx(Nom_Sim,delta_A_Wishart,Nom_Sim,delta_A_X,Nom_Sim,delta_A_U,Nom_Sim,delta_A)
%     semilogx(Nom_Sim,EM_Ainv2_Wishart,Nom_Sim,EM_Ainv2_X,Nom_Sim,EM_Ainv2_U,Nom_Sim,EM_Ainv2)
    
%     pause(0.2)
    
end
figure
loglog(Nom_Sim,1./(Nom_Sim).^0.5,Nom_Sim,Conv_moy_A_Wishart,Nom_Sim,Conv_moy_A_X,Nom_Sim,Conv_moy_A_U,Nom_Sim,Conv_moy_A)
tx = xlabel('Nombre Simulation (log)'); ty = ylabel('Residu (log)');
set(tx,'Fontsize',10); set(ty,'Fontsize',10);
tt=title(['Etude Convergence : esp¨¦rance math¨¦matique vers moyenne']) ;
set(tt,'FontName','Bell MT','Fontsize',12); set(gca,'FontName','Times','FontSize',10) ;
figure
semilogx(Nom_Sim,delta_A_Wishart,Nom_Sim,delta_A_X,Nom_Sim,delta_A_U,Nom_Sim,delta_A)
tx = xlabel('Nombre Simulation (log)'); ty = ylabel('\delta_X');
set(tx,'Fontsize',10); set(ty,'Fontsize',10);
tt=title(['Etude Convergence : coefficient de variation']) ;
set(tt,'FontName','Bell MT','Fontsize',12); set(gca,'FontName','Times','FontSize',10) ;
figure
semilogx(Nom_Sim,EM_Ainv2_Wishart,Nom_Sim,EM_Ainv2_X,Nom_Sim,EM_Ainv2_U,Nom_Sim,EM_Ainv2)
tx = xlabel('Nombre Simulation (log)'); ty = ylabel('E(\mid\mid A^{-1} \mid\mid^2)');
set(tx,'Fontsize',10); set(ty,'Fontsize',10);
tt=title(['Etude Convergence : moment d''ordre 2 de A^{-1}']) ;
set(tt,'FontName','Bell MT','Fontsize',12); set(gca,'FontName','Times','FontSize',10) ;

