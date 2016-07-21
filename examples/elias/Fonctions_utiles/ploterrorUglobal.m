function [resultU] = ploterrorUglobal(w,U,Ufinal,dd,varargin)
%function [resultU] = ploterrorUwUglobal(w,U,Ufinal,dd,varargin)
%chemin = getcharin('chemin',varargin,0)
chemin = getcharin('chemin',varargin,0)


%clr = {'k',':b','--r','--m','--c','--g'}

%mk = ['+','o','<','*','x','s','d','^','>','p','h','v'];

[Upc,resultU] = spectral_decomposition(U,'display',true,'reference',U);


figure(746)
semilogy(resultU.error,'linewidth',dd+1)
hold on
set(gca,'fontsize',16)
xlabel('Rank m')
set(gca,'fontsize',16)
ylabel('Error')
hold on;
xlim([0 35])
ylim([10^(-15) 1])
legend('L=0.5','L=1','L=1.5');

figure(587)
semilogy(length(w)*norm(U)*resultU.error/norm(Ufinal),':b','linewidth',dd+1)
hold on

set(gca,'fontsize',16)
xlabel('Rank m')
set(gca,'fontsize',16)
ylabel('Error')
hold on;
xlim([0 35])
ylim([10^(-15) 1])
legend('L=0.5','L=1','L=1.5');


if (chemin~=0)

    figure(746+dd)
    myprint(chemin,['plot_spectral_decomp_U_L_4patches' num2str(dd)],{'jpeg','epsc2'});

    figure(587+dd)
    myprint(chemin,['plot_spectral_decomp_U_L_4patches_refUfinal' num2str(dd)],{'jpeg','epsc2'});
end
end
