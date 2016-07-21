function [resultU,resultw,resultUfinal] = ploterrorUwU(w,U,Ufinal,dd,varargin)
%function [resultU,resultw,resultUfinal] = ploterrorUwU(w,U,Ufinal,dd,varargin)
%chemin = getcharin('chemin',varargin,0)
chemin = getcharin('chemin',varargin,0)
errref = getcharin('errref',varargin,0)
wpc = cell(1,length(w)+1);
resultw = cell(1,length(w)+1);

clr = {'k',':b','--r','--m','--c','--g'}

mk = ['+','o','<','*','x','s','d','^','>','p','h','v'];
for k=2:length(w)
    [wpc{k},resultw{k}] = spectral_decomposition(w{k},'display',true,'reference',w{k});
end
[Upc,resultU] = spectral_decomposition(U,'display',true,'reference',U);
[Ufinalpc,resultUfinal] = spectral_decomposition(Ufinal,'display',true,'reference',Ufinal);

figure(7546+dd)
if errref==0
    for k=2:length(w)
        set(gca,'fontsize',16)
        semilogy(resultw{k}.error,clr{k+1},'linewidth',k-1)%,'marker',mk(k))
        hold on;
    end
    semilogy(resultU.error,':b','linewidth',3)
    hold on;
    semilogy(resultUfinal.error,'k','linewidth',3)
    hold on;
    set(gca,'fontsize',16)
    xlabel('Rank m')
    set(gca,'fontsize',16)
    ylabel('Error')
    hold on;
%xlim([0 25])
%ylim([10^(-15) 1])
% elseif(errref==1)
%  for k=2:length(w)
% set(gca,'fontsize',16)
% semilogy(norm(w{k})\resultw{k}.error,'linewidth',k-1,'color',rand(1,3))%,'marker',mk(k))
% hold on;
% end
% semilogy(norm(U)\resultU.error,':g','linewidth',3)
% hold on;
% semilogy(norm(Ufinal)\resultUfinal.error,'-.k','linewidth',3)
% hold on;
% set(gca,'fontsize',16)
% xlabel('Rank m')
% set(gca,'fontsize',16)
% ylabel('Error')
% hold on;   
end
if length(w)==5
    legend('w_1','w_2','w_3','w_4','U','u_{ref}');
end
if length(w)==2
    legend('w','U','u_{ref}');
end

if length(w)==3
    legend('w_1','w_2','U','u_{ref}');
end

if length(w)==8
    legend('w_1','w_2','w_3','w_4','w_5','w_6','w_7','U','u_{ref}');
end

if (chemin~=0)
    if length(w)==5
        figure(7546+dd)
        myprint(chemin,['plot_error_superp_4patch_dirichlet' num2str(dd)],{'jpeg','epsc2'});
    end
    if length(w)==2
        figure(7546+dd)
        myprint(chemin,['plot_error_superp_dirichlet' num2str(dd)],{'jpeg','epsc2'});
    end

    if length(w)==8
        figure(7546+dd)
        myprint(chemin,['plot_error_spectrale_diffusion_dirichlet' num2str(dd)],{'jpeg','epsc2'});
    end

    if length(w)==3
        figure(7546+dd)
        myprint(chemin,['plot_error_approximation_dirichlet' num2str(dd)],{'jpeg','epsc2'});
    end

end
