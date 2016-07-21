function [resultU,resultw,resultUfinal] = ploterrorUwU_Uref(w,U,Ufinal,dd,varargin)
%function [resultU,resultw,resultUfinal] = ploterrorUwU(w,U,Ufinal,dd,varargin)
%chemin = getcharin('chemin',varargin,0)
chemin = getcharin('chemin',varargin,0)
wpc = cell(1,length(w)+1);
resultw = cell(1,length(w)+1);

if length(w)==5
    clr = {'k',':b','--r','--m','--c','--g'}
elseif length(w)==2
    clr = {'--r','--r','--r'}%,'--c','--g'}

end

mk = ['+','o','<','*','x','s','d','^','>','p','h','v'];
for k=2:length(w)
    [wpc{k},resultw{k}] = spectral_decomposition(w{k},'display',true,'reference',w{k});
end
[Upc,resultU] = spectral_decomposition(U,'display',true,'reference',U);
[Ufinalpc,resultUfinal] = spectral_decomposition(Ufinal,'display',true,'reference',Ufinal);

figure(1234+dd)

for k=2:length(w)
    set(gca,'fontsize',16)
    semilogy(length(w)*norm(w{k})*(resultw{k}.error)/norm(Ufinal),clr{k+1},'linewidth',k-1)%,'marker',mk(k))
    hold on;
end
semilogy(length(w)*norm(U)*resultU.error/norm(Ufinal),':b','linewidth',3)
hold on;
semilogy(length(w)*resultUfinal.error,'k','linewidth',3)
hold on;
set(gca,'fontsize',16)
xlabel('Rank m')
set(gca,'fontsize',16)
ylabel('Error')
hold on;
xlim([0 25])
ylim([10^(-15) 1])

if length(w)==5
    legend('w_1','w_2','w_3','w_4','U','u_{ref}');
end
if length(w)==2
    legend('w','U','u_{ref}');
end

if (chemin~=0)
    if length(w)==5
        figure(1234+dd)
        myprint(chemin,['plot_error_superp_4patch_refUfinal_dirichlet' num2str(dd)],{'jpeg','epsc2'});
    end
    if length(w)==2
        figure(1234+dd)
        myprint(chemin,['plot_error_superp_refUfinal_dirichlet' num2str(dd)],{'jpeg','epsc2'});
    end
end
