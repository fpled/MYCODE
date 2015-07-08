function carte_stat_patch_nonconf(S,Spart,Spartf,surface,U,w,Uref,Uzref,wref,P,maxv,maxe,varargin)
% function carte_stat_patch_nonconf(S,Spart,Spartf,surface,U,w,Uref,Uzref,wref,P,maxv,maxe,varargin)
% Pour afficher les moyennes, les variances,les ecarts types, lesindices de sensibilites

chemin = getcharin('chemin',varargin,0)

%
figure(75)

plot_sol(Spart{1},P{1}*expect(U),'edgecolor','none')
hold on
for k=2:length(w)
    plot_sol(Spartf{k},expect(w{k}),'edgecolor','none')
end
colorbar
caxexpect=caxis

% figure(775)
% plot_sol(S,expect(U),'edgecolor','none')
% caxis(caxexpect);
% colorbar

for k=2:length(surface)
    plot(surface{k})
    hold on
end

figure(76)
plot_sol(Spart{1},P{1}*mavar(U),'edgecolor','none')
hold on
for k=2:length(w)
    plot_sol(Spartf{k},mavar(w{k}),'edgecolor','none')
end
colorbar
caxvar = caxis;
%
for k=2:length(surface)
    plot(surface{k},'w')
    hold on
end
%
if(isa(Uzref,'PCMATRIX'))
    figure(77)
    plot_sol(Spart{1},expect(Uref),'edgecolor','none')
    hold on
    for k=2:length(w)
        plot_sol(Spartf{k},expect(wref{k}),'edgecolor','none')
    end
    colorbar
    for k=2:length(surface)
        plot(surface{k})
        hold on
    end
    
    figure(78)
    plot_sol(Spart{1},mavar(Uref),'edgecolor','none')
    hold on
    for k=2:length(w)
        plot_sol(Spartf{k},mavar(wref{k}),'edgecolor','none')
    end
    colorbar
    hold on
    for k=2:length(surface)
        plot(surface{k})
        hold on
    end
end
%
figure(79)
plot_sol(S,expect(U),'edgecolor','none')
caxis(caxexpect);
colorbar
%
figure(80)
plot_sol(S,mavar(U),'edgecolor','none')
caxis(caxvar);
colorbar
%
wr = cell(1,length(w));
for k=2:length(w)
    figure(90+k)
    %set(gca,'fontsize',16);
    %title(['w' num2str(k-1)]);
    plot_sol(Spartf{k},expect(w{k}),'edgecolor','none')
    colorbar
    figure(135+k)
    plot_sol(Spartf{k},mavar(w{k}),'edgecolor','none')
    colorbar
end

for sb=2:length(w)
    figure(800+sb)
    plot_sol(Spart{1},sobol_indices(P{1}*U,sb-1),'edgecolor','none')
    hold on
    for k=2:length(w)
        plot_sol(Spartf{k},sobol_indices(w{k},sb-1),'edgecolor','none')
        plot(surface{sb})
        hold on
    end
    colorbar
end


for sb=2:length(w)
    figure(11000+sb)
    plot_sol(Spart{1},var_expect_cond(P{1}*U,sb-1)/maxv,'edgecolor','none')
    hold on
    for k=2:length(w)
        plot_sol(Spartf{k},var_expect_cond(w{k},sb-1)/maxv,'edgecolor','none')
        plot(surface{k})
        hold on
    end
    colorbar
end

for sb=2:length(w)
    figure(14560+sb)
    plot_sol(Spart{1},var_expect_cond(P{1}*U,sb-1)/maxe(sb-1),'edgecolor','none')
    hold on
    for k=2:length(w)
        plot_sol(Spartf{k},var_expect_cond(w{k},sb-1)/maxe(sb-1),'edgecolor','none')
        plot(surface{k})
        hold on
    end
    colorbar
end


if(chemin~=0)
    
    figure(75)
    myprint(chemin,'plot_Ufinal_moyenne_4patchdiffusion' ,{'jpeg','epsc2'});
    figure(76)
    myprint(chemin,'plot_Ufinal_variance_4patchdiffusion' ,{'jpeg','epsc2'});
    if(isa(Uzref,'PCMATRIX'))
        figure(77)
        myprint(chemin,'plot_Uref_moyenne_4patchdiffusion' ,{'jpeg','epsc2'});
        figure(78)
        myprint(chemin,'plot_Uref_variance_4patchdiffusion' ,{'jpeg','epsc2'});
    end
    figure(79)
    myprint(chemin,'plot_U_moyenne_4patchdiffusion' ,{'jpeg','epsc2'});
    figure(80)
    myprint(chemin,'plot_U_variance_4patchdiffusion' ,{'jpeg','epsc2'});
    
    for k=2:length(w)
        figure(90+k)
        myprint(chemin,['plot_w_moyenne_4patchdiffusion' num2str(k-1)] ,{'jpeg','epsc2'});
        figure(135+k)
        myprint(chemin,['plot_w_variance_4patchdiffusion' num2str(k-1)] ,{'jpeg','epsc2'});
        
        figure(800+k)
        myprint(chemin,['plot_indice_sensib_4patchdiffusion' num2str(k-1)] ,{'jpeg','epsc2'});
        figure(11000+k)
        myprint(chemin,['plot_indice_sensib_maxvar_4patchdiffusion' num2str(k-1)] ,{'jpeg','epsc2'});
        figure(14560+k)
        myprint(chemin,['plot_indice_sensib_maxexpect_4patchdiffusion' num2str(k-1)] ,{'jpeg','epsc2'});
    end
end

