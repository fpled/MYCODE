function carte_stat_patch(S,Spart,surface,U,w,Ufinal,Uref,varargin)
% function carte_stat_patch(S,Spart,U,w,Ufinal,Uref,varargin)
% Pour afficher les moyennes, les variances,les ecarts types, les indices de sensibilites

chemin = getcharin('chemin',varargin,0)

%
figure(75)
plot_sol(S,expect(Ufinal),'edgecolor','none')
colorbar


figure(76)
plot_sol(S,mavar(Ufinal),'edgecolor','none')
colorbar
hold on
for k=2:length(surface)
    plot(surface{k})
    hold on
end
%
if(Uref~=0)
    figure(77)
    plot_sol(S,expect(Uref),'edgecolor','none')
    colorbar
    figure(78)
    plot_sol(S,mavar(Uref),'edgecolor','none')
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
colorbar
%
figure(80)
plot_sol(S,mavar(U),'edgecolor','none')
colorbar
%
wr = cell(1,length(w));
for k=2:length(w)
    figure(90+k)
    %set(gca,'fontsize',16);
    %title(['w' num2str(k-1)]);
    plot_sol(Spart{k},expect(w{k}),'edgecolor','none')
    colorbar
    figure(135+k)
    plot_sol(Spart{k},mavar(w{k}),'edgecolor','none')
    colorbar
    figure(800+k)
    plot_sol(S,sobol_indices(Ufinal,k-1),'edgecolor','none')
    hold on
    for kk=2:length(w)
        plot(surface{kk})
        hold on
    end
    colorbar
end


if(chemin~=0)
    figure(75)
    myprint(chemin,'plot_Ufinal_moyenne' ,{'jpeg','epsc2'});
    figure(76)
    myprint(chemin,'plot_Ufinal_variance' ,{'jpeg','epsc2'});
    if(Uref~=0)
        figure(77)
        myprint(chemin,'plot_Uref_moyenne' ,{'jpeg','epsc2'});
        figure(78)
        myprint(chemin,'plot_Uref_variance' ,{'jpeg','epsc2'});
    end
    figure(79)
    myprint(chemin,'plot_U_moyenne' ,{'jpeg','epsc2'});
    figure(80)
    myprint(chemin,'plot_U_variance' ,{'jpeg','epsc2'});
    
    for k=2:length(w)
        figure(90+k)
        myprint(chemin,['plot_w_moyenne' num2str(k-1)] ,{'jpeg','epsc2'});
        figure(135+k)
        myprint(chemin,['plot_w_variance' num2str(k-1)] ,{'jpeg','epsc2'});
        
        figure(800+k)
        myprint(chemin,['plot_indice_sensib' num2str(k-1)] ,{'jpeg','epsc2'});
    end
end

