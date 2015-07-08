function plot_patch(S,Spart,U,w,Ufinal,Uref,varargin)
%function plot_patch(S,Spart,U,w,Ufinal,Uref,varargin)
% ('m','n') m et n sont les param�tres de subplot
% m*n est le nombre de r�alisations
% ('grad',1) pour afficher les gradients
grad = getcharin('grad',varargin,0);
m =  getcharin('m',varargin,3);
n =  getcharin('n',varargin,1);
chemin = getcharin('chemin',varargin,0)
for nb=1:(m*n)
    [ur,r] = random(U);
    figure(45)
    subplot(m,n,nb)
    %title('U')
    plot_sol(S,ur,'edgecolor','none')
    colorbar
    
    wr = cell(1,length(w));
    for k=2:length(w)
        wr{k} = randomeval(w{k},r);
        figure(46+k)
        %set(gca,'fontsize',16);
        %title(['w' num2str(k-1)]);
        subplot(m,n,nb)
        plot_sol(Spart{k},wr{k},'edgecolor','none')
        colorbar
    end
    
    ufr = randomeval(Ufinal,r);
    ufr = unfreevector(S,ufr);
    figure(1000)
    subplot(m,n,nb)
    %title('Ufinal,S')
    plot_sol(S,ufr,'edgecolor','none')
    colorbar
    
    Urefr = randomeval(Uref,r);
    Urefr = unfreevector(S,Urefr);
    
    figure(1035)
    subplot(m,n,nb)
    plot_sol(S,Urefr,'edgecolor','none')
    %title('Uref');
    %cax=caxis;
    %caxis(cax);
    colorbar
    
    
end

if(chemin~=0)
    figure(45)
    myprint(chemin,'plot_sol_U2_4patch' ,{'jpeg','epsc2'});
    %myprint(chemin,'plot_sol_U2' ,{'jpeg','eps2'});
    for k=2:length(w)
        figure(46+k)
        myprint(chemin,['plot_sol_w_4patch' num2str(k-1)] ,{'jpeg','epsc2'});
        %myprint(chemin,['plot_sol_w' num2str(k-1)] ,{'jpeg','eps2'});
    end
    
    figure(1000)
    myprint(chemin,'plot_sol_Ufinal_4patch' ,{'jpeg','epsc2'});
    %myprint(chemin,'plot_sol_Ufinal' ,{'jpeg','eps2'});
    
    
    figure(1035)
    myprint(chemin,'plot_sol_Ureference_4patch' ,{'jpeg','epsc2'});
    %myprint(chemin,'plot_sol_Ureference' ,{'jpeg','eps2'});
end



if(grad)
    for nb=1:m*n
        figure(4500)
        subplot(m,n,nb)
        %title('gradx U,S')
        plot_sol(S,ur,'epsilon',1,'edgecolor','none')
        colorbar
        
        
        for k=2:length(w)
            figure(4501+k)
            %title(['gradx w' num2str(k-1)]);
            subplot(m,n,nb)
            plot_sol(Spart{k},wr{k},'epsilon',1,'edgecolor','none')
            colorbar
        end
        
        
        figure(4600)
        subplot(m,n,nb)
        %title('gradx Ufinal,S')
        plot_sol(S,ufr,'epsilon',1,'edgecolor','none')
        colorbar
        
        figure(3500)
        subplot(m,n,nb)
        %title('grady U,S')
        plot_sol(S,ur,'epsilon',2,'edgecolor','none')
        colorbar
        
        for k=2:length(w)
            figure(3501+k)
            %title(['grady w' num2str(k-1)]);
            subplot(m,n,nb)
            plot_sol(Spart{k},wr{k},'epsilon',2,'edgecolor','none')
            colorbar
        end
        
        figure(3600)
        subplot(m,n,nb)
        %title('grady Ufinal,S')
        plot_sol(S,ufr,'epsilon',2,'edgecolor','none')
        colorbar
    end
    
    
    if(chemin~=0)
        
        figure(4500)
        myprint(chemin,'plot_sol_gradxU_4patch' ,{'jpeg','epsc2'});
        %myprint(chemin,'plot_sol_gradxU' ,{'jpeg','eps2'});
        for k=2:length(w)
            figure(4501+k)
            myprint(chemin,['plot_sol_gradxw_4patch' num2str(k-1)] ,{'jpeg','epsc2'});
            %myprint(chemin,['plot_sol_gradxw' num2str(k-1)] ,{'jpeg','eps2'});
        end
        figure(4600)
        myprint(chemin,'plot_sol_gradxUfinal_4patch' ,{'jpeg','epsc2'});
        %myprint(chemin,'plot_sol_gradxUfinal' ,{'jpeg','eps2'});
        figure(3500)
        myprint(chemin,'plot_sol_gradyU_4patch' ,{'jpeg','epsc2'});
        %myprint(chemin,'plot_sol_gradyU' ,{'jpeg','eps2'});
        for k=2:length(w)
            figure(3501+k)
            myprint(chemin,['plot_sol_gradyw_4patch' num2str(k-1)] ,{'jpeg','epsc2'});
            %myprint(chemin,['plot_sol_gradyw' num2str(k-1)] ,{'jpeg','eps2'});
        end
        
        figure(3600)
        myprint(chemin,'plot_sol_gradyUfinal_4patch' ,{'jpeg','epsc2'});
        %myprint(chemin,'plot_sol_gradyUfinal' ,{'jpeg','eps2'});
        
    end
    
end

