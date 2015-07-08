function error = plot_comparaison(S,Uref,U2,varargin)
%function error = plot_comparaison(S,Uref,U2,varargin)
%chemin = getcharin('chemin',varargin,0)
%nbreal = getcharin('nbreal',varargin,3);

chemin = getcharin('chemin',varargin,0);
nbreal = getcharin('nbreal',varargin,3);

error = zeros(1,nbreal);
for nb=1:nbreal
    [Ur,r] = random(Uref);
    Ur = unfreevector(S,Ur);
    Ur2 = randomeval(U2,r);
    Ur2 = unfreevector(S,Ur2);
    
    figure(35)
    subplot(nbreal,1,nb)
    plot_sol(S,Ur,'edgecolor','none')
    %title('Uref');
    %cax=caxis;
    %caxis(cax);
    colorbar
    
    figure(36)
    subplot(nbreal,1,nb)
    plot_sol(S,Ur2,'edgecolor','none','numnode')
    %title('Ufinal');
    %caxis(cax);
    colorbar
    error(nb) = norm(Ur-Ur2)/norm(Ur)
end
if(chemin~=0)
    figure(35)
    myprint(chemin,'plot_sol_Uref_4patch' ,{'jpeg','eps2'});
    %myprint(chemin,'plot_sol_Uref' ,{'jpeg','eps2'});
    figure(36)
    myprint(chemin,'plot_sol_Ufinal2_4patch' ,{'jpeg','eps2'}');
    %myprint(chemin,'plot_sol_Ufinal2' ,{'jpeg','eps2'});
end

end
