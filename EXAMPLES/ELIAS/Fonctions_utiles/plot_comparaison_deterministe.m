function [error,Udeter] = plot_comparaison_deterministe(S,Spart,U,P,kappa,varargin)
%function error = plot_comparaison(S,Uref,U2,varargin)
%chemin = getcharin('chemin',varargin,0)
%nbreal = getcharin('nbreal',varargin,3);

chemin = getcharin('chemin',varargin,0);
nbreal = getcharin('nbreal',varargin,3);

error = zeros(1,nbreal);
for nb=1:nbreal
    [Ur,r] = random(U);
    Ur = unfreevector(S,Ur);
    Udeter = calcsoldeter_patch(S,Spart,P,kappa,r);
    
    figure(25)
    subplot(1,nbreal,nb)
    plot_sol(S,Ur,'edgecolor','none')
    %title('Ufinal');
    %cax=caxis;
    %caxis(cax);
    colorbar
    
    figure(26)
    subplot(1,nbreal,nb)
    plot_sol(S,Udeter,'edgecolor','none','numnode')
    %title('Udeter');
    %caxis(cax);
    colorbar
    error(nb) = norm(Ur-Udeter)/norm(Udeter)
end
if(chemin~=0)
    figure(25)
    myprint(chemin,'plot_sol_Uref_4patch' ,{'jpeg','eps2'});
    %myprint(chemin,'plot_sol_Uref' ,{'jpeg','eps2'});
    figure(26)
    myprint(chemin,'plot_sol_Udeter_4patch' ,{'jpeg','eps2'});
    %myprint(chemin,'plot_sol_Udeter' ,{'jpeg','eps2'});
end

end
