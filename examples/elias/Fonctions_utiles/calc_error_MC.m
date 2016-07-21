function errL2 = calc_error_MC(y1,y2,Sr,nbreal)
% function errL2 = calc_error_MC(y1,y2,Sr,ls,nbreal)
% errL2 = norm(y2-y1)/norm(y1) avec MC en eloignant un peu de la frontiere
% Sr : Model 
% ls : level-set deterministe
%   
rep = find(getvalue(Sr.ls{1})<-1e-4);
den=0;
num=0;
%for ixi=1:nbreal
%pourcentage(ixi,nbreal)  
num = num + 1/nbreal*norm(y2(rep)-y1(rep))^2;
den = den + 1/nbreal*norm(y1(rep))^2;
% figure(4)
% clf
% plot_sol(Sr,umr,'sedgecolor','none')
% plot(POINT(x(rep,:)),'*')
%pause
%end
errL2 = sqrt(num/den);

return


