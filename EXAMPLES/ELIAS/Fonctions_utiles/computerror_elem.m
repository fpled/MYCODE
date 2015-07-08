function [ke,ke1] = computerror_elem(eleme,xnodee,lse,ue,uex,xi,varargin)
% function [ke,ke1] = computerror_elem(eleme,xnodee,lse,ue,uex,xi,varargin)
% eleme : element fini triangle
% xnode : coordonn�es des noeuds du triangle
%
% varargin : "method" --> 'penal' ou 'nitsche' ou 'carac' ou 'pas de
% methode' pour rien faire


calc_cut = getcharin('calc_cut',varargin);
method = getcharin('method',varargin);
switch method
    case {'penal','pas de methode'}
        
        if all(lse<=0)
            gauss = calc_gauss(eleme,8);
        elseif ~all(lse>=0)
            [gauss,gaussout] = calc_lssubgauss(eleme,lse,8);
        end
        
        if ~calc_cut
            if any(lse>0) && any(lse<=0)
                gauss.w(:)=0;
            end
        end
        N = calc_N(eleme,xnodee,gauss.coord);
        x = calc_x(eleme,xnodee,gauss.coord);
        deltau = N*ue - uex(x,xi);
        detJ = calc_detJ(eleme,xnodee,gauss.coord);
        ke = double(sum(gauss.w*deltau.*deltau*abs(detJ),4));
        ke1 = double(sum(gauss.w*uex(x,xi).*uex(x,xi)*abs(detJ),4));
        
    case 'nitsche'
        
        if all(lse<=0)
            gauss = calc_gauss(eleme,2);
        elseif ~all(lse>=0)
            [gauss,gaussout] = calc_lssubgauss(eleme,lse,2);
        end
        
        if ~calc_cut
            if any(lse>0) && any(lse<=0)
                gauss.w(:)=0;
            end
        end
        
        N = calc_N(eleme,xnodee,gauss.coord);
        x = calc_x(eleme,xnodee,gauss.coord);
        deltau = N*ue - uex(x,xi);
        detJ = calc_detJ(eleme,xnodee,gauss.coord);
        ke = double(sum(gauss.w*deltau.*deltau*abs(detJ),4));
        ke1 = double(sum(gauss.w*uex(x,xi).*uex(x,xi)*abs(detJ),4));
    case 'carac'
        
        if all(lse<=0)
            gauss = calc_gauss(eleme,8);
        elseif ~all(lse>=0)
            [gauss,gaussout] = calc_lssubgauss(eleme,lse,8);
        end
        
        if ~calc_cut
            if any(lse>0) && any(lse<=0)
                gauss.w(:)=0;
            end
        end
        N = calc_N(eleme,xnodee,gauss.coord);
        x = calc_x(eleme,xnodee,gauss.coord);
        
        psi = -N*lse;
        
        %psi = (N*lse).^2;
        %  psi =  exp(-N*lse)-1;
        
        deltau = (N*ue)*psi - uex(x,xi);
        detJ = calc_detJ(eleme,xnodee,gauss.coord);
        ke= double(sum(gauss.w*deltau.*deltau*abs(detJ),4));
        ke1= double(sum(gauss.w*uex(x,xi).*uex(x,xi)*abs(detJ),4));
        %keyboard
        %    DN = getDN()
        
        
end

