function [norme,fieldtype,fieldstorage,fieldddl] = dusurdn(elem,node,ls,u,varargin)
% function [norme,fieldtype,fieldstorage,fieldddl] = dusurdn(elem,node,ls,u,varargin)
% Calcule le double integrale (sur Xi et sur la frontiere de Omega) de (drond u/drond n)^2

fieldtype='scalar';
fieldstorage='center';
fieldddl=DDL('error');
switch getlstype(elem)
    case {'in','out'}
        norme{1} = zerosND(1,1,getnbelem(elem));
        
    case {'cut'}
        PC = getPC(u);
        %
        nbgausssto = getcharin('nbgausssto',varargin,2*(getorder(PC)+1));
        nbsubgausssto =  getcharin('nbsubgausssto',varargin,2*4);
        %
        type=getpoly(PC,1);
        if  strcmp(class(type),'POLYFE')
            gauss_sto = calc_gausspoints(PC,nbgausssto);
        elseif  strcmp(class(type),'POLYLEGENDRE')
            if nbsubgausssto>1
                gauss_sto = calc_subgausspoints(PC,nbgausssto,nbsubgausssto);
            else
                gauss_sto = calc_gausspoints(PC,nbgausssto);
            end
        elseif strcmp(class(type),'POLYFELAGRANGE')
            
            gauss_sto = calc_gausslobattopoints(type,nbgausssto);
        end
        
        norme{1} = zerosND(1,1,getnbelem(elem));
        
        xnode = getcoord(node,elem);
        connec = getconnec(elem);
        ls = randomeval(ls,gauss_sto.coord,RANDVARS(PC));
        ueval = randomeval(u,gauss_sto.coord);
        %xi = transfer(RANDVARS(PC),RV,gauss_sto.coord);
        
        
        
        for e=1:getnbelem(elem);
            pourcentage(e,getnbelem(elem),10)
            eleme= getelem(elem,e);
            connece = getconnec(eleme);
            xnodee = xnode(:,:,e);
            lse=ls(connece);
            ue = localize(eleme,ueval);
            
            
            for i=1:gauss_sto.nbgauss
                
                if ~(all(lse{i}>=0) || all(lse{i}<0))
                    xnodeelocal = nodelocalcoord(eleme);
                    xnodeseg = contour_oneelem(lse{i},xnodeelocal);
                    
                    elemseg = SEG2(xnodeseg,1,1:2);
                    xvectseg = xnodeseg(1,:)-xnodeseg(2,:);
                    
                    xvectnormal(1) = -xvectseg(1,2);
                    xvectnormal(2) = xvectseg(1,1);
                    
                    xvectnormal = normalize(xvectnormal);
                    
                    gauss = calc_gauss(elemseg,2);
                    xi = calc_x(elemseg,xnodeseg,gauss.coord);
                    xnodeseg=double(xnodeseg);
                    x1 = calc_x(eleme,xnodee,xnodeseg(1,:));
                    x2 = calc_x(eleme,xnodee,xnodeseg(2,:));
                    le = norm(x2-x1);
                    DN = calc_DN(eleme,xnodee,xi);
                    du = DN*ue(:,i);
                    
                    ke= double((sum(gauss.w*(le/2)*(xvectnormal*du)*(xvectnormal*du))));
                    
                    
                    norms(i) = ke;
                    
                else
                    norms(i) = 0;
                end
                
            end
            norme{1}(:,:,e) = sum(norms.*gauss_sto.w);
        end
end





