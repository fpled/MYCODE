function [norme,fieldtype,fieldstorage,fieldddl] = int_gradientu(elem,node,ls,u,varargin)
%function [norme,fieldtype,fieldstorage,fieldddl] =int_gradientu(elem,node,ls,u,varargin)
% Calcule le double int�grale du (gradient u)^2
%"calc_cut" --> 1 ou 0, pour prendre en compte l'erreur sur les �l�ments cut (1) ou non (0)

fieldtype='scalar';
fieldstorage='center';
fieldddl=DDL('error');
switch getlstype(elem)
    case {'out'}
        norme{1} = zerosND(1,1,getnbelem(elem));
        
    otherwise
        PC = getPC(u);
        calc_cut = getcharin('calc_cut',varargin,1);
        nbgausssto = getcharin('nbgausssto',varargin,2*(getorder(PC)+1));
        nbsubgausssto =  getcharin('nbsubgausssto',varargin,2*4);
        
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
        
        
        for e=1:getnbelem(elem);
            pourcentage(e,getnbelem(elem),10)
            eleme= getelem(elem,e);
            connece = getconnec(eleme);
            xnodee = xnode(:,:,e);
            lse=ls(connece);
            ue = localize(eleme,ueval);
            for j=1:gauss_sto.nbgauss
                if ~all(lse{j}>0)
                    
                    if all(lse{j}<=0)
                        gauss = calc_gauss(eleme,8);
                    elseif ~all(lse{j}>=0)
                        [gauss,gaussout] = calc_lssubgauss(eleme,lse{j},8);
                    end
                    
                    if ~calc_cut
                        if any(lse{j}>0) && any(lse{j}<=0)
                            gauss.w(:)=0;
                        end
                    end
                    DN = calc_DN(eleme,xnodee,gauss.coord);
                    x = calc_x(eleme,xnodee,gauss.coord);
                    gradu = DN*ue(:,j);
                    detJ = calc_detJ(eleme,xnodee,gauss.coord);
                    ke = double(sum(gauss.w*(gradu'*gradu)*abs(detJ),4));
                    norms(j) = ke;
                else
                    norms(j) = 0;
                end
            end
            norme{1}(:,:,e) = sum(norms.*gauss_sto.w);
            
        end
end

