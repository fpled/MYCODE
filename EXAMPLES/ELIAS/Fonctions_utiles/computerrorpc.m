function [norm,fieldtype,fieldstorage,fieldddl] = computerrorpc(elem,node,ls,RV,u,uex,varargin)
%function [norm,fieldtype,fieldstorage,fieldddl]=computerrorpc(elem,node,ls,RV,u,uex,varargin)
fieldtype='scalar';
fieldstorage='center';
fieldddl=DDL('error');
switch getlstype(elem)
    case {'out'}
        norm{1} = zerosND(1,1,getnbelem(elem));
        norm{2} = zerosND(1,1,getnbelem(elem));
    otherwise
        PC = getPC(u);
        calc_cut = getcharin('calc_cut',varargin,1);
        nbgausssto = getcharin('nbgausssto',varargin,2*(getorder(PC)+1));
        nbsubgausssto =  getcharin('nbsubgausssto',varargin,2*4);
        method = getcharin('method',varargin);
        if nbsubgausssto>1
            gauss_sto = calc_subgausspoints(PC,nbgausssto,nbsubgausssto);
        else
            gauss_sto = calc_gausspoints(PC,nbgausssto);
        end
        norm{1} = zerosND(1,1,getnbelem(elem));
        norm{2} = zerosND(1,1,getnbelem(elem));
        xnode = getcoord(node,elem);
        connec = getconnec(elem);
        ls = randomeval(ls,gauss_sto.coord,RANDVARS(PC));
        ueval = randomeval(u,gauss_sto.coord);
        xi = transfer(RANDVARS(PC),RV,gauss_sto.coord);
        
        for e=1:getnbelem(elem)
            pourcentage(e,getnbelem(elem),10)
            eleme = getelem(elem,e);
            connece = connec(e,:);
            lse = ls(connece);
            xnodee = xnode(:,:,e);
            ue = localize(eleme,ueval);
            for j=1:gauss_sto.nbgauss
                if ~all(lse{j}>0)
                    [norms(j) normsref(j)] = computerror_elem(eleme,xnodee,lse{j},ue(:,j),uex,xi(j,:),'method',method,'calc_cut',calc_cut);
                else
                    norms(j) = 0;
                    normsref(j) = eps;
                end
            end
            norm{1}(:,:,e) = sum(norms.*gauss_sto.w);
            norm{2}(:,:,e) = sum(normsref.*gauss_sto.w);
        end
end

