function S_phase = setphasefieldproperties(S_phase,materials)
% function S_phase = setphasefieldproperties(S_phase,materials)

node_phase = getnode(S_phase);
for m=1:length(materials)
    mat = materials{m};
    if isparam(mat,'delta') && any(getparam(mat,'delta')>0) && isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
        elem = getgroupelem(S_phase,m);
        nbelem = getnbelem(elem);
        xnode = node_phase(elem);
        gauss = calc_gauss(elem,'mass');
        xgauss = gauss.coord;
        x = calc_x(elem,xnode,xgauss);
        x = getcoord(NODE(POINT(x(:,:,:))));
        shinozukaPF = getparam(mat,'shinozuka');
        Xi = shinozukaPF(x); % sample for bivariate Gaussian random field with statistically independent normalized Gaussian components
        k = evalparam(mat,'k',elem,xnode,xgauss);
        r = evalparam(mat,'r',elem,xnode,xgauss);
        gc = sqrt(k.*r); % mean fracture toughness
        l = sqrt(k./r); % mean regularization parameter
        delta = getparam(mat,'delta'); % coefficients of variation for fracture toughness and regularization parameter
        if length(delta)==1
            deltaGc = delta; % 0 <= deltaGc < 1/sqrt(2). coefficient of variation for fracture toughness
            deltaL = delta; % 0 <= deltaL < 1/sqrt(2). coefficient of variation for regularization parameter
        else
            deltaGc = delta(1); % 0 <= deltaGc < 1/sqrt(2). coefficient of variation for fracture toughness
            deltaL = delta(2); % 0 <= deltaL < 1/sqrt(2). coefficient of variation for regularization parameter
        end
        if deltaGc<0 || deltaGc>=1/sqrt(2)
            error('Coefficient of variation delta = %g for fracture toughness should be between 0 and %g',deltaGc,1/sqrt(2))
        end
        if deltaL<0 || deltaL>=1/sqrt(2)
            error('Coefficient of variation delta = %g for regularization parameter should be between 0 and %g',deltaL,1/sqrt(2))
        end
        aGc = 1/deltaGc^2; % aGc > 2
        bGc = gc/aGc; % 0 < bGc = gc/aGc < gc/2 since gc > 0 and aGc > 2
        aL = 1/deltaL^2; % aL > 2
        bL = l/aL; % 0 < bL = l/aL < l/2 since l > 0 and aL > 2
        if deltaGc && deltaL
            rho = 0;
            if isparam(mat,'rcorr')
                rho = getparam(mat,'rcorr'); % correlation coefficient between fracture toughness and regularization parameter
            end
            gc = gaminv(normcdf(Xi(:,1)),aGc,bGc); % sample for fracture toughness [N/m^2]
            l = gaminv(normcdf(rho*Xi(:,1) + sqrt(1-rho^2)*Xi(:,2)),aL,bL); % sample for regularization parameter [m]
        elseif deltaGc
            gc = gaminv(normcdf(Xi(:,1)),aGc,bGc); % sample for fracture toughness [N/m^2]
        else
            l = gaminv(normcdf(Xi(:,1)),aL,bL); % sample for regularization parameter [m]
        end
        if deltaGc
            gc = reshape(gc,1,1,nbelem,gauss.nbgauss);
            gc = MYDOUBLEND(gc);
            gc = FEELEMFIELD({gc},'storage','gauss','type','scalar','ddl',DDL('gc'));
        end
        if deltaL
            l = reshape(l,1,1,nbelem,gauss.nbgauss);
            l = MYDOUBLEND(l);
            l = FEELEMFIELD({l},'storage','gauss','type','scalar','ddl',DDL('l'));
        end
        k = gc.*l;
        r = gc./l;
        mat = setparam(mat,'k',k);
        mat = setparam(mat,'r',r);
    end
    S_phase = setmaterial(S_phase,mat,m);
end
