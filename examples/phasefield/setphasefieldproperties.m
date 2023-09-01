function S_phase = setphasefieldproperties(S_phase,materials)
% function S_phase = setphasefieldproperties(S_phase,materials)

node_phase = getnode(S_phase);
for m=1:length(materials)
    mat = materials{m};
%     if isparam(mat,'delta') && any(getparam(mat,'delta')>0) && isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
%         elem = getgroupelem(S_phase,m);
%         nbelem = getnbelem(elem);
%         xnode = node_phase(elem);
%         gauss = calc_gauss(elem,'mass');
%         xgauss = gauss.coord;
%         x = calc_x(elem,xnode,xgauss);
%         x = getcoord(NODE(POINT(x(:,:,:))));
%         shinozukaPF = getparam(mat,'shinozuka');
%         Xi = shinozukaPF(x); % sample for bivariate Gaussian random field with statistically independent normalized Gaussian components
%         PFregularization = getparam(mat,'PFregularization');
%         k = evalparam(mat,'k',elem,xnode,xgauss);
%         switch lower(PFregularization)
%             case 'at1'
%                 c0 = 8/3;
%                 qn = evalparam(mat,'qn',elem,xnode,xgauss);
%                 gc = c0*sqrt(-k.*qn/2); % mean fracture toughness
%                 l = sqrt(-k./qn/2); % mean regularization parameter
%             case 'at2'
%                 % c0 = 2;
%                 r = evalparam(mat,'r',elem,xnode,xgauss);
%                 gc = sqrt(k.*r); % mean fracture toughness
%                 l = sqrt(k./r); % mean regularization parameter
%             otherwise
%                 error('Wrong regularization model');
%         end
%         delta = getparam(mat,'delta'); % coefficients of variation for fracture toughness and regularization parameter
%         if length(delta)==1
%             deltaGc = delta; % 0 <= deltaGc < 1/sqrt(2). coefficient of variation for fracture toughness
%             deltaL = delta; % 0 <= deltaL < 1/sqrt(2). coefficient of variation for regularization parameter
%         else
%             deltaGc = delta(1); % 0 <= deltaGc < 1/sqrt(2). coefficient of variation for fracture toughness
%             deltaL = delta(2); % 0 <= deltaL < 1/sqrt(2). coefficient of variation for regularization parameter
%         end
%         if deltaGc<0 || deltaGc>=1/sqrt(2)
%             error('Coefficient of variation delta = %g for fracture toughness should be between 0 and %g',deltaGc,1/sqrt(2))
%         end
%         if deltaL<0 || deltaL>=1/sqrt(2)
%             error('Coefficient of variation delta = %g for regularization parameter should be between 0 and %g',deltaL,1/sqrt(2))
%         end
%         aGc = 1/deltaGc^2; % aGc > 2
%         bGc = gc/aGc; % 0 < bGc = gc/aGc < gc/2 since gc > 0 and aGc > 2
%         aL = 1/deltaL^2; % aL > 2
%         bL = l/aL; % 0 < bL = l/aL < l/2 since l > 0 and aL > 2
%         if deltaGc && deltaL
%             rho = 0;
%             if isparam(mat,'rcorr')
%                 rho = getparam(mat,'rcorr'); % correlation coefficient between fracture toughness and regularization parameter
%             end
%             gc = gaminv(normcdf(Xi(:,1)),aGc,bGc); % sample for fracture toughness [N/m]
%             l = gaminv(normcdf(rho*Xi(:,1) + sqrt(1-rho^2)*Xi(:,2)),aL,bL); % sample for regularization parameter [m]
%         elseif deltaGc
%             gc = gaminv(normcdf(Xi(:,1)),aGc,bGc); % sample for fracture toughness [N/m]
%         else
%             l = gaminv(normcdf(Xi(:,1)),aL,bL); % sample for regularization parameter [m]
%         end
%         if deltaGc
%             gc = reshape(gc,1,1,nbelem,gauss.nbgauss);
%             gc = MYDOUBLEND(gc);
%             gc = FEELEMFIELD({gc},'storage','gauss','type','scalar','ddl',DDL('gc'));
%         end
%         if deltaL
%             l = reshape(l,1,1,nbelem,gauss.nbgauss);
%             l = MYDOUBLEND(l);
%             l = FEELEMFIELD({l},'storage','gauss','type','scalar','ddl',DDL('l'));
%         end
%         switch lower(PFregularization)
%             case 'at1'
%                 % c0 = 8/3;
%                 k = 3/4*gc.*l; % k = 2*(gc.*l)/c0;
%                 qn = -3/8*gc./l; % qn = -(gc./l)/c0;
%                 mat = setparam(mat,'k',k);
%                 mat = setparam(mat,'qn',qn);
%             case 'at2'
%                 % c0 = 2;
%                 k = gc.*l; % k = 2*(gc.*l)/c0;
%                 r = gc./l; % r = 2*(gc./l)/c0;
%                 mat = setparam(mat,'k',k);
%                 mat = setparam(mat,'r',r);
%             otherwise
%                 error('Wrong regularization model');
%         end
%     end
    if isparam(mat,'aGc') && isparam(mat,'bGc') && any(getparam(mat,'aGc')>0) && any(getparam(mat,'bGc')>0)...
            && isparam(mat,'lcorr') && ~all(isinf(getparam(mat,'lcorr'))) % random field model
        elem = getgroupelem(S_phase,m);
        nbelem = getnbelem(elem);
        xnode = node_phase(elem);
        gauss = calc_gauss(elem,'mass');
        xgauss = gauss.coord;
        x = calc_x(elem,xnode,xgauss);
        x = getcoord(NODE(POINT(x(:,:,:))));
        shinozukaPF = getparam(mat,'shinozuka');
        Xi = shinozukaPF(x); % sample for univariate Gaussian random field
        PFregularization = getparam(mat,'PFregularization');
        k = evalparam(mat,'k',elem,xnode,xgauss);
        switch lower(PFregularization)
            case 'at1'
                qn = evalparam(mat,'qn',elem,xnode,xgauss);
                % l = sqrt(-k./qn/2); % regularization parameter
                r = -2*qn;
            case 'at2'
                r = evalparam(mat,'r',elem,xnode,xgauss);
            otherwise
                error('Wrong regularization model');
        end
        l = sqrt(k./r); % regularization parameter
        aGc = getparam(mat,'aGc'); % lower bound for fracture toughness aGc > 0
        bGc = getparam(mat,'bGc'); % upper bound for fracture toughness bGc > aGc > 0
        if aGc<0 || bGc<0
            error('Lower bound a = %g and upper bound b = %g for fracture toughness must be positive (superior to 0)',aGc,bGc)
        end
        if aGc>bGc
            error('Lower bound a = %g must be inferior to upper bound b = %g for fracture toughness',aGc,bGc)
        end
        gc = unifinv(normcdf(Xi),aGc,bGc); % sample for fracture toughness [N/m]
        gc = reshape(gc,1,1,nbelem,gauss.nbgauss);
        gc = MYDOUBLEND(gc);
        gc = FEELEMFIELD({gc},'storage','gauss','type','scalar','ddl',DDL('gc'));
        switch lower(PFregularization)
            case 'at1'
                % c0 = 8/3;
                k = 3/4*gc.*l; % k = 2*(gc.*l)/c0;
                qn = -3/8*gc./l; % qn = -(gc./l)/c0;
                mat = setparam(mat,'qn',qn);
            case 'at2'
                % c0 = 2;
                k = gc.*l; % k = 2*(gc.*l)/c0;
                r = gc./l; % r = 2*(gc./l)/c0;
                mat = setparam(mat,'r',r);
            otherwise
                error('Wrong regularization model');
        end
        mat = setparam(mat,'k',k);
    end
    S_phase = setmaterial(S_phase,mat,m);
end
