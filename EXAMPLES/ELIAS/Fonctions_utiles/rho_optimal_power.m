function [lmax,V,lmin,r]=rho_optimal_power(maxit,tol)
%function [lmax,V,lmin,ro]=rho_optimal_power(nbit,tol)

%if (nargin==0)
maxit = 10
tol=1e-3
%elseif (nargin==1)
%tol=1e-20
%end
prog_patch;
V = rand(getnbddlfree(S),1);
Y = rand(getnbddlfree(S),1);

for it =1:maxit
    AV = iterpatch(V);
    AY = 0%iterpatch(Y);
    errV = norm(AV/norm(AV)-V)/norm(V);
    % errY = norm(AY/norm(AY)-Y)/norm(Y);
    if (errV<tol)
        break
    end
    if isa(V,'double')
        lmax = expect(V'*AV)/(expect(V')*expect(V))%(norm(V)^2);
    elseif isa(V,'PCMATRIX')
        lmax = (expect(V')*expect(AV))/(expect(V')*expect(V))%(norm(V)^2);
    end
    
    lmin = lmax - (expect(Y')*expect(lmax*Y-AY))/(expect(Y')*expect(Y))
    
    V = AV/norm(AV);
    Y = (lmax*Y-AY)/norm(lmax*Y-AY);
    
    rooptim = 2/(lmin+lmax);
    
    fprintf('iteration = %d, erreur v :%d, lmax=%d, lmin=%d, ro=%d\n',it,errV,lmax(1),lmin(1),rooptim)
end

%end
