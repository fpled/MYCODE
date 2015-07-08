function lmin = calc_lmin_patch(lmax,S,Spartfd,AS,Atilde,Apatch,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf)
% function lmin = calc_lmin_patch(lmax,S,Spartfd,AS,Atilde,Apatch,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf)

a = @(V) rho_optimal(V,S,Spartfd,AS,Atilde,Apatch,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf)  
OPTS.v0 = rand(getnbddlfree(S),1);
% OPTS.isreal = false;
OPTS.maxit=100;
% OPTS.tol=[1e-2];

lmin = lmax-eigs(@(V) lmax*V-a(V), getnbddlfree(S),1,'lm',OPTS);

end
