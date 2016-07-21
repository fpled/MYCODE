function lmax = calc_lmax_patch(S,Spartfd,AS,Atilde,Apatch,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf)
% function lmax = calc_lmax_patch(S,Spartfd,AS,Atilde,Apatch,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf)

a = @(V) rho_optimal(V,S,Spartfd,AS,Atilde,Apatch,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf)  
OPTS.v0 = rand(getnbddlfree(S),1);
% OPTS.isreal = false;
OPTS.maxit=100;
% OPTS.tol=[1e-2];

lmax = eigs(a,getnbddlfree(S),1,'lm',OPTS)

end
