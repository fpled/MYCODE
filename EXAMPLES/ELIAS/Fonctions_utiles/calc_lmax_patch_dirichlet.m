function lmax = calc_lmax_patch_dirichlet(S,Spartf,Spartfd,AS,Atilde,Apatch,A2patch,psim,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf)
% function lmax = calc_lmax_patch_dirichlet(S,Spartf,Spartfd,AS,Atilde,Apatch,A2patch,psim,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf)

a = @(V) rho_optimal_dirichlet(V,S,Spartf,Spartfd,AS,Atilde,Apatch,A2patch,psim,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf)  
OPTS.v0 = rand(getnbddlfree(S),1);
% OPTS.isreal = false;
OPTS.maxit=100;
OPTS.tol=[1e-20];

lmax = eigs(a,getnbddlfree(S),1,'lm',OPTS)

end
