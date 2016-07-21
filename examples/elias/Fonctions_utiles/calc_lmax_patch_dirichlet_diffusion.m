function lmax = calc_lmax_patch_dirichlet_diffusion(S,Spartf,Spartfd,AS,Atilde,PC,Apatchsep,A2patchsep,psitp,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlfsep,Dlf,numpatchdiff,numpatchdir)
% function lmax = calc_lmax_patch_dirichlet_diffusion(S,Spartf,Spartfd,AS,Atilde,PC,Apatchsep,A2patchsep,psitp,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlfsep,Dlf,numpatchdiff,numpatchdir)

a = @(V) rho_optimal_dirichlet_diffusion(V,S,Spartf,Spartfd,AS,Atilde,PC,Apatchsep,A2patchsep,psitp,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlfsep,Dlf,numpatchdiff,numpatchdir)  
OPTS.v0 = rand(getnbddlfree(S),1);
% OPTS.isreal = false;
OPTS.maxit=10;
OPTS.tol=[1e-3];

lmax = eigs(a,getnbddlfree(S),1,'lm',OPTS)

end
