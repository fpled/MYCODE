function repA = rho_optimal_dirichlet_diffusion(V,S,Spartf,Spartfd,AS,Atilde,PC,Apatchsep,A2patchsep,psitp,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlfsep,Dlf,numpatchdiff,numpatchdir)
%function repA = rho_optimal_dirichlet_diffusion(V,S,Spartf,Spartfd,AS,Atilde,PC,Apatchsep,A2patchsep,psitp,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlfsep,Dlf,numpatchdiff,numpatchdir)
% liens
tolPGD = 5e-3;
PGDrank = 100;
n = getnbddlfree(S);
V = SEPMATRIX([{V}  getphi(one(PC))],1)
tempphi = SEPMATRIX([{zeros(getnbddlfree(S),1)} getphi(one(PC))],1);

for ii=1:length(numpatchdiff)
    k = numpatchdiff(ii)+1;
    Bpatchrsep{k} = - freevector(Apatchsep{k}*mtimes(Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k},V,1),1,Spartfd{k});
    Apatchrsep{k} = freematrix(Apatchsep{k},1,Spartfd{k});
    
    solver = SEPSOLVER(getdim(Apatchrsep{k}),'tol',tolPGD,...
        'update',2,'updatedim',2:getdim(Apatchrsep{k}),...
        'maxorder',PGDrank,'maxiter',5,'reference',[],...
        'errorindicator','none','itercrit',1e-3,...
        'righthandSD',false,'display',true);
    
    [w1{k},resultsep{k}] = solve(Apatchrsep{k},Bpatchrsep{k},solver);
    w1{k} = unfreevector(w1{k},1,Spartfd{k});
    %
    % w1{k} = dw1{k} + unfreevector(Spartfd{k},w1prec{k});
    % w{k} = w1{k} + mtimes(Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k},U,1);
    w{k} = calc_ximasse(PCTPMATRIX(w1{k},PC))+ ...
        Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k}*PCTPMATRIX(V,PC);
    %
    blambda = mtimes(Ppatchfgammaf{k},Apatchsep{k},1)*mtimes(Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k},V,1) + mtimes(Ppatchfgammaf{k},Apatchsep{k}',1)*unfreevector(w1{k},1,Spartfd{k});
    %
    lambda{k} = mldivide(random(Dlf{k}),blambda,1);
end

for ii=1:length(numpatchdir)
    k= numpatchdir(ii)+1;
    Bpatchrsep{k} = - freevector(A2patchsep{k}*mtimes(Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k},V,1),1,Spartfd{k});
    Apatchrsep{k} = freematrix(Apatchsep{k},1,Spartfd{k});
    
    solver = SEPSOLVER(getdim(Apatchrsep{k}),'tol',tolPGD,...
        'update',2,'updatedim',2:getdim(Apatchrsep{k}),...
        'maxorder',PGDrank,'maxiter',5,'reference',[],...
        'errorindicator','none','itercrit',1e-3,...
        'righthandSD',false,'display',true);
    
    [w1{k},resultsep{k}] = solve(Apatchrsep{k},Bpatchrsep{k},solver);
    % w1{k} = dw1{k}+w1prec{k};
    w{k} = psitp{k}.*calc_ximasse((unfreevector(Spartfd{k},PCTPMATRIX(w1{k},PC))))+ ...
        Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k}*PCTPMATRIX(V,PC);
    %
    blambdadir = mtimes(Ppatchfgammaf{k}*calc_matrix(BILINFORM(1,1),Spartf{k})*Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k},V,1)+mtimes(Ppatchfgammaf{k},A2patchsep{k}',1)*unfreevector(w1{k},1,Spartfd{k});
    lambda{k} = mldivide(random(Dlf{k}),blambdadir,1);
end

for k=2:length(Apatchsep)
    tempphi = tempphi-mtimes(Pgamma{k}'*Pgammaoutgammaf{k}',Dlfsep{k},1)*lambda{k} +(mtimes(Atilde{k},V,1));
end

phi = mldivide(AS,tempphi,1);
repA = expect(PCTPMATRIX(V-phi,PC));
end
