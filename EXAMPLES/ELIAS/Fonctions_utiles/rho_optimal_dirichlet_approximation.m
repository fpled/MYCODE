function repA = rho_optimal_dirichlet_approximation(V,S,Spartf,Spartfd,AS,Atilde,Apatch,A2patch,psim,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf,PC,PCphi)
% function repA = rho_optimal_dirichlet_approximation(V,S,Spartf,Spartfd,AS,Atilde,Apatch,A2patch,psim,Ppatchfgammaf,Pgammaoutgammaf,Pgamma,Dlf,PC,PCphi)

n = getnbddlfree(S);
tempphi = zeros(1,n);
for k=2:length(Apatch)
    Bpatchr{k} = - freevector(Spartfd{k},A2patch{k}*(Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k}*project(V,PCphi{k})));
    Apatchr{k} = freematrix(Spartfd{k},Apatch{k});
    %
    w1{k} = solve(Apatchr{k},Bpatchr{k});
    w1{k} = unfreevector(Spartfd{k},w1{k});
    %
    w{k} = project(psim{k},PC).*w1{k} + Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k}*project(V,PCphi{k});
    %
    % blambda = Ppatchfgammaf{k}*Apatch{k}*w{k} - Ppatchfgammaf{k}*calc_vector(LINFORM(0),Spartf{k});
    blambda = project(Ppatchfgammaf{k}*calc_matrix(BILINFORM(1,1),Spartf{k})*Ppatchfgammaf{k}'*Pgammaoutgammaf{k}*Pgamma{k}*project(V,PCphi{k}),PC)+project(Ppatchfgammaf{k}*A2patch{k}'*unfreevector(Spartfd{k},w1{k}),PC);
    %
    lambda{k} = solve((random(Dlf{k})),blambda) ;
    tempphi = tempphi - Pgamma{k}'*(Pgammaoutgammaf{k}'*Dlf{k}*(lambda{k}))+ Atilde{k}*V;
end

% temp = zeros(1,getnbddlfree(S))';
% for k=2:length(myboxes)+1
%     temp = temp - Pgamma{k}'*(Pgammaoutgammaf{k}'*Dlf{k}*(beta{k}));
%     temp = temp +  Atilde{k}*V;
% end

% Uhat = cgs(AS,  P{1}'*bpatch{1} + temp,tolcgs,1000);
phi = solve(AS,tempphi);
repA = expect(V-phi);
end
