function delta = solveCyclicCompression(x,N)
% function delta = solveCyclicCompression(x,N)

delta0 = x(1);
Kinf = x(2);
Nref = x(3);

delta = (delta0+Kinf*N).*(1-exp(-N/Nref));

end
