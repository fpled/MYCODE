function V = shinozuka_scalar(c,z,phi,k,x)
% function V = shinozuka_scalar(c,z,phi,k,x)

[nx,dim] = size(x);
nu = size(k,2);

V = zeros(nx,1);
switch dim
    case 1
        for i=1:nu
            Vi = sqrt(2) * c(i) * z(i) * cos(phi(i) + k(i)*x);
            V = V + Vi;
        end
    case 2
        z = reshape(z,[nu,nu]);
        phi = reshape(phi,[nu,nu]);
        for i=1:nu
            for j=1:nu
                Vij = sqrt(2) * c(1,i)*c(2,j) * z(i,j) * cos(phi(i,j) + k(1,i)*x(:,1) + k(2,j)*x(:,2));
                V = V + Vij;
            end
        end
    case 3
        z = reshape(z,[nu,nu,nu]);
        phi = reshape(phi,[nu,nu,nu]);
        for i=1:nu
            for j=1:nu
                for kk=1:nu
                    Vijk = sqrt(2) * c(1,i)*c(2,j)*c(3,kk) * z(i,j,kk) * cos(phi(i,j,kk) + k(1,i)*x(:,1) + k(2,j)*x(:,2) + k(3,kk)*x(:,3));
                    V = V + Vijk;
                end
            end
        end
end

end