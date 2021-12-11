function V = shinozukaRand_scalar(c,z,phi,k,x)
% function V = shinozukaRand_scalar(c,z,phi,k,x)

[nx,dim] = size(x);
order = size(k,2);

V = zeros(nx,1);
switch dim
    case 1
        for i=1:order
            Vi = sqrt(2) * c(i) * z(i) * cos(phi(i) + k(i)*x);
            V = V + Vi;
        end
    case 2
        for i=1:order
            Vi = sqrt(2) * c(i) * z(i) * cos(phi(i) + k(1,i)*x(:,1) + k(2,i)*x(:,2));
            V = V + Vi;
        end
    case 3
        for i=1:order
            Vi = sqrt(2) * c(i) * z(i) * cos(phi(i) + k(1,i)*x(:,1) + k(2,i)*x(:,2) + k(3,i)*x(:,3));
            V = V + Vi;
        end
end

end