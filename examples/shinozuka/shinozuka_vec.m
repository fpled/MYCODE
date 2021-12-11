function V = shinozuka_vec(c,z,phi,k,x)
% function V = shinozuka_vec(c,z,phi,k,x)

dim = size(x,2);

switch dim
    case 1
        K = k;
        C = c';
    case 2
        [c1,c2] = ndgrid(c(1,:),c(2,:));
        [k1,k2] = ndgrid(k(1,:),k(2,:));
        c1 = c1(:); c2 = c2(:);
        k1 = k1(:); k2 = k2(:);
        K = [k1 k2]';
        C = [c1 c2];
    case 3
        [c1,c2,c3] = ndgrid(c(1,:),c(2,:),c(3,:));
        [k1,k2,k3] = ndgrid(k(1,:),k(2,:),k(3,:));
        c1 = c1(:); c2 = c2(:); c3 = c3(:);
        k1 = k1(:); k2 = k2(:); k3 = k3(:);
        K = [k1 k2 k3]';
        C = [c1 c2 c3];
end
V =  sqrt(2) * cos(phi + x*K) * (prod(C,2) .* z');

end