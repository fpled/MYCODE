function V = shinozuka_sum(c,z,phi,k,x)
% function V = shinozuka_sum(c,z,phi,k,x)

dim  = size(x,2);

switch dim
    case 1
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            V = sqrt(2) * sum(bsxfun(@times,c.*z,cos(phi + bsxfun(@times,k,x))),2);
        else
            V = sqrt(2) * sum(c .* z .* cos(phi + k.*x),2);
        end
    case 2
        [c1,c2] = ndgrid(c(1,:),c(2,:));
        [k1,k2] = ndgrid(k(1,:),k(2,:));
        c1 = c1(:)'; c2 = c2(:)';
        k1 = k1(:)'; k2 = k2(:)';
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            V = sqrt(2) * sum(bsxfun(@times,c1.*c2 .* z,cos(phi + bsxfun(@times,k1,x(:,1)) + bsxfun(@times,k2,x(:,2)))),2);
        else
            V = sqrt(2) * sum(c1.*c2 .* z .* cos(phi + k1.*x(:,1) + k2.*x(:,2)),2);
        end
    case 3
        [c1,c2,c3] = ndgrid(c(1,:),c(2,:),c(3,:));
        [k1,k2,k3] = ndgrid(k(1,:),k(2,:),k(3,:));
        c1 = c1(:)'; c2 = c2(:)'; c3 = c3(:)';
        k1 = k1(:)'; k2 = k2(:)'; k3 = k3(:)';
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            V = sqrt(2) * sum(bsxfun(@times,c1.*c2.*c3 .* z,cos(phi + bsxfun(@times,k1,x(:,1)) + bsxfun(@times,k2,x(:,2)) + bsxfun(@times,k3,x(:,3)))),2);
        else
            V = sqrt(2) * sum(c1.*c2.*c3 .* z .* cos(phi + k1.*x(:,1) + k2.*x(:,2) + k3.*x(:,3)),2);
        end
end

end