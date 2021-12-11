function V = shinozukaRand_sum(c,z,phi,k,x)
% function V = shinozukaRand_sum(c,z,phi,k,x)

dim  = size(x,2);

switch dim
    case 1
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            V = sqrt(2) * sum(bsxfun(@times,c.*z,cos(phi + bsxfun(@times,k,x))),2);
        else
            V = sqrt(2) * sum(c .* z .* cos(phi + k.*x),2);
        end
    case 2
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            V = sqrt(2) * sum(bsxfun(@times,c.*z,cos(phi + bsxfun(@times,k(1,:),x(:,1)) + bsxfun(@times,k(2,:),x(:,2)))),2);
        else
            V = sqrt(2) * sum(c .* z .* cos(phi + k(1,:).*x(:,1) + k(2,:).*x(:,2)),2);
        end
    case 3
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            V = sqrt(2) * sum(bsxfun(@times,c.*z,cos(phi + bsxfun(@times,k(1,:),x(:,1)) + bsxfun(@times,k(2,:),x(:,2)) + bsxfun(@times,k(3,:),x(:,3)))),2);
        else
            V = sqrt(2) * sum(c .* z .* cos(phi + k(1,:).*x(:,1) + k(2,:).*x(:,2) + k(3,:).*x(:,3)),2);
        end
end

end