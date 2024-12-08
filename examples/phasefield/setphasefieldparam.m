function [k,r,qn] = setphasefieldparam(gc,l,PFregularization)
% function [k,r,qn] = setphasefieldparam(gc,l,PFregularization)

switch lower(PFregularization)
    case 'at1'
        % c0 = 8/3;
        k = 3/4*gc.*l; % k = 2*(gc.*l)/c0;
        r = 0;
        qn = -3/8*gc./l; % qn = -(gc./l)/c0;
    case 'at2'
        % c0 = 2;
        k = gc.*l; % k = 2*(gc.*l)/c0;
        r = gc./l; % r = 2*(gc./l)/c0;
        qn = 0;
    otherwise
        error('Wrong regularization model');
end
