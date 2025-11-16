function [k,r,qn] = sethealingfieldparam(hc,l,PFregularization)
% function [k,r,qn] = sethealingfieldparam(hc,l,PFregularization)

switch lower(PFregularization)
    case 'at1'
        % c0 = 8/3;
        k = 3/4*hc.*l; % k = 2*(hc.*l)/c0;
        r = 0;
        qn = -3/8*hc./l; % qn = -(hc./l)/c0;
    case 'at2'
        % c0 = 2;
        k = hc.*l; % k = 2*(hc.*l)/c0;
        r = hc./l; % r = 2*(hc./l)/c0;
        qn = 0;
    otherwise
        error('Wrong regularization model');
end
