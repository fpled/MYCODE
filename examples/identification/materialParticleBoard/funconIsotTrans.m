function [c,ceq,gc,gceq] = funconIsotTrans(lambda)
% function [c,ceq,gc,gceq] = funconIsotTrans(lambda)

useRedParam = (length(lambda)==4); % reduced parameterization

% Compute nonlinear inequality constraints
c = lambda(3)-2*sqrt(lambda(1)*lambda(2)); % c(lambda) < 0
% c(1) = lambda(3)-2*sqrt(lambda(1)*lambda(2)); % c(lambda) < 0
% if useRedParam, c(2) = lambda(4); else, c(2) = lambda(6); end % c(lambda) < 0

% Compute nonlinear equality constraints
ceq = []; % ceq(lambda) = 0

if nargout>2
    % Compute gradient of nonlinear inequality constraints
    if useRedParam
        gc = [-sqrt(lambda(2)./lambda(1)); -sqrt(lambda(1)./lambda(2)); 1; 0];
        % gc = [-sqrt(lambda(2)./lambda(1)), 0;
        %       -sqrt(lambda(1)./lambda(2)), 0;
        %                                 1, 0;
        %                                 0, 1];
    else
        gc = [-sqrt(lambda(2)./lambda(1)); -sqrt(lambda(1)./lambda(2)); 1; 0; 0; 0];
        % gc = [-sqrt(lambda(2)./lambda(1)), 0;
        %       -sqrt(lambda(1)./lambda(2)), 0;
        %                                 1, 0;
        %                                 0, 0;
        %                                 0, 0;
        %                                 0, 1];
    end
    % Compute gradient of nonlinear equality constraints
    gceq = [];
end
