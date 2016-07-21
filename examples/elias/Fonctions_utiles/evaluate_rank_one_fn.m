function [feval]=evaluate_rank_one_fn(PC,lopt,rs)
if (isempty(getphi0(lopt)))
    feval=[];
else
    modelTestBasis=polyval(PC,rs);
    modelapprox=ones(size(rs,1),1)*lopt{0};
    for i=1:getnbgroups(PC)
        modelapprox=modelapprox.*(modelTestBasis{i}*lopt{i});
    end
    feval=modelapprox;
end