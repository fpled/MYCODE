function rep = multiplicationpoint(u,v,PCphi)
% function rep = multiplication(u,v,PCphi)
% pour multiplier de PCTPMATRIXSUM u.*v

rep = PCTPMATRIX(PCphi,zeros(size(getphi0(u{1}),1),1));

for i=1:getm(u)
    for j=1:getm(v)
        rep = rep + u{i}.*v{j};
    end
end

return
