function ind = indicator_diffusion(xi,yi,cx,cy,lc)
%function ind = indicator_diffusion(xi,yi,cx,cy,lc)
ind  = (((xi-cx)<lc).*((xi-cx)>-lc)).*(((yi-cy)<lc).*((yi-cy)>-lc));
end
