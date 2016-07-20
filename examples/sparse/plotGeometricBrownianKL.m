function plotGeometricBrownianKL(Xref,Xpc,varargin)
% function plotGeometricBrownianKL(Xref,Xpc,varargin)
% Display geometric Brownian motion

p = ImprovedInputParser;
addParamValue(p,'legend',true,@islogical);
addParamValue(p,'label',true,@islogical);
addParamValue(p,'grid',true,@islogical);
addParamValue(p,'box',true,@islogical);
addParamValue(p,'FontSize',16,@isscalar);
addParamValue(p,'LineWidth',1,@isscalar);
addParamValue(p,'Interpreter','latex',@ischar);
parse(p,varargin{:})

figure('Name','Evolution of Brownian motion w.r.t. time')
% set(gcf,'Name','Evolution of Brownian motion w.r.t. time')
clf
t = linspace(0,1,size(Xref,1));
plot(t,[Xref,Xpc],'LineWidth',p.Results.LineWidth);
if p.Results.grid
    grid on
end
if p.Results.box
    box on
end
set(gca,'FontSize',p.Results.FontSize)
xlabel('Time')
ylabel('Brownian motion')
if p.Results.legend
    l = legend({'$B_t^{\mathrm{ref}}$','$B_t^{\mathrm{approx}}$'});
    set(l,'Interpreter',p.Results.Interpreter)
end

end
