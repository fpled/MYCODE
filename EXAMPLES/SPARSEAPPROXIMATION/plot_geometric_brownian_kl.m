function plot_geometric_brownian_kl(Xref,Xpc,varargin)
% function geometric_brownian_kl(Xref,Xpc,varargin)
% Display geometric Brownian motion

nolegend = ischarin('nolegend',varargin);
nogrid = ischarin('nogrid',varargin);
nobox = ischarin('nobox',varargin);
fontsize = getcharin('FontSize',varargin,16);
linewidth = getcharin('LineWidth',varargin,1);
interpreter = getcharin('Interpreter',varargin,'latex');

figure('Name','Evolution of Brownian motion w.r.t. time')
% set(gcf,'Name','Evolution of Brownian motion w.r.t. time')
clf
t = linspace(0,1,size(Xref,1));
plot(t,[Xref,Xpc],'LineWidth',linewidth);
if ~nogrid
    grid on
end
if ~nobox
    box on
end
set(gca,'FontSize',fontsize)
xlabel('Time')
ylabel('Brownian motion')
if ~nolegend
    l = legend({'$B_t^{\mathrm{ref}}$','$B_t^{\mathrm{approx}}$'});
    set(l,'Interpreter',interpreter)
end

end
