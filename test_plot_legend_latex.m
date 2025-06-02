x = (1:10)';
y = [95.4843 127.597 82.333 187.835 136.958 198.104 175.525 227.794 250.54 273.702]';
figure
semilogy(x,y,'LineStyle','-','Color',getfacecolor(3),'LineWidth',1);
grid on
box on
set(gca,'FontSize',16)
xlabel('$m$','Interpreter','latex')
legend('$E^2_{\mathrm{CRE}}$','Interpreter','latex')
mymatlab2tikz('.','toto.tex');
