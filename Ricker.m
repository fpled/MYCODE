%% Ricker wavelet and its Fourier transform %%
%%------------------------------------------%%

freq = 1/3; % dominant frequency
td = 1/freq; % time delay 
T = 12; % final time
fs = 1e3; % sampling frequency
dt = 1/fs; % sampling period (time step)
t = 0:dt:T; % time vector
N = length(t)-1; % number of time steps
% N = T/dt;
% t = linspace(0,T,N+1);

X = (2*(pi*freq)^2*(t-td).^2 - 1).*exp(-(pi*freq)^2*(t-td).^2); % Ricker wavelet

figure
plot(t,X,'LineStyle','-','Color','b','LineWidth',1)
grid on
box on
set(gca,'FontSize',16)
xlabel('Time $t$ (s)','Interpreter','latex')
ylabel('Amplitude $r(t)$','Interpreter','latex')
title('Ricker wavelet in Time Domain')

% n = 2^nextpow2(N+1); % number of frequency points
n = 2^20;

Y = dt*fftshift(fft(X,n)); % discrete Fourier transform (DFT)
% Y = dt*fft((-1).^(0:N).*X,n); 

P = abs(Y);

f = fs/n*(0:(n/2)); % frequency points
% f = 1/(n*dt)*(0:(n/2)); 

figure
plot(f,P(n/2:end),'LineStyle','-','Color','r','LineWidth',1)
grid on
box on
xlim([0 2]);
set(gca,'FontSize',16)
xlabel('Frequency $f$ (Hz)','Interpreter','latex')
ylabel('Amplitude $|\hat{r}(f)|$','Interpreter','latex')
title('Ricker wavelet in Frequency Domain')
