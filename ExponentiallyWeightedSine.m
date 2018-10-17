%% Exponentially weighted sine wavelet and its Fourier transform %%
%%---------------------------------------------------------------%%

ampl = 100; % amplitude
fd = 1e6; % dominant frequency
td = 1/fd; % time delay 
T = 5e-6; % final time
fs = 1e9; % sampling frequency
dt = 1/fs; % sampling period (time step)
t = 0:dt:T; % time vector
N = length(t)-1; % number of time steps
% N = T/dt;
% t = linspace(0,T,N+1);

X = ampl*sin(2*pi*fd*t).*exp(-4*(fd*t-1).^2); % wavelet

figure
plot(t*1e6,X,'LineStyle','-','Color','b','LineWidth',0.8)
grid on
box on
xlim([0 T*1e6]);
set(gca,'FontSize',16)
xlabel('Time $t$ [$\mu$s]','Interpreter','latex')
ylabel('Amplitude $F(t)$','Interpreter','latex')
% title('Exponentially weighted sine wavelet in Time Domain')
% mymatlab2tikz('./','source_time.tex');

% n = 2^nextpow2(N+1); % number of frequency points
n = 2^20;

Y = dt*fftshift(fft(X,n)); % discrete Fourier transform (DFT)
% Y = dt*fft((-1).^(0:N).*X,n); 

P = abs(Y);

f = fs/n*(0:(n/2)); % frequency points
% f = 1/(n*dt)*(0:(n/2)); 

I = find(f<=5e6);

figure
% plot(f*1e-6,P(n/2:end),'LineStyle','-','Color','r','LineWidth',0.8)
plot(f(I)*1e-6,P(n/2:(n/2+length(I)-1)),'LineStyle','-','Color','r','LineWidth',0.8)
grid on
box on
xlim([0 5]);
set(gca,'FontSize',16)
xlabel('Frequency $f$ [MHz]','Interpreter','latex')
ylabel('Amplitude $|\hat{F}(f)|$','Interpreter','latex')
% title('Exponentially weighted sine wavelet in Frequency Domain')
% mymatlab2tikz('./','source_freq.tex');
