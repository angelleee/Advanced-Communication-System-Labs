clc; clear; close all;
data = load('Lab8_given_data\S8HW.mat');
x = data.x(:).';
Nx = length(x);

xf = fftshift(fft(x));
N = length(xf);
Nf = (-(N-1)/2:(N-1)/2)/N;
figure;
stem(Nf, abs(xf));

M1 = 4;
M2 = 8;
SRRC = srrcfun(0.3, M1, 5);
L = length(SRRC);

xu = zeros(1, M1*Nx);
xu(1:M1:end) = x;

xuf = fftshift(fft(xu));
N = length(xuf);
Nf = (-(N-1)/2:(N-1)/2)/N;
figure;
hold on;
stem(Nf, abs(xuf)/max(abs(xuf)));
NSRRC = length(SRRC);
NSRRCf = (-(NSRRC-1)/2:(NSRRC-1)/2)/NSRRC;
plot(NSRRCf, abs(SRRC/max(SRRC)));
ylabel('Normalized Magnitude');
legend('upsampling(M1=4)', 'SRRC');
hold off;
grid;

xu_srrc = conv(xu, SRRC);

xu_srrcf = fftshift(fft(xu_srrc));
N = length(xu_srrcf);
Nf = (-(N-1)/2:(N-1)/2)/N;
figure;
stem(Nf, abs(xu_srrcf));

xu_srrc_u = zeros(1, M2*length(xu_srrc));
xu_srrc_u(1:M2:end) = xu_srrc;

delta = [ones(1,1), zeros(1, length(xu_srrc_u)-1)];
h = filter(IIR_hw, delta);
Hf = fftshift(fft(h));
xu_srrc_u_f = fftshift(fft(xu_srrc_u));
N = length(xu_srrc_u_f);
Nf = (-(N-1)/2:(N-1)/2)/N;
figure;
hold on;
stem(Nf, abs(xu_srrc_u_f)/max(abs(xu_srrc_u_f)));
plot(Nf, abs(Hf));
ylabel('Normalized Magnitude');
legend('upsampling(M2=8)', 'IIR');
hold off;

Tx = filter(IIR_hw, xu_srrc_u);

Txf = fftshift(fft(Tx));
N = length(Txf);
Nf = (-(N-1)/2:(N-1)/2)/N;
figure;
stem(Nf, abs(Txf));

Rx = filter(IIR_hw, Tx);
delay = 2*8;
Rx_d2 = Rx(1+delay:M2:end)*M2;

Rx_srrc = conv(Rx_d2, SRRC);
Rx_d1 = Rx_srrc(L:M1:end)*M1;

figure;
hold on;
plot(x);
plot(1:Nx, Rx_d1(1:Nx));
hold off;
grid;