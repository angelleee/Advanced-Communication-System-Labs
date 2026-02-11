clc; clear; close all;

data = load('Lab10_given_data\S9P1.mat');
b = data.b(:).';
bi = b(1:2:end);
bq = b(2:2:end);
b = bi + 1j*bq;
Nb = length(b);

M1 = 16;
M2 = 2;
SRRC = srrcfun(0.3, M1, 5);
L = length(SRRC);

bu1 = zeros(1, Nb*M1);
bu1(1:M1:end) = b;
bu1 = conv(bu1, SRRC);

bu2 = zeros(1, length(bu1)*M2);
bu2(1:M2:end) = bu1;
bu2 = filter(hw_IIR_tx, bu2);

m = 0:length(bu2)-1;
fc = 8e6;
fs = 32e6;
%tx = bu2.*exp(1i*2*pi*(fc/fs)*m);
%tx = real(tx) * 2; % amplitude *2
g = 1.2;
phi = deg2rad(15);
tx = (real(bu2).*cos(2*pi*(fc/fs)*m) - g*imag(bu2).*sin(2*pi*(fc/fs)*m + phi)) * 2;

fif = 4e6;
rx = (tx.*cos(2*pi*((fc-fif)/fs)*m)) * 2; % amplitude *2
rx = filter(hw_IIR_rx, rx);
delay = 1 + 1;
rxd2 = rx(1+delay:M2:end) * M2; % amplitude *M2

m = 0:length(rxd2)-1;
fs = 32e6/M2;
rxd2 = rxd2.*exp(-1i*2*pi*(fif/fs)*m);
delay = 1;
rxd2 = rxd2*exp(-1i*2*pi*(fif/fs)*(delay/M2)); %delayphase

rxd2 = conv(SRRC, rxd2);
delay = L-1;
rxd1 = rxd2(1+delay:M1:end) * M1; % amplitude *M1

figure;
hold on;
plot(b);
plot(rxd1(1:Nb));
legend('b', 'IF demodulation');
title(sprintf('b part (g = %.2f, \\phi = %.1f^\\circ)', g, rad2deg(phi)));
hold off;

figure;
hold on;
stem(bi);
stem(real(rxd1(1:Nb)));
legend('b', 'IF demodulation');
title(sprintf('b real part (g = %.2f, \\phi = %.1f^\\circ)', g, rad2deg(phi)));
hold off;

figure;
hold on;
stem(bq);
stem(imag(rxd1(1:Nb)));
legend('b', 'IF demodulation');
title(sprintf('b imag part (g = %.2f, \\phi = %.1f^\\circ)', g, rad2deg(phi)));
hold off;

theoi = bi-bq*g*sin(phi);
theoq = bq*g*cos(phi);
theo = theoi + 1j*theoq;

figure;
hold on;
plot(theo);
plot(rxd1(1:Nb));
legend('theo', 'IF demodulation');
title(sprintf('theo (g = %.2f, \\phi = %.1f^\\circ)', g, rad2deg(phi)));
hold off;

figure;
hold on;
stem(theoi);
stem(real(rxd1(1:Nb)));
legend('theo', 'IF demodulation');
title(sprintf('theo real part (g = %.2f, \\phi = %.1f^\\circ)', g, rad2deg(phi)));
hold off;

figure;
hold on;
stem(theoq);
stem(imag(rxd1(1:Nb)));
legend('theo', 'IF demodulation');
title(sprintf('theo imag part (g = %.2f, \\phi = %.1f^\\circ)', g, rad2deg(phi)));
hold off;