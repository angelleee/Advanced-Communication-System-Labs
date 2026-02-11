clc; clear; close all;

% find equivalent baseband ZF equalizer
h = [0.3, 1, 0.3];
%r = roots(h);
Ne = 10;
ne = 0:Ne-1;
e1 = -(15/4) * (-3).^(ne-Ne); %noncausal
e2 = (-5/12) * (-1/3).^ne; %causal
e = zeros(1, 2*Ne);
e(1:Ne) = e1;
e(Ne+1:2*Ne) = e2;

figure;
stem(e);
title('delayed e[n]');
grid;
eh = conv(e, h);
figure;
stem(eh);
title('conv(e,h)');
grid;

% Tx-Rx
%Nb = 16;
%b = [ones(1,1), zeros(1,Nb-1)];
data = load('S11P3.mat');
b = data.b(:).';
Nb = length(b);

M1 = 4;
M2 = 8;
SRRC = srrcfun(0.3, M1, 5);
L = length(SRRC);

bu1 = zeros(1, Nb*M1);
bu1(1:M1:end) = b;
bu1 = conv(bu1, SRRC);

bu2 = zeros(1, length(bu1)*M2);
bu2(1:M2:end) = bu1;
bu2 = filter(IIR_tx, bu2);

m = 1:length(bu2);
fc = 8e6;
fs = 32e6;
fm = fc/fs;
carrier = exp(1i*2*pi*fm*m);
tx = bu2.*carrier;
tx = real(tx) * sqrt(2);

channel = zeros(1, M1*M2*length(h));
for i = 1:length(h)
    channel(1 + M1*M2*(i-1)) = h(i);
end
figure;
stem(channel);
title('channel effect in passband');
grid;
tx = conv(tx, channel);

m = 1:length(tx);
carrier = exp(-1i*2*pi*fm*m);
rx = tx.*carrier;
rx = filter(IIR_tx, rx);
delay = 10*2;
rxd2 = rx(1+delay:M2:end) * M2;

rxd2 = conv(SRRC, rxd2);
delay = L-1;
rxd1 = rxd2(1+delay:M1:end) * M1;

figure;
hold on;
plot(b, '*');
plot(rxd1(1:Nb) * sqrt(2), '.');
legend('b', 'demodulated signal');
title('without equalizer');
hold off;
grid;

figure;
hold on;
plot(real(b));
plot(real(rxd1(1:Nb)) * sqrt(2));
legend('b', 'demodulated signal');
title('without equalizer real part');
hold off;
grid;

figure;
hold on;
plot(imag(b));
plot(imag(rxd1(1:Nb)) * sqrt(2));
legend('b', 'demodulated signal');
title('without equalizer imag part');
hold off;
grid;

% equalizer
rxd1 = conv(rxd1, e);
rxd1 = rxd1(1+Ne:end);

figure;
hold on;
plot(b, '*');
plot(rxd1(1:Nb) * sqrt(2), '.');
legend('b', 'demodulated signal');
title('with equalizer');
hold off;
grid;

figure;
hold on;
plot(real(b));
plot(real(rxd1(1:Nb)) * sqrt(2));
legend('b', 'demodulated signal');
title('with equalizer real part');
hold off;
grid;

figure;
hold on;
plot(imag(b));
plot(imag(rxd1(1:Nb)) * sqrt(2));
legend('b', 'demodulated signal');
title('with equalizer imag part');
hold off;
grid;