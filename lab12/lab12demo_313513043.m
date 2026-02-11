clc; clear; close all;

data = load('S12P2.mat');
b = data.b(:).';
Nb = length(b);

M1 = 16;
M2 = 4;
BT = 0.5;
fs = 1e6;
fd = 1.5e5;
fIF = 2e6;

gaufilter = gaufilter(M1, BT, 1);
Lg = length(gaufilter);

SRRC = srrcfun(0.3, 4, 1);
Ls = length(SRRC);

rect = (1/M1)*ones(1,M1);
Lr = length(rect);

bu1 = zeros(1, Nb*M1);
bu1(1:M1:end) = b;
bu1 = conv(bu1, rect);
bu1 = conv(gaufilter, bu1);

sigma = cumsum(bu1);
x = exp(1j*2*pi*(fd/(fs*M1))*sigma);
x = x .* exp(1j*2*pi*(fIF/(fs*M1))*(1:length(x)));
x = real(x) * sqrt(2);

bu2 = zeros(1, length(x)*M2);
bu2(1:M2:end) = x;
bu2 = filter(DMA, bu2);

% ignoring the analog carriers

bd2 = filter(DMA, bu2);
delay = 6 + 6;
bd2 = bd2(1+delay:M2:end) * M2;

bd2 = bd2 .* exp(-1j*2*pi*(fIF/(fs*M1))*(1:length(bd2)));
bd2 = conv(SRRC, bd2);

sigma = angle(bd2);
sigma = sigma / (2*pi*fd*(1/(M1*fs)));
y = [sigma(1), diff(sigma)];

bd1 = conv(rect, y);
figure;
hold on;
plot(bu1(floor((length(bu1)-Nb*M1)/2):floor((length(bu1)+Nb*M1)/2)));
plot(bd1(floor((length(bd1)-Nb*M1)/2):floor((length(bd1)+Nb*M1)/2)));
hold off;
grid;

bd1 = conv(gaufilter, bd1);
figure;
hold on;
plot(bu1(floor((length(bu1)-Nb*M1)/2):floor((length(bu1)+Nb*M1)/2)));
plot(bd1(floor((length(bd1)-Nb*M1)/2):floor((length(bd1)+Nb*M1)/2)));
hold off;
grid;

delay = (Lr-1) + (Lg-1) + (Ls-1)/2;
bd1 = bd1(1+delay:M1:end) * M1;

figure;
hold on;
stem(b);
stem(bd1(1:Nb));
legend('b', 'demod');
hold off;
grid;