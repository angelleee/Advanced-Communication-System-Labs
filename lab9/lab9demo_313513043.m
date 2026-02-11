%% p1
clc; clear; close all;

data = load('Lab9_given_data\S9P1.mat');
b = data.b(:).';
Nb = length(b);

M1 = 16;
M2 = 2;
SRRC = srrcfun(0.3, M1, 5);
L = length(SRRC);

bu1 = zeros(1, Nb*M1);
bu1(1:M1:end) = b;
bu1 = conv(bu1, SRRC);

bu2 = zeros(1, (length(bu1)-(L-1)/2)*M2);
bu2(1:M2:end) = bu1(1+(L-1)/2:end);
%bu2 = zeros(1, length(bu1)*M2);
%bu2(1:M2:end) = bu1(1:end);
bu2 = filter(IIR_tx, bu2);

m = 1:length(bu2);
fc = 8e6;
fs = 32e6;
carrier = exp(1i*2*pi*(fc/fs)*m);
tx = bu2.*carrier;
tx = real(tx);

carrier = exp(-1i*2*pi*(fc/fs)*m);
rx = tx.*carrier;
%rx = real(rx);
rx = filter(IIR_rx, rx);
delay = 1 + 19;
%delay = (L-1)/2*2 + 1 + 19;
rxd = rx(1+delay:M1*M2:end) * M1*M2;

figure;
hold on;
plot(b);
plot(real(rxd(1:Nb))*2);
hold off;
grid;

%% p2

bu1 = zeros(1, Nb*M1);
bu1(1:M1:end) = b;
bu1 = conv(bu1, SRRC);

bu2 = zeros(1, length(bu1)*M2);
bu2(1:M2:end) = bu1;
bu2 = filter(IIR_tx, bu2);

m = 1:length(bu2);
fc = 8e6;
fs = 32e6;
carrier = exp(1i*2*pi*(fc/fs)*m);
tx = bu2.*carrier;
tx = real(tx);

carrier = exp(-1i*2*pi*(fc/fs)*m);
rx = tx.*carrier;
%rx = real(rx);
rx = filter(IIR_tx, rx);
delay = 1 + 1;
rxd2 = rx(1+delay:M2:end) * M2;

rxd2 = conv(SRRC, rxd2);
delay = L-1;
rxd1 = rxd2(1+delay:M1:end) * M1;

figure;
hold on;
plot(b);
plot(real(rxd1(1:Nb))*2);
hold off;
grid;

%% p3

bu1 = zeros(1, Nb*M1);
bu1(1:M1:end) = b;
bu1 = conv(bu1, SRRC);

bu2 = zeros(1, length(bu1)*M2);
bu2(1:M2:end) = bu1;
bu2 = filter(IIR_tx, bu2);

m = 1:length(bu2);
fc = 8e6;
fs = 32e6;
carrier = exp(1i*2*pi*(fc/fs)*m);
tx = bu2.*carrier;
tx = real(tx) * 2; % amplitude *2

txf = fftshift(fft(tx));
N = length(txf);
Nf = (-(N-1)/2:(N-1)/2)/N;
figure;
stem(Nf, abs(txf));
title('tx');

fif = 4e6;
IF1 = cos(2*pi*((fc-fif)/fs)*m) * 2; % amplitude *2
rx = tx.*IF1;

rxf = fftshift(fft(rx));
N = length(rxf);
Nf = (-(N-1)/2:(N-1)/2)/N;
figure;
stem(Nf, abs(rxf));
title('rx*cos');

rx = filter(IIR_rx3, rx);
delay = 1 + 1;
rxd2 = rx(1+delay:M2:end) * M2; % amplitude *M2

rxd2f = fftshift(fft(rxd2));
N = length(rxd2f);
Nf = (-(N-1)/2:(N-1)/2)/N;
figure;
stem(Nf, abs(rxd2f));
title('rx*cos down');

m = 1:length(rxd2);
fs = 32e6/M2;
IF2 = exp(-1i*2*pi*(fif/fs)*m);
delayphase = exp(-1i*2*pi*(fif/fs)*(1/M2)); % delay of Rx DMA filter = 1
%delayphase = exp(-1i*(1/4)*pi);
rxd2 = rxd2.*IF2*delayphase;

rxd2f = fftshift(fft(rxd2));
N = length(rxd2f);
Nf = (-(N-1)/2:(N-1)/2)/N;
figure;
stem(Nf, abs(rxd2f));
title('rx*cos down *IF');

rxd2 = conv(SRRC, rxd2);

rxd2f = fftshift(fft(rxd2));
N = length(rxd2f);
Nf = (-(N-1)/2:(N-1)/2)/N;
figure;
stem(Nf, abs(rxd2f));
title('rx*cos down *IF srrc');

delay = L-1;
rxd1 = rxd2(1+delay:M1:end) * M1; % amplitude *M1
rxd1_draw = real(rxd1) .* (abs(rxd1)./abs(real(rxd1))); % amplitude *abs(rxd1)./abs(real(rxd1))

rxd1f = fftshift(fft(rxd1));
N = length(rxd1f);
Nf = (-(N-1)/2:(N-1)/2)/N;
figure;
stem(Nf, abs(rxd1f));
title('rx*cos down *IF srrc down');

figure;
hold on;
plot(b);
plot(rxd1_draw(1:Nb));
hold off;
grid;