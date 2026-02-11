%% p1-1
clc; clear; close all;

Nb = 16;
b = [ones(1,1), zeros(1,Nb-1)];

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

channel = zeros(1, 1 + M1*M2);
channel(1 + M1*M2*0) = 1;
channel(1 + M1*M2*1) = 0.5;
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
stem(real(rxd1(1:Nb)) * sqrt(2));
stem(imag(rxd1(1:Nb)) * sqrt(2));
legend('real', 'imag');
title('equivalent baseband channel response');
hold off;
grid;

%% p1-2
Nb = 16;
b = [ones(1,1), zeros(1,Nb-1)];

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

channel = zeros(1, 1 + M1*M2*3.5);
channel(1 + M1*M2*0) = 1;
channel(1 + M1*M2*3.5) = 0.5;
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
stem(real(rxd1(1:Nb)) * sqrt(2));
stem(imag(rxd1(1:Nb)) * sqrt(2));
legend('real', 'imag');
title('equivalent baseband channel response');
hold off;
grid;

%% p2
clc; clear; close all;

Ne = 10; % H((Nh+Ne-1)*Ne)e(Ne*1) = delta(Nh+Ne-1)
n = 0:Ne-1;

% baseband channel response 1
h1 = [1, 0.5];
%Nh = length(h1);
%H1 = zeros(Nh+Ne-1, Ne);
%for i = 1:Ne
%    H1(i:i+Nh-1,i) = h1;
%end
%delta = zeros(Nh+Ne-1, 1);
%delta(1) = 1;
%e1 = H1\delta;
e1 = (-0.5).^n;

figure;
stem(e1);
title('equalizer 1');
grid;

ce1 = conv(h1, e1);
figure;
stem(ce1);
grid;

% baseband channel response 2
h2 = [0.5, 1];
%Nh = length(h2);
%H2 = zeros(Nh+Ne-1, Ne);
%for i = 1:Ne
%    H2(i:i+Nh-1,i) = h2;
%end
%delta = zeros(Nh+Ne-1, 1);
%delta(1+Ne) = 1;
%e2 = H2\delta;
e2 = -2 * (-2).^(n-Ne); %e2 = (-0.5).^(-(n-Ne)-1);

figure;
stem(e2);
title('equalizer 2');
grid;

ce2 = conv(h2, e2);
figure;
stem(ce2);
grid;

%% p3
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

channel = zeros(1, 1 + M1*M2);
channel(1 + M1*M2*0) = 1;
channel(1 + M1*M2*1) = 0.5;
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
rxd1 = conv(rxd1, e1);

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