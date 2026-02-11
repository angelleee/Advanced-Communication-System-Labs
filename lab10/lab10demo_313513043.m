%% p1
clc; clear; close all;

data = load('Lab10_given_data\S10P1.mat');
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
%carrier = exp(1i*2*pi*fm*m);
%tx = bu2.*carrier;
%tx = real(tx) * 2;
g = 1.2;
phi = deg2rad(15);
tx = (real(bu2).*cos(2*pi*fm*m) - g*imag(bu2).*sin(2*pi*fm*m + phi)) * 2;

carrier = exp(-1i*2*pi*fm*m);
rx = tx.*carrier;
rx = filter(IIR_tx, rx);
delay = 10*2;
rxd2 = rx(1+delay:M2:end) * M2;

rxd2 = conv(SRRC, rxd2);
delay = L-1;
rxd1 = rxd2(1+delay:M1:end) * M1;

alpha = (1+g*exp(1j*phi))/2;
beta = (1-g*exp(1j*phi))/2;
theo = alpha * b + beta * conj(b);

figure;
hold on;
plot(b);
plot(rxd1(1:Nb));
legend('transmitted', 'received');
hold off;
grid;

figure;
hold on;
plot(theo);
plot(rxd1(1:Nb));
legend('theoretical', 'received');
hold off;
grid;

%% p2
clc; clear; close all;

data = load('Lab10_given_data\S10P2.mat');
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
%carrier = exp(1i*2*pi*fm*m);
%tx = bu2.*carrier;
%tx = real(tx) * 2;
g = 1.2;
phi = deg2rad(15);
tx = (real(bu2).*cos(2*pi*fm*m) - g*imag(bu2).*sin(2*pi*fm*m + phi)) * 2;

carrier = exp(-1i*2*pi*fm*m);
rx = tx.*carrier;
rx = filter(IIR_tx, rx);
delay = 10*2;
rxd2 = rx(1+delay:M2:end) * M2;

rxd2 = conv(SRRC, rxd2);
delay = L-1;
rxd1 = rxd2(1+delay:M1:end) * M1;

alpha = (1+g*exp(1j*phi))/2;
beta = (1-g*exp(1j*phi))/2;
theo = alpha * b + beta * conj(b);

figure;
hold on;
plot(b);
plot(rxd1(1:Nb));
hold off;
grid;

figure;
hold on;
plot(theo);
plot(rxd1(1:Nb));
hold off;
grid;

% IQ imbalance compensation
H = [1, -g*sin(phi);
    0, g*cos(phi)];
v = H \ [real(rxd1); imag(rxd1)];
rxd1_comI = v(1,:);
rxd1_comQ = v(2,:);
rxd1_com = rxd1_comI + 1j*rxd1_comQ;

figure;
hold on;
plot(real(b));
plot(real(rxd1_com(1:Nb)));
legend('transmitted', 'compensation');
title('real part');
hold off;
grid;

figure;
hold on;
plot(imag(b));
plot(imag(rxd1_com(1:Nb)));
legend('transmitted', 'compensation');
title('imag part');
hold off;
grid;

figure;
hold on;
plot(b);
plot(rxd1_com(1:Nb));
legend('transmitted', 'compensation');
hold off;
grid;

%% p3
clc; clear; close all;

data = load('Lab10_given_data\S10P2.mat');
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
tx = real(tx) * 2;
%g = 1.2;
%phi = deg2rad(15);
%tx = (real(bu2).*cos(2*pi*fm*m) - g*imag(bu2).*sin(2*pi*fm*m + phi)) * 2;

dw = 0.01; % CFO
carrier = exp(-1i*(2*pi*fm + dw)*m);
rx = tx.*carrier;
rx = filter(IIR_tx, rx);
delay = 10*2;
rxd2 = rx(1+delay:M2:end) * M2;

rxd2 = conv(SRRC, rxd2);
delay = L-1;
rxd1 = rxd2(1+delay:M1:end) * M1;

figure;
hold on;
plot(b);
scatter(real(rxd1(1:Nb)), imag(rxd1(1:Nb)));
hold off;
grid;

figure;
hold on;
plot(real(b));
plot(real(rxd1(1:Nb)));
title('real part');
hold off;
grid;

figure;
hold on;
plot(imag(b));
plot(imag(rxd1(1:Nb)));
title('imag part');
hold off;
grid;
