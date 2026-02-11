clc; clear; close all;

% p1
data = load('Lab8_given_data\S8P1.mat');
b = data.b(:).';
Nb = length(b);

M = 4;
SRRC = srrcfun(0.3, M, 5);
L = length(SRRC);

bu1 = zeros(1, Nb*M);
bu1(1:M:end) = b;
bu1 = conv(bu1, SRRC);
bu2 = zeros(1, length(bu1)*M);
bu2(1:M:end) = bu1;

bu2f = fftshift(fft(bu2));
fn = length(bu2f);
FN = (-(fn-1)/2:(fn-1)/2)/fn;
figure;
plot(FN, abs(bu2f));

bu2 = conv(bu2, SRRC);

figure;
plot(bu2);
grid;

% p2
bd1_ = conv(bu2, SRRC);
bd1 = bd1_(L:M:end)*M; % 1+delay = 1+(L-1) = L
bd2_ = conv(bd1, SRRC);
bd2 = bd2_(L:M:end)*M; % 1+delay = 1+(L-1) = L

figure;
hold on;
plot(b);
plot(1:Nb, bd2(1:Nb));
hold off;
grid;

% p3
bu1 = zeros(1, Nb*M);
bu1(1:M:end) = b;
bu1 = conv(bu1, SRRC);
bu2 = zeros(1, length(bu1)*M);
bu2(1:M:end) = bu1;
bu2 = filter(IIR, bu2);
bd1_ = filter(IIR, bu2);
delay = 2*6;
bd1 = bd1_(1+delay:M:end)*M;
bd2_ = conv(bd1, SRRC);
bd2 = bd2_(L:M:end)*M; % 1+delay = 1+(L-1) = L

figure;
hold on;
plot(b);
plot(1:Nb, bd2(1:Nb));
hold off;
grid;