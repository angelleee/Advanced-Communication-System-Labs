clc; clear; close all;

% p1
data = load('S14P1.mat');
x1 = data.x1(:).';
x2 = data.x2(:).';
x3 = data.x3(:).';

nb = 4;
[x1q, SQNR1] = quantize(x1, nb, 1);
[x2q, SQNR2] = quantize(x2, nb, 1);
fprintf('SQNR_rampup = %f dB.\n', SQNR1);
fprintf('SQNR_uniform = %f dB.\n', SQNR2);
theoSQNR = (1/3)/((1/3)*(1/2)^(2*nb)); % 1/3 = E[x^2], x~uniform(-1,1)
theoSQNR = 10*log10(theoSQNR);
fprintf('SQNR_theo = %f dB.\n', theoSQNR);

[x3q, SQNR3] = quantize(x3, 7, 3);
fprintf('SQNR_gaussian = %f dB.\n', SQNR3);

figure;
hold on;
plot(x1);
plot(x1q);
hold off;
grid;

figure;
hold on;
plot(x2);
plot(x2q);
hold off;
grid;

figure;
hold on;
plot(x3);
plot(x3q);
hold off;
grid;

% p2
N = 1e4;
x = randn(1, N);
h = [0.82, 3.1, 1.2, 0.53];
out = zeros(1, N);
figure;
histogram(x, 32);
title(sprintf('sin(n)'));
grid on;
for i = 1:4
    si = [zeros(1, 4-i), x(1:end-(4-i))] * h(5-i);
    figure;
    histogram(si, 32);
    title(sprintf('s%d(n)', 5-i));
    grid on;

    out = out + si;
    if(i ~= 1)
        figure;
        histogram(out, 32);
        title(sprintf('s%d(n)', i+3));
        grid on;
    end
end

% p3
nb = 8;
[xq, ~] = quantize(x, nb, 4);
xq1 = [0, xq(1:end-1)];
xq2 = [0, 0, xq(1:end-2)];
xq3 = [0, 0, 0, xq(1:end-3)];
[hq, ~] = quantize(h, nb, 4);
outq = zeros(1, N);

[s1q, ~] = quantize(xq*hq(1), nb, 4);
[s2q, ~] = quantize(xq1*hq(2), nb, 16);
[s3q, ~] = quantize(xq2*hq(3), nb, 8);
[s4q, ~] = quantize(xq3*hq(4), nb, 2);

[s5q, ~] = quantize(s3q+s4q, nb, 4);
[s6q, ~] = quantize(s2q+s5q, nb, 16);
[soutq, soutSQNR] = quantize(s1q+s6q, nb, 16);
fprintf('SQNR_Problem3 = %f dB.\n', soutSQNR);

figure;
plot(out-soutq);

% quantize function
function [xQuan, SQNR] = quantize(x, nb, dr)
    L = 2^nb;
    delta = 2*dr/L;
    xClip = min(max(x, -dr), dr);
    xQuan = delta * min(round(xClip/delta), L/2-1);

    SQNR = mean(x.^2)/mean((x-xQuan).^2);
    SQNR = 10*log10(SQNR);
end