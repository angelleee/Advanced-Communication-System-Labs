clc; clear; close all;

N = 1e4;
x = randn(1, N);
x1 = [0, x(1:end-1)];
x2 = [0, 0, x(1:end-2)];
h = [0.82, 3.1, 1.2, 0.53];

s1 = x * h(1);
figure;
histogram(s1, 32);
title('s1(n)');
grid on;

s2 = x1 * h(2);
figure;
histogram(s2, 32);
title('s2(n)');
grid on;

s3 = x1 * h(3);
figure;
histogram(s3, 32);
title('s3(n)');
grid on;

s4 = x2 * h(4);
figure;
histogram(s4, 32);
title('s4(n)');
grid on;

s5 = s3 + s4;
figure;
histogram(s5, 32);
title('s5(n)');
grid on;

s51 = [0, s5(1:end-1)];
s6 = s2 + s51;
figure;
histogram(s6, 32);
title('s6(n)');
grid on;

sout = s1 + s6;
figure;
histogram(sout, 32);
title('sout(n)');
grid on;

[xq, ~] = quantize(x, 8, 4);
xq1 = [0, xq(1:end-1)];
xq2 = [0, 0, xq(1:end-2)];
[hq, ~] = quantize(h, 8, 4);
[s1q, ~] = quantize(xq*hq(1), 8, 4);
[s2q, ~] = quantize(xq1*hq(2), 9, 16);
[s3q, ~] = quantize(xq1*hq(3), 8, 4);
[s4q, ~] = quantize(xq2*hq(4), 8, 2);
[s5q, ~] = quantize(s3q+s4q, 8, 4);
s5q1 = [0, s5q(1:end-1)];
[s6q, ~] = quantize(s2q+s5q1, 9, 16);
[soutq, ~] = quantize(s1q+s6q, 9, 16);

SQNR = mean(soutq.^2)/mean((soutq-sout).^2);
SQNR = 10*log10(SQNR);
figure;
plot(sout-soutq);
title(sprintf('Sout(n) - SoutQuan(n) --> SQNR = %f dB', SQNR));

% quantize function
function [xQuan, SQNR] = quantize(x, nb, dr)
    L = 2^nb;
    delta = 2*dr/L;
    xClip = min(max(x, -dr), dr);
    xQuan = delta * min(round(xClip/delta), L/2-1);

    SQNR = mean(x.^2)/mean((x-xQuan).^2);
    SQNR = 10*log10(SQNR);
end