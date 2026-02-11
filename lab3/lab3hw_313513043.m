clear; close all; clc;

data = load('./Lab3_given_data/S3HW.mat');
x = data.x(:).';
y = data.y(:).';

figure;
hold on;
Ny = length(y);
Y = fftshift(fft(y));
stem(-Ny/2:(Ny-1)/2,abs(Y));
Nx = length(y);
X = fftshift(fft(x));
stem(-Nx/2:(Nx-1)/2,abs(X));
xlabel('k');
ylabel('Magnitude');
legend('Y','X');
title('X and Y');
hold off;

Nh = 16;
zeros_n = Nh/4;
H_f = [zeros(1,zeros_n), ones(1,Nh-2*zeros_n+1), zeros(1,zeros_n-1)];
h = ifftshift(ifft(H_f));
figure;
stem(1:Nh,H_f);
xlabel('k'); ylabel('H[k]');
title('Frequency Response of the DFT-based Highpass Filter');
figure;
stem(1:Nh,h);
xlabel('n'); ylabel('h[n]');
title('Impulse Response of the DFT-based Highpass Filter');

out = conv(y,h);
out = out(Nh:Ny);
figure;
hold on;
plot(x((Nh/2):(Ny-Nh/2+1)));
plot(out);
legend('x','out');
xlabel('n');
ylabel('Magnitude');
title('x[n] and out[n]');
hold off;

OUT = fftshift(fft(out));
Nout = length(out);
f_out = (-Nout/2:(Nout-1)/2) / Nout;
f_x = (-Nx/2:(Nx-1)/2) / Nx;
figure;
hold on;
stem(f_out,abs(OUT));
stem(f_x,abs(X));
xlabel('Normalized frequency (k/N)');
ylabel('Magnitude');
legend('OUT','X');
hold off;