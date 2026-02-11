clc; close all;
data = load("Lab6_given_data\S6HW.mat");
x = data.x(:).';
N = length(x);

Xf = fftshift(fft(x));
Xfidx = ((-N/2):(N-1)/2)/N;
figure;
plot(Xfidx, abs(Xf));

% Downsample
M = floor((N/2)/60);
Nd = N/M;
xd = x(1:M:end);
Xdf = fftshift(fft(xd));
Xdfidx = ((-Nd/2):(Nd-1)/2)/Nd;
figure;
plot(Xdfidx, abs(Xdf));

% Upsample
xu = zeros(1,N);
xu(1:M:end) = xd;
Xuf = fftshift(fft(xu));
figure;
hold on;
plot(Xfidx, abs(Xuf));
plot(Xfidx, abs(Xf)/M);
legend('Upsampling', 'original*M');
hold off;

% IIR
delta = [ones(1,1) zeros(1,N-1)];
h = filter(IIR_hw, delta);
Hf = fftshift(fft(h));
figure;
hold on;
plot(Xfidx, abs(Xuf));
plot(Xfidx, 12*abs(Hf));
hold off;

z = filter(IIR_hw, xu);
z = z*M;
delay = 6;
figure;
hold on;
plot(1:(N-delay), x(1:end-delay));
plot(1:(N-delay), z(delay+1:end));
legend('x', 'z');
hold off;
e_IIR = 10*log10(mean((z(delay+1:end)-x(1:end-delay)).^2));
fprintf('e_IIR: %.4f dB\n', e_IIR);

% FIR
h = filter(FIR_hw, delta);
Hf = fftshift(fft(h));
figure;
hold on;
plot(Xfidx, abs(Xuf));
plot(Xfidx, 12*abs(Hf));
hold off;

z = filter(FIR_hw, xu);
z = z*M;
delay = 7;
figure;
hold on;
plot(1:(N-delay), x(1:end-delay));
plot(1:(N-delay), z(delay+1:end));
legend('x', 'z');
hold off;
e_FIR = 10*log10(mean((z(delay+1:end)-x(1:end-delay)).^2));
fprintf('e_FIR: %.4f dB\n', e_FIR);