%% p1
N = 128;
n = 1:N;
f = 1/20;
x = cos(2*pi*f*n);
M = 4;

Xf = fftshift(fft(x));
Xf_idx = (-N/2:(N-1)/2)/N;
figure;
stem(Xf_idx, abs(Xf));

% Downsample
xd = x(1:M:end);
Xdf = fftshift(fft(xd));
Nd = N/M;
Xdf_idx = (-Nd/2:(Nd-1)/2)/Nd;
figure;
stem(Xdf_idx, abs(Xdf));

% Upsample
xu = zeros(1,N);
xu(1:M:end) = xd;
figure;
hold on;
%stem(n, x);
stem(n, xu);
%legend('x', 'xu');
hold off;
Xuf = fftshift(fft(xu));
figure;
stem(Xf_idx, abs(Xuf));

%% p2
%FIR
delta = [ones(1,1) zeros(1,N-1)];
h = filter(LPF, delta);
Hf = fftshift(fft(h));
figure;
hold on;
plot(Xf_idx, abs(Xuf));
plot(Xf_idx, 12*abs(Hf));
hold off;

z = filter(LPF, xu);
z = z*M;
delay = 8;
figure;
hold on;
plot(n, x);
plot(1:(N-delay), z(delay+1:end));
legend('x', 'z');
hold off;
e_FIR = 10*log10(mean((z(delay+1:end)-x(1:end-delay)).^2));

%IIR
delta = [ones(1,1) zeros(1,N-1)];
h = filter(IIR, delta);
Hf = fftshift(fft(h));
figure;
hold on;
plot(Xf_idx, abs(Xuf));
plot(Xf_idx, 12*abs(Hf));
hold off;

z = filter(IIR, xu);
z = z*M;
delay = 6;
figure;
hold on;
plot(n, x);
plot(1:(N-delay), z(delay+1:end));
legend('x', 'z');
hold off;
e_IIR = 10*log10(mean((z(delay+1:end)-x(1:end-delay)).^2));

fprintf('MSE_FIR(dB) = %.4f\n', e_FIR);
fprintf('MSE_IIR(dB) = %.4f\n', e_IIR);

%% p3
data = load('./Lab6_given_data/S6P3.mat');
x = data.x(:).';
N = length(x);
figure;
plot(1:N, x);

% Upsample/DAC
M = 32;
NM = N*M;
xu = zeros(1,NM);
xu(1:M:end) = x;

Xuf = fftshift(fft(xu));
figure;
Xuf_idx = (-NM/2:(NM-1)/2)/NM;
plot(Xuf_idx, abs(Xuf));

y = filter(DAC, xu);
delay = 2;
figure;
plot(1:(NM-delay), y(delay+1:end));

% Downsample/ADC
xd = y(1:M:end)*M;
figure;
hold on;
plot(1:(N-delay), xd(delay+1:end));
plot(1:N, x);
legend('ADC', 'original')
hold off;
