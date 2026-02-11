%% p1
clear; close all; clc;

data = load('./Lab3_given_data/S3P1.mat');
x = data.x(:).';
h = [1 2 3 4 -2 -1];
Nx = length(x);
Nh = length(h);
Ny = Nx + Nh - 1;
y = zeros(1,Ny);

for n = 1:Nh
    y = y + h(n)*[zeros(1,(n-1)), x, zeros(1,(Nh-1-(n-1)))];
end
%disp(y);

answer = conv(x,h);
disp(all(abs(answer - y) < 1e-6));

%% p2
clear; close all; clc;

data = load('./Lab3_given_data/S3P2.mat');
x = data.x(:).';
x = [zeros(1,2),x];
N = length(x);
y = zeros(1,N);

for n =3:N
    y(n) = 0.5*y(n-1) + 0.1*y(n-2) + x(n) + 2*x(n-1) + 3*x(n-2);
end
y = y(1,3:N);
%disp(y);

b = [1 2 3];
a = [1 -0.5 -0.1];
answer = filter(b,a,x);
answer = answer(1,3:N);
disp(all(abs(answer - y) < 1e-6));

%% p3
clear; close all; clc;

N = 128;
n = 1:128;
x_low = cos(2*pi*n*0.03);
x_high = cos(2*pi*n*0.3);
x = x_low + x_high;

zeros_pos = [-1, exp(1j*31*pi/32), exp(-1j*31*pi/32), exp(1j*30*pi/32), exp(-1j*30*pi/32)];
poles_pos = [exp(1j*5*pi/16), exp(-1j*5*pi/16), exp(1j*6*pi/16), exp(-1j*6*pi/16)];
poles_pos = 0.5*poles_pos;

b = poly(zeros_pos);
a = poly(poles_pos);
H0 = sum(b)/sum(a);
b = b/H0;

[h,w] = freqz(b,a,128);
f = w/(2*pi);
figure;
plot(f, abs(h));
xlabel('Normalized Frequency (0~0.5)');
ylabel('|H(f)|');
title('Frequency Response of the Filter');
grid on;
figure;
plot(w, angle(h));
xlabel('w');
ylabel('Phase (radians)');
title('Phase Response of the Filter');
grid on;

figure;
zplane(b,a);
title('Pole-Zero Plot of the Filter');

delta = [1, zeros(1,N-1)];
h = filter(b,a,delta);
figure;
stem(1:16,h(1:16)); %impulse response length never mind
title('Impulse Response of the Filter');
xlabel('n');
ylabel('h[n]');

y = filter(b,a,x);
x_draw = x(1:125);
y_draw = y(4:128);
figure;
hold on;
plot(0:124,x_draw);
plot(0:124,y_draw);
xlabel('n');
hold off;

x_low_draw = x_low(1:125);
figure;
hold on;
plot(0:124,x_low_draw);
plot(0:124,y_draw);
xlabel('n');
hold off;

%% p4
clear; close all; clc;

N = 128;
n = 1:128;
x_low = cos(2*pi*n*0.03);
x_high = cos(2*pi*n*0.3);
x = x_low + x_high;

FFT_size = 16;
ones_n = 3;
H_f = [ones(1,ones_n), zeros(1,FFT_size-2*ones_n+1), ones(1,ones_n-1)];
h = ifft(H_f);
h = ifftshift(h); % i dont know
h = h * sqrt(FFT_size);
DC_gain = sum(h);
h = h / DC_gain;
disp("gain: " + num2str(sum(h)));

y = conv(x,h);
figure;
stem(1:16,h);
xlabel('n'); ylabel('h[n]');
title('Impulse Response of the DFT-based Lowpass Filter');

delay=9;
figure;
hold on;
plot(n, x_low);
plot(n, y(delay+1:delay+N));
xlabel('n'); ylabel('Amplitude');
title('Original vs Filtered Signal');
legend('Original','Filtered');
hold off;