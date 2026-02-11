clc; clear; close all;

data = load('Lab7_given_data\S7HW.mat');
b = data.b(:).';
Nb = length(b);
figure;
stem(b);
title('bit stream');

%upsampling
M = 64;
bu = zeros(1,M*Nb);
bu(1:M:end) = b;
figure;
stem(bu);
title('upsampling');
rect_pulse = (1/M) * ones(1,M);
%bu_prac = conv(bu, rect_pulse, 'same');
bu_prac = bu;
figure;
stem(bu_prac * M);
title('rectangular filtering');

%SRRC filter
a = 0.3; M = 64; S = 5;
SRRC = srrcfun(a, M, S);
b_ps_prac = conv(bu_prac, SRRC, 'same');
figure;
plot(b_ps_prac * M);
title('TX SRRC filtering');
b_ps_prac_f = fftshift(fft(b_ps_prac));
figure;
fn = (-(M*Nb)/2:(M*Nb-1)/2) / (M*Nb);
plot(fn, abs(b_ps_prac_f));
title('TX SRRC filtering (frequency spectrum)');

%up-conversion
R = 1e6; fc = 8e6;
fs = fc/(R*M); % fc = sampling_rate*fs = R*M*fs;
n = 1:M*Nb;
b_ps_prac_up = b_ps_prac.*exp(1j*2*pi*fs*n);
b_ps_prac_up = real(b_ps_prac_up);
b_ps_prac_up_f = fftshift(fft(b_ps_prac_up));
figure;
plot(fn, abs(b_ps_prac_up_f));
title('TX up-conversion (frequency spectrum)');

%RX
%down-conversion
b_ps_prac_down = b_ps_prac_up.*exp(-1j*2*pi*fs*n);
b_ps_prac_down_f = fftshift(fft(b_ps_prac_down));
b_ps_prac_down_f = b_ps_prac_down_f.*[zeros(1, (M*Nb-1024)/2), ones(1,1024), zeros(1, (M*Nb-1024)/2)];
figure;
plot(fn, abs(b_ps_prac_down_f));
title('RX down-conversion (frequency spectrum)');
b_ps_prac_down = ifft(ifftshift(b_ps_prac_down_f));
b_ps_prac_down = real(b_ps_prac_down);

%SRRC filter
b_ps_prac_rec = conv(b_ps_prac_down, SRRC, 'same');
figure;
hold on;
plot(b_ps_prac_rec * M * 2);
plot(b_ps_prac * M);
legend('b-ps-rec', 'b-ps');
title('RX SRRC filtering');
hold off;

%downsampling
b_rec = b_ps_prac_rec(1:M:end);
figure;
hold on;
stem(b_rec * M * 2);
stem(b);
legend('b-rec', 'b');
title('RX downsampling');
hold off;