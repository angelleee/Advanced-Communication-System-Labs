clc; clear; close all;

fc = 1/8;
N = 128;
n = 1:N;

%lowpass
Nh = 32;
lowpass = 5;
Hf = [ones(1,lowpass) zeros(1,Nh+1-2*lowpass) ones(1,lowpass-1)];
h = ifftshift(ifft(Hf));
h = h*sqrt(Nh);
%figure;
%stem(1:Nh,h);
DC_gain = sum(h);
h = h/DC_gain;

%m1
tri = [0:64 63:-1:1];
tri = tri/64;
%m2
f = 1/128;
sin_ = sin(2*pi*f*n);
%m
signal = tri + 1i*sin_;

%carrier
carrier = exp(1i*2*pi*fc*n)*sqrt(2);
%mod
Tx = signal.*carrier;
Tx = real(Tx);
%demod
carrier_demod_0 = exp(-1i*(2*pi*fc*n))*sqrt(2);
Rx_0 = Tx.*carrier_demod_0;
Rxf_0 = fftshift(fft(Rx_0));
theta = deg2rad(-90);
carrier_demod = exp(-1i*(2*pi*fc*n-theta))*sqrt(2);
Rx = Tx.*carrier_demod;
Rxf = fftshift(fft(Rx));
figure;
hold on;
plot(-N/2:(N-1)/2,abs(Rxf_0));
plot(-N/2:(N-1)/2,abs(Rxf));
legend('phase shift = 0', 'phase shift = 5');
title('Amplitude Spectrum of demodulated signal');
hold off;
figure;
hold on;
plot(-N/2:(N-1)/2,angle(Rxf_0));
plot(-N/2:(N-1)/2,angle(Rxf));
legend('phase shift = 0', 'phase shift = 5');
title('Phase Spectrum of demodulated signal');
hold off;

%lowpass
out = conv(h,Rx); % lowpass
outf = fftshift(fft(out));
outn = length(outf);
figure;
plot(-outn/2:(outn-1)/2,abs(outf));
outi = real(out);
outq = (out-outi)/1j;

z1 = tri*cos(theta)-sin_*sin(theta);
z2 = tri*sin(theta)+sin_*cos(theta);

delay = 16;
figure;
hold on;
plot(1:N,tri);
plot(1:N,z1);
plot(1:N,outi(delay+1:delay+N));
legend('m_1(t)','z_1(t)','demod signal');
hold off;
figure;
hold on;
plot(1:N,sin_);
plot(1:N,z2);
plot(1:N,outq(delay+1:delay+N));
legend('m_2(t)','z_2(t)','demod signal');
hold off;