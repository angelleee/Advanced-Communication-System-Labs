%% p1
clc; clear; close all;

N = 128;
tri = [0:64 63:-1:1];
n = 1:N;
f = 1/8;
SIN = cos(2*pi*f*n);
Txsignal = SIN.*tri; % mod

figure;
hold on;
plot(1:N,tri);
plot(1:N,Txsignal);
hold off;

figure;
plot(1:N,abs(Txsignal));

Rxsignal = Txsignal.*SIN; % demod

Nh = 32;
lowpass = 5;
Hf = [ones(1,lowpass) zeros(1,Nh+1-2*lowpass) ones(1,lowpass-1)];
h = ifftshift(ifft(Hf));
h = h*sqrt(Nh);
%figure;
%stem(1:Nh,h);
DC_gain = sum(h);
h = h/DC_gain;
out = conv(h,Rxsignal); % lowpass
out = out*2;

delay = 16;
figure;
hold on;
plot(1:N,tri);
plot(1:N,out(delay+1:delay+N));
hold off;

%% p2
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

%I-branch
tri = [0:64 63:-1:1];
tri = tri/64;
SINi = cos(2*pi*fc*n)*sqrt(2);
Txi = SINi.*tri; % mod
Rxi = SINi.*Txi; % demod
outi = conv(h,Rxi); % lowpass

%Q-branch
f = 1/128;
sin_ = sin(2*pi*f*n);
SINq = -sin(2*pi*fc*n)*sqrt(2);
Txq = SINq.*sin_; % mod
Rxq = SINq.*Txq; % demod
outq = conv(h,Rxq); % lowpass

delay = 16;
figure;
hold on;
plot(1:N,tri);
plot(1:N,outi(delay+1:delay+N));
hold off;
figure;
hold on;
plot(1:N,sin_);
plot(1:N,outq(delay+1:delay+N));
hold off;

%% p3
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
carrier_demod = exp(-1i*2*pi*fc*n)*sqrt(2);
Rx = Tx.*carrier_demod;
Rxf = fftshift(fft(Rx));
figure;
plot(-N/2:(N-1)/2,abs(Rxf));

%lowpass
out = conv(h,Rx); % lowpass
outf = fftshift(fft(out));
outn = length(outf);
figure;
plot(-outn/2:(outn-1)/2,abs(outf));
outi = real(out);
outq = (out-outi)/1j;

delay = 16;
figure;
hold on;
plot(1:N,tri);
plot(1:N,outi(delay+1:delay+N));
hold off;
figure;
hold on;
plot(1:N,sin_);
plot(1:N,outq(delay+1:delay+N));
hold off;