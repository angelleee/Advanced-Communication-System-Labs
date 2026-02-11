clc; clear; close all;

data = load('Lab9_given_data\S9HWs.mat');
b = data.b(:).';
Nb = length(b);

M1 = 16;
M2 = 4;
SRRC = srrcfun(0.3, M1, 5);
L = length(SRRC);

bu1 = zeros(1, Nb*M1);
bu1(1:M1:end) = b;
bu1 = conv(SRRC, bu1);

bu2 = zeros(1, length(bu1)*M2);
bu2(1:M2:end) = bu1;
bu2 = filter(LPFTxhw, bu2);

fc = 16e6;
fif = 10e6;
fs = 1e6*M1*M2; % symbol rate * M1 * M2
Tx = bu2.*exp(1i*2*pi*(fc/fs)*(1:length(bu2)));
Tx = real(Tx) * 2; % amplitude *2

fi = 2*fif - fc;
I = cos(2*pi*(fi/fs)*(1:length(Tx))) * max(abs(Tx));
scaler = max(abs(fftshift(fft(Tx)))) / max(abs(fftshift(fft(I))));
Tx = Tx + I*scaler;

figure;
hold on;
bf = fftshift(fft(Tx));
Nf = length(bf);
nf = (-Nf/2:Nf/2-1)/Nf;
plot(nf, abs(bf)/max(abs(bf)));
Hd = BPFhw();
[B,A] = tf(Hd);
[H,w] = freqz(B,A,'whole');  % whole: 0~2pi
H = fftshift(H);
w = (w/(2*pi)) - 0.5; % 0~1 to -0.5~0.5
plot(w, abs(H)/max(abs(H)));
legend('Tx','BPF');
ylabel('normalized magnitude');
xlabel('normalized frequency');
hold off;

Tx = filter(BPFhw, Tx); % image rejection filter

fs = 1e6*M1*M2; % symbol rate * M1 * M2
Rx = Tx.*cos(2*pi*((fc-fif)/fs)*(1:length(Tx))) *2; % amplitude *2

Rx = filter(LPFRxhw, Rx);
delay = 1 + 6 + 2; % LPF(Tx) + BPF + LPF(Rx)
bd2 = Rx(1+delay:M2:end) * M2; % amplitude *M2

fs = 1e6*M1; % symbol rate * M1
bd2 = bd2.*exp(-1i*2*pi*(fif/fs)*(1:length(bd2)));
bd2 = bd2*exp(-1i*2*pi*(fif/fs)*(delay/M2)); % phase delay

bd2 = conv(SRRC, bd2);
delay = L-1; % SRRC(Tx) + SRRC(Rx) = ((L-1)/2)*2
bd1 = bd2(1+delay:M1:end) * M1; % amplitude *M1
bd1 = real(bd1) .* (abs(bd1)./abs(real(bd1))); % amplitude *abs(bd1)./abs(real(bd1))

figure;
hold on;
plot(b);
plot(bd1(1:Nb));
hold off;