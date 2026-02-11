clc; clear; close all;
%% p1
% raised cosine
a = 0.3;
M = 4;
S = 5;
RC = zeros(1,2*M*S+1);
for n = -M*S:M*S
    idx = n+M*S+1;
    RC(idx) = sinc(n/M)*(cos(pi*a*(n/M)) / (1 - 4*(a*n/M)^2));
end

DC_gain = sum(RC);
RC = RC/DC_gain;

figure;
stem(RC);
grid;

%% p2
% squared raised cosine
a = 0.3;
M = 4;
S = 5;
SRRC = zeros(1,2*M*S+1);
for n = -M*S:M*S
    idx = n+M*S+1;
    if n ~= 0
        SRRC(idx) = (4*a/pi) * (cos((1+a)*pi*n/M) + M*sin((1-a)*pi*n/M)/(4*a*n)) / (1-(4*a*n/M)^2);
    else
        SRRC(idx) = (1-a+(4*a/pi));
    end
end

DC_gain = sum(SRRC);
SRRC = SRRC/DC_gain;

figure;
stem(SRRC);
grid;

figure;
hold on;
plot(RC);
plot(SRRC);
legend('RC', 'SRRC');
hold off;
grid;

RC_test = conv(SRRC, SRRC);
L = length(RC_test);
mid = ceil(L/2);
figure;
hold on;
plot(RC);
plot(1:(2*M*S+1), RC_test((mid-M*S):(mid+M*S)));
legend('RC', 'SRRC*SRRC');
hold off;
grid;

%% p3
figure;
for i = 1:4
    RC_down = RC(i:M:end);
    RC_down_zero = [RC_down, zeros(1,256-length(RC_down))];
    RC_downf = fft(RC_down_zero);
    subplot(4,2,i*2-1);
    stem(RC_down);
    grid on;
    subplot(4,2,i*2);
    plot(abs(RC_downf));
    ylim([0, 0.35]);
    yticks(0:0.05:0.35);
    xlim([0, 256]);
    grid on;
end

%% p4
clc; clear; close all;
data = load("Lab7_given_data\S7P4.mat");
b = data.b(:).';
Nb = length(b);

a = 0.3;
M = 32;
S = 5;
RC = rcfun(a, M, S);

bu = zeros(1,Nb*M);
bu(1:M:end) = b;

b_ps = conv(bu, RC, 'same');
figure;
plot(b_ps*M);
xlim([0,2048]);
grid;

b_re = b_ps(1:M:end);

%for i = 1:length(b_re)
%    if b_re(i) > 0
%        b_re(i) = 1;
%    else
%        b_re(i) = -1;
%    end
%end

figure;
hold on;
stem(b);
stem(M * b_re);
legend('b', 'b-recover');
hold off;
grid;

errors = zeros(1,M);
for phase = 0:M-1
    b_re_test = b_ps(phase+1:M:end);
    for i = 1:length(b_re_test)
        if b_re_test(i) > 0
            b_re_test(i) = 1;
        else
            b_re_test(i) = -1;
        end
    end
    errors(phase+1) = norm(b - b_re_test);
    fprintf('phase = %d, error = %.2e\n', phase, errors(phase+1));
end

% p5
rect_pulse = (1/M)*ones(1,M);
bu_practical = conv(bu, rect_pulse, 'same');
b_ps_practical = conv(bu_practical, RC, 'same');
figure;
hold on;
plot(b_ps);
plot(b_ps_practical);
xlim([0,2048]);
legend('ideal', 'practical');
hold off;
grid;

RC_rect = conv(RC, rect_pulse, 'same');
figure;
hold on;
plot(RC);
plot(RC_rect);
legend('RC', 'RC-rect');
hold off;
grid;

b_re_practical = b_ps_practical(1:M:end);

figure;
hold on;
stem(b_re * M);
stem(b_re_practical * M);
legend('ideal', 'practical');
hold off;
grid;

%% p6
%clc; clear; close all;
data = load("Lab7_given_data\S7P6.mat");
b = data.b(:).';
Nb = length(b);
M = 32;

delta = [ones(1,1), zeros(1,Nb*M-1)];
h = filter(IIR, delta);
figure;
plot(h);
xlim([0,256]);

bu = zeros(1,Nb*M);
bu(1:M:end) = b;
bu_IIR = filter(IIR, bu);
b_ps_IIR = filter(IIR, bu_IIR);
b_re = b_ps_IIR(43:M:end);

figure;
hold on;
stem(b);
stem(b_re * M);
legend('b', 'b-recover');
hold off;
grid;