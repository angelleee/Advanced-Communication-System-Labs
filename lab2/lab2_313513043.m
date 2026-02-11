%% homework 1
clc; close all; clear;
f = 1e6;
N = 256;
n = 0:N-1;
fs = [4e6, 1.5e6];

for i = 1:2
    x = cos(2*pi*(f/fs(i))*n);

    x_fft = fft(x);

    figure;
    stem(0:N-1, abs(x_fft));
    title(['DFT Magnitude Spectrum, fs = ' num2str(fs(i)/1e6) ' MHz']);
    xlabel('n');
    ylabel('|X[n]|');
end

%% homework 2
clc; close all; clear;
N = 128;
peak = 64;
n = 0:N-1;
trial = 100;

tri = [0:peak, peak-1:-1:1]';
tri_power = mean(tri.^2);
SNR_dB = 15;
SNR = 10^(SNR_dB/10);
noise_power = tri_power/SNR;

square5 = zeros(N,1);
square5(N/2-2:N/2+2) = 1;
square7 = zeros(N,1);
square7(N/2-3:N/2+3) = 1;
square9 = zeros(N,1);
square9(N/2-4:N/2+4) = 1;
square11 = zeros(N,1);
square11(N/2-5:N/2+5) = 1;

windows = {hanning(N), hamming(N), blackman(N), square5, square7, square9, square11};
window_names = {'Hanning','Hamming','Blackman','Square-5','Square-7','Square-9','Square-11'};

mse_avg = zeros(length(windows),1);
for i = 1:length(windows)
    window = windows{i};

    %figure;
    %hold on;

    mse = zeros(trial,1);
    for j = 1:trial
        noise = sqrt(noise_power) * randn(N,1);
        x = tri + noise;
        x_fft = fft(x);
        x_fft = fftshift(x_fft);
        x_fft = x_fft.*window;

        %stem(-N/2:(N-1)/2, abs(x_fft));
        %plot(-N/2:(N-1)/2, window);
        
        x_fft = ifftshift(x_fft);
        x_ifft = ifft(x_fft, 'symmetric');

        mse(j) = mean((x_ifft-tri).^2);
    end
    %hold off;
    
    mse_avg(i) = mean(mse);
    disp([window_names{i}, ' average MSE = ', num2str(mse_avg(i))]);
end

[~, best_idx] = min(mse_avg);
disp(['Optimum windowing function with min. MSE_avg: ', window_names{best_idx}]);
