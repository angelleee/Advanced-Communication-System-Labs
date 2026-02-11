clc; clear; close all;

data = load('./Lab5_given_data/S5P3.mat');
bs = data.bs(:).';
L = length(bs);

%16QAM
n = 2;
N = 2^n;
levels = -(N-1):2:(N-1);

grayCodes = zeros(N,1);
for i = 0:N-1 %binary
    grayCodes(i+1) = bitxor(i, bitshift(i,-1));
end

k = 2*n;
Nsymbols = floor(L/k);
symbols_I = zeros(Nsymbols,1);
symbols_Q = zeros(Nsymbols,1);
for i = 1:Nsymbols
    seg_I = bs((i-1)*k+1 : (i-1)*k+n);
    segVal_I = seg_I * (2.^(n-1:-1:0)).';
    idx_I = find(segVal_I == grayCodes);
    symbols_I(i) = levels(idx_I);

    seg_Q = bs((i-1)*k+n+1 : i*k);
    segVal_Q = seg_Q * (2.^(n-1:-1:0)).';
    idx_Q = find(segVal_Q == grayCodes);
    symbols_Q(i) = levels(idx_Q);
end
symbols = symbols_I + 1j*symbols_Q;

Es = mean(abs(symbols).^2);
SNRdb_vec = 0:2:18;
SER = zeros(1,length(SNRdb_vec));
theo_SER = zeros(1,length(SNRdb_vec));

for j = 1:length(SNRdb_vec)
    SNRdb = SNRdb_vec(j);
    SNR = 10^(SNRdb/10);
    N0 = Es/SNR;
    noise = sqrt(N0/2)*(randn(size(symbols)) + 1j*randn(size(symbols)));
    rx = symbols + noise;
    
    noise_energy = mean(abs(noise).^2);
    test_SNR = Es/noise_energy;
    test_SNR = 10*log10(test_SNR);
    fprintf('test SNR(dB) = %.4f\n', test_SNR);
    
    figure;
    hold on;
    scatter(real(rx), imag(rx));
    scatter(real(symbols), imag(symbols), "filled");
    xlabel('I');
    ylabel('Q');
    legend('Noisy symbols','Original symbols');
    hold off;
    grid on;
    
    rx_I = real(rx);
    rx_Q = imag(rx);
    det_I = zeros(size(rx_I));
    det_Q = zeros(size(rx_Q));
    for i = 1:length(rx)
        [~, idx_I] = min(abs(rx_I(i) - levels));
        [~, idx_Q] = min(abs(rx_Q(i) - levels));
        det_I(i) = levels(idx_I);
        det_Q(i) = levels(idx_Q);
    end
    rx_detect = det_I + 1j*det_Q;
    
    num_errors = sum(rx_detect ~= symbols);
    SER(j) = num_errors / length(symbols);
    %fprintf('SER = %.4f\n', SER(i));
    
    E0 = 1;
    p = (6/4)*qfunc(sqrt(E0/(N0/2)));
    theo_SER(j) = 1-(1-p)^2;
    %fprintf('theoretical SER = %.4f\n', theo_SER(i);
end

figure;
hold on;
plot(SNRdb_vec, SER);
plot(SNRdb_vec, theo_SER);
legend('SERs', 'theoretical SERs');
xlabel('SNR(dB)');
ylabel('SER');
hold off;
grid on;