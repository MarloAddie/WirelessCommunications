clear all
close all
clc

load('ofdm_map.mat')
%load('benchmark_parameter_174623_1472.mat');
load('benchmark_rece_data_174623_1472.mat');

LFM = rece_data_ofdm_bench;

figure
plot(LFM)

fc = 24e3; % carrier frequency
samplingRate = 256e3; %Hz

%YPB = bandpass(LFM, [-1000+fc, 8000+fc], samplingRate);
YPB = bandpass(LFM, [fc-4000, fc+4000], samplingRate);
%YPB2 = bandpass(LFM, [fc, 8000+fc], samplingRate);
Trx = 8.2687;  % Approximated from spectrogram
Ttx = 8.2695; % sec
%a_hat = Ttx/Trx - 1;
a_hat = 5e-5;

% Resample
Y_PB_re = resample(YPB,round((1+a_hat)*1E5), 1E5);

Ls = 192; 
Ms = 256;
Lp = 24;
N = Lp*Ls-1;
h = Ls*fir1(N,1/Ms,kaiser(N+1,7.8562));
Y_PB_re_tilde = upfirdn(Y_PB_re, h, Ls, Ms);

figure
plot(Y_PB_re)
hold on
plot(Y_PB_re_tilde)
legend('YPBre', 'YPBretilde')

load('pilot_signal_for_synchronization.mat')
%load('test_rece_data_173048_1472.mat')

c = (OFDM_data_pre_old);
%c = (rece_data_ofdm_test);
figure
plot(real(c))

b = xcorr(Y_PB_re_tilde,c);
figure
plot((b))

[peak,peakindex] = max(b);

peakindex = (peakindex)-length(Y_PB_re_tilde);

Y_PB_re_tilde = Y_PB_re_tilde(peakindex:end);
figure
plot(Y_PB_re_tilde)

load('itc_1007_compfilter.mat')
hcomp = h_comp;

yhtilde = conv(Y_PB_re_tilde, hcomp);
yhtilde = yhtilde(50+1:end);
figure
plot(yhtilde)

Fs = 192e3;
ts = 1/Fs;

% Passband to Baseband
for n = 0:length(yhtilde)-1
    Y_BB(n+1,1) = 2*(yhtilde(n+1)*cos(2*pi*fc*n*ts)) + 2*(-1i*yhtilde(n+1)*sin(2*pi*fc*n*ts));
end

beta = 0.125; % roll-off factor
delay = 100; % samples before oversampling
span = 2*delay; % number of zero crossings
lambda = 24; % oversampling rate

R = rcosdesign(beta,span,lambda,'sqrt'); % square-root raised cosine

ybb = filter(R,1,Y_BB); % convolution
ybb = ybb(lambda*delay+1:end);

figure
plot(real(ybb))

L = 200; % # of zeros (padding)
k = 2048; % # of QPSK symbols

% Generate DFT Matrix:
F = zeros(k,k);
for m = 0:k-1
    for i = 0:k+L-1
        F(m+1,i+1) = exp((-1j*2*pi*m*i)/k); % DFT Matrix
    end
end

n0 = zeros(1,21);
eps = zeros(1,21);
n01_idx = 0;
for n01 = [2200:1:2400]
    n01_idx = n01_idx + 1;
    eps1_idx = 0;
    for eps1 = [-2:0.1:2]
        eps1_idx = eps1_idx +1;
        n_idx = transpose(n01+[0:lambda:(k+L)*lambda - 1]);
        YBB_1 = ybb(n_idx).*exp(-1i*2*pi*eps1*n_idx*ts); % CFO Compensation + Downsampled
        
        Z_m = F*YBB_1; % FFT
        P_null(n01_idx, eps1_idx) = sum((abs(Z_m(ofdm_map==0))).^2);
    end
end

[M, idx] = min(P_null(:));
[row, col] = ind2sub(size(P_null),idx);

n01 = [2200:1:2400];
eps1 = [-2:0.1:2];

min_eps = eps1(col);
min_n01 = n01(row);

n0(1) = min_n01 -(k+L)*lambda;

for W = 1:21
n01_idx = 0;
    for n01 = n0(W) +(k+L)*lambda+[-2*lambda:1:2*lambda] 
        n01_idx = n01_idx + 1;
        eps1_idx = 0;
        for eps1 = [-2:0.1:2]
            eps1_idx = eps1_idx +1;
            n_idx = transpose(n01+[0:lambda:(k+L)*lambda - 1]);
            YBB_1 = ybb(n_idx).*exp(-1i*2*pi*eps1*n_idx*ts); % CFO Compensation + Downsampled

            Z_m = F*YBB_1; % FFT
            P_null2(n01_idx, eps1_idx) = sum((abs(Z_m(ofdm_map==0))).^2);
        end
    end
    [M, idx] = min(P_null2(:));
    [row, col] = ind2sub(size(P_null2),idx);

    n01_v = n0(W) +(k+L)*lambda+[-2*lambda:1:2*lambda];
    eps1_v = [-2:0.1:2];

    eps(W+1) = eps1_v(col);
    n0(W+1) = n01_v(row);
end

