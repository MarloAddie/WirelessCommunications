%% ECE 6961 - Project Part II
% Emma Dingman, Collin Griswold, Marlo Esperson

%% 
clear all, clc, close all

load('ofdm_map.mat') % load ofdm map
load('benchmark_rece_data_172648_1.mat')

%load('test_rece_data_173048_1472.mat'); % load received data

%y = rece_data_ofdm_test; % received data
y = rece_data_ofdm_bench; 

fc = 24e3; % carrier frequency
samplingRate = 256e3; % receiver sampling rate

YPB = bandpass(y, [fc-4000, fc+4000], samplingRate); % passband filter

figure
plot(YPB)

% Find Trx from plot:

Trx = 8.2641; % estimated received signal duration
Ttx = 8.2695; % given transmitted signal duration

a_hat = Ttx/Trx - 1; % estimated doppler rate
% a_hat = 6.4856e-4; 

Y_PB_re = resample(YPB,round((1+a_hat)*1E5), 1E5); % resample with a_hat

% Resample from 256 kHz to 192 kHz:
Ls = 192; 
Ms = 256;
Lp = 24;
N = Lp*Ls-1;
h = Ls*fir1(N,1/Ms,kaiser(N+1,7.8562));
Y_PB_re_tilde = upfirdn(Y_PB_re, h, Ls, Ms);

load('pilot_signal_for_synchronization.mat') % load pilot signal
%load('test_rece_data_173048_1472.mat')

c = OFDM_data_pre_old; % pilot signal

% cross correlate the pilot with the signal see where the data begins
b = xcorr(Y_PB_re_tilde,c);

[~,peakindex] = max(b); % peak from xcorr

peakindex = (peakindex)-length(Y_PB_re_tilde); % data's starting index 

Y_PB_re_tilde = Y_PB_re_tilde(peakindex:end); % extract data

load('itc_1007_compfilter.mat') % load vector for match filtering

yhtilde = conv(Y_PB_re_tilde, h_comp); % match filtering (first stage)
yhtilde = yhtilde(50+1:end); % remove the delay caused by the filter

Fs = 192e3; % sampling rate
ts = 1/Fs; % sampling period

% Passband to Baseband:
for n = 1:length(yhtilde)
    Y_BB(n,1) = yhtilde(n)*2*cos(2*pi*fc*n*ts) + -1i*yhtilde(n)*2*sin(2*pi*fc*n*ts);
end

beta = 0.125; % roll-off factor
delay = 100; % samples before oversampling
span = 2*delay; % number of zero crossings
lambda = 24; % oversampling rate

R = rcosdesign(beta,span,lambda,'sqrt'); % square-root raised cosine

ybb = filter(R,1,Y_BB); % match filtering (stage 2)
ybb = ybb(lambda*delay+1:end); % remove the delay caused by the filter

L = 200; % # of zero-padded symbols
k = 2048; % # of subcarriers

% Generate DFT Matrix:
F = zeros(k,k);
for m = 0:k+L-1
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
        eps1_idx = eps1_idx + 1;
        n_idx = transpose(n01+[0:lambda:(k+L)*lambda - 1]); % indices
        YBB_1 = ybb(n_idx).*exp(-1i*2*pi*eps1*n_idx*ts); % CFO Compensation + Downsampled
        
        Z_m = F*YBB_1; % DFT
        P_null1(n01_idx, eps1_idx) = sum((abs(Z_m(ofdm_map==0))).^2); % power over null subcarriers
    end
end

% Find minimum eps1 and n01:
[M, idx] = min(P_null1(:));
[row, col] = ind2sub(size(P_null1),idx);
n01 = [2200:1:2400];
eps1 = [-2:0.1:2];

min_eps = eps1(col); % minimum eps1
min_n01 = n01(row); % minimum n01 = starting index

n0(1) = min_n01 -(k+L)*lambda; % initial starting point

for W = 1:21
n0_idx = 0;
    for n01 = n0(W) +(k+L)*lambda+[-2*lambda:1:2*lambda] 
        n0_idx = n0_idx + 1;
        eps_idx = 0;
        for eps1 = [-2:0.1:2]
            eps_idx = eps_idx +1;
            n_idx = transpose(n01+[0:lambda:(k+L)*lambda - 1]); % indices
            YBB_1 = ybb(n_idx).*exp(-1i*2*pi*eps1*n_idx*ts); % CFO Compensation + Downsampled

            Z_m = F*YBB_1; % DFT
            P_null(n0_idx, eps_idx) = sum((abs(Z_m(ofdm_map==0))).^2); % power over null subcarriers
        end
    end

    % Find minimum eps and n0 values:
    [M, idx] = min(P_null(:));
    [row, col] = ind2sub(size(P_null),idx);
    n01 = n0(W) +(k+L)*lambda+[-2*lambda:1:2*lambda];
    eps1 = [-2:0.1:2];
    
    eps(W+1) = eps1(col); % minimum eps
    n0(W+1) = n01(row); % minimum n0
end

% Erase unecessary minimum values from step 9:
n0 = n0(2:end);
eps = eps(2:end);

% Generate the frequency data using optimal eps and n0 values:
for n = 1:length(eps)
  n_idx = transpose(n0(n) + [0:lambda:(k+L)*lambda - 1]);  
  YBB = ybb(n_idx).*exp(-1i*2*pi*eps(n)*n_idx*ts); % CFO Compensation + Downsampled
  Z_w(:,n) = F*YBB; % FFT
end


%% Validation:

load('benchmark_parameter_174623_1472.mat')
load('benchmark_Zw_174623_1472.mat')

