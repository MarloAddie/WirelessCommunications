clear all
close all
clc

%load('benchmark_parameter_174623_1472.mat');
load('benchmark_rece_data_174623_1472.mat');

LFM = rece_data_ofdm_bench;

figure
plot(LFM)

fc = 24e3; % carrier frequency
samplingRate = 256e3; %Hz

YPB = bandpass(LFM, [-1000+fc, 8000+fc], samplingRate);
%YPB2 = bandpass(LFM, [fc, 8000+fc], samplingRate);
Trx = 8.2687;  % Approximated from spectrogram
Ttx = 8.2695; % sec
a_hat = Ttx/Trx - 1;

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

peaks = zeros(length(b),1);
for i = 1:length(b)-10
    peaks(i) = var(b(i:i+10));
end

figure
plot(peaks)