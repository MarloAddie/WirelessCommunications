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

[peak,peakindex] = max(b);

peakindex = (peakindex)-length(Y_PB_re_tilde);

Y_PB_re_tilde = Y_PB_re_tilde(peakindex:end);
figure
plot(Y_PB_re_tilde)

load('itc_1007_compfilter.mat')
hcomp = h_comp;

yhtilde = conv(Y_PB_re_tilde, hcomp);
yhtilde = yhtilde(50:end);
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

for n01 = [2200:1:2400]
    for eps1 = [-2:0.1:2]
        for n = 0:(k+L)*lambda - 1
            YBB_1(n+1) = ybb(n01+n)*exp(-1i*2*pi*eps1*(n01+n)*ts);
        end
        YBBv1 = downsample(YBB_1,lambda);
    end
end