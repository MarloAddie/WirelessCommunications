%% ECE 6961 - Project Part II
% Emma Dingman, Collin Griswold, Marlo Esperson

clear all, clc, close all

% Parameters:
L = 201; % # of taps
M = 2048; % # of subcarriers
N = 1; % # of OFDM symbols
numBits = M*2; % # of bits to generate M QPSK symbols
Fs=192e3;
lambda = 24; % oversampling rate
beta = 0.125;
delay = 100;
span = 2*delay;

x = [];
d = [];

for n = 1:N
    
    % Generate random binary data
    binaryData = randi([0 1], numBits, 1);
    
    real_comp = binaryData(1:2:end);
    imag_comp = binaryData(2:2:end);
    s = [real_comp imag_comp];
    
    % QPSK mapping
    d_tilde = zeros(length(real_comp),1);
    for i = 1:length(real_comp)
        if (real_comp(i) == imag_comp(i))
            re = 1;
        else
            re = -1;
        end
        if (real_comp(i) < 1)
            im = 1i;
        else
            im = -1i;
        end
        d_tilde(i) = 1/sqrt(2)*(re+im); % QPSK symbols
    end
    
    d_i = ifft(d_tilde); % time-domain symbols
    
    x_i = [zeros(L-1,1);d_i(:)]; % zero padding
    x = cat(1, x,x_i); % x for all OFDM symbols

    d = cat(1, d,d_i); % d for all OFDM symbols
end
x = [x;zeros(L-1,1)]; % add zeros to the back

% x_ovf = zeros((N*(L-1+M)+L-1)*lambda,1);
% i = 1;
% for n = 0:(N*(L-1+M)+L-1)*lambda-1
%     if mod(n,lambda) == 0
%         x_ovf(n+1) = x(i);
%         i = i+1;
%     end
% end

x_ovf = upsample(x,24);

R = rcosdesign(beta,span,lambda,'sqrt');
% x_bb = conv(x_ovf,R);
x_bb = filter(R,1,x_ovf); % convolution


% figure
% plot(abs(fftshift(fft(R))))
fft_x_bb=fftshift(fft(x_bb));
frequency = [-length(fft_x_bb)/2:length(fft_x_bb)/2-1]*Fs/length(fft_x_bb);
figure; plot(frequency, abs(fftshift(fft(x_bb))))


% Step 6:
fc = 24e3; % carrier frequency 
B = 8e3; % bandwidth
Ts = 1/B; % sampling period
ts = Ts/lambda;

for n = 0:length(x_bb)-1
    x_pb(n+1,1) = real(x_bb(n+1)*exp(1i*2*pi*fc*n*ts));
end

t1 = 0.05;
f0 = fc-4000;
f1 = fc+4000;
t= 0:1/Fs:0.05;

chirp_sig = chirp(t,f0,t1,f1);
chirp_sig = [chirp_sig.';zeros((0.25*Fs)-length(chirp_sig),1)];
chirp_sig = [chirp_sig;chirp_sig;chirp_sig;chirp_sig];

t = 0:1/Fs:1-1/Fs;
figure
plot(t, chirp_sig)

LFM = [chirp_sig; x_pb; chirp_sig];

figure
colormap jet
Nfft2 = 2048;
spectrogram(LFM, Nfft2, Nfft2*3/4, Nfft2, Fs, 'yaxis')


% Step 8:
% passband to baseband:
for n = 0:length(x_pb)-1
    X_BB(n+1,1) = (2*x_pb(n+1)*cos(2*pi*fc*n*ts)) + (-2*1i*x_pb(n+1)*sin(2*pi*fc*n*ts));
end

X_OVF = filter(R,1,X_BB); % match filtering 

X = X_OVF((lambda)*(span)+1:lambda:end);

X = X(L:end);
%X_OVF = X_OVF((lambda)*(span):end);
%X = downsample(X_OVF,lambda-1);

for k = 0:M-1
    Y_k = 0;
    for n = 0:M-1
        Y_k = X(n+1)*exp((-1i*2*pi*n*k)/M) + Y_k;
    end
    Y(k+1,1) = Y_k;
end

figure
plot(abs(Y))
hold on
Y2 = fft(x);

% figure
plot(abs(fft(x)))
