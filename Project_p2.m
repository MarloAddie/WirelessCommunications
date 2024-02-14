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

x_ovf = zeros((N*(L-1+M)+L-1)*lambda,1);
i = 1;
for n = 0:(N*(L-1+M)+L-1)*lambda-1
    if mod(n,lambda) == 0
        x_ovf(n+1) = x(i);
        i = i+1;
    end
end

R = rcosdesign(beta,span,lambda,'sqrt');
%x_bb = conv(x_ovf,R);
x_bb = filter(R,1,x_ovf);


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
    x_pb(n+1) = real(x_bb(n+1)*exp(1j*2*pi*fc*n*ts));
end
x_pb = transpose(x_pb);

t1 = 0.05;
f0 = fc-4000;
f1 = fc+4000;
beta = (f1-f0)/t1;
i = 0;

for t = 0:1/Fs:1
    i = i+1;
    if t <= 0.05
        f = f0+beta*t;
        chirp(i) = cos(2*pi*f*t);
    elseif t >=0.25 && t <= 0.3
        f = f0+beta*t;
        chirp(i) = cos(2*pi*f*t);
    elseif t >=0.5 && t <= 0.55
        f = f0+beta*t;
        chirp(i) = cos(2*pi*f*t);
    elseif t >= 0.75 && t <= 0.8
        f = f0+beta*t;
        chirp(i) = cos(2*pi*f*t);
    else
        chirp(i) = 0;
    end
end

t = 0:1/Fs:1;
figure
plot(t, chirp)

LFM = [chirp.'; x_pb; chirp.'];
Nfft2 = 4096;

figure
colormap jet
spectrogram(LFM, Nfft2, Nfft2*3/4, Nfft2, Fs, 'yaxis')
