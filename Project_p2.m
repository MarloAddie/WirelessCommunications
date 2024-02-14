%% ECE 6961 - Project Part II
% Emma Dingman, Collin Griswold, Marlo Esperson

clear all, clc, close all

% Parameters:
L = 5; % # of taps
M = 2048; % # of subcarriers
N = 4; % # of OFDM symbols
numBits = M*2; % # of bits to generate M QPSK symbols
lambda = 24; % oversampling rate
beta = 0.125;
delay = 100;
span = 2*delay;

x = [];
d = [];

for n = 1:N
    
    % Generate random binary data
    binaryData = randi([0 1], numBits, 1);
    
    real = binaryData(1:2:end);
    imag = binaryData(2:2:end);
    s = [real imag];
    
    % QPSK mapping
    d_tilde = zeros(length(real),1);
    for i = 1:length(real)
        if (real(i) == imag(i))
            re = 1;
        else
            re = -1;
        end
        if (real(i) < 1)
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
x = [x;zeros(L-1,1)];

x_ovf = zeros((N*(L-1+M)+L-1)*lambda,1);
i = 1;
for n = 0:(N*(L-1+M)+L-1)*lambda-1
    if mod(n,lambda) == 0
        x_ovf(n+1) = x(i);
        i = i+1;
    end
end

R = rcosdesign(beta,span,lambda,'sqrt');
x_bb = conv(x_ovf,R);


figure
plot(abs(fftshift(fft(R))))

figure
plot(abs(fftshift(fft(x_bb))))





