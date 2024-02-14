%% ECE 6961 - Project Part I
% Emma Dingman, Collin Griswold, Marlo Esperson

%% CP-OFDM

clear all, clc, close all

L = 5; % # of taps
M = 16; % # of subcarriers
N = 4; % # of OFDM symbols

% initialize outputs
y_tilde = [];
d_tilde = [];
d_tilde_rx = [];

for n = 1:N
    [new_y_tilde,new_d_tilde,new_d_tilde_rx] = CP_OFDM(L,M);
    y_tilde = cat(1,new_y_tilde, y_tilde); % received signal for N OFDM symbols
    d_tilde = cat(1,new_d_tilde, d_tilde); % M*N transmitted QPSK symbols
    d_tilde_rx = cat(1,new_d_tilde_rx, d_tilde_rx); % M*N received QPSK symbols (should match d_tilde)
end

%% ZP-OFDM

clear all, clc, close all

L = 5; % # of taps
M = 16; % # of subcarriers
N = 4; % # of OFDM symbols

% initialize outputs
y_tilde = [];
d_tilde = [];
d_tilde_rx = [];

for n = 1:N
    [new_y_tilde,new_d_tilde,new_d_tilde_rx] = ZP_OFDM(L,M);
    y_tilde = cat(1,new_y_tilde, y_tilde); % received signal for N OFDM symbols
    d_tilde = cat(1,new_d_tilde, d_tilde); % M*N transmitted QPSK symbols
    d_tilde_rx = cat(1,new_d_tilde_rx, d_tilde_rx); % M*N received QPSK symbols (should match d_tilde)
end












%% Functions:

% CP-OFDM:
function [y_tilde,d_tilde,d_tilde_rx] = CP_OFDM(L,M)
    cp = L-1; % length of cp
    numBits = M*2; % # of bits to generate M QPSK symbols
    
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
    
    
    d = ifft(d_tilde); % time-domain symbols
    
    x = [d(end-cp+1:end);d(:)]; % add cp
    
    h = [0.227 0.46 0.688 0.46 0.227]; % taps
    
    y = zeros(M,1);
    for m = L:M+L-1 % remove cp
        y_m = 0;
       for l = 0:L-1
            y_m = h(l+1)*x(m-l) + y_m; % convolution
       end
       y(m-L+1) = y_m; % received signals (time-domain)
    end
    
    y_tilde = fft(y); % received signals (freq-domain)
    
    % verify correctness using one-tap equalizer:
    h_tilde = fft(h,M);
    h_tilde = transpose(h_tilde); % channel in freq-domain
    d_tilde_rx = y_tilde./h_tilde; % received symbols are equal to transmitted QPSK symbols

end

% ZP-OFDM:
function [y_tilde,d_tilde,d_tilde_rx] = ZP_OFDM(L,M)
    
    numBits = M*2; % # of bits to generate M QPSK symbols
    
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
    
    d = ifft(d_tilde); % time-domain symbols
    
    x = [zeros(L-1,1);d(:);zeros(L-1,1)]; % pad with zeros
    
    h = [0.227 0.46 0.688 0.46 0.227]; % taps
    
    y = zeros(M,1);
    for m = L:M+2*(L-1) 
        y_m = 0;
       for l = 0:L-1
            y_m = h(l+1)*x(m-l) + y_m; % convolution
       end
       y(m-L+1) = y_m; % received signals (time-domain)
    end
    
    y_tilde = zeros(M,1);
    for i = 0:M-1
        y_i = 0;
        for k = 0:M+L-2
            y_i = y(k+1)*exp((-1j*2*pi*k*i)/M) + y_i; % modified DFT
        end
        y_tilde(i+1) = y_i; % received signals (freq-domain)
    end
    
    
    % verify correctness using one-tap equalizer:
    h_tilde = fft(h,M);
    h_tilde = transpose(h_tilde); % channel in freq-domain
    d_tilde_rx = y_tilde./h_tilde; % received symbols are equal to transmitted QPSK symbols

end

