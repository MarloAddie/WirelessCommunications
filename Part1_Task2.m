%% ECE 6961 - Project Part I Task II
% Emma Dingman, Collin Griswold, Marlo Esperson

clear all, clc, close all

% Parameters:
L = 201; % # of zeros +1
M = 2048; % # of subcarriers
N = 4; % # of OFDM symbols

[d_tilde,Y,LFM] = ZP_OFDM_Tx(L,M,N);



%% ZP-OFDM Transmitter Design & Validation:

function [d_tilde,Y,LFM] = ZP_OFDM_Tx(L,M,N)
    numBits = M*2; % # of bits to generate M QPSK symbols
    x = []; % initialize for concatenation
    
    for n = 1:N
        
        % Generate random binary data
        binaryData = randi([0 1], numBits, 1);
        
        real_comp = binaryData(1:2:end);
        imag_comp = binaryData(2:2:end);
        s = [real_comp imag_comp];
        
        % QPSK mapping
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
            d_tilde(i,n) = 1/sqrt(2)*(re+im); % QPSK symbols
        end
        
        d(:,n) = ifft(d_tilde(:,n)); % time domain symbols
        
        x_i = [zeros(L-1,1);d(:,n)]; % zero padding
        x = cat(1, x, x_i); % x for all OFDM symbols
    
    end
    x = [x;zeros(L-1,1)]; % add zeros to the back
    
    lambda = 24; % oversampling rate
    x_ovf = upsample(x,lambda); % upsample
    
    % Step 5: Filter with Square Root Raised Cosine
    Fs = 192e3; % sampling frequency
    beta = 0.125; % roll-off factor
    delay = 100; % samples before oversampling
    span = 2*delay; % number of zero crossings
    
    R = rcosdesign(beta,span,lambda,'sqrt'); % square-root raised cosine
    x_bb = filter(R,1,x_ovf); % convolution
    
    % Plot:
    fft_x_bb=fftshift(fft(x_bb));
    frequency = [-length(fft_x_bb)/2:length(fft_x_bb)/2-1]*Fs/length(fft_x_bb);
    
    fft_x_ovf=fftshift(fft(x_ovf));
    figure
    plot(frequency, abs(fft_x_ovf))
    title('x_o_v_f')

    figure
    plot(frequency, abs(fft_x_bb))
    title('x_b_b')
    
    % Step 6: 
    fc = 24e3; % carrier frequency 
    B = 8e3; % bandwidth
    Ts = 1/B; % sampling period
    ts = Ts/lambda;
    
    % baseband to passband:
    for n = 0:length(x_bb)-1
        x_pb(n+1,1) = real(x_bb(n+1)*exp(1i*2*pi*fc*n*ts));
    end
    
    % Step 7: Concatenate Linear Frequency Modulation
    t1 = 0.05; % length of each chirp
    f0 = fc-4000; % instantaneous frequency at t=0
    f1 = fc+4000; % instantaneous frequency at t=t1
    t = 0:1/Fs:t1; % time array for each chirp
    
    chirp_sig = chirp(t,f0,t1,f1); % generate 1 chirp
    chirp_sig = [chirp_sig.';zeros((0.25*Fs)-length(chirp_sig),1)]; % add zeros to chirp
    chirp_sig = [chirp_sig;chirp_sig;chirp_sig;chirp_sig]; % create 4 chirps
    
    LFM = [chirp_sig; x_pb; chirp_sig]; % concatenate LFM before and after OFDM symbols
    
    % Plot:
    figure
    colormap jet
    Nfft2 = M;
    spectrogram(LFM, Nfft2, Nfft2*3/4, Nfft2, Fs, 'yaxis')
    title('LFM')
    
    
    % Step 8: Validation
    
    % 1) Passband to Baseband:
    for n = 0:length(x_pb)-1
        X_BB(n+1,1) = 2*(x_pb(n+1)*cos(2*pi*fc*n*ts)) + 2*(-1i*x_pb(n+1)*sin(2*pi*fc*n*ts));
    end
    
    % 2) 
    X_OVF = filter(R,1,X_BB); % match filtering
    X_d = X_OVF(lambda*span+1:lambda:end); % downsampling
    
    % 3)
    X_r = reshape(X_d,[length(X_d)/N,N]); % separate OFDM symbols into columns
    
    X = X_r(L:end,:); % remove zeros for each symbol
    Y = fft(X); % frequency domain symbols
    
    % Plot: 
    figure
    plot(Y)
    title('Recovered Symbols (y)')
    
    figure
    plot(d_tilde)
    title('Original Symbols (s)')
end
