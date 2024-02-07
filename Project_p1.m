% ECE 6961 - Project Part I
% Emma Dingman, Collin Griswold, Marlo Esperson

clear all, clc, close all

% Task 1: Understand CP-OFDM and ZP-OFDM
% Generate a sequence of N*M QPSK information symbols d~ in the frequency
% domain, where M = 16  is the number of subcarriers and N = 4 is the
% number of OFDM symbols.

M = 16;           % # of subcarriers
N = 4;            % # of OFDM symbols
numSymbols = N*M; % # of QPSK info symbols

% create bits and map them to QPSK symbols:
num=32;
b=sign(randn(num,1));
b = (b+1)/2;
b = reshape(b,num/2,2);

for n = 1:length(b)
    if b(n,1) == 0
        imaginary(n,1) = 1;
        if b(n,2) == 0
            real(n,1) = 1;
        else
            real(n,1) = -1;
        end
    else
        imaginary(n,1) = -1;
        if b(n,2) == 0
            real(n,1) = -1;
        else
            real(n,1) = 1;
        end
    end

    s(n,1) = real(n,1) + 1j*imaginary(n,1);
end
s = (1/sqrt(2)).*s; % symbols

% taps:
h0 = 0.227;
h1 = 0.46;
h2 = 0.688;
h3 = 0.46;
h4 = 0.227;

% create channel:
h = [h0 0 0 0 0 0 0 0 0 0 0 0 h4 h3 h2 h1]; % first line (y5)
for n = 2:16
    h(n,:) = circshift(h(n-1,:),1); % channel matrix
end

% generate F_N matrix:
for i = 1:M
    for J = 1:M
        f_i_j(i,J) = exp((-1j*2*pi*(i-1)*(J-1))/M);
    end
end
F = (1/sqrt(M)) .* f_i_j; % F_N matrix

x = F'.*s; % transmitted signal

y = h*x; % received signal (time domain)

Y = F*y; % received signal (freq domain)


% TESTING:
% F_I = F'*F; % identity matrix
% 
% h_f= F*h*F'; % channel in freq domain

% figure; [ch,w]=freqz([h0 h1 h2 h3 h4],1,'whole',2001);
% plot(w/pi/2, 20*log10(abs(ch)));
