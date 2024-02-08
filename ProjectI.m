% ECE 6961 - Project Part I
% Emma Dingman, Collin Griswold, Marlo Esperson

clc
clear all
close all

% Task 1: Understand CP-OFDM and ZP-OFDM
% Generate a sequence of N*M QPSK information symbols d~ in the frequency
% domain, where M = 16  is the number of subcarriers and N = 4 is the
% number of OFDM symbols.

M = 16;           % # of subcarriers
N = 4;            % # of OFDM symbols
numSymbols = N*M; % # of QPSK info symbols

% Generate rand binary data
randSeq = rand(1, numSymbols);
binaryData = (randSeq > 0.5);
s = zeros(1, length(randSeq*4)); % Symbols
cp = 5;

s(1:4) = binaryData(13:16);
s(5:20) = binaryData(1:16);

realComp = s(1:2:end);
imagComp = s(2:2:end);

graycode = zeros(1,length(realComp));

for i = 1:length(realComp)
    if realComp(i) && imagComp(i) == 0
        graycode(i) = (1/sqrt(2))*(1+1i);
    end
    
    if realComp(i) && imagComp(i) == 1
        graycode(i) = (1/sqrt(2))*(1-1i);
    end
    
    if realComp(i) == 0 && imagComp(i) == 1
        graycode(i) = (1/sqrt(2))*(-1+1i);
    end
    
    if realComp(i) == 1 && imagComp(i) == 0
        graycode(i) = (1/sqrt(2))*(-1-1i);
    end
end


