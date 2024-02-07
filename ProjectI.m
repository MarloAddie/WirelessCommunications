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
N = 1;            % # of OFDM symbols
numSymbols = N*M; % # of QPSK info symbols

% Generate rand binary data
randSeq = rand(1, numSymbols);
binaryData = (randSeq > 0.5);
x = zeros(1, length(randSeq*4));

for i = 0:length(binaryData)/2

 if(binaryData(i+1)==binaryData(i+2))
     re = 1;
 else 
     re = -1;
 end
if(x(i+1) == 0)
    im = 1i;
else
    im = -1i;
end
x(i+1) = (1/sqrt(2))*[re+im];
end

% cp = 5;
% 
% x(1:4) = binaryData(13:16);
% x(5:20) = binaryData(1:16);
