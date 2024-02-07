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

% Hello