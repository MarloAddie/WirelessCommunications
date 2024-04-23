% Collin Griswold, Emma Dingman, Marlo Esperson
% One Hydrophone

clc, clear, close all

addpath('C:\Users\marlo\Documents\ECE6961\WirelessCommunications\Project3files\Matlab_Python_LDPC_update\Matlab_Python_LDPC_update\Matlab_Python_LDPC\Matlab_version\mex');


load('OFDM_PILOT.mat')
load('ofdm_map.mat')

load('benchmark_NoiseVar_172648_1.mat')
load('benchmark_intermediate_172648_from_1_single_hydrophones.mat');
load('benchmark_Zw_172648_1.mat')
load('benchmark_Zw_172648_2.mat')
load('benchmark_Zw_172648_3.mat')
load('INTRLVR.mat')
load('CODE.mat')

k = 2048;   % # of subcarriers
L = 199;    % Zero Padding
kp = 512;   % # of pilot subcarriers
null_sub = 112; % # of null subcarriers

load('Z_w1_test.mat');
z_w1 = Z_w; % Frequency data calculated in Part 2 for one hydrophone
[Hw1,noiseVar1] = channelest(z_w1, ofdm_map, OFDM_PILOT); % Task 1
Hw1 = Hw1(ofdm_map==2,:); % Locating the channel for data subcarriers


z_w1 = z_w1(ofdm_map==2,:); % Locating frequency data for data subcarriers

% QPSK Map
X1 = 1/sqrt(2) + (1/sqrt(2))*1i;
X2 = -1/sqrt(2) + (1/sqrt(2))*1i;
X3 = 1/sqrt(2) - (1/sqrt(2))*1i;
X4 = -1/sqrt(2) - (1/sqrt(2))*1i;


for i = 1:21
    LRcount1 = 1;
    LRcount2 = 2;
    for j = 1:1420 % # of data subcarriers
        Hk = [Hw1(j,i)]; % Est. channel for the kth subcarrier of one hydrophone
        Zk = [z_w1(j,i)]; % Received freq domain for the kth subcarrier of one hydrophone
        Hknorm = sqrt(Hk'*Hk);
        qk(j,i) = (Hk'*Zk)/(Hknorm); % Projected signal of kth subcarrier of one hydrophone

        noiseVar(i) = (1/(Hknorm)^2) * ( noiseVar1(i) * Hw1(j,i)' * Hw1(j,i) );

        % Log Likelihood Ratio:

        A = -abs(qk(j,i)-Hknorm*X1)^2/noiseVar(i);
        B = -abs(qk(j,i)-Hknorm*X3)^2/noiseVar(i);
        
        Lb1_num = max(real(A),real(B))+log(1+exp(-abs(real(B)-real(A))));
        
        A = -abs(qk(j,i)-Hknorm*X2)^2/noiseVar(i);
        B = -abs(qk(j,i)-Hknorm*X4)^2/noiseVar(i);
        
        Lb1_den = max(real(A),real(B))+log(1+exp(-abs(real(B)-real(A))));
        
        A = -abs(qk(j,i)-Hknorm*X1)^2/noiseVar(i);
        B = -abs(qk(j,i)-Hknorm*X2)^2/noiseVar(i);
        
        Lb2_num = max(real(A),real(B))+log(1+exp(-abs(real(B)-real(A))));
        
        A = -abs(qk(j,i)-Hknorm*X3)^2/noiseVar(i);
        B = -abs(qk(j,i)-Hknorm*X4)^2/noiseVar(i);
        
        Lb2_den = max(real(A),real(B))+log(1+exp(-abs(real(B)-real(A))));
        
        Lb1 = Lb1_num-Lb1_den;
        Lb2 = Lb2_num-Lb2_den;
        
        % Combining the log likelihood ratios
        LR(LRcount1, i) = Lb1;
        LR(LRcount2, i) = Lb2;
        LRcount1  = LRcount1 + 2;
        LRcount2 = LRcount2 + 2;

    end
end



NAME = './5G_LDPC_M10_N20_Z142_Q2_nonVer.txt';
[address, LDPC_INFOLEN] = ldpc_mex_initial_CAPI([1420,2840,2],NAME);

for i = 2:21

    LR_in_de = zeros(length(LR(:,i)),1);
    LR_in_de(INTRLVR) = LR(:,i);
    APP_code(:,i) = ldpcDecoder_CAPI(address,LR_in_de);
    
    % Hard Decision
    est_code(:,i) = (APP_code(:,i)<0);
    bec(:,i) = sum(abs(est_code(:,i)-CODE(:,i)));


end

BER = sum(bec) / 2840 / 20; % Bit Error Rate
BLER = nnz(bec) / 20;       % Block Error Rate




%% This function defined in task 1 estimates the channel and noise variance
function [Hw, noiseVar] = channelest(z_w, ofdm_map, OFDM_PILOT) 
    
    k = 2048;   % # of subcarriers
    L = 199;    % Zero Padding
    kp = 512;   % # of pilot subcarriers
    null_sub = 112; % # of null subcarriers
    
    idx = find(ofdm_map ==1); % Pilot indices
    
    Z_pw = z_w(ofdm_map==1,:); % Frequency data for pilot subcarriers
    
    for nn = 1:kp
        for n = 1:L+1
            V(nn,n) = exp((-1i*2*pi*(n-1)*(idx(nn)-1))/k); % DFT Matrix
        end
    end
    
    D = OFDM_PILOT(ofdm_map==1);
    D = diag(D); % Transmitted QPSK symbols
    hls = (1/kp)*V'*D'*Z_pw; % Time domain channel est  200 x 1
    
    for nn = 1:k
        for n = 1:L+1
            V2(nn,n) = exp((-1i*2*pi*(n-1)*(nn-1))/k);  % Fourier transform 
        end                                             % Dimension: (k x (L+1))
                                                        % 2048 x (200)
    end
    
    Hw = V2*hls; % Channel Gain. Freq domain. Dimension: k x (L+1)   (2048 x 200)

    % Noise variance estimate:
    for kl =1:21

        znw = z_w(:,kl);
        Z_pw_null = znw((ofdm_map == 0));
        noiseVar(kl) = (1/null_sub)*sum(abs(Z_pw_null).^2);

    end



end