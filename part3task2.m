clc, clear, close all

addpath('C:\Users\marlo\Documents\ECE6961\WirelessCommunications\Project3files\Matlab_Python_LDPC_update\Matlab_Python_LDPC_update\Matlab_Python_LDPC\Matlab_version\mex');


load('OFDM_PILOT.mat')
load('ofdm_map.mat')
%%
load('benchmark_NoiseVar_172648_1.mat')
load('benchmark_intermediate_172648_from_1_single_hydrophones.mat');
load('benchmark_Zw_172648_1.mat')
load('benchmark_Zw_172648_2.mat')
load('benchmark_Zw_172648_3.mat')
load('INTRLVR.mat')
load('CODE.mat')
%%
k = 2048;
L = 199;
kp = 512; % # of Subcarriers
null_sub = 112;

z_w1 = bb_rece_data_172648_1474;
Hw1 = channelest(z_w1, ofdm_map, OFDM_PILOT);
z_w2 = bb_rece_data_172648_1475;
Hw2 = channelest(z_w2, ofdm_map, OFDM_PILOT);
z_w3 = bb_rece_data_172648_1476;
Hw3 = channelest(z_w3, ofdm_map, OFDM_PILOT);

for i = 1:21
    for j = 1:k
        Hk = [Hw1(j,i);Hw2(j,i);Hw3(j,i)];
        Zk = [z_w1(j,i);z_w2(j,i);z_w3(j,i)];
        Hknorm = sqrt(Hk'*Hk);
        %Hknorm(j,i) = sqrt(Hw1(j,i)'*Hw1(j,i) + Hw2(j,i)'*Hw2(j,i) + Hw3(j,i)'*Hw3(j,i));
        qk(j,i) = (Hk'*Zk)/(Hknorm);
    end
end






function [Hw] = channelest(z_w, ofdm_map, OFDM_PILOT) 
    
    k = 2048;
    L = 199;
    kp = 512; % # of Subcarriers
    null_sub = 112;
    
    idx = find(ofdm_map ==1);
    
   % z_w = bb_rece_data_172648_1474; % Symbols after doppler compensation
    Z_pw = z_w(ofdm_map==1,:); 
    
    for nn = 1:kp
        for n = 1:L+1
            V(nn,n) = exp((-1i*2*pi*(n-1)*(idx(nn)-1))/k);
        end
    end
    
    D = OFDM_PILOT(ofdm_map==1);
    
    D = diag(D); % Transmitted QPSK symbols
    hls = (1/kp)*V'*D'*Z_pw; % Time domain channel est  200 x 1
    
    for nn = 1:k
        for n = 1:L+1
            V2(nn,n) = exp((-1i*2*pi*(n-1)*(nn-1))/k); % Fourier transform 
        end                                         % Dimension: (k x (L+1))
                                                    % 2048 x (200)
    end
    
    Hw = V2*hls; % Channel Gain. Freq domain. Dimension: k x (L+1)   (2048 x 200)
    
    %Z_pw_null = z_w((ofdm_map == 0),:);

end