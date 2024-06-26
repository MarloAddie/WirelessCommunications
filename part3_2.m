clc, clear, close all

addpath('C:\Users\marlo\Documents\ECE6961\WirelessCommunications\Project3files\Matlab_Python_LDPC_update\Matlab_Python_LDPC_update\Matlab_Python_LDPC\Matlab_version\mex');


load('OFDM_PILOT.mat')
load('ofdm_map.mat')

load('benchmark_NoiseVar_172648_1.mat')
load('benchmark_intermediate_172648_from_1_single_hydrophones.mat');
load('benchmark_Zw_172648_1.mat')
load('INTRLVR.mat')
load('CODE.mat')


k = 2048;
L = 199;
kp = 512; % # of Subcarriers
null_sub = 112;

idx = find(ofdm_map ==1);

z_w = bb_rece_data_172648_1474; % Symbols after doppler compensation
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

for kl =1:21
    
    znw = z_w(:,kl);
    Z_pw_null = znw((ofdm_map == 0));
    noiseVar(kl) = (1/null_sub)*sum(abs(Z_pw_null).^2);
    
end





sd = find(ofdm_map == 2);

zd = z_w(sd,:);

HD = Hw(sd,:);

X1 = 1/sqrt(2) + (1/sqrt(2))*1i;
X2 = -1/sqrt(2) + (1/sqrt(2))*1i;
X3 = 1/sqrt(2) - (1/sqrt(2))*1i;
X4 = -1/sqrt(2) - (1/sqrt(2))*1i;

LR = [];



for jkl = 1:length(noiseVar)

    LRcount1 = 1;
    LRcount2 = 2;

    for pl = 1:length(HD)
        
        A = -abs(zd(pl,jkl)-HD(pl,jkl)*X1)^2/noiseVar(jkl);
        B = -abs(zd(pl,jkl)-HD(pl,jkl)*X3)^2/noiseVar(jkl);
        
        Lb1_num = max(A,B)+log(1+exp(-abs(B-A)));
        
        A = -abs(zd(pl,jkl)-HD(pl,jkl)*X2)^2/noiseVar(jkl);
        B = -abs(zd(pl,jkl)-HD(pl,jkl)*X4)^2/noiseVar(jkl);
        
        Lb1_den = max(A,B)+log(1+exp(-abs(B-A)));
        
        A = -abs(zd(pl,jkl)-HD(pl,jkl)*X1)^2/noiseVar(jkl);
        B = -abs(zd(pl,jkl)-HD(pl,jkl)*X2)^2/noiseVar(jkl);
        
        Lb2_num = max(A,B)+log(1+exp(-abs(B-A)));
        
        A = -abs(zd(pl,jkl)-HD(pl,jkl)*X3)^2/noiseVar(jkl);
        B = -abs(zd(pl,jkl)-HD(pl,jkl)*X4)^2/noiseVar(jkl);
        
        Lb2_den = max(A,B)+log(1+exp(-abs(B-A)));
        
        Lb1 = Lb1_num-Lb1_den;
        Lb2 = Lb2_num-Lb2_den;
        
        %LR1 = [Lb1; Lb2];
        
        %LR = cat(1, LR, LR1);
        LR(LRcount1, jkl) = Lb1;
        LR(LRcount2, jkl) = Lb2;
        LRcount1  = LRcount1 + 2;
        LRcount2 = LRcount2 + 2;
 
    end

end



NAME = './5G_LDPC_M10_N20_Z142_Q2_nonVer.txt';
[address, LDPC_INFOLEN] = ldpc_mex_initial_CAPI([1420,2840,2],NAME);

for i = 2:21

    LR_in_de = zeros(length(LR(:,i)),1);
    LR_in_de(INTRLVR) = LR(:,i);
    %LR_in_de_2 = zeros(length(LR(:,i)),1);
    %LR_in_de_2(INTRLVR) = Le_OUT(:,i);
    APP_code(:,i) = ldpcDecoder_CAPI(address,LR_in_de);
    %APP_code_2(:,i) = ldpcDecoder_CAPI(address,LR_in_de_2);
    
    % Hard Decision
    
    est_code(:,i) = (APP_code(:,i)<0);
    bec(:,i) = sum(abs(est_code(:,i)-CODE(:,i)));

    % est_code_2 = (APP_code_2(:,i)<0);
    % bec_2(:,i) = sum(abs(est_code_2-CODE(:,i)));

end

BER = sum(bec) / 2840 / 20;

% figure; plot(1:length(LR_in_de_1), LR_in_de_1, '-r',...
%              1:length(LR_in_de_2), LR_in_de_2, '--b');

