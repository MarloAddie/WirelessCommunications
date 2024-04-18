clc, clear, close all

addpath('C:\Users\marlo\Documents\ECE6961\WirelessCommunications\Project3files\Matlab_Python_LDPC_update\Matlab_Python_LDPC_update\Matlab_Python_LDPC\Matlab_version\mex');


load('OFDM_PILOT.mat')
load('ofdm_map.mat')
%%
load('benchmark_NoiseVar_172648_1.mat')
load('benchmark_intermediate_172648_from_1_single_hydrophones.mat');
load('benchmark_Zw_172648_1.mat')
load('INTRLVR.mat')
load('CODE.mat')
%%

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

for i = 1:21

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

% figure; plot(1:length(LR_in_de_1), LR_in_de_1, '-r',...
%              1:length(LR_in_de_2), LR_in_de_2, '--b');


%addpath './mex'
ldpc_checkLen = 128*10;
ldpc_codeLen = 128*20;
ldpc_qnumber = 2.0;
code_name = './forClass/5G_LDPC_M10_N20_Z128_Q2_nonVer.txt';
% Initialization of parameters and H matrix
%[address,ldpc_infoLen] = ldpc_mex_initial_CAPI([ldpc_checkLen ldpc_codeLen ldpc_qnumber], code_name);
[address, ldpc_infoLen] = ldpc_mex_initial_CAPI([1420,2840,2],NAME);


SNR_dB = [0.01 0.05 0.1 0.5 1 2 4];
min_wec = 100;
%mesg = APP_code(:,1);
mesg = est_code(:,2);
for i=1:length(SNR_dB)
    wec = 0; % word error number
    bec = 0; % bit error number
    bec_mesg = 0; % mesg bit error number
    tot = 0; % total block number
    while(wec<min_wec && tot < 10000)
        % generated message
        info_len = ldpc_infoLen;
        mesg = double(rand(info_len,1)>0.5);
        code = ldpcEncoder_CAPI(address, mesg);
        mod_code = 1 - code*2;
        noise_var = 10^(-SNR_dB(i)/10);
        r = mod_code + sqrt(noise_var)*randn(size(mod_code));
        prior = 2*r/noise_var;
        out_APP_c = ldpcDecoder_CAPI(address, prior);
        hard_code = (out_APP_c<0);
        out_APP_m = ldpc_Ext_CAPI(address, out_APP_c);
        hard_mesg = (out_APP_m<0);
        bit_mesg_err = sum(abs(hard_mesg-mesg));
        bit_err = sum(abs(hard_code-code));
        bec = bec + bit_err;
        bec_mesg = bec_mesg + bit_mesg_err;
        if bit_err ~= 0
            wec = wec + 1;
        end
        tot = tot + 1;
    end
    ber = bec / tot / length(code);
    ber_mesg = bec_mesg / tot / length(mesg);
    wer = wec / tot;
    fprintf('SNR(dB) = %f, tot = %f, ber = %f, ber_mesg = %f, wer = %f\n', SNR_dB(i), tot, ber, ber_mesg, wer);
end




