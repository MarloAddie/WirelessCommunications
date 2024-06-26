clc, clear, close all

load('OFDM_PILOT.mat')
load('ofdm_map.mat')
%%
load('benchmark_NoiseVar_172648_1.mat')
load('benchmark_Zw_172648_1.mat')
%%

k = 2048;
L = 199;
kp = 512;
null_sub = 112;

idx = find(ofdm_map ==1);

z_w = bb_rece_data_172648_1474;
Z_pw = z_w(ofdm_map==1);

for nn = 1:kp
    for n = 1:L+1
        V(nn,n) = exp((-1i*2*pi*(n-1)*(idx(nn)-1))/k);
    end
end

D = OFDM_PILOT(ofdm_map==1);

D = diag(D);
hls = (1/kp)*V'*D'*Z_pw;

for nn = 1:k
    for n = 1:L+1
        V2(nn,n) = exp((-1i*2*pi*(n-1)*(nn-1))/k);
    end
end

Hw = V2*hls;

Z_pw_null = z_w((ofdm_map == 0));

for kl =1:21
    
    znw = z_w(:,kl);
    Z_pw_null = znw((ofdm_map == 0));
    noiseVar(kl) = (1/null_sub)*sum(abs(Z_pw_null).^2);
    
end



sd = find(ofdm_map == 2);

zd = z_w(sd-1);

HD = Hw(sd-1);

X1 = 1/sqrt(2) + (1/sqrt(2))*1i;
X2 = -1/sqrt(2) + (1/sqrt(2))*1i;
X3 = 1/sqrt(2) - (1/sqrt(2))*1i;
X4 = -1/sqrt(2) - (1/sqrt(2))*1i;

LR = [];
for pl = 1:length(HD)
    
    A = -abs(zd(pl)-HD(pl)*X1)^2/noiseVar(1);
    B = -abs(zd(pl)-HD(pl)*X3)^2/noiseVar(1);
    
    Lb1_num = max(A,B)+log(1+exp(-abs(B-A)));
    
    A = -abs(zd(pl)-HD(pl)*X2)^2/noiseVar(1);
    B = -abs(zd(pl)-HD(pl)*X4)^2/noiseVar(1);
    
    Lb1_den = max(A,B)+log(1+exp(-abs(B-A)));
    
    A = -abs(zd(pl)-HD(pl)*X1)^2/noiseVar(1);
    B = -abs(zd(pl)-HD(pl)*X2)^2/noiseVar(1);
    
    Lb2_num = max(A,B)+log(1+exp(-abs(B-A)));
    
    A = -abs(zd(pl)-HD(pl)*X3)^2/noiseVar(1);
    B = -abs(zd(pl)-HD(pl)*X4)^2/noiseVar(1);
    
    Lb2_den = max(A,B)+log(1+exp(-abs(B-A)));
    
    Lb1 = Lb1_num-Lb1_den;
    Lb2 = Lb2_num-Lb2_den;
    
    LR1 = [Lb1; Lb2];
    
    LR = cat(1, LR, LR1);
 
    
    
    
    
%     Lb1 = log((exp(-norm(zd(pl)-HD(pl)*X1)^2/noiseVar)+exp(-norm(-abs(zd(pl)-HD(pl)*X3)^2/noiseVar)))...
%         -log((exp(-norm(zd(pl)-HD(pl)*X2)^2/noiseVar)+exp(-norm(zd(pl)-HD(pl)*X4)^2/noiseVar))));
%     
%     Lb2 = log((exp(-norm(zd(pl)-HD(pl)*X1)^2/noiseVar)+exp(-norm(zd(pl)-HD(pl)*X2)^2/noiseVar))...
%         -log((exp(-norm(zd(pl)-HD(pl)*X3)^2/noiseVar)+exp(-norm(zd(pl)-HD(pl)*X4)^2/noiseVar))));
%     
%     LR1 = cat(1, Lb1, Lb2);
%     
%     LR = cat(1, LR, LR1);
%     


end

NAME = './5G_LDPC_M10_N20_Z142_Q2_nonVer.txt';
[address, LDPC_INFOLEN] = ldpc_mex_initial_CAPI([1420,2840,2],NAME);





